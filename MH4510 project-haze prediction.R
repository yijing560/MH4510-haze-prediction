library(tidyverse)
library(lubridate)
library(tidyverse)
library(caret)
library(glmnet)
library(Hmisc) # For Lag() function
# ==============================================================================
# 1. LOAD RAW DATASETS
# ==============================================================================
psi_raw    <- read_csv("Historical1hrPM2.5.csv")
wind_dir   <- read_csv("HistoricalWindDirectionacrossSingapore2019.csv")
wind_spd   <- read_csv("HistoricalWindSpeedacrossSingapore2019.csv")
fire_raw   <- read_csv("fire_archive_M-C61_825464.csv")

# ==============================================================================
# 2. CLEAN & PROCESS NASA FIRMS FIRE HOTSPOTS
# ==============================================================================
fire_daily <- fire_raw %>%
  # Filter satellite detection confidence >= 30% to remove false alarms
  filter(confidence >= 30) %>%
  mutate(
    Date = ymd(acq_date),
    # Segment fires into Sumatra vs. Kalimantan zones using geographic bounding boxes
    source_zone = case_when(
      latitude >= -6 & latitude <= 6  & longitude >= 95  & longitude <= 106 ~ "Sumatra",
      latitude >= -4 & latitude <= 4  & longitude >= 108 & longitude <= 119 ~ "Kalimantan",
      TRUE ~ "Other"
    )
  ) %>%
  # Keep only the target regions
  filter(source_zone %in% c("Sumatra", "Kalimantan")) %>%
  # Aggregate to daily total Fire Radiative Power (FRP) and fire counts
  group_by(Date, source_zone) %>%
  summarise(
    daily_frp   = sum(frp, na.rm = TRUE),
    daily_count = n(),
    .groups     = "drop"
  ) %>%
  # Pivot wider so Sumatra and Kalimantan become separate numeric feature columns
  pivot_wider(
    names_from  = source_zone,
    values_from = c(daily_frp, daily_count),
    values_fill = 0
  )

# ==============================================================================
# 3. CLEAN & PROCESS WIND DATA (VECTOR DECOMPOSITION)
# ==============================================================================
target_stations <- c("S50", "S24", "S109", "S111", "S104") # Representative Singapore stations

# Standardize Wind Direction (Degrees)
clean_dir <- wind_dir %>%
  filter(station_id %in% target_stations) %>%
  mutate(
    Timestamp = floor_date(ymd_hms(timestamp), unit = "hour"),
    wind_dir_deg = as.numeric(reading_value)
  ) %>%
  group_by(Timestamp) %>%
  summarise(mean_wind_dir = mean(wind_dir_deg, na.rm = TRUE), .groups = "drop")

# Standardize Wind Speed (Knots)
clean_spd <- wind_spd %>%
  filter(station_id %in% target_stations) %>%
  mutate(
    Timestamp = floor_date(ymd_hms(timestamp), unit = "hour"),
    wind_spd_knot = as.numeric(reading_value)
  ) %>%
  group_by(Timestamp) %>%
  summarise(mean_wind_spd = mean(wind_spd_knot, na.rm = TRUE), .groups = "drop")

# Merge wind tables and convert circular angles into Cartesian u and w vectors
wind_hourly <- inner_join(clean_dir, clean_spd, by = "Timestamp") %>%
  mutate(
    rad    = mean_wind_dir * pi / 180,
    # u component: East-West wind velocity vector
    wind_u = -mean_wind_spd * sin(rad),
    # w component: North-South wind velocity vector (positive = Southerly winds pushing towards SG)
    wind_w = -mean_wind_spd * cos(rad)
  )

# 4. CLEAN AIR QUALITY DATA & BUILD THE MASTER TABLE
# ==============================================================================
master_haze_2019 <- psi_raw %>%
  rename(datetime_str = 1) %>%
  mutate(
    Timestamp = dmy_hm(datetime_str),
    Date = as_date(Timestamp)
  ) %>%
  # Filter strictly to the 2019 peak Southwest Monsoon haze window (Aug 1 - Oct 31)
  filter(year(Timestamp) == 2019, month(Timestamp) %in% 8:10) %>%
  select(Timestamp, Date, north, south, east, west, central) %>%
  # Join hourly wind vectors
  left_join(wind_hourly, by = "Timestamp") %>%
  # Join daily Indonesian fire data
  left_join(fire_daily, by = "Date") %>%
  # Replace missing fire days with zero
  replace_na(list(
    daily_frp_Sumatra = 0, daily_frp_Kalimantan = 0,
    daily_count_Sumatra = 0, daily_count_Kalimantan = 0
  )) %>%
  # Create physics-based advection interaction term: Southerly wind * Sumatra Fire Power
  mutate(
    smoke_transport_Sumatra = pmax(0, wind_w) * daily_frp_Sumatra
  ) %>%
  arrange(Timestamp)

# Check final structure of your clean master data frame
dim(master_haze_2019)
head(master_haze_2019)



# ==============================================================================
# PART 1: DATA PREPROCESSING (TIMESTAMP CORRECTION) & EDA
# ==============================================================================
# A. Convert wide regional data to long format for EDA time-series plotting
long_data <- master_haze_2019 %>%
  select(Timestamp, north, south, east, west, central) %>%
  pivot_longer(north:central, names_to = "location", values_to = "PSI_Level")

# Inspect a random sample of the long data
long_data %>% sample_n(5)

# B. EDA Multi-Series Line Plot (Regional PSI trends over time)
ggplot(data = long_data, aes(x = Timestamp, y = PSI_Level, group = location, colour = location)) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Exploratory Data Analysis: Regional PM2.5/PSI Trends Across Singapore",
    subtitle = "Comparing multi-station air quality observations over time",
    x = "Timestamp",
    y = "PSI Level",
    colour = "Region"
  ) +
  theme_minimal()


# ==============================================================================
# PART 2: PHASE 2 - ELASTIC NET BASELINE MODEL & LAGGING
# ==============================================================================

# A. Create Lagged Features (Using built-in dplyr::lag)
lagged_obs <- function(df, L) {
  if (L == 0) {
    return(df)
  }
  df %>%
    transmute(across(c(north, south, east, west, central),
                     ~lag(., L),
                     .names = "Lag_{col}_{L}"))
}

# Generate 1 to 5 hour lags across all regional columns
PSI_with_lags <- c(0, 1:5) %>%
  map_dfc(function(L) lagged_obs(master_haze_2019, L)) %>%
  select(west, starts_with("Lag")) %>%
  drop_na()

# B. Custom Train-Test Split (Using runif with set.seed for reproducibility)
set.seed(42)
ind <- runif(nrow(PSI_with_lags)) <= 0.7
train_data <- PSI_with_lags %>% filter(ind)
test_data  <- PSI_with_lags %>% filter(!ind)

cat("Dimensions of training set:", dim(train_data), "\n")
cat("Dimensions of testing set:", dim(test_data), "\n")

# C. Train Elastic Net Model using caret (5-fold CV, scaling, hyperparameter grid)
set.seed(42)
elastic_net_model <- train(
  west ~ ., 
  data = train_data, 
  method = "glmnet",
  tuneGrid = expand.grid(
    alpha = seq(from = 0, to = 1, length = 5), 
    lambda = 10^(seq(from = -2, to = 2, length = 5))
  ),
  trControl = trainControl(method = "cv", number = 5),
  preProcess = c("scale")
)

# Print tuning summary and best hyperparameters found via cross-validation
print(elastic_net_model)
print(elastic_net_model$bestTune)

# D. Evaluate on Out-of-Sample Test Set
predictions <- predict(elastic_net_model, newdata = test_data)

test_rmse <- RMSE(predictions, test_data$west)
test_mae  <- MAE(predictions, test_data$west)
test_r2   <- cor(predictions, test_data$west)^2

cat("--- Out-of-Sample Test Performance ---\n")
cat("Test RMSE:", test_rmse, "\n")
cat("Test MAE:", test_mae, "\n")
cat("Test R-squared:", test_r2, "\n")

# E. Extract Final Model Coefficients
best_alpha <- elastic_net_model$bestTune$alpha
best_lambda <- elastic_net_model$bestTune$lambda

en_coefs <- coef(elastic_net_model$finalModel, s = best_lambda)
print(en_coefs)

# F. Visualisation: Coefficient Paths across Lambda (Lab 4 Style)
lambda_for_plotting <- 10^seq(from = -2, to = 2, length = 100)

en_path_model <- glmnet(
  x = as.matrix(scale(train_data %>% select(-west))),
  y = train_data$west,
  alpha = best_alpha,
  lambda = lambda_for_plotting
)

en_path_df <- as.matrix(coef(en_path_model, lambda_for_plotting)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(lambda = lambda_for_plotting)

en_path_df %>%
  select(-`(Intercept)`) %>%
  pivot_longer(-lambda, names_to = "variable", values_to = "coef") %>%
  ggplot(aes(x = lambda, y = coef, colour = variable)) +
  geom_line() +
  scale_x_log10() +
  labs(
    title = "Elastic Net Regularization Coefficient Paths",
    subtitle = sprintf("Optimal alpha = %.2f", best_alpha),
    x = "Lambda (log scale)",
    y = "Coefficient Estimate",
    colour = "Predictor"
  ) +
  theme_minimal()