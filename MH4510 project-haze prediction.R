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
ggsave("EDA_plot.png", width = 5, height = 3, dpi = 300)
# Create a clean summary table for regional PSI levels
regional_summary <- long_data %>% 
  group_by(location) %>% 
  summarise(
    Mean_PSI = round(mean(PSI_Level, na.rm = TRUE), 1),
    Max_PSI = max(PSI_Level, na.rm = TRUE)
  )
print(regional_summary)
# ==============================================================================
# PART 2: PHASE 2 - ELASTIC NET BASELINE MODEL & LAGGING (TIME-BASED SPLIT)
# ==============================================================================

# A. Create Lagged Features (Using built-in dplyr::lag)
split_date <-ymd("2019-10-21")
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
  select(Timestamp, west, starts_with("Lag")) %>%
  drop_na()

# B. Chronological Time-Based Train-Test Split (No random set.seed needed)
# Training window: August 1 to October 20
# Testing window: October 21 to October 31
train_data <- PSI_with_lags %>% 
  filter(Timestamp < ymd("2019-10-21")) %>% 
  select(-Timestamp) # Remove timestamp column before modeling matrix

test_data <- PSI_with_lags %>% 
  filter(Timestamp >= ymd("2019-10-21")) %>% 
  select(-Timestamp)

cat("Dimensions of training set:", dim(train_data), "\n")
cat("Dimensions of testing set:", dim(test_data), "\n")

# C. Train Elastic Net Model using caret (5-fold CV, scaling, hyperparameter grid)
set.seed(42) # Keep seed for cross-fold partitioning reproducibility only
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

cat("--- Out-of-Sample Test Performance (Time-Based Split) ---\n")
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
    subtitle = sprintf("Optimal alpha = %.2f (Chronological Split)", best_alpha),
    x = "Lambda (log scale)",
    y = "Coefficient Estimate",
    colour = "Predictor"
  ) +
  theme_minimal()

# ==============================================================================
# PHASE 3
# ==============================================================================
# install.packages("randomForest")
library(randomForest)

region_cols <- c("north", "south", "east", "west", "central")
exog_cols   <- c("wind_u", "wind_w", "mean_wind_spd",
                 "daily_frp_Sumatra", "daily_frp_Kalimantan",
                 "daily_count_Sumatra", "daily_count_Kalimantan",
                 "smoke_transport_Sumatra")

build_horizon_data <- function(df, h, n_lags = 5) {
  lag_df <- map_dfc(0:n_lags, function(L) {
    df %>% transmute(across(all_of(region_cols),
                            ~ dplyr::lag(.x, L),
                            .names = paste0("Lag_{col}_", L)))
  })
  bind_cols(
    tibble(Timestamp = df$Timestamp, target = dplyr::lead(df$west, h)),
    lag_df,
    df %>% select(all_of(exog_cols))
  ) %>%
    drop_na()
}

fit_horizon <- function(h, seed = 42) {
  d           <- build_horizon_data(master_haze_2019, h)
  target_time <- d$Timestamp + hours(h)
  train       <- d[target_time <  split_date, ]
  test        <- d[target_time >= split_date, ]
  
  x_train <- train %>% select(-Timestamp, -target) %>% as.data.frame()
  x_test  <- test  %>% select(-Timestamp, -target) %>% as.data.frame()
  y_train <- train$target
  y_test  <- test$target
  p       <- ncol(x_train)
  
  mtry_grid <- sort(unique(pmax(1, round(c(sqrt(p), p / 5, p / 3, p / 2, p)))))
  oob_tbl <- map_dfr(mtry_grid, function(m) {
    set.seed(seed)
    rf_m <- randomForest(x = x_train, y = y_train, ntree = 500, mtry = m)
    tibble(horizon = h, mtry = m, oob_mse = tail(rf_m$mse, 1))
  })
  best_mtry <- oob_tbl$mtry[which.min(oob_tbl$oob_mse)]
  
  set.seed(seed)
  rf_final <- randomForest(x = x_train, y = y_train, ntree = 500,
                           mtry = best_mtry, importance = TRUE)
  rf_pred <- predict(rf_final, newdata = x_test)
  
  set.seed(seed)
  en_fit <- train(
    x = x_train, y = y_train, method = "glmnet",
    tuneGrid = expand.grid(alpha = seq(0, 1, length = 5),
                           lambda = 10^seq(-2, 2, length = 5)),
    trControl = trainControl(method = "cv", number = 5),
    preProcess = c("scale")
  )
  en_pred <- predict(en_fit, newdata = x_test)
  
  naive_pred <- x_test$Lag_west_0
  
  metrics <- function(pred, obs, name) {
    tibble(horizon = h, model = name,
           RMSE = RMSE(pred, obs), MAE = MAE(pred, obs),
           R2 = cor(pred, obs)^2)
  }
  results <- bind_rows(
    metrics(naive_pred, y_test, "Persistence (naive)"),
    metrics(en_pred,    y_test, "Elastic Net"),
    metrics(rf_pred,    y_test, "Random Forest")
  )

  imp <- importance(rf_final, type = 1) %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    rename(perm_importance = `%IncMSE`) %>%
    mutate(horizon = h)
  
  list(results = results, oob = oob_tbl, importance = imp,
       rf = rf_final, best_mtry = best_mtry,
       test_df = tibble(Timestamp = test$Timestamp, actual = y_test,
                        rf = rf_pred, en = en_pred, horizon = h))
}

horizons   <- c(1, 6)
phase3_fit <- map(horizons, fit_horizon)
names(phase3_fit) <- paste0("t+", horizons)

phase3_results <- map_dfr(phase3_fit, "results") %>%
  mutate(horizon = paste0("t+", horizon, " hr"))
print(phase3_results)

phase3_oob <- map_dfr(phase3_fit, "oob")
print(phase3_oob)
print(map_int(phase3_fit, "best_mtry"))

phase3_imp <- map_dfr(phase3_fit, "importance") %>%
  mutate(horizon = paste0("t+", horizon, " hr"))

phase3_imp %>%
  group_by(horizon) %>%
  slice_max(perm_importance, n = 15) %>%
  ungroup() %>%
  mutate(variable = paste(variable, horizon, sep = "___")) %>%
  ggplot(aes(x = perm_importance, y = reorder(variable, perm_importance))) +
  geom_col(fill = "steelblue") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  facet_wrap(~ horizon, scales = "free_y") +
  labs(title = "Random Forest Permutation Variable Importance (West PSI)",
       x = "% Increase in MSE when permuted", y = NULL) +
  theme_minimal()

phase3_group <- phase3_imp %>%
  mutate(
    group = case_when(
      str_detect(variable, "^Lag_west_")                    ~ "Local (west) lags",
      str_detect(variable, "^Lag_")                         ~ "Other-region lags",
      variable %in% c("wind_u", "wind_w", "mean_wind_spd")  ~ "Wind",
      TRUE                                                  ~ "Fire / smoke transport"
    ),
    imp_pos = pmax(perm_importance, 0)
  ) %>%
  group_by(horizon, group) %>%
  summarise(total = sum(imp_pos), .groups = "drop") %>%
  group_by(horizon) %>%
  mutate(share = total / sum(total)) %>%
  ungroup()
print(phase3_group)

ggplot(phase3_group, aes(x = horizon, y = share, fill = group)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Share of Permutation Importance by Predictor Group",
       subtitle = "Does local persistence or upstream burning dominate as horizon lengthens?",
       x = "Forecast horizon", y = "Share of importance", fill = "Group") +
  theme_minimal()

map_dfr(phase3_fit, "test_df") %>%
  mutate(horizon = paste0("t+", horizon, " hr")) %>%
  ggplot(aes(x = Timestamp)) +
  geom_line(aes(y = actual, colour = "Actual"), alpha = 0.7) +
  geom_line(aes(y = rf, colour = "Random Forest"), alpha = 0.7) +
  facet_wrap(~ horizon, ncol = 1) +
  labs(title = "Random Forest Forecasts vs Actual West PSI (Test Period)",
       x = "Timestamp", y = "PSI Level", colour = NULL) +
  theme_minimal()



# ==============================================================================
# PHASE 5: SPATIAL AUTOREGRESSIVE MODEL (SAR)
# ==============================================================================

# install.packages(c("splm", "spdep"))
library(splm)
library(spdep)

# 1. Define spatial weights matrix
regions <- c("north", "south", "east", "west", "central")

W <- matrix(0, 5, 5, dimnames = list(regions, regions))

W["north", c("east", "west", "central")] <- 1
W["south", c("east", "west", "central")] <- 1
W["east", c("north", "south", "central")] <- 1
W["west", c("north", "south", "central")] <- 1
W["central", c("north", "south", "east", "west")] <- 1

W <- W / rowSums(W)
listw <- mat2listw(W, style = "W")

print(W)

# 2. Prepare panel data for each forecast horizon
build_sar_data <- function(h) {

  d <- master_haze_2019 %>%
    arrange(Timestamp) %>%
    mutate(
      future_time = lead(Timestamp, h),
      across(all_of(regions), ~ lead(.x, h),
             .names = "future_{.col}")
    ) %>%
    filter(future_time == Timestamp + hours(h)) %>%
    drop_na()

  map_dfr(regions, function(r) {
    tibble(
      Timestamp = d$Timestamp,
      target_time = d$future_time,
      region = r,
      target = d[[paste0("future_", r)]],
      pm25 = d[[r]],
      wind_u = d$wind_u,
      wind_w = d$wind_w,
      fire_frp = d$daily_frp_Sumatra
    )
  }) %>%
    mutate(region = factor(region, levels = regions)) %>%
    arrange(region, Timestamp)
}

# 3. Estimate OLS and SAR models
fit_sar <- function(h) {

  d <- build_sar_data(h)

  train_sp <- d %>%
    filter(target_time < ymd("2019-10-21"))

  test_sp <- d %>%
    filter(target_time >= ymd("2019-10-21"))

  formula_sp <- target ~ pm25 + wind_u + wind_w +
    fire_frp + region

  # Baseline OLS
  ols <- lm(formula_sp, data = train_sp)

  # Spatial autoregressive model (Maximum Likelihood)
  sar <- spml(
    formula_sp,
    data = train_sp,
    index = c("region", "Timestamp"),
    listw = listw,
    model = "pooling",
    lag = TRUE,
    spatial.error = "none"
  )

  # 4. Predict with OLS
  pred_ols <- predict(ols, newdata = test_sp)

  # 5. Predict with SAR using spatial multiplier
  rho <- as.numeric(sar$arcoef)

  pred_sar <- test_sp %>%
    group_split(Timestamp) %>%
    map_dfr(function(g) {

      g <- arrange(g, region)

      X <- model.matrix(formula_sp, data = g)
      xb <- as.vector(X %*% sar$coefficients)

      g$prediction <- as.vector(
        solve(diag(5) - rho * W, xb)
      )

      g
    })

  # 6. Evaluate forecasts
  results <- tibble(
    Horizon = paste0("t+", h, " hr"),
    Model = c("OLS", "SAR"),
    RMSE = c(
      RMSE(pred_ols, test_sp$target),
      RMSE(pred_sar$prediction, pred_sar$target)
    ),
    MAE = c(
      MAE(pred_ols, test_sp$target),
      MAE(pred_sar$prediction, pred_sar$target)
    )
  )

  list(
    results = results,
    rho = rho,
    sar_model = sar
  )
}

# 7. Run 1-hour and 6-hour models
phase5_fit <- map(c(1, 6), fit_sar)

phase5_results <- map_dfr(phase5_fit, "results")
print(phase5_results)

# 8. Examine spatial autoregressive coefficients
map_dbl(phase5_fit, "rho")

summary(phase5_fit[[1]]$sar_model)
summary(phase5_fit[[2]]$sar_model)

# 9. Visualize forecasting performance
ggplot(phase5_results,
       aes(x = Horizon, y = RMSE, fill = Model)) +
  geom_col(position = "dodge") +
  labs(
    title = "OLS vs Spatial Autoregressive Model",
    subtitle = "Out-of-sample forecasting performance",
    x = "Forecast Horizon",
    y = "Test RMSE"
  ) +
  theme_minimal()
