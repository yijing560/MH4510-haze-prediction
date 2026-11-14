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
# PHASE 3: Random Forest
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
  
  train <- d[
    target_time >= as.POSIXct("2019-08-01 00:00:00",
                              tz = "Asia/Singapore") &
      target_time < as.POSIXct("2019-10-21 00:00:00",
                               tz = "Asia/Singapore"),
  ]
  
  test <- d[
    target_time >= as.POSIXct("2019-10-21 00:00:00",
                              tz = "Asia/Singapore") &
      target_time < as.POSIXct("2019-11-01 00:00:00",
                               tz = "Asia/Singapore"),
  ]
  
  x_train <- train %>% select(-Timestamp, -target) %>% as.data.frame()
  x_test  <- test  %>% select(-Timestamp, -target) %>% as.data.frame()
  y_train <- train$target
  y_test  <- test$target
  cat("\nForecast horizon:", h, "hour(s)\n")
  
  cat("Training period:",
      min(train$target_time), "to",
      max(train$target_time), "\n")
  
  cat("Testing period:",
      min(test$target_time), "to",
      max(test$target_time), "\n")
  
  cat("Training observations:", nrow(train), "\n")
  cat("Testing observations:", nrow(test), "\n")
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
    x = x_train,
    y = y_train,
    method = "glmnet",
    tuneGrid = expand.grid(
      alpha = seq(0, 1, length = 5),
      lambda = 10^seq(-2, 2, length = 5)
    ),
    trControl = trainControl(
      method = "timeslice",
      initialWindow = floor(0.7 * nrow(x_train)),
      horizon = 24,
      fixedWindow = FALSE
    ),
    preProcess = c("scale")
  )
  
  en_pred <- predict(en_fit, newdata = x_test)
  
  naive_pred <- x_test$Lag_west_0
  
  metrics <- function(pred, obs, name) {
    tibble(
      horizon = h,
      model = name,
      RMSE = sqrt(mean((obs - pred)^2)),
      MAE = mean(abs(obs - pred)),
      R2 = 1 - sum((obs - pred)^2) /
        sum((obs - mean(obs))^2)
    )
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
  labs(
    title = "Random Forest Forecasts vs Actual West PSI", subtitle = "Testing period: 21–31 October 2019",
    x = "Target Timestamp", y = "PSI Level", colour = NULL) +
  theme_minimal()
                   
# ==============================================================================
# PHASE 4 -  VECTOR AUTOREGRESSIVE EXOGENOUS (VARX) SYSTEM ----
# ==============================================================================
# install.packages("vars")
# install.packages("tseries")
library(vars)
library(tseries)

# 1. Train/Test sets (Chronological Split) ----

# Training set: August 1 to October 20, 2019
# Test set: October 21 to October 31, 2019
train_df <- master_haze_2019 %>% filter(Timestamp < as.POSIXct("2019-10-21 00:00:00", tz = "UTC"))
test_df  <- master_haze_2019 %>% filter(Timestamp >= as.POSIXct("2019-10-21 00:00:00", tz = "UTC"))
cat("VARX Training observations:", nrow(train_df), "\n")
cat("VARX Testing observations: ", nrow(test_df), "\n")
region_cols <- c("north", "south", "east", "west", "central")
exog_cols   <- c("wind_u", "wind_w", "daily_frp_Sumatra", "smoke_transport_Sumatra")

# Extract Endogenous Matrix Y_t and Exogenous Matrix X_t
Y_train <- as.matrix(train_df[, region_cols])
X_train <- as.matrix(train_df[, exog_cols])
Y_full  <- as.matrix(master_haze_2019[, region_cols])
X_full  <- as.matrix(master_haze_2019[, exog_cols])

# 2. Stationarity Checks (Augmented Dickey-Fuller Test) ----

apply(Y_train, 2, function(x) {
  adf_res <- adf.test(x)
  return(c(p_value = adf_res$p.value))
  })

# 3. Optimal Lag Selection & VARX Fit ----

# Determine optimal lag length (p) using AIC/BIC
lag_selection <- VARselect(Y_train, lag.max = 10, type = "const", exogen = X_train)
print(lag_selection$selection)
optimal_p <- lag_selection$selection["AIC(n)"]
cat("Selected Optimal Lag p (AIC):", optimal_p, "\n")

# Fit the VARX(p) System
varx_fit <- VAR(Y_train, p = optimal_p, type = "const", exogen = X_train)
summary(varx_fit) # Fits linear regression equations for each region

# 4. Granger Causality Tests ----
# (To test directional transmission across regions under prevailing winds, 
#   specifically South and West due to the Southwest Monsoon)

for (reg in region_cols) {
  cat("\nGranger Causality - Cause:", reg, "\n")
  print(causality(varx_fit, cause = reg)$Granger)
}

# 5. Impulse Response Functions (IRF Simulation) ----
# (Traces how an incoming smoke shock at the South/West diffuses across Singapore)

# South
irf_south <- irf(
  varx_fit,
  impulse = "south",
  response = c("north", "south", "east", "west", "central"),
  n.ahead = 24, # 24-hour diffusion window
  boot = TRUE,
  ci = 0.95
)
plot(irf_south, main = "Island-Wide PM2.5 Response to South Shock (24-Hour Horizon)")
# West
irf_west <- irf(
  varx_fit,
  impulse = "west",
  response = c("north", "south", "east", "west", "central"),
  n.ahead = 24, # 24-hour diffusion window
  boot = TRUE,
  ci = 0.95
)
plot(irf_west, main = "Island-Wide PM2.5 Response to West Shock (24-Hour Horizon)")

# 6. Forecasting & Evaluation on test set (t+1 and t+6) ----

# Extract coefficient matrices from fitted VARX(p) model:
# Y_t = c + A_1 Y_{t-1} + ... + A_p Y_{t-p} + B X_t
A_mats <- Acoef(varx_fit) # List of p (5x5) endogenous lag coefficient matrices
B_mat  <- Bcoef(varx_fit) # (5 x (5p + 1 + k_exog)) full coefficient matrix

# Separate intercept and exogenous coefficient matrices
const_vec <- B_mat[, "const"]
exog_mat  <- B_mat[, exog_cols, drop = FALSE]

# Helper function to recursively forecast h steps ahead from every origin t in test window
forecast_varx_horizon <- function(h) {
  # Indices in master_haze_2019 corresponding to target_time in Oct 21 - Oct 31
  target_indices <- which(
    master_haze_2019$Timestamp >= as.POSIXct("2019-10-21 00:00:00", tz = "UTC") &
      master_haze_2019$Timestamp <  as.POSIXct("2019-11-01 00:00:00", tz = "UTC")
  )
  
  preds_west  <- numeric(length(target_indices))
  actual_west <- numeric(length(target_indices))
  target_ts   <- master_haze_2019$Timestamp[target_indices]
  
  for (i in seq_along(target_indices)) {
    idx_target <- target_indices[i]
    idx_origin <- idx_target - h # Information available up to time t
    # History buffer of last p observed Y vectors: Y_{t}, Y_{t-1}, ..., Y_{t-p+1}
    Y_hist <- Y_full[(idx_origin - optimal_p + 1):idx_origin, , drop = FALSE]
    # Recursively step forward from step 1 to step h
    for (step in 1:h) {
      idx_step <- idx_origin + step
      x_step   <- X_full[idx_step, ]
      # Start with intercept + exogenous effect
      y_next <- const_vec + as.vector(exog_mat %*% x_step)
      # Add autoregressive lag contributions: sum_{L=1}^p A_L * Y_{idx_step - L}
      for (L in 1:optimal_p) {
        y_lag  <- Y_hist[nrow(Y_hist) - L + 1, ]
        y_next <- y_next + as.vector(A_mats[[L]] %*% y_lag)
      }
      # Append predicted y_next to history buffer for multi-step recursion
      Y_hist <- rbind(Y_hist, y_next)
    }
    preds_west[i]  <- tail(Y_hist, 1)[, "west"]
    actual_west[i] <- Y_full[idx_target, "west"]
  }
  
  rmse <- sqrt(mean((preds_west - actual_west)^2))
  mae  <- mean(abs(preds_west - actual_west))
  r2   <- 1 - sum((actual_west - preds_west)^2) / sum((actual_west - mean(actual_west))^2)
  
  list(
    metrics = tibble(
      horizon = paste0("t+", h, " hr"),
      model   = "VARX",
      RMSE    = rmse,
      MAE     = mae,
      R2      = r2
    ),
    df = tibble(
      Timestamp = target_ts,
      Actual    = actual_west,
      VARX_Predicted = preds_west,
      horizon   = paste0("t+", h, " hr")
    )
  )
}

# Run t+1 and t+6 forecasts
varx_h1 <- forecast_varx_horizon(1)
varx_h6 <- forecast_varx_horizon(6)
phase4_results <- bind_rows(varx_h1$metrics, varx_h6$metrics)
print(phase4_results) # MAE is 4.45 and 6.86 for t+1 and t+6 respectively

# Plot Actual vs VARX Predicted PM2.5 for West Region (t+1 hr and t+6 hr)
eval_varx_df <- bind_rows(varx_h1$df, varx_h6$df)
ggplot(eval_varx_df, aes(x = Timestamp)) +
  geom_line(aes(y = Actual, colour = "Actual PM2.5"), alpha = 0.8) +
  geom_line(aes(y = VARX_Predicted, colour = "VARX Forecast"), linetype = "longdash") +
  facet_wrap(~ horizon, ncol = 1) +
  labs(
    title = "VARX Model: Out-of-Sample PM2.5 Forecasts (West Region)",
    subtitle = "Evaluation Horizon: October 21 – October 31, 2019 (t+1 hr vs t+6 hr)",
    x = "Target Timestamp",
    y = "PM2.5 Level (µg/m³)",
    colour = "Legend"
  ) +
  scale_colour_manual(values = c("Actual PM2.5" = "black", "VARX Forecast" = "red")) +
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
    }) %>% 
    arrange(region, Timestamp) # <--- ADD THIS LINE TO FIX THE SORT ORDER!
  
  # Ensure test_sp is in the exact same order for fair comparison
  test_sp_ordered <- test_sp %>% arrange(region, Timestamp)
  
  # 6. Evaluate forecasts using aligned data
  results <- tibble(
    Horizon = paste0("t+", h, " hr"),
    Model = c("OLS", "SAR"),
    RMSE = c(
      RMSE(predict(ols, newdata = test_sp_ordered), test_sp_ordered$target),
      RMSE(pred_sar$prediction, test_sp_ordered$target)
    ),
    MAE = c(
      MAE(predict(ols, newdata = test_sp_ordered), test_sp_ordered$target),
      MAE(pred_sar$prediction, test_sp_ordered$target)
    ),
    R2 = c(
      cor(predict(ols, newdata = test_sp_ordered), test_sp_ordered$target)^2,
      cor(pred_sar$prediction, test_sp_ordered$target)^2
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
