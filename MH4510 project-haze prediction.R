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
  filter(confidence >= 30) %>%
  mutate(
    Date = ymd(acq_date),
    source_zone = case_when(
      latitude >= -6 & latitude <= 6  & longitude >= 95  & longitude <= 106 ~ "Sumatra",
      TRUE ~ "Other"
    )
  ) %>%
  filter(source_zone == "Sumatra") %>%
  group_by(Date, source_zone) %>%
  summarise(
    daily_frp   = sum(frp, na.rm = TRUE),
    daily_count = n(),
    .groups     = "drop"
  ) %>%
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
eda_plot <- ggplot(data = long_data, aes(x = Timestamp, y = PSI_Level, group = location, colour = location)) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Exploratory Data Analysis: Regional PM2.5/PSI Trends Across Singapore",
    subtitle = "Comparing multi-station air quality observations over time",
    x = "Timestamp",
    y = "PSI Level",
    colour = "Region"
  ) +
  theme_minimal()
print(eda_plot)
ggsave("EDA_plot.png", plot = eda_plot, width = 5, height = 3, dpi = 300)
# Create a clean summary table for regional PSI levels
regional_summary <- long_data %>% 
  group_by(location) %>% 
  summarise(
    Mean_PSI = round(mean(PSI_Level, na.rm = TRUE), 1),
    Max_PSI = max(PSI_Level, na.rm = TRUE)
  )
print(regional_summary)
# ==============================================================================
# PART 2: PHASE 2 - ELASTIC NET MODEL (t+1 and t+6 Horizons, MAE Metric)
# ==============================================================================

split_date <- ymd("2019-10-21")

# A. Original lagged observation function for regional PSI readings
lagged_obs <- function(df, L) {
  if (L == 0) {
    return(df %>% select(north, south, east, west, central))
  }
  df %>%
    transmute(across(c(north, south, east, west, central),
                     ~lag(., L),
                     .names = "Lag_{col}_{L}"))
}

# Function to run Elastic Net for any forecast horizon h (t+1 and t+6)
fit_elastic_net_horizon <- function(h, seed = 42) {
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Running Elastic Net for Horizon: t+%d\n", h))
  cat(sprintf("========================================\n"))
  
  # Generate 0 to 5 hour regional lags and define target as lead(west, h)
  PSI_with_lags <- c(0, 1:5) %>%
    map_dfc(function(L) lagged_obs(master_haze_2019, L)) %>%
    mutate(
      Timestamp = master_haze_2019$Timestamp,
      target = lead(west, h)
    ) %>%
    select(Timestamp, target, starts_with("Lag")) %>%
    drop_na()
  
  # B. Chronological Time-Based Train-Test Split
  # Training window: August 1 to October 20
  # Testing window: October 21 to October 31
  train_data <- PSI_with_lags %>%  
    filter(Timestamp < split_date) %>%  
    select(-Timestamp) # Remove timestamp column before modeling matrix
  
  test_data <- PSI_with_lags %>%  
    filter(Timestamp >= split_date) %>%  
    select(-Timestamp)
  
  cat("Dimensions of training set:", dim(train_data), "\n")
  cat("Dimensions of testing set:", dim(test_data), "\n")
  
  # C. Train Elastic Net Model using caret (5-fold CV, scaling, hyperparameter grid)
  set.seed(seed) # Keep seed for cross-fold partitioning reproducibility only
  elastic_net_model <- train(
    target ~ .,  
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
  
  # D. Evaluate on Out-of-Sample Test Set (MAE metric standardization)
  predictions <- predict(elastic_net_model, newdata = test_data)
  test_mae  <- MAE(predictions, test_data$target)
  
  cat(sprintf("--- Out-of-Sample Test Performance (t+%d) ---\n", h))
  cat("Test MAE:", test_mae, "\n")
 
   if (h == 1) {
    model_t1 <<- elastic_net_model
    test_mae_enet_t1 <<- test_mae
  } else if (h == 6) {
    model_t6 <<- elastic_net_model
    test_mae_enet_t6 <<- test_mae
  }
  # E. Extract Final Model Coefficients
  best_alpha <- elastic_net_model$bestTune$alpha
  best_lambda <- elastic_net_model$bestTune$lambda
  
  en_coefs <- coef(elastic_net_model$finalModel, s = best_lambda)
  print(en_coefs)
  
  # F. Visualisation: Coefficient Paths across Lambda (Lab 4 Style)
  lambda_for_plotting <- 10^seq(from = -2, to = 2, length = 100)
  
  en_path_model <- glmnet(
    x = as.matrix(scale(train_data %>% select(-target))),
    y = train_data$target,
    alpha = best_alpha,
    lambda = lambda_for_plotting
  )
  
  en_path_df <- as.matrix(coef(en_path_model, lambda_for_plotting)) %>%
    t() %>%
    as.data.frame() %>%
    mutate(lambda = lambda_for_plotting)
  
  path_plot <- en_path_df %>%
    select(-`(Intercept)`) %>%
    pivot_longer(-lambda, names_to = "variable", values_to = "coef") %>%
    ggplot(aes(x = lambda, y = coef, colour = variable)) +
    geom_line() +
    scale_x_log10() +
    labs(
      title = sprintf("Elastic Net Regularization Coefficient Paths (t+%d hr)", h),
      subtitle = sprintf("Optimal alpha = %.2f (Chronological Split)", best_alpha),
      x = "Lambda (log scale)",
      y = "Coefficient Estimate",
      colour = "Predictor"
    ) +
    theme_minimal()
  
  print(path_plot)
  ggsave(
    sprintf("Elastic_Net_Coefficient_Paths_t%d.png", h),
    plot = path_plot, width = 10, height = 6, dpi = 300
  )
  
  tibble(horizon = paste0("t+", h, " hr"), model = "Elastic Net", MAE = test_mae)
}

# Run Elastic Net for both t+1 and t+6 horizons
horizons <- c(1, 6)
elastic_net_results <- map_dfr(horizons, fit_elastic_net_horizon)

print("\n--- Summary Table: Elastic Net Performance Across Horizons ---")
print(elastic_net_results)

# ==============================================================================
# PHASE 3: Random Forest
# ==============================================================================
# install.packages("randomForest")
# Training set: August 1 to October 20, 2019
# Test set: October 21 to October 31, 2019
library(randomForest)

region_cols <- c("north", "south", "east", "west", "central")

exog_cols <- c(
  "wind_u", "wind_w", "mean_wind_spd",
  "daily_frp_Sumatra",
  "daily_count_Sumatra",
  "smoke_transport_Sumatra"
)

phase3_train_start <- as.POSIXct("2019-08-01 00:00:00", tz = "UTC")
phase3_test_start  <- as.POSIXct("2019-10-21 00:00:00", tz = "UTC")
phase3_test_end    <- as.POSIXct("2019-11-01 00:00:00", tz = "UTC")

build_horizon_data <- function(df, h, n_lags = 5) {
  lag_df <- map_dfc(0:n_lags, function(L) {
    df %>%
      transmute(
        across(
          all_of(region_cols),
          ~ dplyr::lag(.x, L),
          .names = paste0("Lag_{col}_", L)
        )
      )
  })
  
  bind_cols(
    tibble(
      Timestamp = df$Timestamp,
      target = dplyr::lead(df$west, h)
    ),
    lag_df,
    df %>% select(all_of(exog_cols))
  ) %>%
    drop_na()
}

fit_horizon <- function(h, seed = 42) {
  d <- build_horizon_data(master_haze_2019, h)
  d$target_time <- d$Timestamp + lubridate::hours(h)

  train <- d %>%
    filter(
      target_time >= phase3_train_start,
      target_time < phase3_test_start
    )

  test <- d %>%
    filter(
      target_time >= phase3_test_start,
      target_time < phase3_test_end
    )
  
  if (nrow(train) == 0 || nrow(test) == 0) {
    stop(
      "Phase 3 produced an empty training or testing set for horizon t+",
      h,
      ". Check the Timestamp range and time zone."
    )
  }
  
  x_train <- train %>%
    select(-Timestamp, -target, -target_time) %>%
    as.data.frame()
  
  x_test <- test %>%
    select(-Timestamp, -target, -target_time) %>%
    as.data.frame()
  
  y_train <- train$target
  y_test  <- test$target
  
  cat("\nForecast horizon:", h, "hour(s)\n")
  cat(
    "Training target period:",
    format(min(train$target_time), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    "to",
    format(max(train$target_time), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    "\n"
  )
  cat(
    "Testing target period:",
    format(min(test$target_time), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    "to",
    format(max(test$target_time), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    "\n"
  )
  cat("Training observations:", nrow(train), "\n")
  cat("Testing observations:", nrow(test), "\n")
  
  p <- ncol(x_train)

  mtry_grid <- sort(unique(as.integer(
    pmax(1, pmin(p, round(c(sqrt(p), p / 5, p / 3, p / 2, p))))
  )))
  
  oob_tbl <- map_dfr(mtry_grid, function(m) {
    set.seed(seed)
    rf_m <- randomForest(
      x = x_train,
      y = y_train,
      ntree = 500,
      mtry = m
    )
    
    tibble(
      horizon = h,
      mtry = m,
      oob_mse = tail(rf_m$mse, 1)
    )
  })
  
  best_mtry <- oob_tbl$mtry[which.min(oob_tbl$oob_mse)]
  
  set.seed(seed)
  rf_final <- randomForest(
    x = x_train,
    y = y_train,
    ntree = 500,
    mtry = best_mtry,
    importance = TRUE
  )
  
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
  
  list(
    results = results,
    oob = oob_tbl,
    importance = imp,
    rf = rf_final,
    best_mtry = best_mtry,
    test_df = tibble(
      Timestamp = test$target_time,
      actual = y_test,
      rf = rf_pred,
      en = en_pred,
      horizon = h
    )
  )
}

horizons <- c(1, 6)
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

importance_plot <- phase3_imp %>%
  group_by(horizon) %>%
  slice_max(perm_importance, n = 15) %>%
  ungroup() %>%
  mutate(variable = paste(variable, horizon, sep = "___")) %>%
  ggplot(aes(x = perm_importance, y = reorder(variable, perm_importance))) +
  geom_col(fill = "steelblue") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  facet_wrap(~ horizon, scales = "free_y") +
  labs(
    title = "Random Forest Permutation Variable Importance (West PSI)",
    x = "% Increase in MSE when permuted",
    y = NULL
  ) +
  theme_minimal()

print(importance_plot)
ggsave("importance_plot.png", importance_plot, width = 7, height = 5, dpi = 300)

phase3_group <- phase3_imp %>%
  mutate(
    group = case_when(
      str_detect(variable, "^Lag_west_") ~ "Local (west) lags",
      str_detect(variable, "^Lag_") ~ "Other-region lags",
      variable %in% c("wind_u", "wind_w", "mean_wind_spd") ~ "Wind",
      variable %in% c(
        "daily_frp_Sumatra",
        "daily_count_Sumatra",
        "smoke_transport_Sumatra"
      ) ~ "Fire / smoke transport",
      TRUE ~ "Other"
    ),
    imp_pos = pmax(perm_importance, 0)
  ) %>%
  group_by(horizon, group) %>%
  summarise(
    total = sum(imp_pos),
    n_variables = n(),
    avg_importance = mean(imp_pos),
    .groups = "drop"
  ) %>%
  group_by(horizon) %>%
  mutate(
    share = if (sum(total) > 0) total / sum(total) else 0
  ) %>%
  ungroup()

print(phase3_group)

group_plot <- ggplot(
  phase3_group,
  aes(x = horizon, y = share, fill = group)
) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Share of Permutation Importance by Predictor Group",
    subtitle = "Comparing local persistence, wind, and Sumatra fire/smoke predictors",
    x = "Forecast horizon",
    y = "Share of importance",
    fill = "Group"
  ) +
  theme_minimal()

print(group_plot)
ggsave("group_plot.png", group_plot, width = 7, height = 4, dpi = 300)

phase3_forecast_df <- map_dfr(phase3_fit, "test_df") %>%
  mutate(horizon = paste0("t+", horizon, " hr"))

forecast_plot <- ggplot(phase3_forecast_df, aes(x = Timestamp)) +
  geom_line(aes(y = actual, colour = "Actual"), alpha = 0.7) +
  geom_line(aes(y = rf, colour = "Random Forest"), alpha = 0.7) +
  facet_wrap(~ horizon, ncol = 1) +
  labs(
    title = "Random Forest Forecasts vs Actual West PSI",
    subtitle = "Testing target period: 21-31 October 2019",
    x = "Target timestamp",
    y = "PSI Level",
    colour = NULL
  ) +
  theme_minimal()

print(forecast_plot)
ggsave("forecast_plot.png", forecast_plot, width = 7, height = 5, dpi = 300)
                   
# ==============================================================================
# SUPPORT VECTOR REGRESSION (SVM): adapted from mh4510_lab_06_svm.Rmd
# ==============================================================================
# The lab compares linear, polynomial and radial kernels with caret/kernlab.
# Here a numeric target invokes regression, with chronological resampling.
library(kernlab)
svm_region_cols <- c("north", "south", "east", "west", "central")
svm_exog_cols <- c("wind_u", "wind_w", "mean_wind_spd", "daily_frp_Sumatra",
                   "daily_count_Sumatra", "smoke_transport_Sumatra")

build_svm_data <- function(h) {
  lag_df <- map_dfc(0:5, function(L) {
    master_haze_2019 %>% transmute(across(all_of(svm_region_cols),
      ~ dplyr::lag(.x, L), .names = paste0("Lag_{col}_", L)))
  })
  bind_cols(tibble(Timestamp = master_haze_2019$Timestamp,
                   target_time = dplyr::lead(master_haze_2019$Timestamp, h),
                   target = dplyr::lead(master_haze_2019$west, h)),
            lag_df, master_haze_2019 %>% dplyr::select(all_of(svm_exog_cols))) %>%
    filter(target_time == Timestamp + hours(h)) %>% drop_na()
}

svm_grids <- list(
  Linear = list(method = "svmLinear", grid = expand.grid(C = 2^(-4:1))),
  Polynomial = list(method = "svmPoly",
    grid = expand.grid(C = c(0.5, 2), scale = c(1/64, 1/16), degree = c(2, 3))),
  Radial = list(method = "svmRadial",
    grid = expand.grid(sigma = c(1/16, 1/8, 1/4), C = c(0.5, 2)))
)

fit_svm_horizon <- function(h, seed = 42) {
  d <- build_svm_data(h)
  train <- d %>% filter(target_time >= phase3_train_start,
                        target_time < phase3_test_start)
  test <- d %>% filter(target_time >= phase3_test_start,
                       target_time < phase3_test_end)
  # Exactly the target observations used for the Phase 3 West forecasts.
  rf_test <- phase3_fit[[match(h, c(1, 6))]]$test_df
  stopifnot(identical(as.numeric(test$target_time), as.numeric(rf_test$Timestamp)),
            isTRUE(all.equal(test$target, rf_test$actual)))
  x_train <- train %>% dplyr::select(-Timestamp, -target_time, -target) %>% as.data.frame()
  x_test <- test %>% dplyr::select(-Timestamp, -target_time, -target) %>% as.data.frame()
  # Three expanding training windows followed by seven-day validation blocks.
  cutoffs <- as.POSIXct(c("2019-09-23", "2019-10-01", "2019-10-09"), tz = "UTC")
  index_out <- lapply(cutoffs, function(cut) which(train$target_time >= cut &
                                                  train$target_time < cut + days(7)))
  # Purge the boundary: every fitted target must precede the first validation origin.
  index <- lapply(index_out, function(ii) which(train$target_time < min(train$Timestamp[ii])))
  names(index) <- names(index_out) <- paste0("TimeBlock", seq_along(index))
  stopifnot(all(lengths(index) > 0), all(lengths(index_out) > 0))
  for (i in seq_along(index)) {
    stopifnot(max(train$target_time[index[[i]]]) < min(train$Timestamp[index_out[[i]]]))
  }
  control <- trainControl(method = "cv", number = 3, index = index, indexOut = index_out,
                          allowParallel = FALSE, returnResamp = "final", verboseIter = TRUE)
  candidates <- map(names(svm_grids), function(kernel) {
    spec <- svm_grids[[kernel]]
    cat("\nSVM horizon", h, "kernel", kernel, "\n")
    set.seed(seed)
    fit <- caret::train(x = x_train, y = train$target, method = spec$method,
      tuneGrid = spec$grid, trControl = control, metric = "MAE", maximize = FALSE,
      preProcess = c("zv", "center", "scale"))
    best <- fit$results
    for (parameter in names(fit$bestTune)) {
      best <- best[best[[parameter]] == fit$bestTune[[parameter]], , drop = FALSE]
    }
    list(kernel = kernel, fit = fit, cv_mae = best$MAE[1])
  })
  names(candidates) <- names(svm_grids)
  best_index <- which.min(map_dbl(candidates, "cv_mae"))
  winner <- candidates[[best_index]]
  # caret refits each tuned candidate on all training rows; only the CV winner
  # is evaluated on the untouched test set (no kernel selection on test errors).
  pred <- as.numeric(predict(winner$fit, newdata = x_test))
  stopifnot(length(pred) == nrow(test), all(is.finite(pred)))
  params <- winner$fit$bestTune
  param_text <- paste(paste(names(params), format(unlist(params), trim = TRUE), sep = "="),
                      collapse = ", ")
  tuning <- map_dfr(candidates, function(candidate) {
    as_tibble(candidate$fit$results) %>%
      mutate(horizon = paste0("t+", h, " hr"), kernel = candidate$kernel)
  })
  list(
    results = tibble(horizon = paste0("t+", h, " hr"), model = "SVM",
      RMSE = sqrt(mean((pred - test$target)^2)), MAE = mean(abs(pred - test$target)),
      R2 = 1 - sum((pred - test$target)^2) / sum((test$target - mean(test$target))^2)),
    selection = tibble(horizon = paste0("t+", h, " hr"), kernel = winner$kernel,
      Parameters = param_text, Validation_MAE = winner$cv_mae,
      Training_n = nrow(train), Testing_n = nrow(test)),
    tuning = tuning, model = winner$fit,
    kernel_summary = tibble(horizon = paste0("t+", h, " hr"),
      kernel = names(candidates), Validation_MAE = map_dbl(candidates, "cv_mae")),
    forecasts = tibble(Timestamp = test$target_time, Actual = test$target,
      SVM = pred, Persistence = test$Lag_west_0, horizon = paste0("t+", h, " hr")),
    resampling = list(index = index, indexOut = index_out)
  )
}

svm_fit <- map(c(1, 6), fit_svm_horizon)
svm_results <- map_dfr(svm_fit, "results")
svm_selection <- map_dfr(svm_fit, "selection")
svm_tuning <- map_dfr(svm_fit, "tuning")
svm_kernel_summary <- map_dfr(svm_fit, "kernel_summary")
svm_forecasts <- map_dfr(svm_fit, "forecasts")
print(svm_results)
print(svm_selection)
print(svm_kernel_summary)

svm_forecast_plot <- svm_forecasts %>%
  pivot_longer(c(Actual, SVM, Persistence), names_to = "Series", values_to = "PM25") %>%
  ggplot(aes(Timestamp, PM25, colour = Series, linetype = Series)) +
  geom_line(linewidth = 0.45) + facet_wrap(~ horizon, ncol = 1) +
  scale_colour_manual(values = c(Actual = "black", SVM = "#0072B2", Persistence = "#D55E00")) +
  scale_linetype_manual(values = c(Actual = "solid", SVM = "solid", Persistence = "dotted")) +
  labs(title = "SVM Forecasts vs Actual PM2.5: West", x = "Target timestamp",
       y = "PM2.5 (ug/m3)", colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
print(svm_forecast_plot)
ggsave("SVM_Forecast_vs_Actual.png", svm_forecast_plot, width = 10, height = 6, dpi = 300)

svm_kernel_plot <- ggplot(svm_kernel_summary, aes(kernel, Validation_MAE, colour = horizon)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  labs(title = "SVM Kernel Selection", subtitle = "Best mean validation MAE per kernel across three chronological blocks",
       x = "Kernel", y = "Validation MAE (ug/m3)", colour = "Horizon") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("SVM_Kernel_Validation.png", svm_kernel_plot, width = 8, height = 4.5, dpi = 300)
write_csv(svm_results, "SVM_Test_Metrics.csv")
write_csv(svm_selection, "SVM_Selected_Models.csv")
write_csv(svm_tuning, "SVM_Tuning_Results.csv")
write_csv(svm_forecasts, "SVM_Test_Forecasts.csv")
saveRDS(list(results = svm_results, selection = svm_selection, tuning = svm_tuning,
             kernel_summary = svm_kernel_summary, forecasts = svm_forecasts, fits = svm_fit),
        "SVM_Analysis_Outputs.rds")

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

# Save the base-R IRF plot for inclusion in the report.
png("South IRF.png", width = 8, height = 10, units = "in", res = 300)
plot(irf_south, main = "Island-Wide PM2.5 Response to South Shock (24-Hour Horizon)")
dev.off()

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

png("West IRF.png", width = 8, height = 10, units = "in", res = 300)
plot(irf_west, main = "Island-Wide PM2.5 Response to West Shock (24-Hour Horizon)")
dev.off()

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
varx_forecast_plot <- ggplot(eval_varx_df, aes(x = Timestamp)) +
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

print(varx_forecast_plot)
ggsave("VARX_Actual_Predicted_West.png", plot = varx_forecast_plot,
       width = 10, height = 7, dpi = 300)

# ==============================================================================
# PHASE 5: SPATIAL AUTOREGRESSIVE MODEL (SAR)
# ==============================================================================
# Fit the original fixed-neighbour SAR model by pooled Gaussian maximum
# likelihood in base R. No spatial-weight or other hyperparameter tuning.
regions <- c("north", "south", "east", "west", "central")
W <- matrix(0, 5, 5, dimnames = list(regions, regions))
W["north", c("east", "west", "central")] <- 1
W["south", c("east", "west", "central")] <- 1
W["east", c("north", "south", "central")] <- 1
W["west", c("north", "south", "central")] <- 1
W["central", c("north", "south", "east", "west")] <- 1
W <- W / rowSums(W)
sar_test_start <- as.POSIXct("2019-10-21", tz = "UTC")
formula_sp <- target ~ pm25 + wind_u + wind_w + fire_frp + region

build_sar_data <- function(h) {
  d <- master_haze_2019 %>%
    arrange(Timestamp) %>%
    mutate(future_time = lead(Timestamp, h),
           across(all_of(regions), ~ lead(.x, h), .names = "future_{.col}")) %>%
    filter(future_time == Timestamp + hours(h)) %>%
    drop_na(all_of(c(regions, paste0("future_", regions), "wind_u", "wind_w",
                     "daily_frp_Sumatra")))
  map_dfr(regions, function(r) {
    tibble(Timestamp = d$Timestamp, target_time = d$future_time,
           region = r, target = d[[paste0("future_", r)]], pm25 = d[[r]],
           wind_u = d$wind_u, wind_w = d$wind_w,
           fire_frp = d$daily_frp_Sumatra)
  }) %>%
    mutate(region = factor(region, levels = regions)) %>%
    arrange(region, Timestamp)
}

fit_pooled_sar <- function(d, weights) {
  X <- model.matrix(formula_sp, data = d)
  y <- d$target
  groups <- split(seq_len(nrow(d)), d$Timestamp)
  stopifnot(all(vapply(groups, length, integer(1)) == length(regions)))
  Wy <- numeric(length(y))
  for (ii in groups) {
    ii <- ii[order(d$region[ii])]
    stopifnot(identical(as.character(d$region[ii]), regions))
    Wy[ii] <- as.vector(weights %*% y[ii])
  }
  qrX <- qr(X)
  if (qrX$rank != ncol(X)) stop("SAR design matrix is rank deficient.")
  ey <- qr.resid(qrX, y)
  ew <- qr.resid(qrX, Wy)
  n <- length(y)
  nt <- length(groups)
  # Profile negative log-likelihood, omitting constants. For fixed rho,
  # beta is OLS of (y-rho*Wy) on X and sigma^2 is the residual MSE.
  nll <- function(rho) {
    rss <- sum((ey - rho * ew)^2)
    ld <- determinant(diag(5) - rho * weights, logarithm = TRUE)
    if (ld$sign <= 0 || rss <= 0) return(Inf)
    n / 2 * log(rss / n) - nt * as.numeric(ld$modulus)
  }
  opt <- optimize(nll, interval = c(-0.995, 0.995), tol = 1e-9)
  rho <- opt$minimum
  beta <- qr.coef(qrX, y - rho * Wy)
  # Numerical check: the optimum must be no worse than a coarse search.
  stopifnot(opt$objective <= min(vapply(seq(-0.99, 0.99, length.out = 101),
                                       nll, numeric(1))) + 1e-6)
  list(arcoef = rho, coefficients = beta, weights = weights,
       sigma2 = sum((ey - rho * ew)^2) / n,
       method = "Pooled Gaussian SAR maximum likelihood with fixed neighbours")
}

predict_pooled_sar <- function(model, d) {
  d %>% group_split(Timestamp) %>% map_dfr(function(g) {
    g <- arrange(g, region)
    stopifnot(identical(as.character(g$region), regions))
    X <- model.matrix(delete.response(terms(formula_sp)), data = g)
    xb <- as.vector(X[, names(model$coefficients), drop = FALSE] %*%
                      model$coefficients)
    g$prediction <- as.vector(solve(diag(5) - model$arcoef * model$weights, xb))
    # Confirm that the joint predictions satisfy the fitted equations.
    stopifnot(max(abs(g$prediction - model$arcoef *
                       as.vector(model$weights %*% g$prediction) - xb)) < 1e-7)
    g
  }) %>% arrange(region, Timestamp)
}

sar_metrics <- function(predictions, h) {
  map_dfr(c("OLS", "SAR"), function(m) {
    pred <- predictions[[m]]
    obs <- predictions$target
    tibble(Horizon = paste0("t+", h, " hr"), Model = m,
           RMSE = sqrt(mean((pred - obs)^2)), MAE = mean(abs(pred - obs)),
           R2 = cor(pred, obs)^2,
           Predictive_R2 = 1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2))
  })
}

fit_sar <- function(h) {
  d <- build_sar_data(h)
  train_sp <- filter(d, target_time < sar_test_start)
  test_sp <- filter(d, target_time >= sar_test_start)
  stopifnot(nrow(train_sp) > 0, nrow(test_sp) > 0)
  # The same original neighbour matrix is used at both horizons.
  # Rho and regression coefficients are estimated, not hyperparameters tuned.
  sar <- fit_pooled_sar(train_sp, W)
  ols <- lm(formula_sp, data = train_sp)
  forecasts <- predict_pooled_sar(sar, test_sp) %>%
    rename(SAR = prediction) %>%
    mutate(OLS = as.numeric(predict(ols, newdata = .)),
           Horizon = paste0("t+", h, " hr"))
  list(results = sar_metrics(filter(forecasts, region == "west"), h),
       pooled_results = sar_metrics(forecasts, h),
       rho = sar$arcoef, sar_model = sar, ols_model = ols,
       forecasts = forecasts)
}

phase5_fit <- map(c(1, 6), fit_sar)
names(phase5_fit) <- c("t+1 hr", "t+6 hr")
# Primary scores use West, matching Phases 3 and 4. Keep pooled scores separate.
phase5_results <- map_dfr(phase5_fit, "results")
phase5_pooled_results <- map_dfr(phase5_fit, "pooled_results")
phase5_forecasts <- map_dfr(phase5_fit, "forecasts")
phase5_spatial <- map_dfr(seq_along(phase5_fit), function(i) {
  fit <- phase5_fit[[i]]
  tibble(Horizon = names(phase5_fit)[i], rho = fit$rho)
})
print(phase5_results)
print(phase5_pooled_results)
print(phase5_spatial)

sar_forecast_plot <- phase5_forecasts %>% filter(region == "west") %>%
  dplyr::select(target_time, Horizon, Actual = target, OLS, SAR) %>%
  pivot_longer(c(Actual, OLS, SAR), names_to = "Series", values_to = "PM25") %>%
  ggplot(aes(x = target_time, y = PM25, colour = Series, linetype = Series)) +
  geom_line(linewidth = 0.45) + facet_wrap(~ Horizon, ncol = 1) +
  scale_colour_manual(values = c(Actual = "black", OLS = "#D55E00", SAR = "#0072B2")) +
  scale_linetype_manual(values = c(Actual = "solid", OLS = "dotted", SAR = "solid")) +
  labs(title = "SAR and OLS Forecasts vs Actual PM2.5: West",
       subtitle = "Held-out target period: 21-31 October 2019",
       x = "Target timestamp", y = "PM2.5 (ug/m3)", colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
print(sar_forecast_plot)
ggsave("SAR_Forecast_vs_Actual.png", sar_forecast_plot, width = 10, height = 6, dpi = 300)

sar_regional_plot <- phase5_forecasts %>%
  dplyr::select(target_time, region, Horizon, Actual = target, SAR) %>%
  pivot_longer(c(Actual, SAR), names_to = "Series", values_to = "PM25") %>%
  ggplot(aes(target_time, PM25, colour = Series)) + geom_line(linewidth = 0.35) +
  facet_grid(region ~ Horizon) +
  scale_colour_manual(values = c(Actual = "black", SAR = "#0072B2")) +
  labs(title = "SAR Forecasts Across All Five Regions", x = "Target timestamp",
       y = "PM2.5 (ug/m3)", colour = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave("SAR_Regional_Forecasts.png", sar_regional_plot, width = 12, height = 10, dpi = 300)

sar_spatial_plot <- ggplot(phase5_spatial, aes(x = Horizon, y = rho)) +
  geom_col(fill = "#0072B2", width = 0.45) +
  geom_text(aes(label = sprintf("rho = %.3f", rho)), vjust = -0.5, size = 4) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Spatial Dependence at Each Forecast Horizon",
       subtitle = "The same fixed neighbour links are used for both models",
       x = "Forecast horizon", y = "Estimated spatial coefficient (rho)") +
  theme_minimal(base_size = 11)
print(sar_spatial_plot)
ggsave("SAR_Spatial_Coefficient.png", sar_spatial_plot, width = 8, height = 4.5, dpi = 300)

write_csv(phase5_results, "SAR_Test_Metrics_West.csv")
write_csv(phase5_pooled_results, "SAR_Test_Metrics_All_Regions.csv")
write_csv(phase5_spatial, "SAR_Spatial_Coefficients.csv")
write_csv(phase5_forecasts, "SAR_Test_Forecasts.csv")
saveRDS(list(results = phase5_results, pooled_results = phase5_pooled_results,
             spatial = phase5_spatial,
             forecasts = phase5_forecasts, fits = phase5_fit), "SAR_Analysis_Outputs.rds")
# Keep compact objects needed to render the report without rerunning training.
saveRDS(mget(c("long_data", "model_t1", "model_t6", "elastic_net_results",
               "phase3_results", "phase3_group", "phase3_oob", "phase3_fit",
               "svm_results", "svm_selection", "svm_kernel_summary", "svm_forecasts",
               "phase4_results", "phase5_results", "phase5_pooled_results",
               "phase5_spatial", "phase5_forecasts")),
        "Report_Analysis_Outputs.rds")
