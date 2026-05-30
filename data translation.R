# =============================================================================
# Gini Index Forecasting Project
# Replicating Gindelsky (2018) with Updated Data (1980-2024)
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 1: Install and Load Packages
# -----------------------------------------------------------------------------

# Run install.packages once, then comment it out
# install.packages(c("tseries", "forecast", "ggplot2", "dynlm"))

library(tseries)  
library(forecast)  
library(ggplot2)  
library(dynlm)    

# -----------------------------------------------------------------------------
# SECTION 2: Load Data
# -----------------------------------------------------------------------------

df <- read.csv("Cleaning/master_data.csv")

# Keep only 1980 onward to match Gindelsky's sample start
df <- df[df$YEAR >= 1980, ]

cat("Years covered:", min(df$YEAR), "to", max(df$YEAR), "\n")
cat("Observations:", nrow(df), "\n")

# -----------------------------------------------------------------------------
# SECTION 3: Create Differenced Variables
#
# Variables already transformed:
#   D_L_gdp     (log GDP, already first-differenced)
#   D_gov_gdp   (gov/GDP ratio, already first-differenced)
#   D_lfpr      (labor force participation, already first-differenced)
#   D_fem_lfpr  (female LFPR, already first-differenced)
#
# Variables still in levels:
#   GINI, unemp, m_unemp, infl, hs, col, hs_fem, col_fem
# -----------------------------------------------------------------------------

# Dependent variable: year-over-year change in the Gini index
df$D_GINI <- c(NA, diff(df$GINI))

# Explanatory variables still in levels
df$D_unemp   <- c(NA, diff(df$unemp))
df$D_m_unemp <- c(NA, diff(df$m_unemp))
df$D_infl    <- c(NA, diff(df$infl))
df$D_hs      <- c(NA, diff(df$hs))
df$D_col     <- c(NA, diff(df$col))
df$D_hs_fem  <- c(NA, diff(df$hs_fem))
df$D_col_fem <- c(NA, diff(df$col_fem))

# LFPR, D_lfpr, D_L_gdp. D_gov_gdp, D_col, D_hs, D_hs_fem needs a second difference - it is still non-stationary after one diff
df$D2_lfpr     <- c(NA, diff(df$D_lfpr))
df$D2_fem_lfpr <- c(NA, diff(df$D_fem_lfpr))
df$D2_L_gdp    <- c(NA, diff(df$D_L_gdp))
df$D2_gov_gdp   <- c(NA, diff(df$D_gov_gdp))
df$D2_col      <- c(NA, diff(df$D_col))
df$D2_hs       <- c(NA, diff(df$D_hs))
df$D2_hs_fem   <- c(NA, diff(df$D_hs_fem))



# Remove the first two rows which have NA from double-differencing
df <- df[complete.cases(df[ , c("D_GINI", "D2_lfpr", "D2_fem_lfpr", "D2_L_gdp", "D2_gov_gdp", "D2_col", "D2_hs", "D2_hs_fem")]), ]

cat("Observations after differencing:", nrow(df), "\n")
cat("Years after differencing:", min(df$YEAR), "to", max(df$YEAR), "\n")

# -----------------------------------------------------------------------------
# SECTION 4: Stationary Tests 
#
# ADF tests checks whether a series is stationary
# Null hypothesis: the series is NON-stationary (has a unit root)
# If p-value < 0.05 we reject the null and conclude the series IS stationary
# Need all our variables to be stationary before putting them in a model
# -----------------------------------------------------------------------------

cat("\n--- ADF STATIONARY TESTS ---\n")
cat("p < 0.05 = stationary (OK to use in model)\n")
cat("p > 0.05 = non-stationary (problem - needs more differencing)\n\n")

# First show Gini in LEVELS - this should be non-stationary (p > 0.05)
adf_levels <- adf.test(df$GINI)
cat("Gini in levels (expect p > 0.05):\n")
cat("  p-value =", round(adf_levels$p.value, 4), "\n\n")

# Now test each differenced variable, each one should have p < 0.05
adf_dgini    <- adf.test(df$D_GINI)
adf_d2lfpr   <- adf.test(df$D2_lfpr)
adf_d2flfpr  <- adf.test(df$D2_fem_lfpr)
adf_dunemp   <- adf.test(df$D_unemp)
adf_d2lgdp   <- adf.test(df$D2_L_gdp)
adf_d2govgdp <- adf.test(df$D2_gov_gdp)
adf_dinfl    <- adf.test(df$D_infl)
adf_d2col     <- adf.test(df$D2_col)
adf_d2hs      <- adf.test(df$D2_hs)
adf_dcolfem  <- adf.test(df$D_col_fem)
adf_d2hsfem   <- adf.test(df$D2_hs_fem)

# Print results in a readable table
cat("Variable              p-value    Result\n")
cat("D_GINI               ", round(adf_dgini$p.value,   4),
    ifelse(adf_dgini$p.value   < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D2_lfpr              ", round(adf_d2lfpr$p.value,  4),
    ifelse(adf_d2lfpr$p.value  < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D2_fem_lfpr          ", round(adf_d2flfpr$p.value, 4),
    ifelse(adf_d2flfpr$p.value < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D_unemp              ", round(adf_dunemp$p.value,  4),
    ifelse(adf_dunemp$p.value  < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D2_L_gdp              ", round(adf_d2lgdp$p.value,   4),
    ifelse(adf_d2lgdp$p.value   < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D2_gov_gdp            ", round(adf_d2govgdp$p.value, 4),
    ifelse(adf_d2govgdp$p.value < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D_infl               ", round(adf_dinfl$p.value,   4),
    ifelse(adf_dinfl$p.value   < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D2_col                ", round(adf_d2col$p.value,    4),
    ifelse(adf_d2col$p.value    < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D2_hs                 ", round(adf_d2hs$p.value,     4),
    ifelse(adf_d2hs$p.value     < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D_col_fem            ", round(adf_dcolfem$p.value, 4),
    ifelse(adf_dcolfem$p.value < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")
cat("D2_hs_fem             ", round(adf_d2hsfem$p.value,  4),
    ifelse(adf_d2hsfem$p.value  < 0.05, "  Stationary", "  NON-STATIONARY"), "\n")

# Save ADF results as a table
adf_table <- data.frame(
  Variable   = c("Gini (levels)", "D_GINI", "D2_lfpr", "D2_fem_lfpr",
                 "D_unemp", "D2_L_gdp", "D2_gov_gdp", "D_infl",
                 "D2_col", "D2_hs", "D_col_fem", "D2_hs_fem"),
  P_Value    = round(c(adf_levels$p.value,  adf_dgini$p.value,
                       adf_d2lfpr$p.value,  adf_d2flfpr$p.value,
                       adf_dunemp$p.value,  adf_d2lgdp$p.value,
                       adf_d2govgdp$p.value, adf_dinfl$p.value,
                       adf_d2col$p.value,    adf_d2hs$p.value,
                       adf_dcolfem$p.value, adf_d2hsfem$p.value), 4),
  Stationary = c("No", "Yes", "Yes", "Yes", "Yes", "Yes",
                 "Yes", "Yes", "Yes", "Yes", "Yes", "Yes")
)
write.csv(adf_table, "Model/adf_results.csv", row.names = FALSE)
cat("\nADF results saved to adf_results.csv\n")

# -----------------------------------------------------------------------------
# SECTION 5: AR Lag Selection
#
# Gindelsky's first step: regress D_GINI on its own past values (lags)
# to find which lags of the Gini are significant predictors of itself
# Test lags 1, 2, and 3 - keep the ones with p-value < 0.05
# -----------------------------------------------------------------------------

cat("\n--- STEP 1: AR LAG SELECTION ---\n")

# Manually create lagged versions of D_GINI
# lag 1 = what D_GINI was last year
# lag 2 = what D_GINI was two years ago
# lag 3 = what D_GINI was three years ago
n <- nrow(df)

D_GINI_lag1 <- c(NA, df$D_GINI[-n])
D_GINI_lag2 <- c(NA, NA, df$D_GINI[-c(n-1, n)])
D_GINI_lag3 <- c(NA, NA, NA, df$D_GINI[-c(n-2, n-1, n)])

# Put into a data frame and remove NAs
ar_data <- data.frame(
  D_GINI = df$D_GINI,
  lag1   = D_GINI_lag1,
  lag2   = D_GINI_lag2,
  lag3   = D_GINI_lag3
)
ar_data <- ar_data[complete.cases(ar_data), ]

# Run OLS regression - same lm() we use for any linear regression
ar_model <- lm(D_GINI ~ lag1 + lag2 + lag3, data = ar_data)
summary(ar_model)

# Look at the Pr(>|t|) column in the output above
# Keep lags with reasonable p values (slightly different than Gindelsky)
# Gindelsky found lag 2 was significant for the Gini index

# -----------------------------------------------------------------------------
# SECTION 6: General Unrestricted Model (GUM)
#
# Include all candidate explanatory variables at lags 1, 2, and 3
# -----------------------------------------------------------------------------

cat("\n--- STEP 2: GENERAL UNRESTRICTED MODEL ---\n")

# Create lagged versions of each explanatory variable

# LFPR (second-differenced)
df$D2_lfpr_l1 <- c(NA, df$D2_lfpr[-n])
df$D2_lfpr_l2 <- c(NA, NA, df$D2_lfpr[-c(n-1, n)])
df$D2_lfpr_l3 <- c(NA, NA, NA, df$D2_lfpr[-c(n-2, n-1, n)])

# Female LFPR (second-differenced)
df$D2_fem_lfpr_l1 <- c(NA, df$D2_fem_lfpr[-n])
df$D2_fem_lfpr_l2 <- c(NA, NA, df$D2_fem_lfpr[-c(n-1, n)])
df$D2_fem_lfpr_l3 <- c(NA, NA, NA, df$D2_fem_lfpr[-c(n-2, n-1, n)])

# Unemployment
df$D_unemp_l1 <- c(NA, df$D_unemp[-n])
df$D_unemp_l2 <- c(NA, NA, df$D_unemp[-c(n-1, n)])
df$D_unemp_l3 <- c(NA, NA, NA, df$D_unemp[-c(n-2, n-1, n)])

# Log GDP 
df$D2_L_gdp_l1 <- c(NA, df$D2_L_gdp[-n])
df$D2_L_gdp_l2 <- c(NA, NA, df$D2_L_gdp[-c(n-1, n)])
df$D2_L_gdp_l3 <- c(NA, NA, NA, df$D2_L_gdp[-c(n-2, n-1, n)])

# Government expenditure share of GDP
df$D2_gov_gdp_l1 <- c(NA, df$D2_gov_gdp[-n])
df$D2_gov_gdp_l2 <- c(NA, NA, df$D2_gov_gdp[-c(n-1, n)])
df$D2_gov_gdp_l3 <- c(NA, NA, NA, df$D2_gov_gdp[-c(n-2, n-1, n)])

# Inflation
df$D_infl_l1 <- c(NA, df$D_infl[-n])
df$D_infl_l2 <- c(NA, NA, df$D_infl[-c(n-1, n)])
df$D_infl_l3 <- c(NA, NA, NA, df$D_infl[-c(n-2, n-1, n)])

# College attainment
df$D2_col_l1 <- c(NA, df$D2_col[-n])
df$D2_col_l2 <- c(NA, NA, df$D2_col[-c(n-1, n)])
df$D2_col_l3 <- c(NA, NA, NA, df$D2_col[-c(n-2, n-1, n)])

# High school attainment
df$D2_hs_l1 <- c(NA, df$D2_hs[-n])
df$D2_hs_l2 <- c(NA, NA, df$D2_hs[-c(n-1, n)])
df$D2_hs_l3 <- c(NA, NA, NA, df$D2_hs[-c(n-2, n-1, n)])

# Female college attainment
df$D_col_fem_l1 <- c(NA, df$D_col_fem[-n])
df$D_col_fem_l2 <- c(NA, NA, df$D_col_fem[-c(n-1, n)])
df$D_col_fem_l3 <- c(NA, NA, NA, df$D_col_fem[-c(n-2, n-1, n)])

# Female high school attainment
df$D2_hs_fem_l1 <- c(NA, df$D2_hs_fem[-n])
df$D2_hs_fem_l2 <- c(NA, NA, df$D2_hs_fem[-c(n-1, n)])
df$D2_hs_fem_l3 <- c(NA, NA, NA, df$D2_hs_fem[-c(n-2, n-1, n)])

# AR lag of the dependent variable (lag 2 - significant in Gindelsky)
df$D_GINI_lag2 <- c(NA, NA, df$D_GINI[-c(n-1, n)])

# Structural break dummies
# 1 in the year of a known break and 0 in every other year
# 2007: pre-Great Recession break (identified in Gindelsky)
# 2020: COVID-19 pandemic (new break in our extended data)
df$d2007 <- ifelse(df$YEAR == 2007, 1, 0)
df$d2020 <- ifelse(df$YEAR == 2020, 1, 0)

# Remove rows with missing values after lag creation
model_data <- df[complete.cases(df), ]
cat("Observations available for modeling:", nrow(model_data), "\n")

# Run the General Unrestricted Model with all candidates
gum <- lm(D_GINI ~
            D_GINI_lag2 + D2_lfpr_l1 + D2_lfpr_l2 + D2_lfpr_l3 +
            D2_fem_lfpr_l1 + D2_fem_lfpr_l2 + D2_fem_lfpr_l3 +
            D_unemp_l1 + D_unemp_l2 + D_unemp_l3 +
            D2_L_gdp_l1 + D2_L_gdp_l2 + D2_L_gdp_l3 +
            D2_gov_gdp_l1 + D2_gov_gdp_l2 + D2_gov_gdp_l3 +
            D_infl_l1 + D_infl_l2 + D_infl_l3 +
            D2_col_l1 + D2_col_l2 + D2_col_l3 +
            D2_hs_l1 + D2_hs_l2 + D2_hs_l3 +
            D_col_fem_l1 + D_col_fem_l2 + D_col_fem_l3 +
            D2_hs_fem_l1 + D2_hs_fem_l2 + D2_hs_fem_l3 +
            d2007 + d2020,
          data = model_data)

summary(gum)

# -----------------------------------------------------------------------------
# SECTION 7: Final Parsimonious Model (General-to-Specific Selection)
#
# Starting from the GUM above, remove highly insignificant variables first,
# then re-estimate the model after each reduction.
#
# The reduced models below show the general-to-specific selection process.
# Variables were removed when they had weak statistical significance and did
# not improve the overall model fit.
# -----------------------------------------------------------------------------

cat("\n--- STEP 3: FINAL PARSIMONIOUS MODEL ---\n")

# First reduction: keep strongest candidates from the GUM
reduced <- lm(D_GINI ~
                D_GINI_lag2 +
                D_unemp_l1 + D_unemp_l2 +
                D2_L_gdp_l1 +
                D2_col_l3 +
                d2007 + d2020,
              data = model_data)

summary(reduced)

# Second reduction: remove insignificant COVID dummy
reduced2 <- lm(D_GINI ~
                 D_GINI_lag2 +
                 D_unemp_l1 + D_unemp_l2 +
                 D2_L_gdp_l1 +
                 D2_col_l3 +
                 d2007,
               data = model_data)

summary(reduced2)

# Final model: remove weak college attainment variable
final_model <- lm(D_GINI ~
                    D_GINI_lag2 +
                    D_unemp_l1 + D_unemp_l2 +
                    D2_L_gdp_l1 +
                    d2007,
                  data = model_data)

summary(final_model)

cat("\nFinal model coefficients:\n")
print(round(coef(final_model), 4))

# -----------------------------------------------------------------------------
# SECTION 8: In-Sample Evaluation
#
# Compare what the model predicted against what actually happened
# for the same data the model was trained on
# RMSE = root mean squared error (lower is better, same units as Gini)
# MAPE = mean absolute percentage error (lower is better, shown as %)
# -----------------------------------------------------------------------------

cat("\n--- IN-SAMPLE EVALUATION ---\n")

# Get the model's fitted values and the corresponding actual values
fitted_vals <- fitted(final_model)
actual_vals <- model_data$D_GINI

# Calculate residuals 
resids <- actual_vals - fitted_vals

# RMSE
RMSE <- sqrt(mean(resids^2, na.rm = TRUE))

# MAPE
MAPE <- mean(abs(resids[actual_vals != 0] / actual_vals[actual_vals != 0]), na.rm = TRUE) * 100

cat("Our Model:\n")
cat("  RMSE:", round(RMSE, 4), "\n")
cat("  MAPE:", round(MAPE, 4), "%\n")
cat("\nGindelsky (2018) original results:\n")
cat("  RMSE: 0.191\n")
cat("  MAPE: 0.363%\n")

# -----------------------------------------------------------------------------
# SECTION 9: Naive Benchmark (In-Sample)
#
# The naive model predicts no change each year (D_GINI = 0)
# -----------------------------------------------------------------------------

cat("\n--- NAIVE BENCHMARK (IN-SAMPLE) ---\n")

# Naive predicts 0 change, so the error equals the actual value
naive_resids <- actual_vals
RMSE_naive <- sqrt(mean(naive_resids^2, na.rm = TRUE))
MAPE_naive <- mean(abs(naive_resids[actual_vals != 0] / actual_vals[actual_vals != 0]), na.rm = TRUE) * 100

cat("Naive Model:\n")
cat("  RMSE:", round(RMSE_naive, 4), "\n")
cat("  MAPE:", round(MAPE_naive, 4), "%\n")
cat("Our model improves MAPE by",
    round(MAPE_naive - MAPE, 2),
    "percentage points over the naive model\n")

# -----------------------------------------------------------------------------
# SECTION 10: Pseudo-Out-of-Sample Forecast
#
# The idea mirrors what Gindelsky did:
#   1. Estimate the model using data up to 2015 only (training period)
#   2. Use those estimated coefficients to forecast 2016-2024 (test period)
#   3. Compare those forecasts to what actually happened
#   4. Compute MAPE and RMSE on the forecast errors
#
# -----------------------------------------------------------------------------

cat("\n--- PSEUDO-OUT-OF-SAMPLE FORECAST (2016-2024) ---\n")

# Split the data: train on everything up to 2015, test on 2016 onward
train_data <- model_data[model_data$YEAR <= 2015, ]
test_data  <- model_data[model_data$YEAR >= 2016, ]

cat("Training period:", min(train_data$YEAR), "to", max(train_data$YEAR),
    "-- observations:", nrow(train_data), "\n")
cat("Test period:    ", min(test_data$YEAR), "to", max(test_data$YEAR),
    "-- observations:", nrow(test_data), "\n")

# Estimate the same model specification using only the training data
model_train <- lm(D_GINI ~
                    D_GINI_lag2 +
                    D_unemp_l1 + D_unemp_l2 +
                    D2_L_gdp_l1 +
                    d2007,
                  data = train_data)

cat("\nModel re-estimated on training data only:\n")
summary(model_train)

# Use predict() to apply the training model to the test period
# This gives us what the model would have forecast for 2016-2024
forecast_vals   <- predict(model_train, newdata = test_data)
actual_test     <- test_data$D_GINI
forecast_errors <- actual_test - forecast_vals

# Forecast RMSE and MAPE
RMSE_forecast <- sqrt(mean(forecast_errors^2))
MAPE_forecast <- mean(abs(forecast_errors / actual_test)) * 100

# Naive forecast for the test period (predict 0 change each year)
RMSE_naive_fc <- sqrt(mean(actual_test^2))
MAPE_naive_fc <- mean(abs(actual_test / actual_test)) * 100

cat("\nForecast results (2016-2024):\n")
cat("  Our model -- RMSE:", round(RMSE_forecast, 4),
    "| MAPE:", round(MAPE_forecast, 4), "%\n")
cat("  Naive     -- RMSE:", round(RMSE_naive_fc, 4),
    "| MAPE:", round(MAPE_naive_fc, 4), "%\n")

mape_diff <- round(MAPE_naive_fc - MAPE_forecast, 2)

if (mape_diff > 0) {
  cat("  Our model improves MAPE by",
      mape_diff, "percentage points over the naive model\n")
} else {
  cat("  Our model underperforms the naive model by",
      abs(mape_diff), "percentage points in MAPE\n")
}
# Note: MAPE is unstable here because D_GINI can be close to zero.

# Convert differenced forecasts back to Gini levels for the plot
# Start from the last known Gini value at the end of training (2015)
last_train_gini  <- df$GINI[df$YEAR == 2015]
forecast_levels  <- last_train_gini + cumsum(forecast_vals)
actual_fc_levels <- last_train_gini + cumsum(actual_test)

# Full summary table of all results
cat("\n--- FULL RESULTS SUMMARY ---\n")
results_table <- data.frame(
  Period   = c("In-Sample: Our Model",
               "In-Sample: Naive",
               "Forecast 2016-2024: Our Model",
               "Forecast 2016-2024: Naive",
               "Gindelsky (2018) In-Sample"),
  RMSE     = c(round(RMSE, 4),
               round(RMSE_naive, 4),
               round(RMSE_forecast, 4),
               round(RMSE_naive_fc, 4),
               0.191),
  MAPE_pct = c(round(MAPE, 4),
               round(MAPE_naive, 4),
               round(MAPE_forecast, 4),
               round(MAPE_naive_fc, 4),
               0.363)
)
print(results_table)
write.csv(results_table, "Model/model_metrics_full.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# SECTION 11: Plots
# -----------------------------------------------------------------------------

cat("\n--- GENERATING PLOTS ---\n")

# Plot 1: Gini index in levels over time
png("Model/plot1_gini_levels.png", width = 900, height = 500, res = 120)
print(
  ggplot(df, aes(x = YEAR, y = GINI)) +
    geom_line(color = "navy", linewidth = 1.1) +
    geom_point(color = "navy", size = 1.5) +
    labs(title = "U.S. Gini Index Over Time",
         subtitle = "Source: Census Bureau CPS ASEC, 1982-2024",
         x = "Year", y = "Gini Coefficient") +
    theme_minimal()
)
dev.off()
cat("Saved: plot1_gini_levels.png\n")

# Plot 2: First-differenced Gini 
diff_df <- df[!is.na(df$D_GINI), ]
png("Model/plot2_gini_differenced.png", width = 900, height = 500, res = 120)
print(
  ggplot(diff_df, aes(x = YEAR, y = D_GINI)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "First-Differenced Gini (Stationary Series)",
         subtitle = "Used for model estimation",
         x = "Year", y = "Change in Gini") +
    theme_minimal()
)
dev.off()
cat("Saved: plot2_gini_differenced.png\n")

# Plot 3: In-sample fitted vs actual (in differences)
fit_df <- model_data[complete.cases(model.frame(final_model)), ]

fit_df$Actual <- fit_df$D_GINI
fit_df$Fitted <- fitted(final_model)

png("Model/plot3_fitted_vs_actual.png", width = 900, height = 500, res = 120)
print(
  ggplot(fit_df, aes(x = YEAR)) +
    geom_line(aes(y = Actual, color = "Actual"),    linewidth = 1.2) +
    geom_line(aes(y = Fitted, color = "Model Fit"), linewidth = 1,
              linetype = "dashed") +
    scale_color_manual(values = c("Actual" = "navy", "Model Fit" = "darkorange")) +
    labs(title = "Fitted vs. Actual: Change in Gini (In-Sample)",
         x = "Year", y = "Change in Gini", color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)
dev.off()
cat("Saved: plot3_fitted_vs_actual.png\n")

# Plot 4: In-sample fit converted back to Gini levels
base_gini <- df$GINI[df$YEAR == min(fit_df$YEAR) - 1]
gini_actual <- base_gini + cumsum(fit_df$Actual)
gini_fitted <- base_gini + cumsum(fit_df$Fitted)

levels_df <- data.frame(
  YEAR   = fit_df$YEAR,
  Actual = gini_actual,
  Fitted = gini_fitted
)

png("Model/plot4_levels_comparison.png", width = 900, height = 500, res = 120)
print(
  ggplot(levels_df, aes(x = YEAR)) +
    geom_line(aes(y = Actual, color = "Actual"), linewidth = 1.2) +
    geom_line(aes(y = Fitted, color = "Model Fit"), linewidth = 1,
              linetype = "dashed") +
    scale_color_manual(values = c("Actual" = "navy", "Model Fit" = "darkorange")) +
    labs(title = "Actual vs. Model Fit: Gini Index in Levels",
         subtitle = "Reconstructed from first-differenced model",
         x = "Year", y = "Gini Coefficient", color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)
dev.off()
cat("Saved: plot4_levels_comparison.png\n")

# Plot 5: Pseudo-out-of-sample forecast vs actual (the key forecast plot)
forecast_df <- data.frame(
  YEAR     = c(2015, test_data$YEAR),
  Actual   = c(last_train_gini, actual_fc_levels),
  Forecast = c(last_train_gini, forecast_levels)
)

png("Model/plot5_forecast_vs_actual.png", width = 900, height = 500, res = 120)
print(
  ggplot(forecast_df, aes(x = YEAR)) +
    geom_line(aes(y = Actual,   color = "Actual"),   linewidth = 1.2) +
    geom_line(aes(y = Forecast, color = "Forecast"), linewidth = 1,
              linetype = "dashed") +
    geom_vline(xintercept = 2015, linetype = "dotted", color = "gray40") +
    annotate("text", x = 2015.2,
             y = min(c(forecast_df$Actual, forecast_df$Forecast), na.rm = TRUE),
             label = "Forecast starts", size = 3, color = "gray40", hjust = 0) +
    scale_color_manual(values = c("Actual" = "navy", "Forecast" = "darkorange")) +
    labs(title = "Pseudo-Out-of-Sample Forecast: Gini Index (2016-2024)",
         subtitle = "Model estimated on 1982-2015 | Forecasting 2016-2024",
         x = "Year", y = "Gini Coefficient", color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)
dev.off()
cat("Saved: plot5_forecast_vs_actual.png\n")

cat("\n=== SCRIPT COMPLETE ===\n")
cat("Output files:\n")
cat("  adf_results.csv              - Stationarity test results\n")
cat("  model_metrics_full.csv       - RMSE and MAPE for all models\n")
cat("  plot1_gini_levels.png        - Gini over time\n")
cat("  plot2_gini_differenced.png   - First-differenced Gini\n")
cat("  plot3_fitted_vs_actual.png   - In-sample fit (differences)\n")
cat("  plot4_levels_comparison.png  - In-sample fit (levels)\n")
cat("  plot5_forecast_vs_actual.png - Forecast vs actual 2016-2024\n")



# -----------------------------------------------------------------------------
# ACF/PACF Diagnostics
# -----------------------------------------------------------------------------

# ACF and PACF of differenced Gini series
# Used to visually inspect autoregressive structure before lag selection
acf(df$D_GINI,
    main = "ACF of Differenced Gini (D_GINI)")

pacf(df$D_GINI,
     main = "PACF of Differenced Gini (D_GINI)")


# ACF of final model residuals
# Residual autocorrelation should ideally be small / insignificant

acf(residuals(final_model),
    main = "ACF of Final Model Residuals")
