# ------------------------------ Code Info -------------------------------------
# title: "3.explain_heterogeneity.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, 
# CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-05"

# -------------- (Optional) Clear the environment and free memory --------------
rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# ------------------------------ Libraries -------------------------------------
library("data.table")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("rgl")
library("randomForest")
library("caret")
library("e1071")
library("lme4")
library("glmnet")
library("leaps")
library("pdp")
library("ggplot2")
library("DALEX")
library("tuneRanger")
library("mlr")
library("OpenML")

# ------ Define working dir as the directory where the script is located -------
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# --------------------------- Import useful functions --------------------------
source("functions_lidar.R")
source("../functions_plots.R")

# ----------------------------- Plotting Parameters ----------------------------
default_par <- par(no.readonly = TRUE)
par(mar=c(5, 5, 4, 2) + 0.1) # Adjust margins as per your preference
par(oma=c(0, 0, 0, 0)) # Adjust outer margins as per your preference
par(cex=1.6, cex.axis=1.6, cex.names=1.6)
# par(default_par)
# ---------------------------- Directories & Setup  ---------------------------- 
data <- "../../01_DATA"
site <- "Mormal"
data_dir <- file.path(data, site)

results_dir <- '../../03_RESULTS'
results_path <- file.path("../../03_RESULTS", site)
lidar_dir <- "LiDAR"
lidr_dir <- file.path(results_path, lidar_dir, "PAI/lidR")
masks_dir <- file.path(results_path, lidar_dir, "Heterogeneity_Masks")
dir.create(path = masks_dir, showWarnings = F, recursive = T)
profiles_dir <- file.path(results_path, lidar_dir, "Profiles")
deciles_dir <- "Deciles"
equal_intervals_dir <- "Equal_Intervals"
deciles_dir <- file.path(masks_dir, "Deciles")
equal_intervals_dir <- file.path(masks_dir, "Equal_Intervals")
class10_dir <- file.path(deciles_dir, "10_Quantiles")
class3_dir <- file.path(deciles_dir, "3_Quantiles")

# ----------------------------- Localisation -----------------------------------
# training_points_file <- paste(data_dir,
#                               "Shape/shp_mormal_onf/foret_domaniale_mormal.shp", 
#                               sep = "/")
# training_points_bbox <- paste(data_dir, 
#                               "Shape/foret_mormal_shape_reprojected_bbox.shp", 
#                               sep = "/")
# training_points <- vect(training_points_file)
# bbox <- vect(training_points_bbox)
# crs_bbox <- crs(bbox)
# training_points_projected <- project(training_points, crs_bbox)

# ---------------------------- Data Preparation --------------------------------
# Call rasters
# Mask
mask_10m <- terra::rast(
  file.path(masks_dir, "artifacts_low_vegetation_majority_90_p_res_10_m.envi"))
chm <- terra::rast(file.path(masks_dir, "mnc_masked_with_average_mnc_thresholded_90_p_kept_res_10_m.png.envi"))

# We want to predict LiDAR LAI
lai_lidar_raster <- terra::rast(file.path(masks_dir, "lai_lidar_masked_res_10_m.envi")) 
# Predictors
lai_s2_raster <- terra::rast(file.path(masks_dir, "lai_s2_masked_res_10_m.envi"))
mean_heights_raster <- terra::rast(file.path(masks_dir,
                                             "mnc_mean_heights_res_10_m.envi"))
max_h_raster <- terra::rast(file.path(masks_dir,
                                      "max_res_10_m.envi"))
std_raster <- terra::rast(file.path(masks_dir,
                                    "mnc_std_res_10_m.envi"))
cv_raster <- terra::rast(file.path(masks_dir,
                                   "mnc_coeff_variation_res_10_m.envi"))
variance_raster <- terra::rast(file.path(masks_dir,
                                         "mnc_variance_res_10_m.envi"))
rumple_raster <- terra::rast(file.path(lidr_dir, "rumple_res_10_m_non_norm.tif"))
vci_raster <- terra::rast(file.path(lidr_dir, "vci_res_10_m_non_norm.tif"))
lcv_raster <- terra::rast(file.path(lidr_dir, "lcv_res_10_m_non_norm.tif"))
lskew_raster <- terra::rast(file.path(lidr_dir, "lskew_res_10_m_non_norm.tif"))
vdr_raster <- terra::rast(file.path(lidr_dir, "vdr_res_10_m_non_norm.tif"))

# Define indices
indices <- list(
  na_indices_lmoments = countNA(lcv_raster),
  na_indices_std = countNA(std_raster)
)

# Loop through indices
for (index_name in names(indices)) {
  na_index <- indices[[index_name]]
  lai_lidar_raster[na_index == 1] <- NA
  lai_s2_raster[na_index == 1] <- NA
  mean_heights_raster[na_index == 1] <- NA
  max_h_raster[na_index == 1] <- NA
  std_raster[na_index == 1] <- NA
  cv_raster[na_index == 1] <- NA
  variance_raster[na_index == 1] <- NA
  rumple_raster[na_index == 1] <- NA
  vci_raster[na_index == 1] <- NA
  lcv_raster[na_index == 1] <- NA
  lskew_raster[na_index == 1] <- NA
  vdr_raster[na_index == 1] <- NA
}

# Values 
lai_lidar <- values(lai_lidar_raster)
lai_s2 <- values(lai_s2_raster)
mean_h <- values(mean_heights_raster)
max_h <- values(max_h_raster)
std <- values(std_raster)
cv <- values(cv_raster)
variance <- values(variance_raster)
rumple <- values(rumple_raster)
vci <- values(vci_raster)
lcv <- values(lcv_raster)
lskew <- values(lskew_raster)
vdr <- values(vdr_raster)

# No NA
lai_lidar <- lai_lidar[complete.cases(lai_lidar), ]
lai_s2 <- lai_s2[complete.cases(lai_s2), ]
mean_h <- mean_h[complete.cases(mean_h), ]
max_h <- max_h[complete.cases(max_h), ]
std <- std[complete.cases(std), ]
cv <- cv[complete.cases(cv), ]
variance <- variance[complete.cases(variance), ]
rumple <- rumple[complete.cases(rumple), ]
vci <- vci[complete.cases(vci), ]
lcv <- lcv[complete.cases(lcv), ]
lskew <- lskew[complete.cases(lskew), ]
vdr <- vdr[complete.cases(vdr), ]

# Initial Correlation
plot_density_scatterplot(lai_s2,
                         lai_lidar,
                         "Sentinel-2 LAI",
                         "LiDAR LAI",
                         title = "Sentinel-2 LAI vs LiDAR LAI Values",
                         xlimits = c(0,10),
                         ylimits = c(0,10))

# Prepare predicted raster
lai_lidar_matrix <- as.matrix(lai_lidar_raster)
na_indices <- which(is.na(lai_lidar_matrix), arr.ind = TRUE)
non_na_indices <- which(!is.na(lai_lidar_matrix), arr.ind = TRUE)
na_indices <- na_indices[,1]
non_na_indices <- non_na_indices[,1]
new_raster <- rast(lai_lidar_raster)

# -------------------------------- Training Data -------------------------------
# Create Training DF
training_data_prep <- as.data.frame(cbind(
  lai_lidar, # Predict
  lai_s2,
  std,
  mean_h,
  max_h,
  cv,
  variance,
  rumple,
  vci,
  lcv,
  lskew,
  vdr
))

# Visualize correlation matrix
correlation_matrix <- cor(training_data_prep[, -which(names(training_data_prep) == "lai_lidar")])
corrplot::corrplot(correlation_matrix, method = "number", type = "upper")

# Step 1: Remove Variables with Low Variance
# min_variance <- 0.1  # Set the minimum variance threshold
# low_variance_vars <- names(training_dataa_df)[apply(training_dataa_df, 2, stats::var) < min_variance]
# training_data_cleaned <- training_dataa_df[, !names(training_dataa_df) %in% low_variance_vars]

# Step 2: Select Pertinent Variables (manually)
# pertinent_vars <- c("lai_lidar_val", "mean_h_val", "std_val")  # Example of pertinent variables
# training_data_cleaned <- training_data_cleaned[, c("lai_s2_val", pertinent_vars)]

# Step 3: Remove Highly Correlated Variables
# correlation_threshold <- 0.7  # Set the correlation coefficient threshold
# # correlation_matrix <- cor(training_data_cleaned)
# highly_correlated_pairs <- which(correlation_matrix > correlation_threshold & correlation_matrix < 1, arr.ind = TRUE)
# vars_to_remove <- unique(c(highly_correlated_pairs[, 1], highly_correlated_pairs[, 2]))
# training_data_cleaned <- training_dataa_df[, !names(training_dataa_df) %in% vars_to_remove]
# 
# # Final cleaned dataset
# training_data_cleaned

# Train-test split
standardize_data <- function(data) {
  standardized_data <- data
  for (id_var in length(data)) {
    mean_val <- mean(data[, id_var])
    sd_val <- sd(data[, id_var])
    standardized_data[, id_var] <- (data[, id_var] - mean_val) / sd_val
  }
  return(standardized_data)
}

# Standardize
standardized_data <- standardize_data(training_data_prep)

# Sampling and subset standardized data
set.seed(0) # Reproducibility
sample_size <- 150000
sample_indices <- sample(nrow(standardized_data), sample_size)
sampled_data <- standardized_data[sample_indices, ]

# Subset training_data using the row indices
# train_prop <- 0.7
# train_indices <- createDataPartition(sampled_data$lai_lidar, p = train_prop, list = FALSE)
# training_data <- sampled_data[train_indices, ]
# test_data <- sampled_data[-train_indices, ]
# test_data <- standardized_data
training_data <- head(sampled_data, 5000)
test_data <- tail(sampled_data, 100000)

# Check dimensions
cat("Dimensions of training_data:", dim(training_data),
    "\nDimensions of test_data:", dim(test_data), "\n")

# --------------------------- RF (Cross Validation) ----------------------------

# Initialization
ctrl <- trainControl(method = "LOOCV",
                     number = 5)
# Create formula
formula <- lai_lidar ~ 
  lai_s2 + 
  mean_h + 
  std + 
  cv +
  variance + 
  rumple +
  vci + 
  lcv +
  lskew +
  vdr

# Train the Random Forest model with cross-validation
mtry <- sqrt(ncol(predictors))
tunegrid <- expand.grid(.mtry=c(2:8))
metric <- "RMSE"
predictors = training_data[, -which(names(training_data) == "lai_lidar")]
outcome = training_data$lai_lidar
rf_default <- train(x = predictors,
                    y = outcome, 
                    method="rf", 
                    ntree = 15,
                    metric=metric, 
                    tuneGrid=tunegrid, 
                    trControl=ctrl)
# Check results
print(rf_model)
importance <- varImp(rf_model)
print(importance)
plot(importance)

# ------------------------------- RF tuneRanger --------------------------------
# Define your regression task
regression_task <- makeRegrTask(data = training_data, target = "lai_lidar")

# Estimate runtime
estimateTimeTuneRanger(regression_task)

# Tuning parameters
tuned_ranger <- tuneRanger(
  regression_task,
  measure = list(mse),                # Specify the evaluation measure (e.g., RMSE)
  num.trees = 1000,                   # Number of trees
  num.threads = 4,                    # Number of threads for parallelization
  iters = 100,                         # Number of tuning iterations
  iters.warmup = 30                   # Number of warm-up iterations
)

# Display the results
tuned_ranger
# Print the recommended parameter settings
tuned_ranger$results

# Get the tuned Ranger model
tuned_ranger_model <- tuned_ranger$model

# Make predictions on the test data
predictions_ranger <- predict(tuned_ranger_model, newdata = test_data)
comparison_ranger <- data.frame(Actual = test_data$lai_lidar,
                                Predicted = predictions_ranger)

# Compute metrics
correlation <- cor(comparison_ranger$Actual, comparison_ranger$Predicted.response)
rmse <- Metrics::rmse(comparison_ranger$Actual, comparison_ranger$Predicted.response)
rsquared <- R2(comparison_ranger$Actual, comparison_ranger$Predicted.response)
bias <- Metrics::bias(comparison_ranger$Actual, comparison_ranger$Predicted.response)
mae <- Metrics::mae(comparison_ranger$Actual, comparison_ranger$Predicted.response)

# Print metrics
cat("Root Mean Squared Error:", rmse, 
    "\nR-squared:", rsquared, 
    "\nBias:", bias, 
    "\nCorrelation:", correlation, 
    "\nMAE:", mae)

plot_density_scatterplot(comparison_ranger$Actual,
                         comparison_ranger$Predicted.response,
                         "Actual LiDAR LAI",
                         "Predicted LiDAR LAI",
                         title = "Predicted LiDAR LAI vs Actual LiDAR LAI Values on the Test Set",
                         xlimits = c(0,10),
                         ylimits = c(0,10))

# -------------------------------- Simple RF -----------------------------------

# Set hyperparameters
ntree <- 15
mtry <- 7

# Create formula
formula <- lai_lidar ~ 
  lai_s2 +
  mean_h +
  max_h +
  std +
  cv +
  variance +
  rumple +
  vci +
  lcv +
  lskew +
  vdr

rf_model <- randomForest(formula,
                         data = training_data,
                         ntree = ntree,
                         mtry = mtry,
                         importance = TRUE,
                         do.trace = TRUE
)

# Print the number of trees used
num_trees <- length(rf_model$mse)
print(paste("Number of trees:", num_trees))
print(summary(rf_model))

importance <- varImp(rf_model)
normalized_importance <- importance$Overall / sum(importance$Overall) * 100
names(normalized_importance) <- rownames(importance)
normalized_importance <- sort(normalized_importance, decreasing = TRUE)
df_importance <- data.frame(variable = names(normalized_importance), 
                            importance = normalized_importance)
barplot(df_importance$importance, 
        names.arg = df_importance$variable, 
        main = "Normalized Importance Score",
        xlab = "Variable", 
        ylab = "Importance",
        cex.lab = 2,
        cex.axis = 2,
        cex.main = 2,
        col = "#e6e6e6"
)
print(normalized_importance)

# Predict using the random forest model
predictions <- predict(rf_model, newdata = test_data)
comparison <- data.frame(Actual = test_data$lai_lidar,
                         Predicted = predictions)

# predictions <- predict(rf_model, newdata = lai_lidar)
# comparison <- data.frame(Actual = lai_lidar,
#                          Predicted = predictions)

# Compute metrics
correlation <- cor(comparison$Actual, comparison$Predicted)
rmse <- Metrics::rmse(comparison$Actual, comparison$Predicted)
rsquared <- R2(comparison$Actual, comparison$Predicted)
bias <- Metrics::bias(comparison$Actual, comparison$Predicted)
mae <- Metrics::mae(comparison$Actual, comparison$Predicted)

# Print metrics
cat("Root Mean Squared Error:", rmse, 
    "\nR-squared:", rsquared, 
    "\nBias:", bias, 
    "\nCorrelation:", correlation, 
    "\nMAE:", mae)

plot_density_scatterplot(comparison$Actual,
                         comparison$Predicted,
                         "Actual LiDAR LAI",
                         "Predicted LiDAR LAI",
                         title = "Predicted LiDAR LAI vs Actual LiDAR LAI Values on the Test Set",
                         xlimits = c(0,10),
                         ylimits = c(0,10))

# Visualize predicted raster
new_raster[non_na_indices] <- predictions
plot_histogram(list(predictions, predictions_ranger$data$response, test_data$lai_lidar, test_data$lai_s2),
               var_labs = c("Predicted", "Predicted max", "Actual", "S2"))
# plot_histogram(list(predictions, predictions_ranger$data$response),
#                var_labs = c("Predicted", "Predicted max"))

# Calculate permutation importance without variable names
calculate_permutation_importance <- function(model, data, target_column_index) {
  perm_importance <- numeric(ncol(data))
  baseline_rmse <- Metrics::rmse(predict(model, data), data[[target_column_index]])
  
  for (i in seq_along(data)) {
    perm_data <- data
    perm_data[[i]] <- sample(perm_data[[i]])
    perm_rmse <- Metrics::rmse(predict(model, perm_data), data[[target_column_index]])
    perm_importance[i] <- baseline_rmse - perm_rmse
  }
  names(perm_importance) <- names(validation_data)
  perm_importance <- sort(perm_importance, decreasing = TRUE)
  return(data.frame(variable = names(perm_importance), 
                    importance = perm_importance))
}

# Calculate permutation importance
perm_importance <- calculate_permutation_importance(rf_model, validation_data, 1)
barplot(perm_importance$importance[-1], 
        names.arg = perm_importance$variable[-1],
        main = "Permutation Importance",
        xlab = "Predictor Variable Index",
        ylab = "Importance")
print(perm_importance)

# Combine rankings into a single dataframe
combined_rankings <- data.frame(
  variable = names(normalized_importance),
  norm_rank = 1:length(normalized_importance),
  perm_rank = match(names(normalized_importance), perm_importance$variable)
)
combined_rankings$perm_rank <- nrow(combined_rankings) +2 - combined_rankings$perm_rank

# Set up colors for the bars
colors <- c("blue", "red")

# Plot the rankings for each variable
barplot(
  rbind(combined_rankings$norm_rank, combined_rankings$perm_rank),
  beside = TRUE,
  names.arg = combined_rankings$variable,
  col = colors,
  main = "Variable Rankings",
  xlab = "Variable",
  ylab = "Ranking"
)
legend("topright", 
       legend = c("Normalized Importance", "Permutation Importance"), 
       fill = colors)

# Bonus Plots
explainer <- DALEX::explain(model = rf_model,
                            data = as.matrix(training_data[-1]),
                            y = training_data$lai_lidar_val,
                            verbose = FALSE)
pdp_plot <- predict_profile(explainer, pred.var = predictor_vars, new_observation = training_data[1, ])
plot(pdp_plot)
complexity_plot <- model_profile(explainer)
plot(complexity_plot)

# 1. Calculate RMSE of the trained model
train_predictions <- predict(rf_model, newdata = training_data)
train_rmse <- sqrt(mean((training_data$lai_lidar - train_predictions)^2))
print(paste("RMSE on training data:", train_rmse))

# 2. Make predictions on the validation data
predictions <- predict(rf_model, newdata = test_data)

# 3. Calculate RMSE on the predicted data
validation_rmse <- sqrt(mean((test_data$lai_lidar - predictions)^2))
print(paste("RMSE on validation data:", validation_rmse))

# ---------------------------- RF Train Test  ----------------------------------
# Define the range of number of trees to evaluate
num_trees <- seq(1, 1000, by = 1)

predictors = training_data[, -which(names(training_data) == "lai_lidar")]
outcome = training_data$lai_lidar

# Initialize vectors to store RMSE values
train_rmse_values <- numeric(length(num_trees))
val_rmse_values <- numeric(length(num_trees))
bias_values <- numeric(length(num_trees))
variance_values <- numeric(length(num_trees))

# CV
ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 5,
                     search = "random")
mtry <- sqrt(ncol(predictors))
tunegrid <- expand.grid(.mtry=c(1:15))
metric <- "RMSE"
# rf_default <- train(x = predictors,
#                     y = outcome, 
#                     method="rf", 
#                     ntree = 15,
#                     metric=metric, 
#                     tuneGrid=tunegrid, 
#                     trControl=ctrl)
# print(rf_default)
# plot(rf_default)

# Train models and store RMSE values
ctrl <- trainControl(method = "LOOCV",
                     number = 3)
# Define random forest model
# rf_model2 <- train(
#   x = predictors,
#   y = outcome,
#   method = "rf",
#   ntree = num_trees[i], 
#   trControl = ctrl,
#   tuneGrid = data.frame(mtry = 4),
#   do.trace = TRUE
# )

rf_model <- randomForest(formula,
                         data = training_data,
                         ntree = 1000,  # or any other number of trees
                         mtry = mtry,
                         importance = TRUE,
                         do.trace = TRUE)

# Initialize an empty list to store subset models
subset_models <- list()

# Iterate through different numbers of trees
for (i in seq_along(num_trees)) {
  print(i)
  # Extract a subset of trees from the forest
  rf_model_sub <- randomForest::combine(rf_model, k = num_trees[i])
  subset_models[[i]] <- rf_model_sub
  
  # Make predictions using the subset of trees
  predictions_train <- predict(rf_model_sub, newdata = training_data)
  predictions_val <- predict(rf_model_sub, newdata = test_data)
  
  # Calculate performance metrics
  train_rmse_values[i] <- RMSE(predictions_train, training_data$lai_lidar)
  val_rmse_values[i] <- RMSE(predictions_val, test_data$lai_lidar)
  bias_values[i] <- Metrics::bias(predictions_train, training_data$lai_lidar)
  variance_values[i] <- var(predictions_train)
}

# Print the performance metrics for each tree
for (i in seq_along(num_trees)) {
  cat("Number of trees:", num_trees[i], "\n")
  cat("Train RMSE:", train_rmse_values[i], "\n")
  cat("Validation RMSE:", val_rmse_values[i], "\n")
  cat("Bias:", bias_values[i], "\n")
  cat("Variance:", variance_values[i], "\n\n")
}

# Plot RMSE values against number of trees
plot(num_trees, train_rmse_values, type = "l", col = "#51a343", 
     ylim = c(0,1),
     xlab = "Number of Trees", ylab = "RMSE",
     main = paste("RMSE vs. Number of Trees\n",
                  "K-fold Repeated Cross Validation\n",
                  "5 folds / 5 repeats"),
     cex.lab = 1.5,    # Normal font size for axis labels
     cex.axis = 1.5,
     cex.main = 1.5, 
     lwd = 2
     # main = "RMSE vs. Number of Trees"
)
lines(num_trees, val_rmse_values, type = "l", col = "#9a4c00", lwd = 2)
abline(v = 15, col = "red", lwd = 2)
legend("topright", 
       legend = c("Train", "Test", "NTree Choice"), 
       col = c("#51a343", "#9a4c00", "red"), 
       lty = c(1, 1, 1),
       lwd = c(2, 2, 2),
       cex = 1.4
)
# title(main = paste("RMSE vs. Number of Trees\n",
#                    "K-fold Repeated Cross Validation\n",
#                    "5 folds / 5 repeats"))
# legend("topright", 
#        legend = c("Train", "Test"), 
#        col = c("#51a343", "#9a4c00"), 
#        lty = 1)

# ------------------------------- Simple RF (ES) -------------------------------

max_trees <- 100
max_iter_without_improvement <- 10
min_delta_error <- 0.0001

# Initialize variables
best_error <- Inf
iter_without_improvement <- 0

# Training loop
for (i in 1:max_trees) {
  # Train random forest model
  rf_model <- randomForest(
    lai_lidar_val ~ ., 
    data = training_data, 
    ntree = i,  
    mtry = sqrt(ncol(training_data) - 1),  
    importance = TRUE,
    do.trace = TRUE
  )
  
  # Compute out-of-bag error rate
  oob_error <- rf_model$mse[length(rf_model$mse)]
  
  # Check if the change in error rate is below threshold
  delta_error <- best_error - oob_error
  
  if (delta_error < min_delta_error) {
    iter_without_improvement <- iter_without_improvement + 1
  } else {
    best_error <- oob_error
    iter_without_improvement <- 0
  }
  
  # Check early stopping condition
  if (iter_without_improvement >= max_iter_without_improvement) {
    print("Early stopping due to lack of improvement.")
    break
  }
}

# Print the number of trees used
num_trees <- length(rf_model$mse)
print(paste("Number of trees:", num_trees))
print(summary(rf_model))

# Normalize importance scores
importance <- varImp(rf_model)
normalized_importance <- importance$Overall / sum(importance$Overall) * 100
names(normalized_importance) <- rownames(importance)
df_importance <- data.frame(variable = names(normalized_importance), 
                            importance = normalized_importance)
barplot(df_importance$importance, 
        names.arg = df_importance$variable, 
        main = "Normalized Importance Scores", 
        xlab = "Variable", 
        ylab = "Importance"
)
print(normalized_importance)

# Extract MSE and Variance for each tree
mse_values <- rf_model$mse
# variance_values <- 100 - (rf_model$rsq * 100)  # Correct calculation for Variance

# Plot MSE and Variance for each tree
tree_numbers <- 1:length(mse_values)

# Plot MSE values according to tree number
plot(tree_numbers, mse_values, type = "l", col = "blue", xlab = "Tree", ylab = "MSE", main = "MSE by Tree")

# -------------------------------- RF Tuning -----------------------------------
# Define grid of hyperparameters
param_grid <- expand.grid(
  ntree = c(10, 15, 20, 30, 100, 200, 500, 1000),
  mtry = c(3, 4, 5, 6, 7, 8),
  nodesize_range <- c(1, 5, 10),
  sampsize_range <- c(50, 100, 150),
  splitrule_range <- c("gini", "extratrees"),
  stringsAsFactors = FALSE
)

# Perform grid search
best_rmse <- Inf
best_model <- NULL
train_control <- trainControl(method = "cv", number = 5)

# Formula
formula <- as.formula(lai_lidar ~ 
                        lai_s2 +
                        mean_h +
                        max_h +
                        std +
                        cv +
                        variance +
                        rumple +
                        vci +
                        lcv +
                        lskew +
                        vdr)

# Perform grid search
for (i in 1:nrow(param_grid)) {
  # Train random forest model
  rf_model <- train(
    formula,
    data = training_data,
    method = "rf",
    trControl = train_control,
    tuneGrid = data.frame(
      ntree = param_grid$ntree[i],
      mtry = param_grid$mtry[i],
      nodesize = param_grid$nodesize[i],
      sampsize = param_grid$sampsize[i],
      splitrule = param_grid$splitrule[i]
    )
  )
  
  # Evaluate model
  predictions <- predict(rf_model, newdata = test_data)
  rmse <- Metrics::rmse(test_data$lai_lidar, predictions)
  
  # Check if current model is the best
  if (rmse < best_rmse) {
    best_rmse <- rmse
    best_model <- rf_model
  }
}

# Evaluate performance of the best model
cat("Best mtry:", best_mtry, 
    "\nBest ntree:", best_ntree,
    "\nBest RMSE:", best_rmse, "\n")
best_model_performance <- predict(best_model, test_data)

# ---------------------------- RF Tuning Custom --------------------------------
customRF_regression <- list(type = "Regression", 
                            library = "randomForest", 
                            loop = NULL)
customRF_regression$parameters <- data.frame(parameter = c("mtry", 
                                                           "ntree",
                                                           "nodesize",
                                                           "sampsize",
                                                           "splitrule"), 
                                             class = rep("numeric", 5), 
                                             label = c("mtry", 
                                                       "ntree",
                                                       "nodesize",
                                                       "sampsize",
                                                       "splitrule"
                                             )
)
customRF_regression$grid <- function(x, 
                                     y,
                                     len = NULL, 
                                     search = "grid") {}
customRF_regression$fit <- function(x, 
                                    y, 
                                    wts,
                                    param, 
                                    lev, 
                                    last,
                                    weights, 
                                    classProbs, 
                                    ...) {
  randomForest(x, 
               y, 
               mtry = param$mtry,
               ntree=param$ntree, 
               ...)
}
customRF_regression$predict <- function(modelFit, 
                                        newdata, 
                                        preProc = NULL, 
                                        submodels = NULL)
  predict(modelFit, newdata)
customRF_regression$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) NULL
customRF_regression$sort <- function(x) x[order(x[,1]),]
customRF_regression$levels <- function(x) x$classes

# Train model
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(.mtry=c(3, 4, 5, 6, 7, 8), 
                        .ntree=c(10, 15, 20, 30, 100, 200, 500, 1000),
                        .nodesize=c(1, 5, 10),
                        .sampsize=c(50, 100, 150),
                        .splitrule=c("gini", "extratrees")
)
start_time <- Sys.time()
end_time <- Sys.time()
execution_time <- end_time - start_time
execution_time
custom_regression <- caret::train(x = predictors,
                                  y = outcome,
                                  method = customRF_regression,
                                  metric = "RMSE",
                                  trControl = control,
                                  tuneGrid = tunegrid,
                                  do.trace = TRUE
)
summary(custom_regression)
plot(custom_regression)

# ------------------------------- Regression -----------------------------------
# Create formula
formula <- lai_lidar_val ~ 
  lai_s2_val +
  mean_h_val + 
  max_h_val +
  std_val +
  cv_val +
  variance_val +
  rumple_val +
  vci_val +
  lcv_val +
  lskew_val +
  vdr_val

# Fit the multivariate regression model
lm_model <- lm(formula, data = training_data)

# Summary of the model
summary(lm_model)

# Make predictions on the test dataset
test_predictions <- predict(lm_model, newdata = validation_data)

# Calculate RMSE (Root Mean Squared Error) to evaluate model performance
rmse <- Metrics::rmse(validation_data$lai_lidar_val, test_predictions)
print(paste("RMSE:", rmse))

# Plot actual vs predicted values
plot(validation_data$lai_lidar_val, test_predictions, 
     xlab = "Actual LAI", ylab = "Predicted LAI",
     main = "Actual vs Predicted LAI")

# Add a diagonal line to represent perfect prediction
abline(0, 1, col = "red")

# Perform stepwise regression
stepwise_model <- step(lm_model, direction = "both")

# Summary of the stepwise model (AIC)
summary(stepwise_model)
# --------------------------------- GLM ----------------------------------------
# Initialization
training_data <- as.data.frame(training_data)

# Create formula
formula <- lai_lidar_val ~ 
  lai_s2_val + 
  mean_h_val + 
  std_val + 
  cv_val +
  variance_val + 
  rumple_val +
  vci_val + 
  lcv_val +
  lskew_val
vdr_val

# Fit the GLM
glm_model <- glm(formula, 
                 data = training_data, 
                 family = gaussian(link = "identity"))

# Evaluate model performance
model_performance <- caret::trainControl(method = "cv", number = 5) # 5-fold cross-validation
model_results <- caret::train(formula, 
                              data = training_data, 
                              method = "glm", 
                              trControl = model_performance)

# Print summary of the model
summary(glm_model)

# Plot diagnostic plots
plot(glm_model)

# Plot feature importance (if available)
feature_importance <- varImp(glm_model)
print(feature_importance)
