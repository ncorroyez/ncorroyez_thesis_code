# ---
# title: "functions_explain_heterogeneity"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-14"
# ---

standardize_train_data <- function(data) {
  standardized_data <- data
  means <- vector("numeric", length = ncol(data))
  sds <- vector("numeric", length = ncol(data))
  
  for (id_var in seq_along(data)) {
    train_mean <- mean(data[, id_var], na.rm = TRUE)
    train_sd <- sd(data[, id_var], na.rm = TRUE)
    
    means[id_var] <- train_mean
    sds[id_var] <- train_sd
    
    standardized_data[, id_var] <- (data[, id_var] - train_mean) / train_sd
  }
  
  return(list(standardized_data = standardized_data,
              means = means,
              sds = sds)
  )
}

standardize_test_data <- function(data, means, sds) {
  standardized_data <- data
  
  for (id_var in seq_along(data)) {
    standardized_data[, id_var] <- (data[, id_var] - means[id_var]) / sds[id_var]
  }
  
  return(standardized_data)
}

unstandardize_variable <- function(variable, mean_val, sd_val) {
  unstandardized_variable <- (variable * sd_val) + mean_val
  return(unstandardized_variable)
}

initialize_predicted_raster <- function(masks_dir){
  # Predicted raster will be LiDAR LAI
  lai_lidar_raster <- terra::rast(file.path(masks_dir,
                                            "lai_lidar_masked_res_10_m.envi")) 
  lai_lidar_matrix <- as.matrix(lai_lidar_raster)
  na_indices <- which(is.na(lai_lidar_matrix), arr.ind = TRUE)
  non_na_indices <- which(!is.na(lai_lidar_matrix), arr.ind = TRUE)
  na_indices <- na_indices[,1]
  non_na_indices <- non_na_indices[,1]
  predicted_raster <- rast(lai_lidar_raster)
  return(predicted_raster)
}


create_train_test_datasets <- function(predictors_values, 
                                       total_sample_size,
                                       training_sample_size,
                                       test_sample_size){
  
  # Sample size verification
  if(training_sample_size + test_sample_size > total_sample_size){
    stop(paste("Error: Total sample size is:", total_sample_size, ",",
               "it can't be smaller than the sum of the training sample size",
               "and the testing sample size (respectively", 
               training_sample_size, "and", test_sample_size, ".\n"))
  }
  # Visualize correlation matrix
  df_predictors <- as.data.frame(predictors_values)
  correlation_matrix <- cor(df_predictors[, 
                                          -which(
                                            names(df_predictors)
                                            == "lai_lidar")])
  corrplot::corrplot(correlation_matrix, method = "number", type = "upper")
  
  # Train-test split
  
  # Subset training_data using the row indices
  # train_prop <- 0.7
  # train_indices <- createDataPartition(sampled_data$lai_lidar, 
  #                                      p = train_prop, 
  #                                      list = FALSE)
  # training_data <- sampled_data[train_indices, ]
  # test_data <- sampled_data[-train_indices, ]
  # test_data <- standardized_data
  
  # Sampling
  sample_indices <- sample(nrow(df_predictors), total_sample_size)
  sampled_data <- df_predictors[sample_indices, ]
  
  # Split into training and testing sets
  indices <- sample(1:nrow(sampled_data), training_sample_size)
  training_data <- sampled_data[indices, ]
  test_data <- sampled_data[-indices, ]
  
  # Standardize training data
  train_data_info <- standardize_train_data(training_data)
  standardized_training_data <- train_data_info$standardized_data
  means <- train_data_info$means
  sds <- train_data_info$sds
  
  # Standardize test data using mean and sd from training data
  standardized_test_data <- standardize_test_data(test_data, means, sds)
  
  # Check dimensions
  cat("Dimensions of training_data:", dim(standardized_training_data),
      "\nDimensions of test_data:", dim(standardized_test_data), "\n")
  
  return(list(training_data = standardized_training_data,
              test_data = standardized_test_data,
              lai_lidar_train_mean = train_data_info$means[1],
              lai_lidar_train_sd = train_data_info$sds[1]
  )
  )
}