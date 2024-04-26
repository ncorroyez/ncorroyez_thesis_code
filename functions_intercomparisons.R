# ---
# title: "functions_intercomparisons.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-06"
# ---

library(terra)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(hrbrthemes)
library(raster)
library(rasterVis)
library(latticeExtra)
library(RColorBrewer)

create_dataframe <- function(bands, 
                             noises,
                             row_names,
                             col_names){
  
  data_matrix <- matrix(NA, nrow = length(row_names), ncol = length(col_names))
  df <- data.frame(data_matrix)
  
  rownames(df) <- paste(row_names, col_names, sep = "_")
  colnames(df) <- paste(row_names, col_names, sep = "_")
  return(df)
}

create_atbd_bands_noises_vars_list <- function(bands,
                                               noises,
                                               mask,
                                               plots_dir,
                                               s2_output_dir,
                                               image,
                                               estimated_var){
  
  # Create a list to store data matrices
  data_list <- list()
  
  # Load your data into the list, assuming you have them named as: 
  # '10m_addmult', '20m_addmult', etc.
  for (noise in noises) {
    for (band in bands) {
      data_path <- file.path(s2_output_dir,
                             'PRO4SAIL_INVERSION_atbd',
                             noise,
                             band,
                             paste0(image, 
                                    '_', 
                                    'Refl',
                                    '_',
                                    estimated_var, 
                                    ".envi"))
      data_name <- paste0(noise, "_", band)
      data <- terra::rast(data_path) 
      data <- terra::project(data, mask)
      data <- data * mask
      data <- values(data)
      data_list[[data_name]] <- data
    }
  }
  return(data_list)
}

create_atbd_bands_noises_rasters_list <- function(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var){
  
  # Create a list to store data matrices
  data_list <- list()
  
  # Load your data into the list, assuming you have them named as: 
  # '10m_addmult', '20m_addmult', etc.
  for (noise in noises) {
    for (band in bands) {
      data_path <- file.path(s2_output_dir,
                             'PRO4SAIL_INVERSION_atbd',
                             noise,
                             band,
                             paste0(image, 
                                    '_', 
                                    'Refl',
                                    '_',
                                    estimated_var, 
                                    ".envi"))
      data_name <- paste0(noise, "_", band)
      data <- terra::rast(data_path)
      data <- terra::project(data, mask)
      data <- data * mask
      data_list[[data_name]] <- data
    }
  }
  return(data_list)
}

save_df <- function(df, df_path_to_save, df_name_to_save){
  if (is.null(df_path_to_save) || is.null(df_name_to_save)) {
    stop("DF Path or Name cannot be NULL.")
  }
  write.csv(format(df, digits = 3), 
            file.path(df_path_to_save, 
                      paste0(df_name_to_save, ".csv")),
            row.names = TRUE)
  cat("DF has been successfully saved at", 
      paste0(df_name_to_save, ".csv"), "\n")
}

correlation_atbd_bands_noise <- function(bands,
                                         noises,
                                         row_names,
                                         col_names,
                                         mask,
                                         plots_dir,
                                         s2_output_dir,
                                         image,
                                         estimated_var,
                                         comparison_path,
                                         csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  cor_matrix <- cor(do.call(cbind, data_list), use = "pairwise.complete.obs")
  
  for (i in 1:length(row_names)) {
    for (j in 1:length(col_names)) {
      df[i, j] <- cor_matrix[i, j]
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

rmse_atbd_bands_noise <- function(bands,
                                  noises,
                                  row_names,
                                  col_names,
                                  mask,
                                  plots_dir,
                                  s2_output_dir,
                                  image,
                                  estimated_var,
                                  comparison_path,
                                  csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      rmse <- Metrics::rmse(var1, var2)
      df[i, j] <- rmse
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

bias_atbd_bands_noise <- function(bands,
                                  noises,
                                  row_names,
                                  col_names,
                                  mask,
                                  plots_dir,
                                  s2_output_dir,
                                  image,
                                  estimated_var,
                                  comparison_path,
                                  csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      bias <- Metrics::bias(var1, var2)
      df[i, j] <- bias
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

residuals_atbd_bands_noise <- function(bands,
                                       noises,
                                       row_names,
                                       col_names,
                                       mask,
                                       plots_dir,
                                       s2_output_dir,
                                       image,
                                       estimated_var,
                                       comparison_path){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_rasters_list(bands,
                                                     noises,
                                                     mask,
                                                     plots_dir,
                                                     s2_output_dir,
                                                     image,
                                                     estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      if(!identical(data_list[[i]], data_list[[j]])){
        var1 <- data_list[[i]]  # Use double brackets to access list elements
        var2 <- data_list[[j]]  # Use double brackets to access list elements
        residuals <- var1 - var2
        save_envi_file(residuals,
                       paste0(rownames(df)[i], 
                              "_", 
                              colnames(df)[j], 
                              "_residuals"),
                       file.path(comparison_path, 
                                 "residuals",
                                 paste0(rownames(df)[i], 
                                        "_", 
                                        colnames(df)[j])))
      }
    }
  }
}

error_atbd_bands_noise <- function(bands,
                                   noises,
                                   row_names,
                                   col_names,
                                   mask,
                                   plots_dir,
                                   s2_output_dir,
                                   image,
                                   estimated_var,
                                   comparison_path,
                                   csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      residuals <- var1 - var2
      error <- mean(abs(residuals)) 
      df[i, j] <- error
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

r_squared_atbd_bands_noise <- function(bands,
                                       noises,
                                       row_names,
                                       col_names,
                                       mask,
                                       plots_dir,
                                       s2_output_dir,
                                       image,
                                       estimated_var,
                                       comparison_path,
                                       csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      r_squared <- cor(var1, var2)^2
      df[i, j] <- r_squared
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

std_atbd_bands_noise <- function(bands,
                                 noises,
                                 row_names,
                                 col_names,
                                 mask,
                                 plots_dir,
                                 s2_output_dir,
                                 image,
                                 estimated_var,
                                 comparison_path,
                                 csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      residuals <- var1 - var2
      std <- sd(residuals)
      df[i, j] <- std
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

histograms_scatterplots_atbd_bands_noise <- function(bands,
                                                     noises,
                                                     row_names,
                                                     col_names,
                                                     mask,
                                                     plots_dir,
                                                     s2_output_dir,
                                                     image,
                                                     estimated_var,
                                                     comparison_path){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      if(!identical(data_list[[i]], data_list[[j]])){
        if((names(data_list)[i] == "addmult_10m" 
            && names(data_list)[j] == "addmult_20m")
           || (names(data_list)[i] == "mult_10m" 
               && names(data_list)[j] == "mult_20m")
           || (names(data_list)[i] == "addmult_10m" 
               && names(data_list)[j] == "mult_10m")
           || (names(data_list)[i] == "addmult_20m" 
               && names(data_list)[j] == "mult_20m")
        ) {
          residuals <- data_list[[i]] - data_list[[j]]
          vars_list <- list(data_list[[i]], data_list[[j]], residuals)
          var_labs <- c(names(data_list)[i], names(data_list)[j], "residuals")
          plot_histogram(vars_list, 
                         "LAI Distributions", 
                         "LAI",
                         var_labs, 
                         comparison_path, 
                         paste0(
                           "hist",
                           "_",
                           names(data_list)[i],
                           "_",
                           names(data_list)[j],
                           "_",
                           "residuals"
                         ))
          plot_density_scatterplot(var_x = data_list[[i]],
                                   var_y = data_list[[j]],
                                   xlab = names(data_list)[j],
                                   ylab = names(data_list)[i],
                                   dirname = comparison_path, 
                                   filename = paste0(
                                     "scatterplot",
                                     "_",
                                     names(data_list)[i],
                                     "_",
                                     names(data_list)[j],
                                     "_"))
        }
      }
    }
  }
}

calculate_metrics <- function(actual, 
                              predicted, 
                              save_path,
                              actual_lab = NULL,
                              predicted_lab = NULL
) {
  if (is.null(actual_lab)){
    actual_lab <- "actual"
  }
  if (is.null(predicted_lab)){
    predicted_lab <- "predicted"
  }
  
  # Take values
  actual_values <- values(actual)
  predicted_values <- values(predicted)
  
  # Residuals
  residuals_raster <- actual - predicted
  residuals <- values(residuals_raster)
  
  # Check if the lengths of actual and predicted are the same
  if (length(actual) != length(predicted)) {
    stop("Lengths of actual and predicted values must be the same.")
  }
  
  hist_list <- list(actual_values, predicted_values, residuals)
  
  # Plot histogram: actual, predicted, residuals
  plot_histogram(vars_list = hist_list, 
                 title = "LAI Distributions", 
                 xlab = "LAI", 
                 var_labs = sprintf(c(actual_lab, predicted_lab, "Residuals")), 
                 dirname = save_path, 
                 filename = sprintf("hist_%s_%s_residuals",
                                    actual_lab,
                                    predicted_lab))
  
  plot_density_scatterplot(var_x = predicted_values,
                           var_y = actual_values,
                           xlab = sprintf("%s", predicted_lab), 
                           ylab = sprintf("%s", actual_lab),
                           dirname = save_path, 
                           filename = sprintf("scatterplot_%s_%s",
                                              actual_lab,
                                              predicted_lab),
                           xlimits = c(0, max(predicted_values)),
                           ylimits = c(0, max(actual_values)))
  
  # Remove NAs from both actual and predicted
  valid_indices <- complete.cases(actual_values, predicted_values)
  actual_values <- actual_values[valid_indices]
  predicted_values <- predicted_values[valid_indices]
  residuals <- residuals[valid_indices]
  
  # Check if there are remaining values after removing NAs
  if (length(actual) == 0) {
    stop("No valid data points remaining after removing NAs.")
  }
  
  # Calculate metrics
  correlation <- cor(actual_values, predicted_values)
  ssres <- sum((residuals)^2)
  sstot <- sum((actual_values - mean(actual_values))^2)
  print(ssres)
  print(sstot)
  r_squared <- 1 - (ssres / sstot)
  rmse <- sqrt(mean(residuals^2))
  bias <- mean(residuals)
  error <- mean(abs(residuals))
  std_dev <- sd(residuals)
  
  # Store metrics in a data frame
  metrics_df <- data.frame(
    R = correlation,
    R_Squared = r_squared,
    RMSE = rmse,
    Bias = bias,
    Error = error,
    Std_Dev = std_dev
  )
  
  # Save residuals map
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  writeRaster(residuals_raster, 
              file.path(save_path, sprintf("residuals_%s_%s.envi",
                                           actual_lab,
                                           predicted_lab)),
              overwrite = TRUE)
  
  # Save metrics to a CSV file
  write.csv(format(metrics_df, digits = 3), 
            file = file.path(save_path, sprintf("metrics_%s_%s.csv",
                                                actual_lab,
                                                predicted_lab)),
            row.names = FALSE)
  
  # Print the metrics data frame
  print(metrics_df)
  
  # Return the metrics data frame
  return(metrics_df)
}

calculate_metrics_from_a_reference <- function(reference_raster_path,
                                               rasters_dir,
                                               image,
                                               noises,
                                               bands,
                                               mask,
                                               valid_indices,
                                               estimated_var,
                                               reference_type,
                                               comparison_path,
                                               csv_name) {
  # Create an empty list to store the metrics data frames
  metrics_list <- list()
  
  reference <- open_raster_file(reference_raster_path)
  reference <- terra::project(reference, mask)
  reference <- reference * mask
  reference_values <- values(reference)
  
  for (noise in noises){
    for (band in bands){
      # Raster to compare with reference
      raster <- open_raster_file(
        file.path(rasters_dir, 
                  "PRO4SAIL_INVERSION_atbd",
                  noise,
                  band,
                  paste0(image, 
                         '_', 
                         'Refl',
                         '_',
                         estimated_var, 
                         ".envi")))
      raster <- terra::project(raster, mask)
      raster <- raster * mask
      raster_values <- values(raster)
      
      print(length(reference_values))
      print(length(raster_values))
      
      # Reference
      residuals_raster <- reference - raster
      residuals <- values(residuals_raster)
      
      # Remove NAs
      # valid_indices <- complete.cases(reference_values, raster_values)
      # valid_indices <- complete.cases(raster_values)
      reference_valuess <- reference_values[valid_indices]
      raster_values <- raster_values[valid_indices]
      residuals <- residuals[valid_indices]
      
      print(length(reference_valuess))
      print(length(raster_values))
      
      # Calculate metrics
      correlation <- cor(reference_valuess, raster_values)
      ssres <- sum((residuals)^2)
      sstot <- sum((reference_valuess - mean(reference_valuess))^2)
      # print(ssres)
      # print(sstot)
      r_squared <- 1 - (ssres / sstot)
      rmse <- sqrt(mean(residuals^2))
      bias <- mean(residuals)
      error <- mean(abs(residuals))
      std_dev <- sd(residuals)
      
      # Plots
      plot_histogram(list(reference_valuess, raster_values, residuals),
                     "LAI Distributions",
                     "LAI",
                     c(reference_type, 
                       sprintf("atbd_%s_%s", noise, band), 
                       "Residuals"),
                     comparison_path,
                     paste0(
                       "hist",
                       "_",
                       reference_type,
                       "_",
                       sprintf("atbd_%s_%s", noise, band),
                       "_",
                       "residuals"
                     ),
                     c(-5,10))
      plot_density_scatterplot(var_x = raster_values,
                               var_y = reference_valuess,
                               xlab = sprintf("atbd_%s_%s", noise, band),
                               ylab = reference_type,
                               dirname = comparison_path,
                               filename = paste0(
                                 "scatterplot",
                                 "_",
                                 reference_type,
                                 "_",
                                 sprintf("atbd_%s_%s", noise, band),
                                 "_"),
                               xlimits = c(0,10),
                               ylimits = c(0,10))
      save_envi_file(residuals_raster,
                     paste0(reference_type, 
                            "_", 
                            sprintf("atbd_%s_%s", noise, band), 
                            "_residuals"),
                     file.path(comparison_path, 
                               "residuals",
                               paste0(reference_type, 
                                      "_", 
                                      sprintf("atbd_%s_%s", noise, band))))
      
      # Create a data frame for the metrics
      metrics_df <- data.frame(
        Ref_Metrics = paste0(reference_type, "_atbd_", noise, "_", band),
        R = correlation,
        R_squared = r_squared,
        RMSE = rmse,
        Bias = bias,
        Error = error,
        Std_Dev = std_dev,
        stringsAsFactors = FALSE
      )
      
      # Append the data frame to the list
      metrics_list[[noise]] <- rbind(metrics_list[[noise]], metrics_df)
    }
  }
  
  # Combine all data frames into one
  metrics_df <- do.call(rbind, metrics_list)
  
  # Print the metrics data frame
  print(metrics_df)
  
  # Export the data frame to a tabular format (e.g., CSV)
  save_df(metrics_df, 
          comparison_path,
          csv_name)
}
