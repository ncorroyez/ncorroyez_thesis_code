# ---
# title: "1.create_masks.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-27"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("../Sentinel_2/functions_sentinel_2.R")
source("../functions_plots.R")
source("functions_lidar.R")

# Pre-processing Parameters & Directories
data_dir <- '../../01_DATA'
results_dir <- '../../03_RESULTS'
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual

choice_equal_deciles <- "deciles" #equal_intervals #deciles
# choice_equal_deciles <- c("deciles", "equal_intervals")

for (site in sites){
  # Which mask
  composition_masks <- list("Full_Composition",
                            "Deciduous_Flex",
                            "Deciduous_Only")
  for (composition_mask in composition_masks){
    
    # Setup
    results_path <- file.path(results_dir, site)
    masks_dir <- file.path(results_path, 
                           "LiDAR/Heterogeneity_Masks/Quantiles",
                           composition_mask)
    if (!dir.exists(masks_dir)) {
      dir.create(masks_dir, showWarnings = FALSE, recursive = TRUE)
    }
    metrics_dir <- file.path(results_path, "Metrics")
    composition_mask_dir <- file.path(metrics_dir, composition_mask)
    
    # Open heterogeneity metrics 
    mean <- terra::rast(file.path(composition_mask_dir, "mean_res_10_m.tif"))
    max <- terra::rast(file.path(composition_mask_dir, "max_res_10_m.tif"))
    std <- terra::rast(file.path(composition_mask_dir, "std_res_10_m.tif"))
    cv <- terra::rast(file.path(composition_mask_dir, "cv_res_10_m.tif"))
    rumple <- terra::rast(file.path(composition_mask_dir, 
                                    "rumple_res_10_m.tif"))
    heterogeneity_metrics <- list(mean, max, std, cv, rumple)
    
    save_x_y_plot(xvar = mean,
                  yvar = std,
                  dirname = masks_dir,
                  filename = "std_vs_mean_heights.png",
                  xlab = "Mean Heights",
                  ylab = "Standard Deviation",
                  title = "Standard Deviation vs Mean Heights")
    save_x_y_plot(xvar = mean,
                  yvar = cv,
                  dirname = masks_dir,
                  filename = "cv_vs_mean_heights.png",
                  xlab = "Mean Heights",
                  ylab = "Coefficient of Variation",
                  title = "Coefficient of Variation vs Mean Heights")
    save_x_y_plot(xvar = std,
                  yvar = cv,
                  dirname = masks_dir,
                  filename = "cv_vs_std.png",
                  xlab = "Standard Deviation",
                  ylab = "Coefficient of Variation",
                  title = "Coefficient of Variation vs Standard Deviation")
    save_x_y_plot(xvar = std,
                  yvar = rumple,
                  dirname = masks_dir,
                  filename = "rumple_vs_std.png",
                  xlab = "Standard Deviation",
                  ylab = "Rumple Index",
                  title = "Rumple Index vs Standard Deviation")
    save_x_y_plot(xvar = mean,
                  yvar = rumple,
                  dirname = masks_dir,
                  filename = "rumple_vs_mean.png",
                  xlab = "Standard Deviation",
                  ylab = "Rumple Index",
                  title = "Rumple Index vs Standard Deviation")
    
    for (heterogeneity_metric_raster in heterogeneity_metrics){
      if (identical(heterogeneity_metric_raster, std)) {
        cat("Error: Standard Deviation\n")
        heterogeneity_metric_name <- "std"
      }
      else if (identical(heterogeneity_metric_raster, cv)) {
        cat("Error: Coefficient of Variation\n")
        heterogeneity_metric_name <- "cv"
      }
      else if (identical(heterogeneity_metric_raster, mean)) {
        cat("Error: Mean Heights\n")
        heterogeneity_metric_name <- "mean_heights"
      }
      else if (identical(heterogeneity_metric_raster, rumple)) {
        cat("Error: Rumple\n")
        heterogeneity_metric_name <- "rumple"
      }
      else if (identical(heterogeneity_metric_raster, max)) {
        cat("Error: Max\n")
        heterogeneity_metric_name <- "max"
      }
      else{
        stop("Error is not Standard Deviation or Coefficient of Variation
           or Mean or Max or Rumple")
      }
      for (choice_equal_decile in choice_equal_deciles){
        create_heterogeneity_quantiles(heterogeneity_metric_raster,
                                       heterogeneity_metric_name,
                                       site,
                                       masks_dir,
                                       composition_mask_dir,
                                       choice_equal_decile,
                                       resolution = 10)
      }
    }
  }
}




# Heterogeneity Mask

# 
# plot(mean_heights, standard_deviation_proj, main = "Standard Deviation vs Mean Heights")
# plot(mean_heights, coeff_variation, main = "Coefficient of Variation vs Mean Heights")
# plot(standard_deviation_proj, coeff_variation, main = "Coefficient of Variation vs Standard Deviation")
# 
# 
# 
# # Calculate deciles on heterogeneity raster
# errors <- list(standard_deviation_proj, coeff_variation, mean_heights)
# # errors <- list(mean_heights)
# 
# 
# # Test Deciles
# correlation_values <- c()
# for (low_value in seq(range[1], range[2], by = inc)) {
#   low_value <- round(low_value, digits = 2)
#   high_value <- round(low_value + inc, digits = 2)
#   lidar_intervals <- terra::rast(file.path(class_dir, 
#                                            paste0(sprintf("%s_pai_masked_res_%s_m_", 
#                                                           err, 
#                                                           resolution),
#                                                   low_value,  "_", 
#                                                   high_value, ".tif")))
#   lidar_values <- values(lidar_intervals)
#   
#   s2_intervals <- terra::rast(file.path(class_dir, 
#                                         paste0(sprintf("%s_lai_s2_masked_res_%s_m_", 
#                                                        err, 
#                                                        resolution),
#                                                low_value, "_", 
#                                                high_value, ".tif")))
#   s2_values <- values(s2_intervals)
#   
#   correlation_value <- correlation_test_function(lidar_values, s2_values, method = "pearson")
#   correlation_values <- c(correlation_values, correlation_value)
#   
#   plot_density_scatterplot(var_x = s2_values,
#                            var_y = lidar_values,
#                            xlab = "LAI PROSAIL",
#                            ylab = "LAI LiDAR",
#                            dirname = class_dir,
#                            filename = paste0(sprintf("%s_scatter_lai_lidar_prosail_res_%s_m_", 
#                                                      err, resolution),
#                                              low_value, "_", 
#                                              high_value, ".tif"),
#                            xlimits = c(0, 10),
#                            ylimits = c(0, 10)
#   )
# }
# first_range <- sapply(quantiles_list, function(x) x[1])
# quantile_labels <- c(first_range, max(values(error), na.rm=T))
# quantile_labels <- round(quantile_labels, 1)
# png(filename = file.path(class_dir, 
#                          sprintf("z_%s_correlation_values_for_%s_res_%s_m.png", 
#                                  err, choice, resolution)),
#     width = 1920, height = 1080, 
#     units = "px", res = 200
# )
# plot(1, type = "l", 
#      xlim = c(0, 1), ylim = c(0, 0.65),
#      xlab = choice_equal_deciles_for_plot, 
#      ylab = "Correlation Value", 
#      # main = "Correlation Values between LiDAR LAI and Sentinel-2 LAI for 10 Mean Heights Classes",
#      # main = paste("Correlation Values between LiDAR LAI and\n",
#      #              "Sentinel-2 LAI for 10 Mean Heights Classes"),
#      lwd = 2,
#      xaxt = "n",
#      cex.main = 1.2,   # Increase title font size
#      cex.lab = 1.2,    # Normal font size for axis labels
#      cex.axis = 1.2
# )
# title(main = paste("Pearson Correlation between LiDAR LAI and Sentinel-2 LAI\n",
#                    "for 10 Height Classes (Defined as the CHM Mean Heights\n",
#                    "(1m resolution)) over Sentinel-2 pixels (10m resolution)"))
# lines(midpoints, correlation_values, type = "l")
# axis(side = 1, at = seq(0, 1, by = 0.1), labels = quantile_labels)
# points(midpoints, correlation_values, pch = 16, col = "red")
# # axis(side = 1, at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
# dev.off()
