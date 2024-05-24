# ---
# title: "1.create_masks.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-26"
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
site <- "Aigoual" # Mormal Blois Aigoual
results_path <- file.path("../../03_RESULTS", site)
masks_dir <- file.path(results_path, "LiDAR/Heterogeneity_Masks")
dir.create(path = masks_dir, showWarnings = F, recursive = T)

if (site == "Mormal"){
  image <- "L2A_T31UER_A031222_20210614T105443"
  maskname <- "Shape/mormal_utm.shp"
  dateAcq <- '2021-06-14' # yyyy-mm-dd format mandatory
} else if (site == "Blois"){
  image <- "L2A_T31TCN_A031222_20210614T105443"
  maskname <- "Shape/blois_utm.shp"
  dateAcq <- '2021-06-14' # yyyy-mm-dd format mandatory
} else if (site == "Aigoual"){
  image <- "L2A_T31TEJ_A031608_20210711T104217"
  maskname <- "Shape/aigoual_utm.shp"
  dateAcq <- '2021-07-11' # yyyy-mm-dd format mandatory
} else{
  stop("Error: Site must be Mormal, Blois, or Aigoual.\n")
}

# Calculate CHM (MNC) from DTS and DTM
chm <- terra::rast(file.path("../../01_DATA/",
                             site,
                             "LiDAR/2-las_utm/chm/chm.tif"))

site_edges_path <- paste(data_dir, 
                         site, 
                         maskname, 
                         sep = "/")
s2_creation_directory <- paste(data_dir, site, 'S2_Images', sep = "/")
dir.create(path = s2_creation_directory, showWarnings = FALSE, recursive = TRUE)
result_path <- paste(results_dir, site, sep = '/')

lai_lidar <- terra::rast(file.path(results_path,
                                   "LiDAR/PAI/lidR/pai.tif"))

resolution <- 10

choice_equal_deciles <- "deciles" #equal_intervals #deciles
# choice_equal_deciles <- c("deciles", "equal_intervals")

s2_creation_directory <- paste(data_dir, 
                               site, 
                               'S2_Images', 
                               "res_10m",
                               sep = "/")

dir.create(path = s2_creation_directory, showWarnings = FALSE, recursive = TRUE)
result_path <- paste(results_dir, site, sep = '/')
S2source <- 'SAFE'
saveRaw <- TRUE

# S-2 Pre-Processing: Cloud Mask, Reflectance
results <- preprocess_S2(dateAcq, 
                         site_edges_path,
                         s2_creation_directory,
                         result_path, 
                         resolution = resolution, 
                         S2source = 'SAFE',
                         saveRaw = TRUE)

cloud_path <- results$Cloud_File
refl_path <- results$refl_path # Reflectance

reflectance <- terra::rast(refl_path)
reflectance <- subset(reflectance, 1)
reflectance[] <- 1

# Mask

site_edges_mask <- mask_site_edges(reflectance,
                                   site_edges_path,
                                   resolution,
                                   masks_dir
)

# Mask

cloud_mask <- mask_clouds(reflectance,
                          cloud_path,
                          site_edges_mask,
                          resolution,
                          masks_dir
)

# MNC

mnc <- create_mnc(chm, 
                  masks_dir
)

# Threshold

mnc_thresholded <- mask_mnc_threshold(mnc,
                                      threshold = 2,
                                      resolution,
                                      masks_dir
                                      
)

# Project

average_mnc_thresholded <- project_mnc_threshold(mnc_thresholded,
                                                 reflectance,
                                                 resolution,
                                                 masks_dir
)

# Majority

# percs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
percentage <- 0.9
average_mnc_thresholded_perc <- mask_majority_project_mnc_threshold(mnc,
                                                                    average_mnc_thresholded,
                                                                    percentage,
                                                                    resolution,
                                                                    masks_dir
)

# Final Mask

artifacts <- cloud_mask * site_edges_mask
save_envi_file(artifacts, 
               sprintf("artifacts_res_%s_m_viz",
                       resolution), 
               masks_dir)
save_basic_plot(plot_to_save = artifacts,
                dirname = masks_dir,
                filename = sprintf("artifacts_res_%s_m_viz.png",
                                   resolution),
                title = sprintf("artifacts_res_%s_m_viz",
                                resolution))

final_mask <- cloud_mask * site_edges_mask * average_mnc_thresholded_perc
save_envi_file(final_mask, 
               sprintf("artifacts_low_vegetation_majority_%s_p_res_%s_m_viz", 
                       percentage*100,
                       resolution), 
               masks_dir)
save_basic_plot(plot_to_save = final_mask,
                dirname = masks_dir,
                filename = sprintf("artifacts_low_vegetation_majority_%s_p_res_%s_m_viz.png", 
                                   percentage*100,
                                   resolution),
                title = sprintf("artifacts_low_vegetation_majority_%s_p_res_%s_m_viz", 
                                percentage*100,
                                resolution))
final_mask[final_mask == 0] <- NA
plot(final_mask)
save_envi_file(final_mask, 
               sprintf("artifacts_low_vegetation_majority_%s_p_res_%s_m", 
                       percentage*100,
                       resolution), 
               masks_dir)
save_basic_plot(plot_to_save = final_mask,
                dirname = masks_dir,
                filename = sprintf("artifacts_low_vegetation_majority_%s_p_res_%s_m.png", 
                                   percentage*100,
                                   resolution),
                title = sprintf("artifacts_low_vegetation_majority_%s_p_res_%s_m", 
                                percentage*100,
                                resolution))

# LiDAR
lai_lidar <- terra::project(lai_lidar,
                            final_mask)
lai_lidar_masked <- lai_lidar * final_mask
plot(lai_lidar_masked, main="lai lidar masked")
save_envi_file(lai_lidar_masked, 
               sprintf("lai_lidar_masked_res_%s_m",
                       resolution), 
               masks_dir)
save_basic_plot(plot_to_save = lai_lidar_masked,
                dirname = masks_dir,
                filename = sprintf("lai_lidar_masked_res_%s_m.png",
                                   resolution),
                title = sprintf("lai_lidar_masked_res_%s_m",
                                resolution))

# S-2
# 10m
if (site == "Mormal"){
  lai_s2 <- open_raster_file("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/res_10m/PRO4SAIL_INVERSION_atbd/addmult/bands_3_4_8/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi")
}
if (site == "Blois"){
  lai_s2 <- open_raster_file("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/L2A_T31TCN_A031222_20210614T105443/res_10m/PRO4SAIL_INVERSION_modifbandatbd/addmult/bands_3_4_8/L2A_T31TCN_A031222_20210614T105443_Refl_lai.envi")
}
if (site == "Aigoual"){
  lai_s2 <- open_raster_file("/home/corroyez/Documents/NC_Full/03_RESULTS/Aigoual/L2A_T31TEJ_A031608_20210711T104217/res_10m/PRO4SAIL_INVERSION_atbd/addmult/bands_3_4_8/L2A_T31TEJ_A031608_20210711T104217_Refl_lai.envi")
}
lai_s2 <- terra::project(lai_s2,
                         final_mask)
lai_s2_masked <- lai_s2 * final_mask
plot(lai_lidar_masked, main="lai s2 masked")
save_envi_file(lai_s2_masked, 
               sprintf("lai_s2_masked_res_%s_m",
                       resolution), 
               masks_dir)
save_basic_plot(plot_to_save = lai_s2_masked,
                dirname = masks_dir,
                filename = sprintf("lai_s2_masked_res_%s_m.png",
                                   resolution),
                title = sprintf("lai_s2_masked_res_%s_m",
                                resolution))

# Heterogeneity Mask

# MNC
# mnc_masked <- mask(mnc, terra::project(final_mask, mnc))
mnc_proj <- terra::project(mnc, final_mask)
mnc_masked_proj <- mask(mnc_proj, final_mask)
save_envi_file(mnc_masked_proj,
               sprintf("mnc_masked_proj_res_%s_m", resolution),
               masks_dir)

# Max
max <- terra::project(mnc, final_mask, method = 'max')
max <- mask(max, final_mask)
save_envi_file(max,
               sprintf("mnc_max_res_%s_m", resolution),
               masks_dir)

# Standard Deviation
standard_deviation <- aggregate(mnc,
                                fact=resolution,
                                fun="std")
standard_deviation_proj <- terra::project(standard_deviation, final_mask,
                                          method = 'average')

standard_deviation_proj <- mask(standard_deviation_proj, final_mask)
save_envi_file(standard_deviation_proj,
               sprintf("mnc_std_res_%s_m", resolution),
               masks_dir)

# Variance
variance <- standard_deviation_proj * standard_deviation_proj
save_envi_file(variance,
               sprintf("mnc_variance_res_%s_m", resolution),
               masks_dir)

# Coefficient of Variation
mean_heights <- terra::project(mnc, final_mask, method = 'average')
mean_heights <- mask(mean_heights, final_mask)
coeff_variation <- (standard_deviation_proj / mean_heights) * 100

# coeff_variation <- terra::project(coeff_variation, final_mask)

save_envi_file(coeff_variation,
               sprintf("mnc_coeff_variation_res_%s_m", resolution),
               masks_dir)

# Mean Heights
save_envi_file(mean_heights,
               sprintf("mnc_mean_heights_res_%s_m", resolution),
               masks_dir)

# Rumple index
# rumple <- terra::rast(file.path(masks_dir, "rumple.tif"))
# rumple <- terra::project(rumple, final_mask)
# rumple <- mask(rumple, final_mask)

# Plot the standard deviation as a function of mean heights
plot(mean_heights, standard_deviation_proj, main = "Standard Deviation vs Mean Heights")

# Plot the coefficient of variation as a function of mean heights
plot(mean_heights, coeff_variation, main = "Coefficient of Variation vs Mean Heights")

# Plot the coefficient of variation as a function of the standard deviation
plot(standard_deviation_proj, coeff_variation, main = "Coefficient of Variation vs Standard Deviation")

save_x_y_plot(xvar = mean_heights,
              yvar = standard_deviation_proj,
              dirname = masks_dir,
              filename = "std_vs_mean_heights.png",
              xlab = "Mean Heights",
              ylab = "Standard Deviation",
              title = "Standard Deviation vs Mean Heights")
save_x_y_plot(xvar = mean_heights,
              yvar = coeff_variation,
              dirname = masks_dir,
              filename = "cv_vs_mean_heights.png",
              xlab = "Mean Heights",
              ylab = "Coefficient of Variation",
              title = "Coefficient of Variation vs Mean Heights")
save_x_y_plot(xvar = standard_deviation_proj,
              yvar = coeff_variation,
              dirname = masks_dir,
              filename = "cv_vs_std.png",
              xlab = "Standard Deviation",
              ylab = "Coefficient of Variation",
              title = "Coefficient of Variation vs Standard Deviation")
# save_x_y_plot(xvar = standard_deviation_proj,
#               yvar = rumple,
#               dirname = masks_dir,
#               filename = "rumple_vs_std.png",
#               xlab = "Standard Deviation",
#               ylab = "Rumple Index",
#               title = "Rumple Index vs Standard Deviation")
# save_x_y_plot(xvar = mean_heights,
#               yvar = rumple,
#               dirname = masks_dir,
#               filename = "rumple_vs_std.png",
#               xlab = "Standard Deviation",
#               ylab = "Rumple Index",
#               title = "Rumple Index vs Standard Deviation")

# Calculate deciles on heterogeneity raster
errors <- list(standard_deviation_proj, coeff_variation, mean_heights)
# errors <- list(mean_heights)
for (choice in choice_equal_deciles) {
  for (error in errors){
    if (identical(error, standard_deviation_proj)) {
      cat("Error: Standard Deviation\n")
      err <- "std"
    }
    else if (identical(error, coeff_variation)) {
      cat("Error: Coefficient of Variation\n")
      err <- "cv"
    }
    else if (identical(error, mean_heights)) {
      cat("Error: Mean Heights\n")
      err <- "mean_heights"
    }
    else if (identical(error, rumple)) {
      cat("Error: Rumple\n")
      err <- "rumple"
    }
    else if (identical(error, max)) {
      cat("Error: Max\n")
      err <- "max"
    }
    else{
      stop("Error is not Standard Deviation or Coefficient of Variation
           or Mean Heights")
    }
    
    if (choice == "equal_intervals"){
      class_filename <- "Equal_Intervals"
      class_dir <- file.path(masks_dir, class_filename)
      dir.create(path = class_dir, showWarnings = F, recursive = T)
      choice_equal_deciles_for_plot <- "Interval"
      
      # Calculate the range of error values
      error_min <- min(values(error), na.rm = TRUE)
      error_max <- max(values(error), na.rm = TRUE)
      
      # Determine the number of intervals
      num_intervals <- 5
      inc <- 0.2
      range <- c(0, (num_intervals - 1) * inc)
      midpoints <- seq(0.1, 0.9, by = 0.2)
      
      # Adjust dir
      nb_filename <- sprintf("%s_Intervals/", num_intervals)
      class_dir <- file.path(class_dir, nb_filename)
      dir.create(path = class_dir, showWarnings = F, recursive = T)
      
      # Calculate interval width
      interval_width <- (error_max - error_min) / num_intervals
      
      for (i in 1:num_intervals) {
        low_value <- error_min + interval_width * (i - 1)
        high_value <- error_min + interval_width * i
        
        low_quantile <- inc * (i-1)
        
        cat("Interval", i, "\n")
        cat("Low Value", low_value, "\n")
        cat("High Value", high_value, "\n")
        
        heter_raster <- mask(error,
                             mask = error >= high_value | error < low_value, 
                             maskvalue = 1)
        
        plot(heter_raster, 
             main = paste(sprintf("%s heter raster res %s m", 
                                  err, resolution), 
                          low_quantile, "to", low_quantile + inc))
        
        
        writeRaster(heter_raster, 
                    paste(class_dir, sprintf("%s_heter_raster_res_%s_m_", 
                                             err, resolution),
                          low_quantile, "_",
                          low_quantile + inc, ".tif", sep = ''),
                    overwrite = T)
        
        heter_binary_mask <- error >= high_value | error < low_value
        
        plot(heter_binary_mask,
             main = paste(sprintf("%s heter binary mask res %s m", 
                                  err, resolution), low_quantile,
                          "to", low_quantile + inc))
        
        writeRaster(heter_raster, 
                    paste(class_dir, sprintf("%s_heter_binary_mask_res_%s_m_", 
                                             err, resolution), 
                          low_quantile, "_", low_quantile + inc, 
                          ".tif", sep = ''),
                    overwrite = T)
      }
    }
    else if (choice == "deciles"){
      class_filename <- "Deciles"
      class_dir <- file.path(masks_dir, class_filename)
      dir.create(path = class_dir, showWarnings = F, recursive = T)
      choice_equal_deciles_for_plot <- "Quantile"
      
      # Determine the number of quantiles
      inc <- 1/3 # 1/3 0.1
      range <- c(0, 1-inc)
      midpoints <- seq(range[1] + inc/2, range[2] + inc/2, by = inc)
      quantile_range <- seq(range[1], range[2], by = inc)
      quantiles_list <- list()
      
      # Adjust dir
      nb_filename <- sprintf("%s_Quantiles/", length(quantile_range))
      class_dir <- file.path(class_dir, nb_filename)
      dir.create(path = class_dir, showWarnings = F, recursive = T)
      
      for (low_quantile in quantile_range) {
        quantiles <- quantile(values(error), 
                              probs = c(low_quantile, low_quantile + inc),
                              na.rm = TRUE)
        quantiles_list[[as.character(low_quantile)]] <- quantiles
        print(quantiles)
        print(summary(quantiles))
        
        heter_raster <- mask(error,
                             mask = error >= quantiles[2]
                             |error < quantiles[1], maskvalue = 1)
        
        plot(heter_raster, 
             main = paste(sprintf("%s heter raster res %s m", 
                                  err, resolution), 
                          low_quantile, "to", low_quantile + inc))
        
        heter_binary_mask <- error >= quantiles[2] | error < quantiles[1]
        
        plot(heter_binary_mask,
             main = paste(sprintf("%s heter binary mask res %s m", 
                                  err, resolution), low_quantile,
                          "to", low_quantile + inc))
        
        low_quantile <- round(low_quantile, digits = 2)
        high_quantile <- round(low_quantile + inc, digits = 2)
        
        writeRaster(heter_raster, 
                    paste(class_dir, sprintf("%s_heter_raster_res_%s_m_", 
                                             err, resolution),
                          low_quantile, "_",
                          high_quantile, ".tif", sep = ''),
                    overwrite = T)
        
        writeRaster(heter_raster, 
                    paste(class_dir, sprintf("%s_heter_binary_mask_res_%s_m_", 
                                             err, resolution), 
                          low_quantile, "_", 
                          high_quantile, ".tif", sep = ''),
                    overwrite = T)
      }
    }
    else {
      stop("choice_equal_deciles is not Equal Interval or Decile")
    }
    
    for (low_value in seq(range[1], range[2], by = inc)) {
      # Apply deciles or equal intervals
      low_value <- round(low_value, digits = 2)
      high_value <- round(low_value + inc, digits = 2)
      intervals <- terra::rast(file.path(class_dir, 
                                         paste0(sprintf("%s_heter_binary_mask_res_%s_m_", 
                                                        err, resolution),
                                                low_value, "_", 
                                                high_value, ".tif")))
      # LiDAR
      pai_resample <- terra::project(lai_lidar_masked, intervals)
      pai_masked <- mask(pai_resample, intervals)
      plot(pai_masked, main = paste("pai lidar", low_value, "to", high_value))
      writeRaster(pai_masked, 
                  paste(class_dir, sprintf("%s_pai_masked_res_%s_m_", 
                                           err,
                                           resolution), 
                        low_value, "_", high_value, 
                        ".tif", sep = ''),
                  overwrite = T)
      
      # S-2
      lai_s2_resample <- terra::project(lai_s2_masked, intervals)
      lai_s2_mask <- mask(lai_s2_resample, intervals)
      plot(lai_s2_mask, main = paste("lai s2", low_value, "to", high_value))
      writeRaster(lai_s2_mask, 
                  paste(class_dir, sprintf("%s_lai_s2_masked_res_%s_m_", 
                                           err,
                                           resolution), 
                        low_value, "_", high_value, 
                        ".tif", sep = ''),
                  overwrite = T)
    }
    
    # Test Deciles
    correlation_values <- c()
    for (low_value in seq(range[1], range[2], by = inc)) {
      low_value <- round(low_value, digits = 2)
      high_value <- round(low_value + inc, digits = 2)
      lidar_intervals <- terra::rast(file.path(class_dir, 
                                               paste0(sprintf("%s_pai_masked_res_%s_m_", 
                                                              err, 
                                                              resolution),
                                                      low_value,  "_", 
                                                      high_value, ".tif")))
      lidar_values <- values(lidar_intervals)
      
      s2_intervals <- terra::rast(file.path(class_dir, 
                                            paste0(sprintf("%s_lai_s2_masked_res_%s_m_", 
                                                           err, 
                                                           resolution),
                                                   low_value, "_", 
                                                   high_value, ".tif")))
      s2_values <- values(s2_intervals)
      
      correlation_value <- correlation_test_function(lidar_values, s2_values, method = "pearson")
      correlation_values <- c(correlation_values, correlation_value)
      
      plot_density_scatterplot(var_x = s2_values,
                               var_y = lidar_values,
                               xlab = "LAI PROSAIL",
                               ylab = "LAI LiDAR",
                               dirname = class_dir,
                               filename = paste0(sprintf("%s_scatter_lai_lidar_prosail_res_%s_m_", 
                                                         err, resolution),
                                                 low_value, "_", 
                                                 high_value, ".tif"),
                               xlimits = c(0, 10),
                               ylimits = c(0, 10)
      )
    }
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
  }
}

# mnc_val <- values(mnc)
# mnc_mean_value <- mean(mnc_val, na.rm=T)
# diffs <- mnc - mnc_mean_value
# squared_diffs <- diffs^2
# mean_squared_diffs <- terra::project(squared_diffs, final_mask, method='average')
# standard_deviation_proj <- sqrt(mean_squared_diffs)
