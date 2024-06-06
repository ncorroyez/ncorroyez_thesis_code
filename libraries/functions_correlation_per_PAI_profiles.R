# ---
# title: "functions_lidar"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-27"
# ---

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("sf")
library("dplyr")

mask_forest_composition_area <- function(site_edges, forest_composition, masks_dir,
                                         resolution){
  # Read shapefiles
  site_edges_shp <- st_read(site_edges)
  forest_composition_shp <- st_read(forest_composition)
  
  # Convert forest_composition from Lambert-93 to UTM 31N
  if (st_crs(site_edges_shp) != st_crs(forest_composition_shp)) {
    forest_composition_shp <- st_transform(forest_composition_shp, 
                                           st_crs(site_edges_shp))
  }
  
  # Perform intersection: forest composition on the site only
  merged_shp <- st_intersection(site_edges_shp, forest_composition_shp)
  st_write(merged_shp, file.path(masks_dir,"site_full_composition.shp"), 
           append=FALSE)
  
  # Define deciduous and coniferous terms
  deciduous_terms <- c("chêne", "chênes", 
                       "hêtre", "hêtres", 
                       "feuillu", "feuillus")
  coniferous_terms <- c("conifère", "conifères")
  
  # Create a logical vector for inclusion and exclusion
  deciduous_logical <- grepl(paste(deciduous_terms, collapse = "|"), 
                             merged_shp$TFV, 
                             ignore.case = TRUE)
  coniferous_logical <- grepl(paste(coniferous_terms, collapse = "|"), 
                              merged_shp$TFV, 
                              ignore.case = TRUE)
  
  # Filter to create deciduous-only and deciduous+mixed coniferous mask
  deciduous_only_shp <- merged_shp[deciduous_logical & !coniferous_logical, ]
  deciduous_flex_shp <- merged_shp[deciduous_logical, ]
  
  st_write(deciduous_only_shp, file.path(masks_dir,"site_deciduous_only.shp"), 
           append=FALSE)
  st_write(deciduous_flex_shp, file.path(masks_dir,"site_deciduous_flex.shp"), 
           append=FALSE)
  
  # Rasterize full composition, deciduous-flex and deciduous-only mask
  raster_mask <- raster(extent(deciduous_only_shp), resolution = resolution)
  raster_mask[] <- 1  # Set all cells to a default value, e.g., 1
  
  full_composition_mask <- rasterize(merged_shp, raster_mask)
  full_composition_mask <- terra::rast(full_composition_mask)
  full_composition_mask[!is.na(full_composition_mask)] <- 1
  
  deciduous_only_mask <- rasterize(deciduous_only_shp, raster_mask)
  deciduous_only_mask <- terra::rast(deciduous_only_mask)
  deciduous_only_mask[!is.na(deciduous_only_mask)] <- 1
  
  deciduous_flex_mask <- rasterize(deciduous_flex_shp, raster_mask)
  deciduous_flex_mask <- terra::rast(deciduous_flex_mask)
  deciduous_flex_mask[!is.na(deciduous_flex_mask)] <- 1
  
  # Plot and save
  # Full composition
  save_basic_plot(plot_to_save = full_composition_mask,
                  dirname = masks_dir,
                  filename = sprintf("full_composition_mask_res_%s_m.png", 
                                     resolution),
                  title = sprintf("full_composition_mask_res_%s_m", 
                                  resolution))
  plot(full_composition_mask, main=sprintf("full_composition_mask_res_%s_m", 
                                           resolution))
  save_envi_file(full_composition_mask, 
                 sprintf("full_composition_mask_res_%s_m", 
                         resolution), 
                 masks_dir)
  
  # Deciduous Flex
  save_basic_plot(plot_to_save = deciduous_flex_mask,
                  dirname = masks_dir,
                  filename = sprintf("deciduous_flex_mask_res_%s_m.png", 
                                     resolution),
                  title = sprintf("deciduous_flex_mask_res_%s_m", 
                                  resolution))
  plot(deciduous_flex_mask, main=sprintf("deciduous_flex_mask_res_%s_m", 
                                         resolution))
  save_envi_file(deciduous_flex_mask, sprintf("deciduous_flex_mask_res_%s_m", 
                                              resolution), 
                 masks_dir)
  
  # Deciduous Only
  save_basic_plot(plot_to_save = deciduous_only_mask,
                  dirname = masks_dir,
                  filename = sprintf("deciduous_only_mask_res_%s_m.png", 
                                     resolution),
                  title = sprintf("deciduous_only_mask_res_%s_m", 
                                  resolution))
  plot(deciduous_only_mask, main=sprintf("deciduous_only_mask_res_%s_m", 
                                         resolution))
  save_envi_file(deciduous_only_mask, sprintf("deciduous_only_mask_res_%s_m", 
                                              resolution), 
                 masks_dir)
  
  cat(paste0("Full composition, deciduous only and flex masks", 
             "have been successfully saved at :", 
             masks_dir, "\n"))
  
  return(list(full_composition_mask = full_composition_mask,
              deciduous_only_mask = deciduous_only_mask,
              deciduous_flex_mask = deciduous_flex_mask
  )
  )
}

mask_site_edges <- function(reflectance,
                            site_edges_path, # utm format
                            masks_dir,
                            resolution){
  site_edges_mask <- terra::vect(site_edges_path)
  site_edges_mask <- mask(reflectance, site_edges_mask)
  
  viz <- site_edges_mask
  viz[is.na(viz)] <- 0
  viz[viz == 1] <- NA
  
  save_basic_plot(plot_to_save = site_edges_mask,
                  dirname = masks_dir,
                  filename = sprintf("site_edges_mask_res_%s_m.png", 
                                     resolution),
                  title = sprintf("site_edges_mask_res_%s_m", resolution))
  plot(site_edges_mask, main=sprintf("site_edges_mask_res_%s_m", resolution))
  
  save_envi_file(viz, sprintf("site_edges_mask_res_%s_m_viz", resolution), 
                 masks_dir)
  save_envi_file(site_edges_mask, sprintf("site_edges_mask_res_%s_m", 
                                          resolution), 
                 masks_dir)
  
  cat("Site edges mask has been successfully saved at :", 
      masks_dir, "\n")
  return(site_edges_mask)
}

mask_clouds <- function(reflectance,
                        cloud_path,
                        site_edges,
                        masks_dir,
                        resolution){
  cloud_mask <- terra::rast(cloud_path)
  
  
  cloud_mask <- terra::project(cloud_mask, 
                               site_edges)
  
  cloud_mask <- reflectance * cloud_mask
  
  save_basic_plot(plot_to_save = cloud_mask,
                  dirname = masks_dir,
                  filename = sprintf("cloud_mask_res_%s_m.png", resolution),
                  title = sprintf("cloud_mask_res_%s_m", resolution))
  plot(cloud_mask, main=sprintf("cloud_mask_res_%s_m", resolution))
  
  save_envi_file(cloud_mask, sprintf("cloud_mask_res_%s_m", resolution), 
                 masks_dir)
  
  cat("Clouds mask has been successfully saved at :", 
      masks_dir, "\n")
  return(cloud_mask)
}

save_chm <- function(chm,
                     masks_dir){
  
  # chm <- dts_charged - dtm_charged
  save_basic_plot(plot_to_save = chm,
                  dirname = masks_dir,
                  filename = "chm.png",
                  title = "chm")
  save_envi_file(chm, "chm", masks_dir)
  cat("CHM has been successfully saved at :", 
      masks_dir, "\n")
  return(chm)
}

mask_chm_thresh <- function(chm,
                            threshold,
                            masks_dir,
                            resolution){
  chm_thresh <- chm > threshold
  
  save_basic_plot(plot_to_save = chm_thresh,
                  dirname = masks_dir,
                  filename = sprintf("chm_thresholded_%d_m_res_%s_m.png", 
                                     threshold, 
                                     resolution),
                  title = sprintf("chm_thresholded_%d_m_res_%s_m", 
                                  threshold, 
                                  resolution))
  plot(chm_thresh, 
       main=sprintf("chm_thresholded_%d_m_res_%s_m", 
                    threshold, 
                    resolution))
  
  save_envi_file(chm_thresh, 
                 sprintf("chm_thresholded_%d_m_res_%s_m", 
                         threshold, 
                         resolution), 
                 masks_dir)
  cat(sprintf("CHM thresholded %d m has been successfully saved at :", 
              threshold), 
      masks_dir, "\n")
  return(chm_thresh)
}

project_chm_thresh <- function(chm_thresh,
                               raster_with_right_coordinates,
                               masks_dir,
                               resolution){
  average_chm_thresh <- terra::project(chm_thresh, 
                                            raster_with_right_coordinates, 
                                            method = 'average')
  save_basic_plot(plot_to_save = average_chm_thresh,
                  dirname = masks_dir,
                  filename = sprintf(
                    "chm_1_m_thresholded_average_to_res_%s_m.png", resolution),
                  title = sprintf("chm_1_m_thresholded_average_to_res_%s_m", 
                                  resolution))
  plot(average_chm_thresh, 
       main=sprintf("chm_1_m_thresholded_average_to_res_%s_m", resolution))
  
  save_envi_file(average_chm_thresh, 
                 sprintf("chm_1_m_thresholded_average_to_res_%s_m", resolution),
                 masks_dir)
  cat("CHM thresholded projected has been successfully saved at :", 
      masks_dir, "\n")
  return(average_chm_thresh)
}

mask_majority_project_chm_thresh <- function(chm,
                                             average_chm_thresh,
                                             percentage,
                                             masks_dir,
                                             resolution){
  
  # Majority
  average_chm_thresh_percentage <- average_chm_thresh > percentage
  save_basic_plot(plot_to_save = average_chm_thresh_percentage,
                  dirname = masks_dir,
                  filename = sprintf(
                    "average_chm_thresholded_%s_p_kept_res_%s_m.png", 
                    percentage*100, resolution),
                  title = sprintf("average_chm_thresholded_%s_p_kept_res_%s_m", 
                                  percentage*100, resolution))
  plot(average_chm_thresh_percentage, 
       main = sprintf("average_chm_thresholded_%s_p_kept_res_%s_m", 
                      percentage*100, resolution))
  
  save_envi_file(average_chm_thresh_percentage, 
                 sprintf("average_chm_thresholded_%s_p_kept_res_%s_m", 
                         percentage*100, resolution), 
                 masks_dir)
  
  
  # Project CHM
  average_chm <- terra::project(chm, 
                                average_chm_thresh, 
                                method = 'average')
  chm_val <- values(average_chm)
  chm_val[chm_val < 0] <- 0
  
  index <- average_chm_thresh_percentage == 1
  chm_final_mask <- average_chm
  chm_final_mask[!index] <- NA
  
  save_basic_plot(plot_to_save = chm_final_mask,
                  dirname = masks_dir,
                  filename = sprintf(
                    "chm_masked_with_average_chm_thresholded_%s_p_kept_res_%s_m.png", 
                    percentage*100, resolution),
                  title = sprintf(
                    "chm_masked_with_average_chm_thresholded_%s_p_kept_res_%s_m", 
                    percentage*100, resolution))
  plot(chm_final_mask,
       main = sprintf(
         "chm_masked_with_average_chm_thresholded_%s_p_kept_res_%s_m", 
         percentage*100, resolution))
  
  save_envi_file(chm_final_mask, 
                 sprintf(
                   "chm_masked_with_average_chm_thresholded_%s_p_kept_res_%s_m", 
                   percentage*100, resolution), 
                 masks_dir)
  
  # Comparison to assess the mask
  chm_final_mask_val <- values(chm_final_mask)
  chm_final_mask_val[chm_final_mask_val < 0] <- 0
  
  chm_list <- list(chm_val, chm_final_mask_val)
  labs <- c("chm", "chm_final_mask")
  plot_histogram(vars_list = chm_list, 
                 title = sprintf(
                   "CHM vs CHM masked with threshold majorated %s p res %s m", 
                   percentage*100, resolution),
                 xlab = "heights",
                 var_labs = labs,
                 dirname = masks_dir,
                 filename = sprintf(
                   "hist_chm_vs_chm_masked_threshold_majorated_%s_p_kept_res_%s_m", 
                   percentage*100, resolution))
  return(average_chm_thresh_percentage)
}

create_vegetation_forest_mask <- function(data_dir,
                                          site,
                                          shapefiles_dir,
                                          masks_dir,
                                          results_path,
                                          resolution){
  
  utm_init <- file.path(shapefiles_dir, "utm_init.shp")
  bdforet_2 <- file.path(shapefiles_dir, "bdforet_2.shp")
  
  if (site == "Mormal"){
    dateAcq <- '2021-06-14' # yyyy-mm-dd format mandatory
  } else if (site == "Blois"){
    dateAcq <- '2021-06-14' # yyyy-mm-dd format mandatory
  } else if (site == "Aigoual"){
    dateAcq <- '2021-07-11' # yyyy-mm-dd format mandatory
  } else{
    stop("Error: Site must be Mormal, Blois, or Aigoual.\n")
  }
  
  result_path <- paste(results_dir, site, sep = '/')
  
  # Calculate CHM from DTS and DTM
  chm <- terra::rast(file.path(results_path,
                               "LiDAR/chm/chm.tif"))
  
  # Keep deciduous-only and deciduous-flex forest
  deciduous_masks <- mask_forest_composition_area(utm_init,
                                                  bdforet_2,
                                                  masks_dir,
                                                  resolution = 10)
  full_composition_mask <- deciduous_masks$full_composition_mask
  deciduous_only_mask <- deciduous_masks$deciduous_only_mask
  deciduous_flex_mask <- deciduous_masks$deciduous_flex_mask
  
  s2_creation_directory <- paste(data_dir, site, 'S2_Images', sep = "/")
  dir.create(path = s2_creation_directory, 
             showWarnings = FALSE, 
             recursive = TRUE)
  
  # S-2 Pre-Processing: Cloud Mask, Reflectance
  results <- preprocess_S2(dateAcq, 
                           utm_init,
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
  
  # Site edges Mask
  site_edges_mask <- mask_site_edges(reflectance,
                                     utm_init,
                                     masks_dir,
                                     resolution = 10)
  
  # Cloud Mask
  
  cloud_mask <- mask_clouds(reflectance,
                            cloud_path,
                            site_edges_mask,
                            masks_dir,
                            resolution = 10)
  
  # CHM
  
  chm <- save_chm(chm, masks_dir)
  
  # Threshold
  
  chm_thresh <- mask_chm_thresh(chm,
                                threshold = 2,
                                masks_dir,
                                resolution = 10)
  
  # Project
  
  average_chm_thresh <- project_chm_thresh(chm_thresh,
                                           reflectance,
                                           masks_dir,
                                           resolution = 10)
  
  # Majority
  
  # percs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  percentage <- 0.9
  average_chm_thresh_perc <- mask_majority_project_chm_thresh(chm,
                                                              average_chm_thresh,
                                                              percentage,
                                                              masks_dir,
                                                              resolution = 10)
  
  # Final Mask
  
  # Initial artifacts
  artifacts <- cloud_mask * site_edges_mask
  save_envi_file(artifacts, 
                 sprintf(paste("artifacts",
                               "res_%s_m_viz",
                               sep = "_"),
                         resolution), 
                 masks_dir)
  save_basic_plot(plot_to_save = artifacts,
                  dirname = masks_dir,
                  filename = sprintf(paste("artifacts",
                                           "res_%s_m_viz.png",
                                           sep = "_"),
                                     resolution),
                  title = sprintf(paste("artifacts",
                                        "res_%s_m_viz",
                                        sep = "_"),
                                  resolution))
  
  # Assign unique names to each SpatRaster object
  names(full_composition_mask) <- "full_composition"
  names(deciduous_only_mask) <- "deciduous_only"
  names(deciduous_flex_mask) <- "deciduous_flex"
  
  # Composition list
  compositions <- list(full_composition_mask,
                       deciduous_only_mask, 
                       deciduous_flex_mask)
  
  for (composition in compositions) {
    if (identical(values(composition), values(full_composition_mask))
        && names(composition) == "full_composition") {
      type <- "full_composition"
    } else if (identical(values(composition), values(deciduous_only_mask))
               && names(composition) == "deciduous_only") {
      type <- "deciduous_only"
    } else if (identical(values(composition), values(deciduous_flex_mask))
               && names(composition) == "deciduous_flex") {
      type <- "deciduous_flex"
    } else {
      stop("Error: Composition choice is not amongst right ones.\n")
    }
    
    # Print the type determined for the current composition
    print(paste("Composition type:", type))
    composition_resampled <- terra::resample(composition, 
                                             site_edges_mask, 
                                             method="near")
    
    # Artifacts + Composition
    artifacts_composition <- cloud_mask * site_edges_mask * composition_resampled
    save_envi_file(artifacts_composition, 
                   sprintf(paste("artifacts",
                                 type,
                                 "res_%s_m_viz",
                                 sep = "_"),
                           resolution), 
                   masks_dir)
    save_basic_plot(plot_to_save = artifacts_composition,
                    dirname = masks_dir,
                    filename = sprintf(paste("artifacts",
                                             type,
                                             "res_%s_m_viz.png",
                                             sep = "_"),
                                       resolution),
                    title = sprintf(paste("artifacts",
                                          type,
                                          "res_%s_m_viz",
                                          sep = "_"),
                                    resolution))
    
    final_mask <- (cloud_mask * site_edges_mask 
                   * composition_resampled * average_chm_thresh_perc)
    save_envi_file(final_mask, 
                   sprintf(paste("artifacts",
                                 type,
                                 "low_vegetation_majority_%s_p_res_%s_m_viz",
                                 sep = "_"),
                           percentage*100,
                           resolution), 
                   masks_dir)
    save_basic_plot(plot_to_save = final_mask,
                    dirname = masks_dir,
                    filename = sprintf(paste("artifacts",
                                             type,
                                             "low_vegetation_majority_%s_p_res_%s_m_viz.png",
                                             sep = "_"),
                                       percentage*100,
                                       resolution),
                    title = sprintf(paste("artifacts",
                                          type,
                                          "low_vegetation_majority_%s_p_res_%s_m_viz",
                                          sep = "_"),
                                    percentage*100,
                                    resolution))
    final_mask[final_mask == 0] <- NA
    plot(final_mask)
    save_envi_file(final_mask, 
                   sprintf(paste("artifacts",
                                 type,
                                 "low_vegetation_majority_%s_p_res_%s_m",
                                 sep = "_"),
                           percentage*100,
                           resolution), 
                   masks_dir)
    save_basic_plot(plot_to_save = final_mask,
                    dirname = masks_dir,
                    filename = sprintf(paste("artifacts",
                                             type,
                                             "low_vegetation_majority_%s_p_res_%s_m.png",
                                             sep = "_"),
                                       percentage*100,
                                       resolution),
                    title = sprintf(paste("artifacts",
                                          type,
                                          "low_vegetation_majority_%s_p_res_%s_m",
                                          sep = "_"),
                                    percentage*100,
                                    resolution))
  }
}

create_heterogeneity_quantiles <- function(error,
                                           err,
                                           site,
                                           masks_dir,
                                           composition_mask_dir,
                                           choice,
                                           resolution = 10){
  
  # Equal Intervals
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
  # Deciles
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
  
  # Apply quantiles to LiDAR and Sentinel-2 LAI
  lai_lidar <- terra::rast(file.path(composition_mask_dir, 
                                     "lai_lidar_res_10_m.tif"))
  lai_s2 <- terra::rast(file.path(composition_mask_dir, 
                                  "lai_s2_res_10_m.tif"))
  
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
    pai_resample <- terra::project(lai_lidar, intervals)
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
    lai_s2_resample <- terra::project(lai_s2, intervals)
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
}

apply_and_save_masks <- function(raster,
                                 raster_basename,
                                 masks_dir,
                                 metrics_dir){
  
  if (!dir.exists(file.path(metrics_dir, "Full_Composition"))) {
    dir.create(file.path(metrics_dir, "Full_Composition"), 
               showWarnings = FALSE, 
               recursive = TRUE)
  }
  if (!dir.exists(file.path(metrics_dir, "Deciduous_Flex"))) {
    dir.create(file.path(metrics_dir, "Deciduous_Flex"), 
               showWarnings = FALSE, 
               recursive = TRUE)
  }
  if (!dir.exists(file.path(metrics_dir, "Deciduous_Only"))) {
    dir.create(file.path(metrics_dir, "Deciduous_Only"), 
               showWarnings = FALSE, 
               recursive = TRUE)
  }
  if (!dir.exists(file.path(metrics_dir, "Not_Masked"))) {
    dir.create(file.path(metrics_dir, "Not_Masked"), 
               showWarnings = FALSE, 
               recursive = TRUE)
  }
  # Open masks
  full_comp_mask <- terra::rast(
    file.path(masks_dir, "artifacts_full_composition_low_vegetation_majority_90_p_res_10_m.envi"))
  deciduous_flex_mask <- terra::rast(
    file.path(masks_dir, "artifacts_deciduous_flex_low_vegetation_majority_90_p_res_10_m.envi"))
  deciduous_only_mask <- terra::rast(
    file.path(masks_dir, "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi"))
  
  # Project the raster to right coordinates and resolution
  projected_raster <- terra::project(raster, full_comp_mask)
  
  # Apply and save
  writeRaster(mask(projected_raster, full_comp_mask),
              filename = file.path(metrics_dir, 
                                   "Full_Composition", 
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
  writeRaster(mask(projected_raster, deciduous_flex_mask),
              filename = file.path(metrics_dir, 
                                   "Deciduous_Flex", 
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
  writeRaster(mask(projected_raster, deciduous_only_mask),
              filename = file.path(metrics_dir, 
                                   "Deciduous_Only", 
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
  writeRaster(projected_raster,
              filename = file.path(metrics_dir, 
                                   "Not_Masked", 
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
}

# Metrics 
# VCI
VCI_local <- function(z){
  zmax <- max(z)
  vci <- VCI(z, zmax=zmax)
}

# VDR
VDR <- function(z){
  zmax <- max(z)
  vdr <- (zmax - median(z)) / zmax
}

# Rumple
rumple_index_surface = function(las, res)
{
  las <- filter_surfacepoints(las, 1)
  rumple <- pixel_metrics(las, ~rumple_index(X,Y,Z), res)
  return(rumple)
}

# PAI
myPAI <- function(z, zmin, k=0.5){
  Nout = sum(z<zmin)
  Nin = length(z)
  PAI = -log(Nout/Nin)/k
}

# PAD Layers
myPAD <- function(z, zmin=2, zmax) {
  
  # Translation
  z <- z + zmax - max(z)
  z00 <- zmin + 0.5
  ladd <- LAD(z, z0=min(z)+zmin)
  
  # Fill with 0s
  new_z <- seq(z00, zmax)
  missing_z <- setdiff(new_z, ladd$z)
  missing_data <- data.frame(z = missing_z, lad = rep(0, length(missing_z)))
  ladd <- rbind(ladd, missing_data)
  ladd <- ladd[order(ladd$z), ]
  ladd$lad[is.na(ladd$lad)] <- 0
  row.names(ladd) <- NULL
  
  PAIs <- numeric(length(new_z))
  k <- 1
  for (i in new_z) {
    # Calculate PAI for the specified range of strata
    PAD_sum <- sum(ladd$lad[ladd$z >= i])
    PAIs[k] <- PAD_sum
    k <- k + 1
  }
  PAIs_list <- setNames(as.list(PAIs), new_z)
}
