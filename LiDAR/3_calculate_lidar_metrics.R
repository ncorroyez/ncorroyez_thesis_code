# ---
# title: "1.calculate.pai_lidar.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-14"
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
library("rgl")
library("lmom")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("../libraries/functions_lidar.R")
source("../libraries/functions_create_masks.R")
source("../libraries/functions_plots.R")

# Directories 
data <- "../../01_DATA"
# sites <- "Aigoual" # Mormal Blois Aigoual
sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Mormal")

for (site in sites){
  data_site <- file.path(data, site)
  results_path <- file.path("../../03_RESULTS", site)
  masks_dir <- file.path(results_path, "LiDAR/Heterogeneity_Masks")
  # LAS directories
  # UTM: for all metrics except LAI
  las_utm <- "LiDAR/2-las_utm/"
  las_utm_dir <- file.path(data_site, las_utm)
  las_utm_files <- list.files(las_utm_dir, pattern = "\\.las$", full.names = TRUE)
  
  # Normalized UTM: for LAI only
  las_norm_utm <- "LiDAR/3-las_normalized_utm/"
  las_norm_utm_dir <- file.path(data_site, las_norm_utm)
  las_norm_utm_files <- list.files(las_norm_utm_dir, pattern = "\\.las$", full.names = TRUE)
  
  # Output
  results_path <- file.path("../../03_RESULTS", site)
  metrics_dir <- file.path(results_path, "Metrics")
  raw_dir <- file.path(metrics_dir, "Raw")
  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Catalogs
  normalized_ctg <- readLAScatalog(las_norm_utm_dir)
  ctg <- readLAScatalog(las_utm_dir)
  # Resolution
  res <- 10
  
  # LiDR optimization
  opt_chunk_size(normalized_ctg) <- 0 # Processing by files
  opt_chunk_buffer(normalized_ctg) <- 5
  opt_chunk_size(ctg) <- 0 # Processing by files
  opt_chunk_buffer(ctg) <- 5
  set_lidr_threads(1)
  
  # Load a raster to project metrics at its resolution and coordinates
  mask_10m <- terra::rast(file.path(results_path,
                                    "LiDAR/Heterogeneity_Masks/artifacts_deciduous_flex_low_vegetation_majority_90_p_res_10_m.envi"))
  
  # Point clouds-related metrics
  # Lmoms (Lcv and Lskew)
  lmoms_rasters <- lidR::pixel_metrics(ctg, 
                                       ~as.list(lmom::samlmu(Z, nmom=3, ratios=F)),
                                       res = res)
  
  lmom_1 <- lmoms_rasters[[1]]
  lmom_2 <- lmoms_rasters[[2]]
  lmom_3 <- lmoms_rasters[[3]]
  lmom_1 <- terra::project(lmom_1, mask_10m)
  lmom_2 <- terra::project(lmom_2, mask_10m)
  lmom_3 <- terra::project(lmom_3, mask_10m)
  
  # Lcv
  lcv <- lmom_2 / lmom_1
  lcv_file <- "lcv_res_10_m.tif"
  writeRaster(lcv, 
              filename = file.path(raw_dir, lcv_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(lcv,
                       lcv_file,
                       masks_dir,
                       metrics_dir)
  # Lskew
  lskew <- lmom_3 / lmom_2
  lskew_file <- "lskew_res_10_m.tif"
  writeRaster(lskew, 
              filename = file.path(raw_dir, lskew_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(lskew,
                       lskew_file,
                       masks_dir,
                       metrics_dir)
  
  # VCI
  vci <- lidR::pixel_metrics(ctg, ~VCI_local(Z), res = res)
  vci <- terra::project(vci, mask_10m)
  vci_file <- "vci_res_10_m.tif"
  writeRaster(vci, 
              filename = file.path(raw_dir, vci_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(vci,
                       vci_file,
                       masks_dir,
                       metrics_dir)
  
  # VDR
  vdr <- lidR::pixel_metrics(ctg, ~VDR(Z), res = res)
  vdr <- terra::project(vdr, mask_10m)
  vdr_file <- "vdr_res_10_m.tif"
  writeRaster(vdr, 
              filename = file.path(raw_dir, vdr_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(vdr,
                       vdr_file,
                       masks_dir,
                       metrics_dir)
  
  # Rumple index
  rumple <- catalog_map(ctg, 
                        rumple_index_surface, 
                        res = res, 
                        .options = list(raster_alignment = res))
  rumple <- terra::project(rumple, mask_10m)
  rumple_file <- "rumple_res_10_m.tif"
  writeRaster(rumple, 
              filename = file.path(raw_dir, rumple_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(rumple,
                       rumple_file,
                       masks_dir,
                       metrics_dir)
  
  # PAI calculation
  pai <- lidR::pixel_metrics(normalized_ctg, ~myPAI(Z, zmin=2), res=res)
  pai[pai == Inf] <- NA
  pai <- terra::project(pai, mask_10m)
  pai_file <- "lidarlai_res_10_m.tif"
  writeRaster(pai, 
              filename = file.path(raw_dir, pai_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(pai,
                       pai_file,
                       masks_dir,
                       metrics_dir)
  
  # PAD calculation
  z0 <- 2
  # Manual modification: remove outliers (>~40/45)
  zmax <- max(normalized_ctg@data$Max.Z[normalized_ctg@data$Max.Z < 45])
  pad_rasters <- lidR::pixel_metrics(normalized_ctg, 
                                     ~myPAD(Z, zmin=2, zmax=zmax), 
                                     res=res)
  pad_rasters <- terra::project(pad_rasters, mask_10m)
  if (!dir.exists(file.path(raw_dir, "PAD_Profiles"))) {
    dir.create(file.path(raw_dir, "PAD_Profiles"),
               showWarnings = FALSE, 
               recursive = TRUE)
  }
  for (i in 1:nlyr(pad_rasters)) {
    layer <- subset(pad_rasters, i)
    filename <- paste0("PAD_", 
                       names(pad_rasters)[i], 
                       "_", 
                       zmax,
                       ".tif")
    writeRaster(layer, 
                filename = file.path(raw_dir, "PAD_Profiles", filename),
                overwrite = TRUE)
    apply_and_save_masks(layer,
                         filename,
                         masks_dir,
                         file.path(metrics_dir, "PAD_Profiles"))
  }
  
  # Shadows
  shade <- perform_shadows_analysis(data, site)
  shade_file <- "shade_res_10_m.tif"
  writeRaster(shade, 
              filename = file.path(raw_dir, shade_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(shade,
                       shade_file,
                       masks_dir,
                       metrics_dir)
  
  # CHM-related metrics
  chm <- terra::rast(file.path(results_path, "LiDAR/chm/chm.tif"))
  
  # Max
  max <- terra::project(chm, mask_10m, method = 'max')
  max_file <- "max_res_10_m.tif"
  writeRaster(max, 
              filename = file.path(raw_dir, max_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(max,
                       max_file,
                       masks_dir,
                       metrics_dir)
  
  
  # Mean
  mean <- terra::project(chm, mask_10m, method = 'average')
  mean_file <- "mean_res_10_m.tif"
  writeRaster(mean, 
              filename = file.path(raw_dir, mean_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(mean,
                       mean_file,
                       masks_dir,
                       metrics_dir)
  
  # Standard deviation (std)
  std <- aggregate(chm, fact=res, fun="std")
  std <- terra::project(std, 
                        mask_10m, 
                        method = 'bilinear')
  std_file <- "std_res_10_m.tif"
  writeRaster(std, 
              filename = file.path(raw_dir, std_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(std,
                       std_file,
                       masks_dir,
                       metrics_dir)
  
  # Coefficient of variation (cv)
  cv <- (std / mean) * 100
  cv_file <- "cv_res_10_m.tif"
  writeRaster(cv, 
              filename = file.path(raw_dir, cv_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(cv,
                       cv_file,
                       masks_dir,
                       metrics_dir)
  
  # Variance
  variance <- std * std
  variance_file <- "variance_res_10_m.tif"
  writeRaster(variance, 
              filename = file.path(raw_dir, variance_file), 
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(variance,
                       variance_file,
                       masks_dir,
                       metrics_dir)
}
