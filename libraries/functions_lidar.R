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

#' Vertical Complexity Index (VCI)
#'
#' Calculates the Vertical Complexity Index (VCI) based on input data.
#'
#' @param z Numeric vector representing vertical-related data.
#' @return VCI value.
#' @export
VCI_local <- function(z){
  zmax <- max(z)
  vci <- VCI(z, zmax=zmax)
}

#' Vertical Density Ratio (VDR)
#'
#' Calculates the Vertical Density Ratio (VDR) based on input data.
#'
#' @param z Numeric vector representing vertical-related data.
#' @return VDR value.
#' @export
VDR <- function(z){
  zmax <- max(z)
  vdr <- (zmax - median(z)) / zmax
}

#' Rumple Index Surface
#'
#' Calculates the Rumple Index for a given LiDAR point cloud.
#'
#' @param las LiDAR point cloud data.
#' @param res Resolution of the output raster in meters.
#' @return Rumple Index raster.
#' @export
rumple_index_surface = function(las, res)
{
  las <- filter_surfacepoints(las, 1)
  rumple <- pixel_metrics(las, ~rumple_index(X,Y,Z), res)
  return(rumple)
}

#' Plant Area Index (PAI)
#'
#' Calculates the Plant Area Index (PAI) based on input data.
#'
#' @param z Numeric vector representing vegetation-related data.
#' @param zmin Minimum threshold value.
#' @param k Coefficient for PAI calculation.
#' @return PAI value.
#' @export
myPAI <- function(z, zmin, k=0.5){
  Nout = sum(z<zmin)
  Nin = length(z)
  PAI = -log(Nout/Nin)/k
}

#' Plant Area Density (PAD) Layers
#'
#' Calculates the Plant Area Density (PAD) for different layers
#' based on input data.
#'
#' @param z Numeric vector representing vegetation-related data.
#' @param zmin Minimum threshold value.
#' @param zmax Maximum threshold value.
#' @return List of PAD values for each layer.
#' @export
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

# ---
# title: "functions_shadows_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, 
# CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-27"
# ---

# Rotate the matrix 90 degrees clockwise
rotate90_clockwise <- function(mat) {
  t(apply(mat, 2, rev))
}

# Function to compute shadows analysis

#' Perform Shadows Analysis
#'
#' This function performs shadows analysis on a DSM (Digital Surface Model) 
#' raster to determine shaded areas based on sun angle information.
#'
#' @param datapath Path to the data directory.
#' @param site Site identifier.
#'
#' @details This function opens a DSM raster, calculates sun angle based on 
#' site location and time, and computes shadows using ray tracing techniques.
#' The resulting shaded areas are saved as raster files.
#'
#' @return 10-meters shadows raster + save files.
#'
#' @examples
#' perform_shadows_analysis(datapath = "/path/to/data", site = "Mormal")
#' perform_shadows_analysis(datapath = "/path/to/data", site = "Blois")
#' perform_shadows_analysis(datapath = "/path/to/data", site = "Aigoual")
#'
#' @export
perform_shadows_analysis <- function(datapath,
                                     site){
  # Open DSM
  dsm <- terra::rast(file.path(datapath, 
                               site, 
                               "LiDAR/2-las_utm/dsm/rasterize_canopy.vrt"))
  
  # Results storage
  results <- "../../03_RESULTS"
  metrics_dir <- file.path(results,
                           site,
                           "Metrics")
  results_dir <- file.path(metrics_dir,
                           "Raw")
  
  # Open another raster to get the right dimensions
  tmp_raster <- terra::rast(file.path(results_dir,
                                      "lidarlai_res_10_m.tif"))
  
  if (site == "Mormal"){
    t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
    lat <- 50.20
    lon <- 3.74
  } else if (site == "Blois"){
    t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
    lat <- 47.57
    lon <- 1.29
  } else if (site == "Aigoual"){
    t <- as.POSIXct("2021-07-11 10:40:31", tz = "UTC")
    lat <- 44.12
    lon <- 3.52
  } else{
    stop("Error: Site must be Mormal, Blois, or Aigoual.\n")
  }
  
  # Get sun angle and altitude
  df_sun <- oce::sunAngle(t, lat, lon)
  
  dsm_mat <- terra::as.matrix(dsm, wide=TRUE)
  dsm_mat_transpose <- t(dsm_mat)
  
  shade_mat_transpose <- ray_shade(dsm_mat_transpose,
                                   sunaltitude = df_sun$altitude,
                                   sunangle = df_sun$azimuth,
                                   lambert=FALSE)
  shade_mat_transpose <- rotate90_clockwise(shade_mat_transpose)
  
  # Shade mat transpose
  # 1m
  shade_mat_transpose_rast <- terra::rast(shade_mat_transpose, 
                                          ext = terra::ext(dsm),
                                          crs = terra::crs(dsm))
  shade_mat_transpose_rast <- terra::project(shade_mat_transpose_rast, 
                                             dsm, 
                                             method = 'average')
  writeRaster(shade_mat_transpose_rast, 
              file.path(metrics_dir, "shade_res_1_m.tif"),
              overwrite=TRUE)
  
  # 10m
  shade_mat_transpose_rast <- terra::rast(shade_mat_transpose, 
                                          ext = terra::ext(tmp_raster),
                                          crs = terra::crs(tmp_raster))
  shade_mat_transpose_rast <- terra::project(shade_mat_transpose_rast,
                                             tmp_raster,
                                             method = 'average')
  
  return(shade_mat_transpose_rast)
}
