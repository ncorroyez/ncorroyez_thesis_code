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
#' @return This function does not return anything but saves shadow raster files.
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
  results_dir <- file.path(results,
                           site,
                           "Metrics/Raw")
  
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
  shade_mat_transpose_rast <- terra::rast(shade_mat_transpose, ext = terra::ext(dsm),
                                          crs = terra::crs(dsm))
  shade_mat_transpose_rast <- terra::project(shade_mat_transpose_rast, dsm, method = 'average')
  writeRaster(shade_mat_transpose_rast, 
              file.path(results_dir, "shade_res_1_m.tif"),
              overwrite=TRUE)
  
  # 10m
  shade_mat_transpose_rast <- terra::rast(shade_mat_transpose, ext = terra::ext(tmp_raster),
                                          crs = terra::crs(tmp_raster))
  shade_mat_transpose_rast <- terra::project(shade_mat_transpose_rast, tmp_raster, method = 'average')
  writeRaster(shade_mat_transpose_rast, 
              file.path(results_dir, "shade_res_10_m.tif"),
              overwrite=TRUE)
}
