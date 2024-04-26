# ---
# title: "0.normalize_heights.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-15"
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

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("functions_lidar.R")
source("../functions_plots.R")

# Directories 
data <- "../../01_DATA"
site <- "Mormal"
data_site <- file.path(data, site)

# LAS
las <- "LiDAR"
las_utm <- file.path(las,"3-las_normalize_utm/")
las_dir <- file.path(data_site, las_utm)
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)
las_utm_translated <- file.path(las,"4-las_normalize_utm_translated/")
las_dir_translated <- file.path(data_site, las_utm_translated)

# Output
output_dir <- "../../03_RESULTS/Mormal/LiDAR/PAI/lidR/"

# Catalog
cat("lasDpath: ", las_dir, "\n")
ctg <- readLAScatalog(las_dir)
plot(ctg)
# plot(readLAS(las_files[100]))
# Lire le fichier LAS
# las <- readLAS(las_files[100])

# DTM 
dtm_path <- "../../01_DATA/Mormal/LiDAR/3-las_normalize_utm/dtm/res_1m/rasterize_terrain.vrt"
dtm_charged <- terra::rast(dtm_path)
plot(dtm_charged)

# DTS 
dts_path <- "../../01_DATA/Mormal/LiDAR/3-las_normalize_utm/dts/res_1m/rasterize_canopy.vrt"
dts_charged <- terra::rast(dts_path)
plot(dts_charged)

# CHM
chm <- dts_charged - dtm_charged

# LiDR optimization
opt_chunk_size(ctg) <- 0 # Processing by files
opt_chunk_buffer(ctg) <- 0

set_lidr_threads(4)
# get_lidr_threads()

zmax <- max(ctg@data$Max.Z)
opt_output_files(ctg) <- paste0(file.path(las_dir_translated, "{*}_norm.las"))

myNorm <- function(las, zmax) {
  las$Z <- las$Z + zmax - max(las$Z)
  # nlas <- normalize_height(las, tin())
  return(las)
}

laslist <- list()
for (i in 1:length(las_files)) {
  
  las_i <- readLAS(las_files[i])
  las_i <- myNorm(las_i, zmax)
  
  cat(paste0("File: ", i, "\n"))
  flush.console()
  laslist[i] <- las_i
}
