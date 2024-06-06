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
source("../libraries/functions_lidar.R")
source("../libraries/functions_plots.R")

# Directories 
data <- "../../01_DATA"
site <- "Mormal" # Mormal Blois Aigoual
data_site <- file.path(data, site)
lidar <- "LiDAR"

# Directories setup

# UTM
las_utm <- file.path(lidar,"2-las_utm/")
las_dir <- file.path(data_site, las_utm)
dir.create(path = las_dir, showWarnings = FALSE, recursive = TRUE)

# Normalized UTM
las_norm_utm <- file.path(lidar,"3-las_normalized_utm/")
las_norm_dir <- file.path(data_site, las_norm_utm)
dir.create(path = las_norm_dir, showWarnings = FALSE, recursive = TRUE)

# Catalog
cat("lasDpath: ", las_dir, "\n")
ctg <- readLAScatalog(las_dir)
plot(ctg)

# LiDR optimization
opt_chunk_size(ctg) <- 0 # Processing by files
opt_chunk_buffer(ctg) <- 10
set_lidr_threads(1)

# Normalization
opt_output_files(ctg) <- paste0(file.path(las_norm_dir, "{*}_norm"))
normalized_ctg <- normalize_height(ctg, tin())
