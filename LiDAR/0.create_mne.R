# ---
# title: "0_create_chm.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-15"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("data.table")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

datapath = "../../01_DATA"

# las_dir <- file.path(datapath, "Mormal/LiDAR/2-las_utm")
las_dir <- file.path(datapath, "Mormal/LiDAR/2-las_utm")
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)

# Define the filenames to mask and remove them for las_files variable
filenames_to_mask <- c("LAS_754000_7016000.las", "LAS_754000_7016500.las", 
                       "LAS_754500_7016000.las", "LAS_754500_7016500.las")
if (any(basename(las_files) %in% filenames_to_mask)) {
  las_files <- las_files[!basename(las_files) %in% filenames_to_mask]
}

# Catalog
cat("lasDpath: ", las_dir, "\n")
ctg <- readLAScatalog(las_dir)
plot(ctg)
# lastest <- readLAS(las_files[1])

opt_chunk_size(ctg) <- 0 # Processing by files
opt_chunk_buffer(ctg) <- 10

set_lidr_threads(2)
get_lidr_threads()

# DTM
opt_output_files(ctg) <- paste0(file.path(las_dir, "dtm/res_1m/{XLEFT}_{YBOTTOM}"))
dtm <- rasterize_terrain(ctg, res = 1, algorithm = tin())
plot_dtm3d(dtm, bg = "white") 

# DSM
opt_output_files(ctg) <- paste0(file.path(las_dir, "dsm/res_1m/{XLEFT}_{YBOTTOM}"))
dsm <- rasterize_canopy(ctg, res = 1, p2r())
plot(dsm)

# CHM
chm <- dsm - dtm
plot(chm)

chm_file <- file.path(las_dir, basename('chm'))
writeRaster(chm, filename = file.path(chm_file, basename("chm.tif")), overwrite=TRUE)
