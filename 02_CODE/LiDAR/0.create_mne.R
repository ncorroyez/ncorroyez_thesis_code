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
site <- "Aigoual" # Mormal Blois Aigoual
las_dir <- file.path(datapath, site, "LiDAR/2-las_utm")
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)

# Define the filenames to mask and remove them for las_files variable
if (site == "Mormal"){
  filenames_to_mask <- c("LAS_754000_7016000.las", "LAS_754000_7016500.las", 
                         "LAS_754500_7016000.las", "LAS_754500_7016500.las")
  if (any(basename(las_files) %in% filenames_to_mask)) {
    las_files <- las_files[!basename(las_files) %in% filenames_to_mask]
  } 
}

# Function to format the computation time
format_time <- function(start, end) {
  elapsed <- end - start
  sprintf("Elapsed time: %.2f seconds", elapsed[3])
}

# Function to print raster size
print_raster_size <- function(raster) {
  size <- dim(raster)
  sprintf("Raster size: %dx%d (rows x cols)", size[1], size[2])
}

# Catalog
ctg <- readLAScatalog(las_dir)
plot(ctg)

opt_chunk_size(ctg) <- 0 # Processing by files
opt_chunk_buffer(ctg) <- 10
set_lidr_threads(1)

# DTM
opt_output_files(ctg) <- paste0(file.path(las_dir, "dtm/{XLEFT}_{YBOTTOM}"))
start_time <- proc.time()
dtm <- rasterize_terrain(ctg, res = 1, algorithm = tin()) #knnidw() kriging(k = 40)
end_time <- proc.time()
cat("DTM computation time:\n", format_time(start_time, end_time), "\n")
cat("DTM raster size:\n", print_raster_size(dtm), "\n")
cat("DTM algorithm: kriging\n\n")

# Define a common CRS for consistency
common_crs <- st_crs(ctg)$wkt

# Ensure the CRS of DTM matches the common CRS
crs(dtm) <- common_crs

# DSM
opt_output_files(ctg) <- paste0(file.path(las_dir, "dsm/{XLEFT}_{YBOTTOM}"))
start_time <- proc.time()
dsm <- rasterize_canopy(ctg,
                        res = dtm, # 1,
                        algorithm = p2r()
                        # pitfree(thresholds = c(0, 2, 5, 10, 15),
                        #         max_edge = c(0, 1),
                        #         subcircle = 0.5)
                        ) #p2r()

# CHM
opt_output_files(ctg) <- paste0(file.path(las_dir, "chm/{XLEFT}_{YBOTTOM}"))
chm <- dsm - dtm
end_time <- proc.time()
cat("DSM computation time:\n", format_time(start_time, end_time), "\n")
cat("DSM raster size:\n", print_raster_size(dsm), "\n")
cat("DSM algorithm: pitfree\n\n")

chm_file <- file.path(las_dir, basename('chm'))
dir.create(path = chm_file, showWarnings = FALSE, recursive = TRUE)
writeRaster(chm, filename = file.path(chm_file, basename("chm.tif")), overwrite=TRUE)
