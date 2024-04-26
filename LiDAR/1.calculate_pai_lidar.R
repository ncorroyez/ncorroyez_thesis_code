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
source("functions_lidar.R")
source("../functions_plots.R")

# Directories 
data <- "../../01_DATA"
site <- "Mormal"
data_site <- file.path(data, site)

# LAS
# las_utm <- "LiDAR/3-las_normalize_utm/"
las_utm <- "LiDAR/2-las_utm/"
las_dir <- file.path(data_site, las_utm)
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)

# Output
output_dir <- "../../03_RESULTS/Mormal/LiDAR/PAI/lidR/"

# Catalog
cat("lasDpath: ", las_dir, "\n")
ctg <- readLAScatalog(las_dir)
plot(ctg)

# LiDR optimization
opt_chunk_size(ctg) <- 0 # Processing by files
opt_chunk_buffer(ctg) <- 1
# opt_select(ctg) <- "xyzicr"
set_lidr_threads(1)

# ctg_2 <- readLAScatalog(las_files[72:74])
# opt_chunk_size(ctg_2) <- 0 # Processing by files
# opt_chunk_buffer(ctg_2) <- 1
# set_lidr_threads(1)
# las_files_small <- las_files[72:74]

# Masks
mask_10m <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LiDAR/Heterogeneity_Masks/artifacts_low_vegetation_majority_90_p_res_10_m.envi")
# mask_20m <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LiDAR/Heterogeneity_Masks/artifacts_low_vegetation_majority_90_p_res_20_m.envi")

# Lmoms
lmoms_rasters <- lidR::pixel_metrics(ctg, 
                                     ~as.list(lmom::samlmu(Z, nmom=3, ratios=F)),
                                     10)

lmom_1 <- lmoms_rasters[[1]]
lmom_2 <- lmoms_rasters[[2]]
lmom_3 <- lmoms_rasters[[3]]

lmom_1 <- terra::project(lmom_1, mask_10m)
lmom_1 <- mask(lmom_1, mask_10m)
lmom_2 <- terra::project(lmom_2, mask_10m)
lmom_2 <- mask(lmom_2, mask_10m)
lmom_3 <- terra::project(lmom_3, mask_10m)
lmom_3 <- mask(lmom_3, mask_10m)

lcv <- lmom_2 / lmom_1
lskew <- lmom_3 / lmom_2

writeRaster(lcv, 
            filename = file.path(output_dir, "lcv_res_10_m_non_norm.tif"), 
            overwrite = TRUE)
writeRaster(lskew, 
            filename = file.path(output_dir, "lskew_res_10_m_non_norm.tif"), 
            overwrite = TRUE)

VCI_local <- function(z){
  vci <- VCI(z, zmax=max(z))
}

vci <- lidR::pixel_metrics(ctg, ~VCI_local(Z), res = 10)
vci <- terra::project(vci, mask_10m)
vci <- mask(vci, mask_10m)
writeRaster(vci, 
            filename = file.path(output_dir, "vci_res_10_m_non_norm.tif"), 
            overwrite = TRUE)

# Rumple
rumple_index_surface = function(las, res)
{
  las <- filter_surfacepoints(las, 1)
  rumple <- pixel_metrics(las, ~rumple_index(X,Y,Z), res)
  return(rumple)
}

# opt <- list(raster_alignment = 10)
# rumple <- catalog_map(ctg, rumple_index_surface, res = 10, .options = opt)
# zmax <- max(ctg@data$Max.Z)
# vci <- lidR::pixel_metrics(ctg, ~VCI(Z, zmax), res = 10)
# # plot(rumple)
# rumple <- terra::project(rumple, mask_10m)
# rumple <- mask(rumple, mask_10m)

# rumplee <- lidR::pixel_metrics(ctg_2, ~rumple_index(X,Y,Z), res = 10)
# # plot(rumplee)
# writeRaster(rumple, 
#             filename = file.path(output_dir, "rumple_res_10_m_non_norm.tif"), 
#             overwrite = TRUE)

# VDR
VDR <- function(z){
  zmax <- max(z)
  vdr <- (zmax - median(z)) / zmax
}

vdr <- lidR::pixel_metrics(ctg, ~VDR(Z), res = 10)
vdr <- terra::project(vdr, mask_10m)
vdr <- mask(vdr, mask_10m)
writeRaster(vdr, 
            filename = file.path(output_dir, "vdr_res_10_m_non_norm.tif"), 
            overwrite = TRUE)

# PAI calculation
myPAI <- function(z, zmin, k=0.5){
  Nout = sum(z<zmin)
  Nin = length(z)
  PAI = -log(Nout/Nin)/k
  # print(PAI)
}

# PAD calculation
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


# Test on a small catalog
# n = 2
# ctg_2 = head(ctg, n = n)
ctg_2 <- readLAScatalog(las_files[4])

# On the whole catalog
# pai_ctg <- lidR::pixel_metrics(ctg_2, ~myPAI(Z, zmin=2))
# plot(pai_ctg, col = height.colors(50))



# PAD
z0 <- 2
zmax <- max(ctg@data$Max.Z)
res = 1
pad_rasters <- lidR::pixel_metrics(ctg, ~myPAD(Z, zmin=2, zmax=zmax), res=res)
# pad_rasters <- lidR::pixel_metrics(ctg_2, ~myPADD(Z, 
#                                                zmin = 2, 
#                                                zmax = zmax, 
#                                                input_res = 1, 
#                                                output_res = res), 
#                                    res = res)
# pad_rasters_res <- terra::aggregate(pad_rasters, fact=10) 
pad_rasters_res <- terra::project(pad_rasters, mask_10m, method = 'average')
for (i in 1:nlyr(pad_rasters)) {
  layer <- subset(pad_rasters, i)
  filename <- paste0("PAD_", 
                     names(pad_rasters)[i], 
                     "_", 
                     zmax,
                     # "_1over10_",
                     ".tif")
  writeRaster(layer, 
              filename = file.path(output_dir, filename),
              overwrite = TRUE)
}
# Other metrics
vci <- lidR::pixel_metrics(ctg, ~VCI(Z, zmax))
cv_lad <- lidR::pixel_metrics(ctg, ~cv(LAD(Z)$lad))
rumple <- lidR::pixel_metrics(ctg, ~rumple_index(X,Y,Z))

# Delete inf values
# pai_ctg[pai_ctg == Inf] <- NA
# 
# # Get min and max values 
# palette_min <- min(values(pai_ctg), na.rm = TRUE)
# palette_max <- max(values(pai_ctg), na.rm = TRUE)
# 
# # Plot the graph with the viridis colorbar
# num_colors <- 50
# custom_palette <- viridis(num_colors)
# plot(pai_ctg, col = custom_palette)
# legend("topright", legend = c(palette_min, palette_max), fill = custom_palette)
# writeRaster(pai_ctg, 
#             filename = file.path(output_dir, "pai_las_cleaned.tif"), 
#             overwrite = TRUE)

# Export raster to tif
writeRaster(vci, 
            filename = file.path(output_dir, "vci.tif"), 
            overwrite = TRUE)
writeRaster(cv_lad, 
            filename = file.path(output_dir, "cv_lad.tif"), 
            overwrite = TRUE)
writeRaster(rumple, 
            filename = file.path(output_dir, "rumple.tif"), 
            overwrite = TRUE)
