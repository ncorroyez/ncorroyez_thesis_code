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
site <- "Aigoual" # Mormal Blois Aigoual
data_site <- file.path(data, site)

# LAS
las_utm <- "LiDAR/2-las_utm/"
las_utm_dir <- file.path(data_site, las_utm)
las_utm_files <- list.files(las_utm_dir, pattern = "\\.las$", full.names = TRUE)

las_norm_utm <- "LiDAR/3-las_normalized_utm/"
las_norm_utm_dir <- file.path(data_site, las_norm_utm)
las_norm_utm_files <- list.files(las_norm_utm_dir, pattern = "\\.las$", full.names = TRUE)

# Output
output_dir <- file.path("../../03_RESULTS", site, "LiDAR/PAI/lidR")
dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)

# Test
# las_test <- readLAS(las_files[30])

# Catalog
ctg <- readLAScatalog(las_norm_utm_dir)
plot(ctg)

ctg_non_norm <- readLAScatalog(las_utm_dir)
plot(ctg_non_norm)

# LiDR optimization
opt_chunk_size(ctg) <- 0 # Processing by files
opt_chunk_buffer(ctg) <- 1
opt_chunk_size(ctg_non_norm) <- 0 # Processing by files
opt_chunk_buffer(ctg_non_norm) <- 1
# opt_select(ctg) <- "xyzicr"
set_lidr_threads(1)

# ctg_2 <- readLAScatalog(las_files[72:74])
# opt_chunk_size(ctg_2) <- 0 # Processing by files
# opt_chunk_buffer(ctg_2) <- 1
# set_lidr_threads(1)
# las_files_small <- las_files[72:74]

# Masks
mask_10m <- terra::rast(file.path("/home/corroyez/Documents/NC_Full/03_RESULTS",
                                  site,
                                  "LiDAR/Heterogeneity_Masks/artifacts_low_vegetation_majority_90_p_res_10_m.envi"))
# mask_20m <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LiDAR/Heterogeneity_Masks/artifacts_low_vegetation_majority_90_p_res_20_m.envi")

# Lmoms
lmoms_rasters <- lidR::pixel_metrics(ctg_non_norm, 
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
  zmax <- max(z)
  vci <- VCI(z, zmax=zmax)
}
VDR <- function(z){
  zmax <- max(z)
  vdr <- (zmax - median(z)) / zmax
}
rumple_index_surface = function(las, res)
{
  las <- filter_surfacepoints(las, 1)
  rumple <- pixel_metrics(las, ~rumple_index(X,Y,Z), res)
  return(rumple)
}

vci <- lidR::pixel_metrics(ctg_non_norm, ~VCI_local(Z), res = 10)
vci <- terra::project(vci, mask_10m)
vci <- mask(vci, mask_10m)
vdr <- lidR::pixel_metrics(ctg_non_norm, ~VDR(Z), res = 10)
vdr <- terra::project(vdr, mask_10m)
vdr <- mask(vdr, mask_10m)
opt <- list(raster_alignment = 10)
rumple_upgrade <- catalog_map(ctg_non_norm, rumple_index_surface, res = 10, .options = opt)
rumple_upgrade <- terra::project(rumple_upgrade, mask_10m)
rumple_upgrade <- mask(rumple_upgrade, mask_10m)
writeRaster(vci, 
            filename = file.path(output_dir, "vci_res_10_m_non_norm.tif"), 
            overwrite = TRUE)
writeRaster(vdr, 
            filename = file.path(output_dir, "vdr_res_10_m_non_norm.tif"), 
            overwrite = TRUE)
writeRaster(rumple_upgrade,
            filename = file.path(output_dir, "rumple_res_10_m_non_norm.tif"),
            overwrite = TRUE)
# rumple <- lidR::pixel_metrics(ctg, ~rumple_index(X,Y,Z), res = 10)
# rumple <- terra::project(rumple, mask_10m)
# rumple <- mask(rumple, mask_10m)
# writeRaster(rumple,
#             filename = file.path(output_dir, "rumple_res_10_m_non_norm.tif"),
#             overwrite = TRUE)


# rumple <- lidR::pixel_metrics(ctg_2, ~rumple_index(X,Y,Z), res = 10)
# writeRaster(rumple, 
#             filename = file.path(output_dir, "rumple_res_10_m_non_norm.tif"), 
#             overwrite = TRUE)


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
ctg_2 <- readLAScatalog(las_norm_utm_files[4])

# On the whole catalog
# pai_ctg <- lidR::pixel_metrics(ctg_2, ~myPAI(Z, zmin=2))
# plot(pai_ctg, col = height.colors(50))



# PAD
z0 <- 2
zmax <- max(ctg@data$Max.Z[ctg@data$Max.Z < 50])
res = 10
pad_rasters <- lidR::pixel_metrics(ctg, ~myPAD(Z, zmin=2, zmax=zmax), res=res)
pad_rasters_res <- terra::project(pad_rasters, mask_10m, method = 'average')
for (i in 1:nlyr(pad_rasters_res)) {
  layer <- subset(pad_rasters_res, i)
  filename <- paste0("PAD_", 
                     names(pad_rasters_res)[i], 
                     "_", 
                     zmax,
                     ".tif")
  writeRaster(layer, 
              filename = file.path(output_dir, filename),
              overwrite = TRUE)
}
# Other metrics
vci <- lidR::pixel_metrics(ctg, ~VCI(Z, zmax), res=res)
cv_lad <- lidR::pixel_metrics(ctg, ~cv(LAD(Z)$lad), res=res)
rumple <- lidR::pixel_metrics(ctg, ~rumple_index(X,Y,Z), res=res)
pai <- lidR::pixel_metrics(ctg, ~myPAI(Z, zmin=2), res=res)
pai[pai == Inf] <- NA
writeRaster(pai, 
            filename = file.path(output_dir, "pai.tif"), 
            overwrite = TRUE)
writeRaster(cv_lad, 
            filename = file.path(output_dir, "cv_lad.tif"), 
            overwrite = TRUE)









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
