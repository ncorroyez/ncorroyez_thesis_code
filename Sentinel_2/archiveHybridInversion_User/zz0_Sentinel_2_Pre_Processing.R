# ----------------------------------------------------------------------------------------------------------------------
# title: "0_Sentinel_2_Pre_Processing"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-10-24"
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------- (Optional) Clear the environment and free memory ------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# ----------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------- Libraries --------------------------------------------------------

library("lidR")
library("data.table")
library("raster")
library("sf")
library("devtools")
library("sen2r")
library("preprocS2")
library("sf")

# ----------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

dateAcq <- '2021-06-14' # define date of S2 acquisition

path_vector <- '../../01_DATA/Mormal/Shape/foret_mormal_shape.shp'

# define output directory where SAFE zipfile is stored
DirWrite <- '../../01_DATA/S2_Images'
if (!dir.exists(DirWrite)) {
  dir.create(DirWrite, showWarnings = FALSE, recursive = TRUE)
}

Path_S2 <- get_S2_L2A_Image(l2a_path = DirWrite, 
                            spatial_extent = path_vector, 
                            dateAcq = dateAcq,
                            DeleteL1C = TRUE, 
                            Sen2Cor = TRUE,
                            GoogleCloud = TRUE)

##____________________________________________________________________##
##        Define where data is stored and where to write results      ##
##--------------------------------------------------------------------##
# Result directory
result_path <- '../../03_RESULTS'
dir.create(path = result_path,showWarnings = FALSE,recursive = TRUE)

##____________________________________________________________________##
##                  Extract, resample & stack data                    ##
##--------------------------------------------------------------------##
# define resolution
resolution <- 10
# define source of data
S2source <- 'SAFE'
S2obj <- preprocS2::extract_from_S2_L2A(Path_dir_S2 = Path_S2,
                                        path_vector = path_vector,
                                        S2source = S2source,
                                        resolution = resolution)

# update shapefile if needed (reprojection)
path_vector <- S2obj$path_vector

# create specific result directory corresponding to granule name
results_site_path <- file.path(result_path,basename(S2obj$S2_Bands$GRANULE))
dir.create(path = results_site_path,showWarnings = FALSE,recursive = TRUE)
##____________________________________________________________________##
##                        Write CLOUD MASK                            ##
##--------------------------------------------------------------------##
# directory for cloud mask
Cloud_path <- file.path(results_site_path,'CloudMask')
dir.create(path = Cloud_path,showWarnings = FALSE,recursive = TRUE)
# Filename for cloud mask
cloudmasks <- preprocS2::save_cloud_s2(S2_stars = S2obj$S2_Stack,
                                       Cloud_path = Cloud_path,
                                       S2source = S2source, SaveRaw = T)
##____________________________________________________________________##
##                        Write REFLECTANCE                           ##
##--------------------------------------------------------------------##
# directory for Reflectance
Refl_dir <- file.path(results_site_path,'Reflectance')
dir.create(path = Refl_dir,showWarnings = FALSE,recursive = TRUE)
# filename for Reflectance
Refl_path <- file.path(Refl_dir,paste(basename(S2obj$S2_Bands$GRANULE),'_Refl',sep = ''))

# Save Reflectance file as ENVI image with BIL interleaves
# metadata files are important to account for offset applied on S2 L2A products 
tile_S2 <- get_tile(S2obj$S2_Bands$GRANULE)
dateAcq_S2 <- get_date(S2obj$S2_Bands$GRANULE)
preprocS2::save_reflectance_s2(S2_stars = S2obj$S2_Stack, 
                               Refl_path = Refl_path,
                               S2Sat = NULL, 
                               tile_S2 = tile_S2, 
                               dateAcq_S2 = dateAcq_S2,
                               Format = 'ENVI', 
                               datatype = 'Int16', 
                               MTD = S2obj$S2_Bands$metadata, 
                               MTD_MSI = S2obj$S2_Bands$metadata_MSI)