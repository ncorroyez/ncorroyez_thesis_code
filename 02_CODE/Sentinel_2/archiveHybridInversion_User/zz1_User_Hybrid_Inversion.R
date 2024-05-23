# ----------------------------------------------------------------------------------------------------------------------
# title: "1_User_Hybrid_Inversion"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-10-25"
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------- (Optional) Clear the environment and free memory ------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# ----------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------- Libraries --------------------------------------------------------

library("prosail")
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(pracma)

# ----------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# ---------------------------------- Simulate BRF LUT corresponding to InputPROSAIL ------------------------------------

# get sensor response for Sentinel-2
SensorName <- 'Sentinel_2'
SRF <- GetRadiometry(SensorName)


# define parameters to estimate
Parms2Estimate <- c('lai')

# define spectral bands required to train SVR model for each variable
Bands2Select <- list()
for (parm in Parms2Estimate){
  S2BandSelect <- c('B3','B4','B5','B6','B7','B8','B11','B12')
  Bands2Select[[parm]] <- match(S2BandSelect,SRF$Spectral_Bands)
}

# define output directory where LUTs will be saved
PROSAIL_ResPath <- 'HybridInversion_User'
dir.create(path = PROSAIL_ResPath, showWarnings = FALSE,recursive = TRUE)

# define min and max values for PROSAIL input variables
minval <- data.frame('CHL'=10,'CAR'=0,'EWT' = 0.001,'ANT' = 0,'LMA' = 0.001,'N' = 1.0, 'psoil' = 0.0,
                     'BROWN'=0.01, 'LIDFa' = 30, 'lai' = 0.1,'q'=0.1,'tto' = 0,'tts' = 20, 'psi' = 100)
maxval <- data.frame('CHL'=60,'CAR'=15,'EWT' = 0.045,'ANT' = 5,'LMA' = 0.040,'N' = 2.0, 'psoil' = 1.0,
                     'BROWN'=0.5, 'LIDFa' = 70, 'lai' = 5.0,'q'=0.5,'tto' = 5,'tts' = 30, 'psi' = 150)

# Define distribution for PROSAIL input variables: choose between 'Gaussian' and 'Uniform'
TypeDistrib <- data.frame('CHL'='Gaussian', 'CAR'='Uniform', 'EWT' = 'Uniform', 
                          'ANT' = 'Uniform', 'LMA' = 'Uniform', 'N' = 'Uniform',
                          'psoil' = 'Uniform', 'BROWN'='Uniform', 'LIDFa' = 'Uniform', 'lai' = 'Gaussian',
                          'q'='Uniform', 'tto' = 'Uniform', 'tts' = 'Uniform', 'psi' = 'Uniform')

# define mean and STD for gaussian distributions
Mean <- data.frame('CHL'=45,'lai' = 2.5)
Std <- Mean/2
GaussianDistrib <- list('Mean' = Mean, 'Std' = Std)

# define noise level for each variable
NoiseLevel <- list()
NoiseLevel$EWT <- NoiseLevel$CHL <- NoiseLevel$LMA <- NoiseLevel$lai <- 
  NoiseLevel$fCover <- NoiseLevel$fAPAR <- NoiseLevel$albedo <- 0.02

# ----------------------------------------------------------------------------------------------------------------------

# --------------------------------------------- Train SVM Regression ---------------------------------------------------

InputPROSAIL <- get_InputPROSAIL(minval = minval, maxval = maxval, 
                                 TypeDistrib = TypeDistrib, GaussianDistrib = GaussianDistrib)

res <- Generate_LUT_PROSAIL(SAILversion = '4SAIL',
                            InputPROSAIL = InputPROSAIL,
                            SpecPROSPECT = SpecPROSPECT_FullRange,
                            SpecSOIL = SpecSOIL,
                            SpecATM = SpecATM)
BRF_LUT_1nm <- res$BRF
InputPROSAIL$fCover <- res$fCover
InputPROSAIL$fAPAR <- res$fAPAR
InputPROSAIL$albedo <- res$albedo

BRF_LUT <- applySensorCharacteristics(wvl = SpecPROSPECT_FullRange$lambda, 
                                      SRF = SRF, 
                                      InRefl = BRF_LUT_1nm)
# identify spectral bands in LUT
rownames(BRF_LUT) <- SRF$Spectral_Bands

BRF_LUT_Noise <- modelSVR <- list()
for (parm in Parms2Estimate){
  subsetRefl <- BRF_LUT[Bands2Select[[parm]],]
  BRF_LUT_Noise[[parm]] <- subsetRefl + subsetRefl*matrix(rnorm(nrow(subsetRefl)*ncol(subsetRefl),0,NoiseLevel[[parm]]),
                                                          nrow = nrow(subsetRefl))
  modelSVR[[parm]] <- PROSAIL_Hybrid_Train(BRF_LUT = BRF_LUT_Noise[[parm]],
                                           InputVar = InputPROSAIL[[parm]])
}

# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------- Apply Hybrid Inversion --------------------------------------------------

# perform prediction based on models in previous steps
# the prediction returns mean value obtained form the ensemble of regression models for each sample, 
# as well as corresponding standard deviation
MeanEstimate <- StdEstimate <- list()
for (parm in Parms2Estimate){
  HybridRes <- PROSAIL_Hybrid_Apply(RegressionModels = modelSVR[[parm]],
                                    Refl = BRF_LUT_Noise[[parm]])
  MeanEstimate[[parm]] <- HybridRes$MeanEstimate
  StdEstimate[[parm]] <- HybridRes$StdEstimate
}

xylim <- list()
xylim$lai <- c(minval$lai, maxval$lai)
xylim$CHL <- c(minval$CHL, maxval$CHL)
xylim$EWT <- c(minval$EWT, maxval$EWT)
xylim$LMA <- c(minval$LMA, maxval$LMA)
xylim$fCover <- c(0,1)
xylim$fAPAR <- c(0,1)
xylim$albedo <- c(0,0.4)

for (parm in Parms2Estimate){
  # create dataframe and plot results
  Results <- data.frame('measured' = InputPROSAIL[[parm]], 
                        'estimated' = MeanEstimate[[parm]])
  statsReg1 <- cor.test(Results$measured, Results$estimated)$estimate
  statsReg2 <- rmserr(Results$measured, Results$estimated)$rmse
  statsReg <- paste0("r = ", round(statsReg1,2), ", RMSE = ", round(statsReg2,3))
  
  h <- ggplot(data = Results, aes(x=estimated, y=measured)) +
    geom_point(alpha=0.6) +
    geom_smooth(method=lm, aes(group = 1)) +
    coord_fixed(ratio = 1,xlim = xylim[[parm]], ylim = xylim[[parm]]) +
    geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.5,color='gray')+
    ggtitle(statsReg) +
    ylab(paste('Measured', parm)) +
    xlab(paste('Estimated', parm))
  filename <- file.path(PROSAIL_ResPath, paste(parm,'_user.png',sep = ''))
  print(h)
  # ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
  #        scale = 1, width = 5.5, height = 5, units = "in",
  #        dpi = 600, limitsize = TRUE)
}

# ----------------------------------------------------------------------------------------------------------------------

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

########################################################################
##                      COMPUTE SPECTRAL INDEX                        ##
########################################################################
library(prosail)
library(spinR)
library(raster)
library(stars)
# Read raster
Refl <- brick(Refl_path)
# get raster band name and clean format. Expecting band name and wavelength to be documented in image
HDR_Refl <- read_ENVI_header(get_HDR_name(Refl_path))
SensorBands <- HDR_Refl$wavelength
# compute a set of spectral indices defined by IndexList from S2 data
IndexList <- c('NDVI')
# ReflFactor = 10000 when reflectance is coded as INT16
Indices <- spinR::compute_S2SI_Raster(Refl = Refl, SensorBands = SensorBands,
                                      Sel_Indices = IndexList,
                                      ReflFactor = 10000, StackOut=F)

# create directory for Spectral indices
SI_path <- file.path(results_site_path,'SpectralIndices')
dir.create(path = SI_path,showWarnings = FALSE,recursive = TRUE)
# Save spectral indices
for (SpIndx in names(Indices$SpectralIndices)){
  Index_Path <- file.path(SI_path,paste(basename(S2obj$S2_Bands$GRANULE),'_',SpIndx,sep = ''))
  stars::write_stars(st_as_stars(Indices$SpectralIndices[[SpIndx]]), dsn=Index_Path, driver =  "ENVI",type='Float32')
  # write band name in HDR
  HDR <- read_ENVI_header(get_HDR_name(Index_Path))
  HDR$`band names` <- SpIndx
  write_ENVI_header(HDR = HDR,HDRpath = get_HDR_name(Index_Path))
}

# Update Cloud mask based on radiometric filtering
# eliminate pixels with NDVI < NDVI_Thresh because not enough vegetation
NDVI_Thresh <- 0.5
Elim <- which(values(Indices$SpectralIndices[['NDVI']])<NDVI_Thresh)
CloudInit <- stars::read_stars(cloudmasks$BinaryMask)
CloudInit$CloudMask_Binary[Elim] <- 0
# save updated cloud mask
Cloud_File <- file.path(Cloud_path,'CloudMask_Binary_Update')
stars::write_stars(CloudInit, dsn=Cloud_File,driver = "ENVI",type='Byte')


########################################################################
##      COMPUTE BIOPHYSICAL VARIABLES BASED ON PROSAIL INVERSION      ##
########################################################################
# get S2 geometry
# read metadata file from S2 image
xmlfile <- file.path(dirname(Refl_path),'MTD_TL.xml')
S2Geom <- get_S2geometry(MTD_TL_xml = xmlfile)

# Train PROSAIL inversion
minval <- data.frame('CHL'=10,'CAR'=0,'EWT' = 0.005,'ANT' = 0,'LMA' = 0.005,'N' = 1.0,'psoil' = 0.0, 'BROWN'=0.0,
                     'LIDFa' = 30, 'lai' = 0.5,'q'=0.1,'tto' = 0,'tts' = min(S2Geom$SZA), 'psi' = 5)
maxval <- data.frame('CHL'=90,'CAR'=20,'EWT' = 0.04,'ANT' = 3,'LMA' = 0.04,'N' = 2.0, 'psoil' = 1.0, 'BROWN'=0.5,
                     'LIDFa' = 70, 'lai' = 7,'q'=0.25,'tto' = 7,'tts' = max(S2Geom$SZA), 'psi' = 355)

# get sensor response for Sentinel-2
SensorName <- HDR_Refl$`sensor type`
SRF <- GetRadiometry(SensorName,Path_SensorResponse = NULL)
# adjust optical constants from 1nm sampling into spectral S2 spectral sampling
wvl <- SpecPROSPECT$lambda
SpecSensor <- PrepareSensorSimulation(SpecPROSPECT,SpecSOIL,SpecATM,SRF)
SpecPROSPECT_Sensor <- SpecSensor$SpecPROSPECT_Sensor
SpecSOIL_Sensor <- SpecSensor$SpecSOIL_Sensor
SpecATM_Sensor <- SpecSensor$SpecATM_Sensor

# define spectral bands required to train SVR model for each variable
S2BandSelect <- list()
S2BandSelect$CHL <- S2BandSelect$lai <- S2BandSelect$EWT <- S2BandSelect$LMA <- c('B03','B04','B05','B06','B07','B08','B11','B12')
ImgBandNames <- strsplit(HDR_Refl$`band names`,split = ',')[[1]]
# get variable ID for train_prosail_inversion
Bands2Select <- list()
for (bpvar in names(S2BandSelect)){
  Bands2Select[[bpvar]] <- match(S2BandSelect[[bpvar]],ImgBandNames)
}

# define noise level for each variable
NoiseLevel <- list()
NoiseLevel$EWT <- 0.025
NoiseLevel$CHL <- 0.01
NoiseLevel$LMA <- NoiseLevel$lai <- 0.05

# where results will be stored
PROSAIL_ResPath <- file.path(results_site_path,'PRO4SAIL_INVERSION')
dir.create(path = PROSAIL_ResPath,showWarnings = FALSE,recursive = TRUE)

method <- 'liquidSVM' # use 'svmRadial' or 'svmLinear' if you have difficulties with installing liquidSVM
Parms2Estimate <- c('CHL','EWT','LMA','lai')
modelSVR <- train_prosail_inversion(minval = minval, maxval = maxval, TypeDistrib = TypeDistrib, 
                                    GaussianDistrib = GaussianDistrib, Parms2Estimate = Parms2Estimate, 
                                    Bands2Select = Bands2Select, NoiseLevel = NoiseLevel, 
                                    SAILversion = '4SAIL', SpecPROSPECT = SpecPROSPECT_Sensor, 
                                    SpecSOIL = SpecSOIL_Sensor, SpecATM = SpecATM_Sensor,
                                    Path_Results = PROSAIL_ResPath, nbModels = 10, nbSamples = 1000,
                                    FigPlot = FALSE, method = method)

# Apply SVR model on Sentinel-2 data
Apply_prosail_inversion(raster_path = Refl_path, HybridModel = modelSVR, PathOut = PROSAIL_ResPath,
                        SelectedBands = S2BandSelect, bandname = ImgBandNames, MaskRaster = Cloud_File, 
                        MultiplyingFactor = 10000, method = method)