# ----------------------------------------------------------------------------------------------------------------------
# title: "1_Basic_Hybrid_Inversion"
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

# --------------------------------------------- Train SVM Regression ---------------------------------------------------

# get sensor response for Sentinel-2
SensorName <- 'Sentinel_2'
SRF <- GetRadiometry(SensorName)

# define parameters to estimate
Parms2Estimate <- c('lai', 'CHL', 'EWT', 'LMA', 'fCover', 'fAPAR', 'albedo')

# define spectral bands required to train SVR model for each variable
Bands2Select <- list()
for (parm in Parms2Estimate){
  S2BandSelect <- c('B3','B4','B5','B6','B7','B8','B11','B12')
  Bands2Select[[parm]] <- match(S2BandSelect,SRF$Spectral_Bands)
}

# define output directory where LUTs will be saved
PROSAIL_ResPath <- 'HybridInversion'
dir.create(path = PROSAIL_ResPath, showWarnings = FALSE,recursive = TRUE)

GeomAcq <- list()
GeomAcq$min <- GeomAcq$max <- list()
GeomAcq$min$tto <- 0
GeomAcq$max$tto <- 10
GeomAcq$min$tts <- 20
GeomAcq$max$tts <- 30
GeomAcq$min$psi <- 0
GeomAcq$max$psi <- 360

InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                 GeomAcq = GeomAcq)

res <- Generate_LUT_PROSAIL(SAILversion = '4SAIL',
                            InputPROSAIL = InputPROSAIL,
                            SpecPROSPECT = SpecPROSPECT,
                            SpecSOIL = SpecSOIL,
                            SpecATM = SpecATM)
BRF_LUT_1nm <- res$BRF
InputPROSAIL$fCover <- res$fCover
InputPROSAIL$fAPAR <- res$fAPAR
InputPROSAIL$albedo <- res$albedo

BRF_LUT <- applySensorCharacteristics(wvl = SpecPROSPECT$lambda, 
                                      SRF = SRF, 
                                      InRefl = BRF_LUT_1nm)
# identify spectral bands in LUT
rownames(BRF_LUT) <- SRF$Spectral_Bands

BRF_LUT <- applySensorCharacteristics(wvl = SpecPROSPECT$lambda, 
                                      SRF = SRF, 
                                      InRefl = BRF_LUT_1nm)
# identify spectral bands in LUT
rownames(BRF_LUT) <- SRF$Spectral_Bands

BRF_LUT_Noise <- list()
for (parm in Parms2Estimate){
  BRF_LUT_Noise[[parm]] <- apply_noise_atbd(BRF_LUT)
}

modelSVR <- list()
for (parm in Parms2Estimate){
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
xylim$lai <- c(0.0,12)
xylim$CHL <- c(20,90)
xylim$EWT <- c(0.0,0.04)
xylim$LMA <- c(0.003,0.011)
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
  
  ggplot(data = Results, aes(x=estimated, y=measured)) +
    geom_point(alpha=0.6) +
    geom_smooth(method=lm, aes(group = 1)) +
    coord_fixed(ratio = 1,xlim = xylim[[parm]], ylim = xylim[[parm]]) +
    geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.5,color='gray')+
    ggtitle(statsReg) +
    ylab(paste('Measured', parm)) +
    xlab(paste('Estimated', parm))
  filename <- file.path(PROSAIL_ResPath, paste(parm,'_atbd.png',sep = ''))
  ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
         scale = 1, width = 5.5, height = 5, units = "in",
         dpi = 600, limitsize = TRUE)
}