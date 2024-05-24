# ----------------------------------------------------------------------------------------------------------------------
# title: "1_PROSAIL_Forward"
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

# ----------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# define input variables for PROSPECT. 
# refer to prospect tutorial for default values corresponding to undefined PROSPECT variables
CHL <- 40; CAR <- 8; ANT <- 0.0; EWT <- 0.01; LMA <- 0.009; N = 1.5;

# define input variables for SAIL. 
lai <- 5;       # LAI
q <- 0.01;      # Hot spot parameter
TypeLidf <- 2;  LIDFa <- 30;    LIDFb <- NULL;  # leaf inclination distribution function parameters
tts <- 30;      tto <- 10;      psi <- 90;      # geometry of acquisition
rsoil <- SpecSOIL$Dry_Soil                      # soil reflectance (SpecSOIL includes Dry_Soil and Wet_Soil properties)
# run PROSAIL with 4SAIL
Ref_4SAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,
                      CHL = CHL, CAR = CAR, ANT = ANT, EWT = EWT, LMA = LMA, N = N,
                      TypeLidf = TypeLidf,LIDFa = LIDFa,LIDFb = LIDFb,lai = lai,
                      q = q,tts = tts,tto = tto,psi = psi,rsoil = rsoil)

# run PROSAIL with 4SAIL2
fraction_brown <- 0.5
diss <- 0.5
Cv <- 1
Zeta <- 1
# define a couple of leaf chemical constituents corresponding to green and brown leaves
CHL2 <- c(40,5)
CAR2 <- c(8,5)
ANT2 <- c(0,1)
EWT2 <- c(0.01,0.005)
LMA2 <- c(0.009,0.008)
N2 <- c(1.5,2)

# ------------------------------------------ Run PROSAIL in direct mode ------------------------------------------------

Ref_4SAIL2 <- PRO4SAIL(SAILversion = '4SAIL2',Spec_Sensor = SpecPROSPECT,
                       CHL = CHL2, CAR = CAR2, ANT = ANT2, EWT = EWT2, LMA = LMA2, N = N2,
                       TypeLidf = TypeLidf,LIDFa = LIDFa,LIDFb = LIDFb,lai = lai,
                       q = q,tts = tts,tto = tto,psi = psi,rsoil = rsoil,
                       fraction_brown = fraction_brown, diss = diss, Cv = Cv, Zeta = Zeta)

# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------- Compute simplified BRF --------------------------------------------------

# Ref_4SAIL is the variable obtained when running PRO4SAIL as in the previous illustration
# SpecATM corresponds to the direct and diffuse radiation solar spectra
BRF_4SAIL <- Compute_BRF(rdot = Ref_4SAIL$rdot,
                         rsot = Ref_4SAIL$rsot,
                         tts = tts,
                         SpecATM_Sensor = SpecATM)
BRF_4SAIL2 <- Compute_BRF(rdot = Ref_4SAIL2$rdot,
                          rsot = Ref_4SAIL2$rsot,
                          tts = tts,
                          SpecATM_Sensor = SpecATM)

# ----------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------ Compute fAPAR -------------------------------------------------------

# Ref_4SAIL is the variable obtained when running PRO4SAIL as in the previous illustration
# SpecATM corresponds to the direct and diffuse radiation solar spectra
fAPAR_4SAIL <- Compute_fAPAR(abs_dir = Ref_4SAIL$abs_dir,
                             abs_hem = Ref_4SAIL$abs_hem,
                             tts = tts,
                             SpecATM_Sensor = SpecATM)

fAPAR_4SAIL2 <- Compute_fAPAR(abs_dir = Ref_4SAIL2$abs_dir,
                              abs_hem = Ref_4SAIL2$abs_hem,
                              tts = tts,
                              SpecATM_Sensor = SpecATM)
# ----------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------ Compute Albedo ------------------------------------------------------

# Ref_4SAIL is the variable obtained when running PRO4SAIL as in the previous illustration
# SpecATM corresponds to the direct and diffuse radiation solar spectra
albedo_4SAIL <- Compute_albedo(rsdstar = Ref_4SAIL$rsdstar,
                               rddstar = Ref_4SAIL$rddstar,
                               tts = tts,
                               SpecATM_Sensor = SpecATM)

albedo_4SAIL2 <- Compute_albedo(rsdstar = Ref_4SAIL2$rsdstar,
                                rddstar = Ref_4SAIL2$rddstar,
                                tts = tts,
                                SpecATM_Sensor = SpecATM)

# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------- Sensor BRF Simulation ---------------------------------------------------

# get the spectral response for Sentinel-2A
SensorName = 'Sentinel_2A'
# if interested in a different satellite, please use Path_SensorResponse to locate the SRF file expected to be named 'SensorName_Spectral_Response.csv' (separator = tabulations)
SRF <- GetRadiometry(SensorName,Path_SensorResponse = NULL)

# apply sensor characteristics to simulated reflectance
wvl <- SpecPROSPECT$lambda
BRF_4SAIL_S2 <- applySensorCharacteristics(wvl = wvl,
                                           SRF = SRF,
                                           InRefl = BRF_4SAIL)

BRF_4SAIL2_S2 <- applySensorCharacteristics(wvl = wvl,
                                            SRF = SRF,
                                            InRefl = BRF_4SAIL2)

# ----------------------------------------------------------------------------------------------------------------------

