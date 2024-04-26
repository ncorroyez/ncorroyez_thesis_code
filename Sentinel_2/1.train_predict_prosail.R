# ---
# title: "1.train_predict_prosail.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-13"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("functions_sentinel_2.R")
source("../functions_plots.R")

# Pre-processing Parameters
data_dir <- '../../01_DATA'
results_dir <- '../../03_RESULTS'
dateAcq <- '2021-06-14' # yyyy-mm-dd format mandatory
site <- "Mormal"

path_vector <- paste(data_dir, site, "Shape/shp_mormal_onf/foret_domaniale_mormal.shp", sep = "/")

# resolution <- 10 # 10 20
# resolutions <- c(10,20)
resolutions <- c(20) # 10 20

for (resolution in resolutions){
  
  if (resolution == 10){
    s2_creation_directory <- paste(data_dir, 
                                   site, 
                                   'S2_Images', 
                                   "res_10m",
                                   sep = "/")
  } else if (resolution == 20){
    s2_creation_directory <- paste(data_dir, 
                                   site, 
                                   'S2_Images', 
                                   "res_20m",
                                   sep = "/")
  } else{
    stop("Error: Resolution must be 10 or 20.\n")
  }
  dir.create(path = s2_creation_directory, showWarnings = FALSE, recursive = TRUE)
  results_path <- paste(results_dir, site, sep = '/')
  S2source <- 'SAFE'
  saveRaw <- TRUE
  
  # S-2 Pre-Processing: Cloud Mask, Reflectance
  results <- preprocess_S2(dateAcq, 
                           path_vector,
                           s2_creation_directory,
                           results_path, 
                           resolution = resolution, 
                           S2source = 'SAFE',
                           saveRaw = TRUE)
  results_site_path <- results$results_site_path
  refl_path <- results$refl_path # Reflectance
  HDR_Refl <- results$HDR_Refl
  Cloud_File <- results$Cloud_File
  
  # Refl: 10m, 10m no B08, 10m -> 20m, 10m -> 20m no B08
  # 10m
  # refl <- terra::rast(refl_path)
  # 
  # # Modifications
  # res_20m <- "_20m"
  # no_b08 <- "_no_B08"
  # 
  # # 10m no B08
  # layer_index <- which(names(refl) == "B08 (835.1 nanometers)")
  # refl_10m_no_b08 <- refl[[setdiff(1:nlyr(refl), layer_index)]]
  # refl_10m_no_b08_path <- paste0(refl_path, no_b08, ".envi")
  # writeRaster(refl_10m_no_b08, refl_10m_no_b08_path, overwrite=T)
  # 
  # # 20m
  # refl_20m_path <- paste0(refl_path, res_20m, ".envi")
  # writeRaster(terra::aggregate(refl, 
  #                              fact=2, 
  #                              fun='mean'), 
  #             refl_20m_path, overwrite=T)
  # 
  # # 20m no B08
  # refl_20m_no_b08_path <- paste0(refl_path, res_20m, no_b08, ".envi")
  # writeRaster(terra::aggregate(refl[[setdiff(1:nlyr(refl), layer_index)]], 
  #                              fact=2, 
  #                              fun='mean'), 
  #             refl_20m_no_b08_path, overwrite=T)
  # 
  # # List with all reflectances
  # reflectances_list <- list(refl_path, 
  #                           refl_10m_no_b08,
  #                           refl_20m_path,
  #                           refl_20m_no_b08_path)
  
  # Get Sensor Response
  SRF <- get_sensor_response(HDR_Refl)
  
  # Choose Parameters
  SAILversion = '4SAIL'
  nbModels = 10
  nbSamples = 1000
  FigPlot = FALSE
  MultiplyingFactor = 10000
  
  # Define Spectral Bands and Variables to Estimate
  # Parameters
  Parms2Estimate <- c('lai') #lai, CHL, EWT, LMA, fCover, fAPAR, albedo
  
  # Bands
  bands_10m <- c('B3','B4','B8')
  bands_20m <- c('B3','B4','B5','B6','B7','B8A','B11','B12')
  bands_select <- list(bands_10m, bands_20m)
  # bands_select <- list(bands_20m)
  
  # Define Noise Level
  # NULL is AD/MD noise, NoiseLevel is Gaussian
  NoiseLevel <- define_noise_levels()
  noises <- list(NULL_value = NULL, NoiseLevel = NoiseLevel)
  noises <- list(NULL_value = NULL)
  
  # Modify one or more variable distributions
  
  # distribs <- c("atbd", "q_zhang_et_al_2005", "brede_et_al_2020", 
  #               "hauser_et_al_2021", "sinha_et_al_2020",
  #               "verhoef_and_bach_2007", "shiklomanov_et_al_2016")
  # distribs <- c("sinha_et_al_2020",
  # #               "verhoef_and_bach_2007", "shiklomanov_et_al_2016")
  # distribs <- c(
  #   "atbd", "atbd_S2Geom", "atbd_JBGeom",
  #   "atbd_fixed_psi", "atbd_psi_low", "atbd_psi_high",
  #   "atbd_fixed_tto", "atbd_tto_low", "atbd_tto_high",
  #   "atbd_fixed_tts", "atbd_tts_low", "atbd_tts_high",
  #   "q_zhang_et_al_2005", "brede_et_al_2020",
  #   "hauser_et_al_2021", "sinha_et_al_2020",
  #   "verhoef_and_bach_2007", "shiklomanov_et_al_2016",
  #   "atbd_n_high", "atbd_n_low", "atbd_n_full",
  #   "atbd_chl_high", "atbd_chl_low", "atbd_chl_full",
  #   "atbd_brown_high", "atbd_brown_low", "atbd_brown_full",
  #   "atbd_ewt_high", "atbd_ewt_low", "atbd_ewt_full",
  #   "atbd_lma_high", "atbd_lma_low", "atbd_lma_full",
  #   "atbd_lidfa_high", "atbd_lidfa_low", "atbd_lidfa_full",
  #   "atbd_lai_high", "atbd_lai_low", "atbd_lai_full",
  #   "atbd_q_high", "atbd_q_low", "atbd_q_full",
  #   "atbd_psoil_high", "atbd_psoil_low", "atbd_psoil_full"
  # )
  
  distribs <- "atbd"
  
  for (bands in bands_select){
    S2BandSelect <- Bands2Select <- list()
    for (parm in Parms2Estimate){
      S2BandSelect[[parm]] <- bands
      Bands2Select[[parm]] <- match(S2BandSelect[[parm]],SRF$Spectral_Bands)
    }
    
    if (all(bands_20m %in% bands)) {
      res <- "bands_3_4_5_6_7_8A_11_12"
      if (resolution == 10){
        ImgBandNames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')
      } else if (resolution == 20){
        ImgBandNames <- c('B2','B3','B4','B5','B6','B7','B8A','B11','B12')
      } else{
        stop("Error: Resolution must be 10 or 20.\n")
      }
    } 
    else if (all(bands_10m %in% bands)) {
      res <- "bands_3_4_8"
      ImgBandNames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')
    }
    else{
      stop(paste("Error in bands selection. \n"))
    }
    for (noise in noises) {
      if (is.null(noise)) {
        noise_type <- "addmult"
      } 
      else if (is.list(noise)) {
        noise_type <- "mult"
      } 
      else {
        stop("Error: Unknown noise type.\n")
      }
      
      for (distrib in distribs) {
        
        InputPROSAIL <- set_prosail_distribution(distrib, refl_path)
        # Define where Results are Saved: one directory per experiment
        
        # Chosen S-2 Image Resolution
        cat("S-2 Image Resolution", resolution, "\n")
        cat("Reflectance Bands", ImgBandNames, "\n")
        PROSAIL_ResPath <- file.path(results_site_path,
                                     paste0('res_', resolution, 'm'))
        dir.create(path = PROSAIL_ResPath, showWarnings = FALSE, recursive = TRUE)
        
        # Chosen Distribution
        cat("Input Distribution", distrib, "\n")
        PROSAIL_ResPath <- file.path(PROSAIL_ResPath, 
                                     paste0('PRO4SAIL_INVERSION_modifband', distrib))
        dir.create(path = PROSAIL_ResPath, showWarnings = FALSE, recursive = TRUE)
        
        # Chosen Noise
        cat("Noise Type", noise_type, "\n")
        PROSAIL_ResPath <- file.path(PROSAIL_ResPath, noise_type)
        dir.create(path = PROSAIL_ResPath, showWarnings = FALSE, recursive = TRUE)
        
        # Chosen Inversion Bands Resolution
        cat("Inversion Bands Resolution", res, "\n")
        PROSAIL_ResPath <- file.path(PROSAIL_ResPath, res)
        dir.create(path = PROSAIL_ResPath, showWarnings = FALSE, recursive = TRUE)
        
        # Extract and Write Statistics
        summary_df <- extract_summary_stats(InputPROSAIL, PROSAIL_ResPath)
        plot_distributions(InputPROSAIL,
                           save_plots = TRUE,
                           dirname = PROSAIL_ResPath,
                           filename = paste0("hist_", distrib))
        
        # Train
        modelSVR <- train_prosail_inversion(InputPROSAIL = InputPROSAIL,
                                            Parms2Estimate = Parms2Estimate,
                                            Bands2Select = Bands2Select,
                                            NoiseLevel = noise,
                                            SAILversion = '4SAIL',
                                            SRF = SRF,
                                            SpecPROSPECT = NULL,
                                            SpecSOIL = NULL,
                                            SpecATM = NULL,
                                            Path_Results = PROSAIL_ResPath,
                                            nbModels = nbModels,
                                            nbSamples = nbSamples,
                                            FigPlot = FigPlot)
        
        # Predict
        Apply_prosail_inversion(raster_path = refl_path, #reflectance_path
                                HybridModel = modelSVR,
                                PathOut = PROSAIL_ResPath,
                                SelectedBands = S2BandSelect,
                                bandname = ImgBandNames,
                                MultiplyingFactor = MultiplyingFactor)
      }
    }
  }
}