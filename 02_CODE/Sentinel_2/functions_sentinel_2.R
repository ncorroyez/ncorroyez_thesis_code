# ---
# title: "function_sentinel_2.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, 
# CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-11-07"
# ---

library(preprocS2)
library(prosail)
library(spinR)
library(raster)
library(stars)
library(terra)
library(ggplot2)
library(truncnorm)
library(bigRaster)

download_S2_images <- function(dateAcq, path_vector, output_directory) {
  
  # Check if already downloaded
  date_without_hyphens <- gsub("-", "", dateAcq)
  subdirectories <- list.dirs(output_directory, # list all S2_Images
                              full.names = FALSE, recursive = FALSE)
  
  # Skip if yes
  if (any(grepl(date_without_hyphens, subdirectories))) {
    cat("Skipping download for date", dateAcq, 
        "because it has been already downloaded in", output_directory, ".\n")
  } else {
    Path_S2 <- get_S2_L2A_Image(l2a_path = output_directory, 
                                spatial_extent = path_vector, 
                                dateAcq = dateAcq,
                                DeleteL1C = TRUE, 
                                Sen2Cor = TRUE) # Fastest way to download
    cat("S2 image acquired at", dateAcq, 
        "has been successfully created and saved at :", output_directory, ".\n")
  }
}

process_S2_data <- function(output_directory, 
                            path_vector, 
                            result_path, 
                            resolution, 
                            S2source) {
  # Define raster path
  Path_S2 <- file.path(output_directory, list.files(output_directory, 
                                                    pattern = '.SAFE'))
  
  # Create the result directory if it doesn't exist
  dir.create(path = result_path, showWarnings = FALSE, recursive = TRUE)
  
  # Extract, resample, and stack data
  S2obj <- preprocS2::extract_from_S2_L2A(
    Path_dir_S2 = Path_S2,
    path_vector = path_vector,
    S2source = S2source,
    resolution = resolution
  )
  cat("This image has been successfully processed 
      and the results are saved at :", Path_S2, "\n")
  # Update shapefile path if needed (reprojection)
  # forest_mask <- S2obj$forest_mask
  return(list(S2obj = S2obj,
              path_vector = S2obj$path_vector))
}

generate_cloud_mask <- function(S2_Stack, Cloud_path, S2source, saveRaw = TRUE) {
  # Write cloud mask and get filenames
  cloudmasks <- preprocS2::save_cloud_s2(
    S2_stars = S2_Stack,
    Cloud_path = Cloud_path,
    S2source = S2source,
    SaveRaw = saveRaw
  )
  cat("Cloud mask has been successfully created and saved at :", 
      Cloud_path, "\n")
  
  return(cloudmasks)
}

write_reflectance <- function(S2_Stack, S2_Bands, refl_path) {
  # Save Reflectance file as ENVI image with BIL interleaves
  tile_S2 <- preprocS2::get_tile(S2_Bands$GRANULE)
  dateAcq_S2 <- preprocS2::get_date(S2_Bands$GRANULE)
  
  preprocS2::save_reflectance_s2(
    S2_stars = S2_Stack, 
    Refl_path = refl_path,
    S2Sat = NULL, 
    tile_S2 = tile_S2, 
    dateAcq_S2 = dateAcq_S2,
    Format = 'ENVI', 
    datatype = 'Int16', 
    MTD = S2_Bands$metadata, 
    MTD_MSI = S2_Bands$metadata_MSI
  )
  cat("Reflectance has been successfully created and saved at :", 
      refl_path, "\n")
}

read_raster <- function(refl_path) {
  Refl <- brick(refl_path)
  return(Refl)
}

compute_spectral_indices <- function(Refl, SensorBands, IndexList, ReflFactor) {
  Indices <- spinR::compute_S2SI_Raster(
    Refl = Refl,
    SensorBands = SensorBands,
    Sel_Indices = IndexList,
    ReflFactor = ReflFactor,
    StackOut = FALSE
  )
  return(Indices)
}

save_spectral_indices <- function(Indices, 
                                  SI_path, 
                                  S2_Bands) {
  for (SpIndx in names(Indices$SpectralIndices)) {
    Index_Path <- file.path(SI_path, paste(basename(S2_Bands$GRANULE), 
                                           '_', SpIndx, sep = ''))
    stars::write_stars(st_as_stars(Indices$SpectralIndices[[SpIndx]]), 
                       dsn = Index_Path, driver =  "ENVI", type = 'Float32')
    # Write band name in HDR
    HDR <- read_ENVI_header(get_HDR_name(Index_Path))
    HDR$`band names` <- SpIndx
    write_ENVI_header(HDR = HDR, HDRpath = get_HDR_name(Index_Path))
  }
}

update_cloud_mask <- function(CloudPath, NDVI_Thresh, Indices, cloudmasks) {
  Elim <- which(values(Indices$SpectralIndices[['NDVI']]) < NDVI_Thresh)
  CloudInit <- stars::read_stars(cloudmasks$BinaryMask)
  CloudInit$CloudMask_Binary[Elim] <- 0
  # Save the updated cloud mask
  Cloud_File <- file.path(CloudPath, 'CloudMask_Binary_Update')
  stars::write_stars(CloudInit, dsn = Cloud_File, 
                     driver = "ENVI", type = 'Byte')
  return(Cloud_File)
}

preprocess_S2 <- function(dateAcq, 
                          path_vector,
                          output_directory,
                          result_path, 
                          resolution = resolution, 
                          S2source = 'SAFE',
                          saveRaw = TRUE) {
  download_S2_images(dateAcq, path_vector, output_directory)
  process_S2_data <- process_S2_data(output_directory, path_vector, result_path, 
                                     resolution, S2source)
  
  S2obj <- process_S2_data$S2obj
  path_vector <- process_S2_data$path_vector
  
  results_site_path <- file.path(result_path,basename(S2obj$S2_Bands$GRANULE))
  dir.create(path = results_site_path,showWarnings = FALSE,recursive = TRUE)
  
  cloud_path <- file.path(results_site_path,
                          'CloudMask',
                          sprintf("res_%d_m", resolution))
  dir.create(path = cloud_path, showWarnings = FALSE, recursive = TRUE)
  
  refl_dir <- file.path(results_site_path, 
                        'Reflectance',
                        sprintf("res_%d_m", resolution))
  dir.create(path = refl_dir, showWarnings = FALSE, recursive = TRUE)
  refl_path <- file.path(refl_dir, paste(basename(S2obj$S2_Bands$GRANULE), 
                                         '_Refl', sep = ''))
  
  cloudmasks <- generate_cloud_mask(S2obj$S2_Stack, cloud_path, 
                                    S2source, saveRaw)
  write_reflectance(S2obj$S2_Stack, S2obj$S2_Bands, refl_path)
  
  # Compute spectral indices
  IndexList <- c('NDVI')
  ReflFactor <- 10000
  HDR_Refl <- read_ENVI_header(get_HDR_name(refl_path))
  SensorBands <- HDR_Refl$wavelength
  Refl <- read_raster(refl_path)
  Indices <- compute_spectral_indices(Refl, SensorBands, IndexList, ReflFactor)
  
  # Save spectral indices
  SI_path <- file.path(results_site_path, 
                       'SpectralIndices',
                       sprintf("res_%d_m", resolution))
  dir.create(path = SI_path,showWarnings = FALSE, recursive = TRUE)
  save_spectral_indices(Indices, SI_path, S2obj$S2_Bands)
  
  # Update Cloud mask based on radiometric filtering
  NDVI_Thresh <- 0.5
  Cloud_File <- update_cloud_mask(cloud_path, NDVI_Thresh, Indices, cloudmasks)
  
  return(list(results_site_path = results_site_path,
              refl_path = refl_path, 
              HDR_Refl = HDR_Refl, 
              Cloud_File = Cloud_File))
}

calculate_angles <- function(sza, saa, vza, vaa) {
  # Convert angles to radians
  sza_rad <- sza * pi / 180
  saa_rad <- saa * pi / 180
  vza_rad <- vza * pi / 180
  vaa_rad <- vaa * pi / 180
  
  tto <- vza
  tts <- sza
  
  cos_psi <- -cos(sza_rad) * cos(vza_rad) + sin(sza_rad) * sin(vza_rad) * cos(saa_rad - vaa_rad)
  psi <- acos(cos_psi) * 180 / pi  # Convert back to degrees
  
  return(list(tto = tto, tts = tts, psi = psi))
}

get_s2_geometry <- function(refl_path) {
  xmlfile <- file.path(dirname(refl_path), 'MTD_TL.xml')
  S2Geom <- get_S2geometry(MTD_TL_xml = xmlfile) 
  
  # Go from sza, saa, vza, vaa to tto, tts, psi
  S2Geom <- calculate_angles(S2Geom$SZA, S2Geom$SAA, S2Geom$VZA, S2Geom$VAA)
  
  # Create min and max lists
  GeomAcq <- list()
  GeomAcq$min <- GeomAcq$max <- list()
  GeomAcq$min$tto <- min(S2Geom$tto)
  GeomAcq$max$tto <- max(S2Geom$tto)
  GeomAcq$min$tts <- min(S2Geom$tts)
  GeomAcq$max$tts <- max(S2Geom$tts)
  GeomAcq$min$psi <- min(S2Geom$psi)
  GeomAcq$max$psi <- max(S2Geom$psi)
  return(GeomAcq)
}

get_sensor_response <- function(HDR_Refl, Path_SensorResponse = NULL) {
  SensorName <- HDR_Refl$`sensor type` # Sentinel_2A
  SRF <- GetRadiometry(SensorName, Path_SensorResponse = Path_SensorResponse)
  return(SRF)
}

adjust_optical_constants <- function(SRF) {
  if (is.null(SpecPROSPECT)){
    SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
  }
  if (is.null(SpecSOIL)){
    SpecSOIL <- prosail::SpecSOIL
  }
  if (is.null(SpecPROSPECT)){
    SpecATM <- prosail::SpecATM
  }
  wvl <- SpecPROSPECT$lambda
  SpecSensor <- PrepareSensorSimulation(SpecPROSPECT, SpecSOIL, SpecATM, SRF)
  SpecPROSPECT_Sensor <- SpecSensor$SpecPROSPECT_Sensor
  SpecSOIL_Sensor <- SpecSensor$SpecSOIL_Sensor
  SpecATM_Sensor <- SpecSensor$SpecATM_Sensor
  check_SpectralSampling(SpecPROSPECT, SpecSOIL, SpecATM)
  
  return(list(SpecSensor = SpecSensor,
              SpecPROSPECT_Sensor = SpecPROSPECT_Sensor, 
              SpecSOIL_Sensor = SpecSOIL_Sensor, 
              SpecATM_Sensor = SpecATM_Sensor))
}

define_spectral_bands_and_variables <- function(HDR_Refl, S2BandSel) {
  ImgBandNames <- strsplit(HDR_Refl$`band names`, split = ',')[[1]]
  Bands2Select <- list()
  for (bpvar in names(S2BandSel)) {
    Bands2Select[[bpvar]] <- match(S2BandSel[[bpvar]], ImgBandNames)
  }
  for (parm in Parms2Estimate){
    S2BandSelect[[parm]] <- c('B3','B4','B5','B6','B7','B8','B11','B12')
    Bands2Select[[parm]] <- match(S2BandSelect[[parm]], SRF$Spectral_Bands)
  }
  return(list(S2BandSel = S2BandSel, 
              ImgBandNames = ImgBandNames,
              Bands2Select = Bands2Select))
}

define_noise_levels <- function() {
  NoiseLevel <- list()
  NoiseLevel$lai <- 0.05
  return(NoiseLevel)
}         

define_GeomAcq <- function(){
  GeomAcq <- list()
  GeomAcq$min <- GeomAcq$max <- list()
  GeomAcq$min$tto <- 0
  GeomAcq$max$tto <- 10
  GeomAcq$min$tts <- 20
  GeomAcq$max$tts <- 30
  GeomAcq$min$psi <- 0
  GeomAcq$max$psi <- 360
  return(GeomAcq)
}

define_atbd_GeomAcq <- function(){
  GeomAcq <- list()
  GeomAcq$min <- GeomAcq$max <- list()
  GeomAcq$min$tto <- 5
  GeomAcq$max$tto <- 15
  GeomAcq$min$tts <- 25
  GeomAcq$max$tts <- 35
  GeomAcq$min$psi <- 180
  GeomAcq$max$psi <- 210
  return(GeomAcq)
}

# Function to extract min, max, mean, and std for all variables
extract_summary_stats <- function(data, filepath) {
  # Initialize an empty data frame to store the statistics
  statistics_df <- data.frame(Variable = character(0), 
                              Min = numeric(0), 
                              Max = numeric(0), 
                              Mean = numeric(0), 
                              StdDev = numeric(0))
  
  # Loop through each variable in the data frame and compute summary statistics
  for (variable_name in names(data)) {
    variable <- data[[variable_name]]
    min_value <- min(variable, na.rm = TRUE)
    max_value <- max(variable, na.rm = TRUE)
    mean_value <- mean(variable, na.rm = TRUE)
    std_dev_value <- sd(variable, na.rm = TRUE)
    
    # Append the statistics to the data frame
    statistics_df <- rbind(statistics_df, 
                           data.frame(Variable = variable_name, 
                                      Min = min_value,
                                      Max = max_value,
                                      Mean = mean_value,
                                      StdDev = std_dev_value))
  }
  write.table(statistics_df,
              file = file.path(filepath, "PROSAIL_Stats.txt"), 
              sep = "\t", row.names = FALSE)
  
  return(statistics_df)
}

generate_samples <- function(mean, std, min_val, max_val, size) {
  samples <- rtruncnorm(n = size, a = min_val, b = max_val, mean = mean, sd = std)
  return(samples)
}

modify_parameter_distribution <- function(InputPROSAIL,
                                          parameter_name,
                                          desired_distribution,
                                          minval = NULL,
                                          maxval = NULL,
                                          mean = NULL,
                                          std = NULL,
                                          unique_val = NULL) {
  set.seed(0)
  if (parameter_name %in% names(InputPROSAIL)) {
    cat("Modify", parameter_name, "\n")
    if (desired_distribution == 'Gaussian') {
      cat("Desired distribution", desired_distribution, "\n")
      mean_val <- ifelse(!is.null(mean), unname(mean), GaussianDistrib$Mean[[parameter_name]])
      std_val <- ifelse(!is.null(std), unname(std), GaussianDistrib$Std[[parameter_name]])
      min_val <- ifelse(!is.null(minval), unname(minval), -Inf)
      max_val <- ifelse(!is.null(maxval), unname(maxval), Inf)
      
      new_values <- generate_samples(mean_val, std_val, min_val, max_val, 
                                     length(InputPROSAIL[[parameter_name]]))
    } else if (desired_distribution == 'Uniform') {
      cat("Desired distribution", desired_distribution, "\n")
      min_val <- ifelse(!is.null(minval), unname(minval), 0)
      max_val <- ifelse(!is.null(maxval), unname(maxval), 1)
      
      new_values <- runif(length(InputPROSAIL[[parameter_name]]), 
                          min_val, max_val)
    } else if (desired_distribution == 'Unique') {
      cat("Desired distribution", desired_distribution, "\n")
      # Set a unique value for the parameter for all samples
      unique_val <- ifelse(!is.null(unique_val), unname(unique_val), unique_val[[parameter_name]])
      new_values <- rep(unique_val, length(InputPROSAIL[[parameter_name]]))
    } else {
      warning("Unknown distribution type. No changes made.")
      new_values <- InputPROSAIL[[parameter_name]]
    }
    
    # Update the specified parameter with new values
    InputPROSAIL[[parameter_name]] <- new_values
  } else {
    warning(paste("Parameter", parameter_name, "not found in InputPROSAIL. 
                  No changes made. Choose one parameter in this list:", 
                  str(InputPROSAIL)))
  }
  
  return(InputPROSAIL)
}

# ------------------------------------------- Distribution Functions -------------------------------------------------

set_prosail_distribution <- function(distrib_option, refl_path) {
  switch(
    distrib_option,
    "q_zhang_et_al_2005" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 7.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 80)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 4.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.15)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.04)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.00001,
                                                    maxval = 8)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "TypeLidf",
                                                    "Unique",
                                                    unique_val = 2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 10,
                                                    maxval = 89)
    },
    "brede_et_al_2020" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 8)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 80)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 2.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.002,
                                                    maxval = 0.025)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.025)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 1)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tto",
                                                    "Unique",
                                                    unique_val = 0)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Uniform",
                                                    minval = 27.5,
                                                    maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psi",
                                                    "Unique",
                                                    unique_val = 0)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Unique",
                                                    unique_val = 0)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Unique",
                                                    unique_val = 0)
    },
    "hauser_et_al_2021" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Gaussian",
                                                    mean = 2,
                                                    std = 1,
                                                    minval = 0.01,
                                                    maxval = 3.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Gaussian",
                                                    mean = 30,
                                                    std = 20,
                                                    minval = 10,
                                                    maxval = 60)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CAR",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 15)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "ANT",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 10)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1.4,
                                                    maxval = 1.7)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.045)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.040)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Unique",
                                                    unique_val = 0.01)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "TypeLidf",
                                                    "Unique",
                                                    unique_val = 2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 30,
                                                    maxval = 70)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Unique",
                                                    unique_val = 0.01)
    },
    "sinha_et_al_2020" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Gaussian",
                                                    mean = 1.5,
                                                    std = 0.7,
                                                    minval = 0,
                                                    maxval = 2.8)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Gaussian",
                                                    mean = 30,
                                                    std = 20,
                                                    minval = 0,
                                                    maxval = 60)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1.5,
                                                    maxval = 2.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.030,
                                                    maxval = 0.050)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.012,
                                                    maxval = 0.030)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 20,
                                                    maxval = 50)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Unique",
                                                    unique_val = 0.01)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Unique",
                                                    unique_val = 47)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0.2,
                                                    maxval = 1)
    },
    "verhoef_and_bach_2007" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- -0.3
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Unique",
                                                    unique_val = 2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Unique",
                                                    unique_val = 60)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Unique",
                                                    unique_val = 0.009)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Unique",
                                                    unique_val = 0.005)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Unique",
                                                    unique_val = -0.2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Unique",
                                                    unique_val = 0.05)
    },
    "shiklomanov_et_al_2016" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- -0.3
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1.09,
                                                    maxval = 3)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0.78,
                                                    maxval = 106.72)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CAR",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 25.3)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.0043,
                                                    maxval = 0.0439)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.0017,
                                                    maxval = 0.0152)
    },
    "atbd" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE,
                                       GeomAcq = define_atbd_GeomAcq())
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_JBGeom" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE,
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_n_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 3)
    },
    "atbd_n_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 1.2)
    },
    "atbd_n_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 5)
    },
    "atbd_chl_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 40,
                                                    maxval = 100)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
    },
    "atbd_chl_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 60)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
    },
    "atbd_chl_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 100)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
    },
    "atbd_brown_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.00001,
                                                    maxval = 3)
    },
    "atbd_brown_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.01,
                                                    maxval = 1)
    },
    "atbd_brown_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.00001,
                                                    maxval = 8)
    },
    "atbd_ewt_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.02,
                                                    maxval = 0.05)
    },
    "atbd_ewt_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.02)
    },
    "atbd_ewt_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.05)
    },
    "atbd_lma_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.02,
                                                    maxval = 0.04)
    },
    "atbd_lma_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.02)
    },
    "atbd_lma_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.04)
    },
    "atbd_lidfa_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$TypeLidf <- 2
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 50,
                                                    maxval = 90)
    },
    "atbd_lidfa_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$TypeLidf <- 2
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 10,
                                                    maxval = 50)
    },
    "atbd_lidfa_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$TypeLidf <- 2
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 10,
                                                    maxval = 90)
    },
    "atbd_lai_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 7)
    },
    "atbd_lai_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 3)
    },
    "atbd_lai_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 10)
    },
    "atbd_q_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Uniform",
                                                    minval = 0.25,
                                                    maxval = 0.5)
    },
    "atbd_q_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0.25)
    },
    "atbd_q_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0.5)
    },
    "atbd_psoil_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0.5,
                                                    maxval = 1)
    },
    "atbd_psoil_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0.5)
    },
    "atbd_psoil_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 1)
    },
    "atbd_S2Geom" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = get_s2_geometry(refl_path))
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_fixed_tto" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$tto <- 5
    },
    "atbd_tto_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tto",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 5)
    },
    "atbd_tto_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tto",
                                                    "Uniform",
                                                    minval = 5,
                                                    maxval = 10)
    },
    "atbd_fixed_tts" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$tts <- 28
    },
    "atbd_tts_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Uniform",
                                                    minval = 20,
                                                    maxval = 25)
    },
    "atbd_tts_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Uniform",
                                                    minval = 25,
                                                    maxval = 30)
    },
    "atbd_fixed_psi" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$psi <- 150
    },
    "atbd_psi_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psi",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 180)
    },
    "atbd_psi_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psi",
                                                    "Uniform",
                                                    minval = 180,
                                                    maxval = 360)
    },
    "atbd_tts_surprising" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Uniform",
                                                    minval = 27.5,
                                                    maxval = 80)
    },
    {
      stop(paste("Error:", 
                 distrib_option,
                 "is an invalid distribution option.",
                 "Please refer to the function code to see the possibilities"))
    }
  )
  
  return(InputPROSAIL)
}
