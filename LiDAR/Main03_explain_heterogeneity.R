# ------------------------------ Code Info -------------------------------------
# title: "Main03_explain_heterogeneity.R"
# authors: Nathan CORROYEZ & Jean-Baptiste FERET, 
# UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-24"

# -------------- (Optional) Clear the environment and free memory --------------
rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# ------------------------------ Libraries -------------------------------------
library(foreach)
# ------ Define working dir as the directory where the script is located -------
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# --------------------------- Import useful functions --------------------------
source("../libraries/functions_JBF.R")
# ----------------------------- RF hyperparameters ----------------------------
ntree <- 50
mtry <- 3
# ----------------------------- Plotting Parameters ----------------------------
default_par <- par(no.readonly = TRUE)
par(mar=c(5, 5, 4, 2) + 0.1) # Adjust margins as per your preference
par(oma=c(0, 0, 0, 0)) # Adjust outer margins as per your preference
par(cex=1.6, cex.axis=1.6, cex.names=1.6)
# --------------------------- Results Reproducibility --------------------------
training_sample_size <- 5000
test_sample_size <- 50000
# ---------------------------- Directories & Setup  ---------------------------- 
data <- "../../01_DATA"
lidar_dir <- "LiDAR"
figure_path <- '../../04_FIGURES'
dir.create(path = figure_path, showWarnings = F, recursive = T)
sites <- c('Mormal', 'Blois', 'Aigoual')

# ---------------------------------- Process  ---------------------------------- 
# nb workers in parallel for SFS
nbWorkers <- 4
# create site directories
results_path <- lidr_dir <- masks_dir <- list()
for (site in sites){
  results_path[[site]] <- file.path('../../03_RESULTS', site)
  lidr_dir[[site]] <- file.path(results_path[[site]], lidar_dir, "PAI/lidR")
  masks_dir[[site]] <- file.path(results_path[[site]], lidar_dir, "Heterogeneity_Masks")
  dir.create(path = lidr_dir[[site]], showWarnings = F, recursive = T)
  dir.create(path = masks_dir[[site]], showWarnings = F, recursive = T)
}

# --------------------------- Data processing ----------------------------------
regression_data <- train_test_data <- list()
SelectedVars <- EvolCorr <- SelectedVars_RF <- EvolCorr_RF <- list()
# process each site independently
for (site in sites){
  figure_path_site <- file.path(figure_path, site)
  dir.create(path = figure_path_site, showWarnings = F, recursive = T)

  # define path for LAI file derived from S2 and LiDAR
  S2_LAI <- file.path(results_path[[site]], "lai_s2_masked_res_10_m.envi")
  lidar_LAI <- file.path(results_path[[site]], "lai_lidar_masked_res_10_m.envi")
  
  # define path for additional lidar metrics
  lidar_metrics_path <- file.path(results_path[[site]],'lidar_metrics')
  lidar_metrics <- list.files(lidar_metrics_path)
  #-- truc degeulasse a virer quand tu auras harmonise les sorties
  lidar_metrics <- lidar_metrics[stringr::str_detect(string = lidar_metrics,'aux',negate = T)]
  lidar_metrics <- lidar_metrics[stringr::str_detect(string = lidar_metrics,'hdr',negate = T)]
  #-- 
  
  # create a dataframe for path of predicted value
  predicted_path <- data.frame('lidar_LAI' = lidar_LAI)
  
  # create a dataframe for path of predictors
  names_lidar_metrics <- unlist(lapply(stringr::str_split(string = lidar_metrics, 
                                                          pattern = '_'), '[[',1))
  predictors_path <- data.frame('name' = names_lidar_metrics, 
                           'file' = file.path(results_path[[site]],
                                              'lidar_metrics', lidar_metrics))
  predictors_path <- rbind(predictors_path,data.frame('name' = 'S2_LAI', 
                                            'file' = S2_LAI))
  
  # extract values to be used during regression
  regression_data[[site]] <- extract_raster_info(predicted_path = predicted_path, 
                                                 predictors_path = predictors_path)
  
  # scatterplot for direct comparison between S2 and LiDAR LAI
  filename <- paste0('S2_vs_LiDAR_LAI_', site)
  plot_density_scatterplot(regression_data[[site]]$predictors_val$S2_LAI,
                           regression_data[[site]]$predicted_val$lidar_LAI,
                           "Sentinel-2 LAI", "LiDAR LAI",
                           title = paste("S2 LAI vs LiDAR LAI:", site),
                           xlimits = c(0,15), ylimits = c(0,15), 
                           dirname = figure_path_site, filename = filename)
  
  # ----------- apply RF regression model using S2 & lidar data ----------------
  # create training and test data
  train_test_data[[site]] <- train_test(predicted = regression_data[[site]]$predicted_val, 
                                        predictors = regression_data[[site]]$predictors_val, 
                                        training_sample_size = training_sample_size,
                                        test_sample_size = test_sample_size)
  
  # train random forest using full set of variables
  TrainingData <- data.frame(train_test_data[[site]]$trainingSet)
  names(TrainingData) <- c('target', names(train_test_data[[site]]$trainingSet$predictors))
  formula <- target ~  S2_LAI + meanh + cv + rumple + vci + lcv + lskew
  rf_model <- randomForest::randomForest(formula, data = TrainingData, # mormal_ blois_ aigoual_ mix_
                                         ntree = ntree, mtry = mtry,
                                         importance = TRUE, do.trace = TRUE)
  
  # apply random forest on test dataset
  TestData <- data.frame(train_test_data[[site]]$testSet)
  names(TestData) <- c('target', names(train_test_data[[site]]$testSet$predictors))
  predictions <- predict(rf_model, newdata = TestData)
  cor.test(predictions,TestData$target)

  filename <- paste0('LiDAR_LAI_model_', site)
  plot_density_scatterplot(c(predictions), TestData$target,
                           "Estimated LAI", "Measured LiDAR LAI",
                           title = paste("Estimated LAI vs LiDAR LAI:", site),
                           xlimits = c(0,15), ylimits = c(0,15), 
                           dirname = figure_path_site, filename = filename)
  
  # perform SFS to assess LAI using optimal feature selection
  # ---------------------------------- SFS RF  ---------------------------------- 
  message('perform SFS RF')
  Vars_SFS <- names(train_test_data[[site]]$training$predictors)
  outSFS <- SFS_Regression(algo = 'RF', TrainingData = TrainingData, TestData = TestData, 
                           Vars_SFS = Vars_SFS, ntree = ntree, mtry = mtry, site = site, 
                           figure_path_site = figure_path_site, 
                           nbWorkers = nbWorkers)
  SelectedVars_RF[[site]] <- outSFS$SelectedVars
  EvolCorr_RF[[site]] <- outSFS$EvolCorr
}
