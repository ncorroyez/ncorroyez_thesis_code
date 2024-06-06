# ---
# title: "functions_explain_heterogeneity"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-14"
# ---

#' extract raster information in preparation for regression tasks: 
#' - requires path corresponding to rasters used as value to predict 
#' (predicted_path) and predictors (predictors_path)
#' - extract raster data and eliminates NA from extracted data
#' @param predicted_path data.frame path for individual variable to be predicted
#' @param predictors_path data.frame path for set of predictors
#' @return list including data.frames for predicted and predictors
#' @importFrom terra rast values countNA
#' @export

extract_raster_info <- function(predicted_path, predictors_path){
  # read predicted value and predictors included in raster data
  predicted <- terra::rast(predicted_path[[1]]) 
  predictors <- terra::rast(predictors_path$file)
  names(predictors) <- predictors_path$name
  # eliminate NA from all data when one layer has NA
  na_index_predictors <- terra::countNA(predictors)
  na_index_predicted <- terra::countNA(predicted)
  na_index <- na_index_predicted + na_index_predictors
  predictors[na_index > 0] <- NA
  predicted[na_index > 0] <- NA
  # extract values and eliminate NA
  predicted_val <- terra::values(predicted)
  predicted_val <- data.frame(predicted_val[complete.cases(predicted_val), ])
  names(predicted_val) <- names(predicted_path)
  predictors_val <- terra::values(predictors)
  predictors_val <- data.frame(predictors_val[complete.cases(predictors_val), ])
  names(predictors_val) <- predictors_path$name
  return(list('predictors_val' = predictors_val, 
              'predicted_val' = predicted_val))
}

#' produce training & test sets in preparation for ML 
#' initial data includes predicted and predictors
#' - requires predicted and predictors values, 
#' - extract raster data and eliminates NA from extracted data
#' @param predicted data.frame predicted values
#' @param predictors data.frame predictors values
#' @param training_sample_size numeric number of training samples
#' @param test_sample_size numeric number of test samples
#' @return list 
#' @importFrom terra rast values countNA
#' @export

train_test <- function(predicted, predictors, 
                       training_sample_size = 10000, test_sample_size = NULL){
  
  nbsamples <- nrow(predicted)
  # check sample size
  if (nbsamples < training_sample_size)
    stop('Error: training sample size smaller than number of samples available')
  if (nbsamples == training_sample_size)
    stop('Error: training sample size equals number of samples available')
  if (is.null(test_sample_size)) 
    test_sample_size <- nbsamples - training_sample_size
  if (training_sample_size + test_sample_size > nbsamples){
    test_sample_size <- nbsamples - training_sample_size
    warning('Warning: test sample size too large for nb of samples available')
    warning(paste('test sample size reduced to', test_sample_size))
    warning('adjust training sample size tpo increase test sample size')
  }
  # Train-test split
  train_ind <- sample(nbsamples, training_sample_size)
  trainingSet <- testSet <- list()
  trainingSet[['predicted']] <- data.frame(predicted[train_ind,])
  trainingSet[['predictors']] <- predictors[train_ind,]
  names(trainingSet[['predicted']]) <- names(predicted)
  
  testSet[['predicted']] <- data.frame(predicted[-train_ind,])
  testSet[['predictors']] <- predictors[-train_ind,]
  # adjust to desired number of samples
  testSet[['predicted']] <- data.frame(testSet[['predicted']][seq_len(test_sample_size), ])
  names(testSet[['predicted']]) <- names(predicted)
  testSet[['predictors']] <- testSet[['predictors']][seq_len(test_sample_size), ]
  return(list('trainingSet' = trainingSet,
              'testSet' = testSet))
}


#' Forward sequential feature selection (SFS) using a machine learning 
#' algorithm for regression task: rank features optimizing a criterion 
#' (here correlation) for regression task, using training and test sets
#' expects training and test dataset, 
#' can also save density scatterplot for each combination of features tested   
#' @param algo character. currently RF or liquidSVM
#' @param TrainingData numeric. training dataset
#' @param TestData numeric. test dataset
#' @param Vars_SFS character. variable names to be tested from training and test sets
#' @param ntree numeric. nb of trees for RF algo
#' @param site character. name of site
#' @param figure_path_site character. path where to save scatterplot
#' @param nbWorkers numeric. number of workers for parallel computing
#' @return list 
#' @importFrom terra rast values countNA
#' @export

SFS_Regression <- function(algo = 'RF', 
                           TrainingData,
                           TestData, 
                           Vars_SFS,
                           force_s2lai = FALSE,
                           ntree = NULL,
                           site = site, 
                           figure_path_site,
                           nbWorkers = 1,
                           corr_threshold = 0.01){
  NbPCs_To_Keep <- length(Vars_SFS)
  CorrSFS <- RMSESFS <- AssessSFS <- list()
  SelectedVars <- EvolCorr <- EvolRMSE <- c()
  last_corr <- -Inf  # Initialize last correlation with a very low value
  
  # Ensure s2lai is first in the list
  if ("s2lai_res_10_m" %in% Vars_SFS && force_s2lai == TRUE) {
    Vars_SFS <- c("s2lai_res_10_m",
                  Vars_SFS[-which(Vars_SFS == "s2lai_res_10_m")])
  }
  
  # multi-thread feature selection (SFS)
  doFuture::registerDoFuture()
  cl <- parallel::makeCluster(nbWorkers)
  plan("cluster", workers = cl)
  # progressbar
  pb <- progress::progress_bar$new(
    format = "Perform feature selection [:bar] :percent in :elapsedfull",
    total = NbPCs_To_Keep, clear = FALSE, width= 100)
  # explore each variable
  for (nbvars2select in seq_len(NbPCs_To_Keep)){
    NumVar_list <- as.list(seq_len(length(Vars_SFS)))
    subSFS <- SFS_Regression_step(algo = algo, 
                                  NumVar_list = NumVar_list, 
                                  SelectedVars = SelectedVars,
                                  Vars_SFS = Vars_SFS, 
                                  TrainingData = TrainingData,
                                  TestData = TestData, 
                                  ntree = ntree, 
                                  site = site, 
                                  figure_path_site = figure_path_site)
    pb$tick()
    
    CorrSFS[[nbvars2select]] <- data.frame(
      matrix(unlist(lapply(subSFS,'[[','coeffcorr')),ncol = 1))
    RMSESFS[[nbvars2select]] <- data.frame(
      matrix(unlist(lapply(subSFS, '[[', 'rmse')), ncol = 1))
    
    rownames(CorrSFS[[nbvars2select]]) <- Vars_SFS
    rownames(RMSESFS[[nbvars2select]]) <- Vars_SFS
    
    SelVar <- which(CorrSFS[[nbvars2select]] == max(
      CorrSFS[[nbvars2select]],na.rm = T))
    AssessSFS[[nbvars2select]] <- c(
      unlist((lapply(subSFS,'[[','AssessedVal'))[SelVar[1]]))
    
    current_corr <- CorrSFS[[nbvars2select]][SelVar[1], ]
    
    if (nbvars2select == 1 && force_s2lai == TRUE){
      WhichVar <- "s2lai_res_10_m"
      # add selected component to selected vars
      SelectedVars <- c(SelectedVars,WhichVar)
      
      # delete selected component from Vars_SFS
      Vars_SFS <- Vars_SFS[-which(Vars_SFS==WhichVar)]
      EvolCorr <- c(EvolCorr,CorrSFS[[nbvars2select]][SelVar[1],])
      EvolRMSE <- c(EvolRMSE, RMSESFS[[nbvars2select]][SelVar[1], ])
      last_corr <- current_corr  # Update last correlation
    }
    else{
      if ((current_corr - last_corr) < corr_threshold 
          && !is.null(corr_threshold)) {
        break  # Stop if the improvement is less than the threshold
      }
      WhichVar <- rownames(CorrSFS[[nbvars2select]])[SelVar[1]]
      SelectedVars <- c(SelectedVars, WhichVar)
      Vars_SFS <- Vars_SFS[-which(Vars_SFS == WhichVar)]
      EvolCorr <- c(EvolCorr, current_corr)
      EvolRMSE <- c(EvolRMSE, RMSESFS[[nbvars2select]][SelVar[1], ])
      last_corr <- current_corr  # Update last correlation
    }
  }
  parallel::stopCluster(cl)
  plan(sequential)
  return(list('SelectedVars' = SelectedVars, 
              'EvolCorr' = round(EvolCorr, 2),
              'EvolRMSE' = round(EvolRMSE, 2)))
}


#' subprocess for forward sequential feature selection using a machine learning 
#' algorithm for regression task: tests all features for a given number of 
#' features already extracted
#' expects training and test dataset, 
#' can also save density scatterpolot for each combination of features tested   
#' @param algo character. currently RF or liquidSVM
#' @param NumVar_list list. list of variables
#' @param SelectedVars character. 
#' @param Vars_SFS character. 
#' @param TrainingData numeric. training dataset
#' @param TestData numeric. test dataset
#' @param TestData numeric. test dataset
#' @param ntree numeric. nb of trees for RF algo
#' @param site character. name of site
#' @param figure_path_site character. path where to save scatterplot
#' @return list 
#' @importFrom terra rast values countNA
#' @export

SFS_Regression_step <- function(algo = 'RF',
                                NumVar_list, 
                                SelectedVars,
                                Vars_SFS, 
                                TrainingData,
                                TestData,
                                ntree = NULL,
                                site = NULL,
                                figure_path_site = NULL) {
  foreach(numvar = NumVar_list) %dopar% {
    SelFeat_tmp <- c(SelectedVars,Vars_SFS[[numvar]])
    Trainingsubset <- TrainingData[c('target',SelFeat_tmp)]
    Testsubset <- TestData[c('target',SelFeat_tmp)]
    # train random forest using full set of variables
    if (algo == 'RF'){
      formula <- target ~  .
      rf_model <- randomForest::randomForest(formula,
                                             data = Trainingsubset,
                                             ntree = ntree,
                                             importance = TRUE,
                                             do.trace = FALSE)
      
      # apply random forest on test dataset
      predictions <- predict(rf_model, newdata = Testsubset)
    }
    
    else if (algo == 'liquidSVM'){
      svmModel <- liquidSVM::svmRegression(as.matrix(TrainingData[SelFeat_tmp]),
                                           as.matrix(TrainingData['target']))
      # apply SVR on test dataset
      predictions <- predict(svmModel, TestData[SelFeat_tmp])
    } else {
      stop("Unsupported algorithm")
    }
    # produce density scatterplot
    if (!is.na(figure_path_site)){
      filename <- paste0('LiDAR_LAI_model_', 
                         algo, 
                         '_vars_',
                         paste(SelFeat_tmp, collapse = '_'))
      plot_density_scatterplot(c(predictions), 
                               Testsubset$target,
                               "Estimated LAI", 
                               "Measured LiDAR LAI",
                               title = paste("Estimated LAI vs LiDAR LAI:", 
                                             site),
                               xlimits = c(0,15), 
                               ylimits = c(0,15), 
                               dirname = figure_path_site, 
                               filename = filename)
    }
    CorrVal <- cor.test(c(predictions),Testsubset$target)$estimate
    RMSE <- Metrics::rmse(c(predictions), Testsubset$target)
    return(list('coeffcorr' = CorrVal, 
                'rmse' = RMSE,
                'AssessedVal' = c(predictions)))
  }
}

#' Save SFS Results to a Text File
#'
#' This function saves the results of Sequential Feature Selection (SFS), 
#' including the selected variables and their corresponding correlation values,
#' to a text file.
#'
#' @param selected_vars A character vector of selected variables.
#' @param evol_corr A numeric vector of correlation values corresponding to 
#' the selected variables.
#' @param evol_rmse A numeric vector of RMSE values corresponding to 
#' the selected variables.
#' @param figure_path_site A character string specifying the directory where 
#' the results file should be saved.
#' @param site A character string specifying the name of the site.
#' @param mask_name A character string specifying the name of the mask.
#' 
#' @return None. This function writes a text file to the specified location.
#' 
#' @details
#' The function combines the selected variables and their correlation values
#' into a data frame and writes this data frame to a text file. The file is 
#' named using the site and mask name, and is saved in the specified directory.
#' 
#' @examples
#' \dontrun{
#' save_results_to_file(selected_vars = c("var1", "var2"), 
#'                      evol_corr = c(0.85, 0.75), 
#'                      figure_path_site = "path/to/directory", 
#'                      site = "Site1", 
#'                      mask_name = "MaskA")
#' }
#'
#' @export
save_results_to_file <- function(selected_vars,
                                 evol_corr,
                                 evol_rmse,
                                 figure_path_site,
                                 site,
                                 mask_name) {
  results_file <- file.path(figure_path_site,
                            paste0(site, "_", mask_name, "_SFS_results.txt"))
  
  # Combine SelectedVars, EvolCorr, and EvolRMSE into a data frame
  results_df <- data.frame(
    SelectedVars = selected_vars,
    EvolCorr = evol_corr,
    EvolRMSE = evol_rmse,
    stringsAsFactors = FALSE
  )
  
  # Write the results to a text file
  write.table(results_df, 
              file = results_file,
              row.names = FALSE, 
              col.names = TRUE, 
              sep = "\t", 
              quote = FALSE)
}
