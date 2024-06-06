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
#' @param mtry numeric. nb of attempts for RF algo
#' @param site character. name of site
#' @param figure_path_site character. path where to save scatterplot
#' @param nbWorkers numeric. number of workers for parallel computing
#' @return list 
#' @importFrom terra rast values countNA
#' @export

SFS_Regression <- function(algo = 'RF', TrainingData, TestData, 
                           Vars_SFS, ntree = NULL, mtry = NULL, site = site, 
                           figure_path_site, nbWorkers = 1){
  NbPCs_To_Keep <- length(Vars_SFS)
  CorrSFS <- AssessSFS <- list()
  SelectedVars <- EvolCorr <- c()
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
    subSFS <- SFS_Regression_step(algo = algo, NumVar_list = NumVar_list, 
                                  SelectedVars = SelectedVars, Vars_SFS = Vars_SFS, 
                                  TrainingData = TrainingData, TestData = TestData, 
                                  ntree = ntree, mtry = mtry,
                                  site = site, figure_path_site = figure_path_site)
    pb$tick()
    CorrSFS[[nbvars2select]] <- data.frame(matrix(unlist(lapply(subSFS,'[[','coeffcorr')),ncol = 1))
    rownames(CorrSFS[[nbvars2select]]) <- Vars_SFS
    SelVar <- which(CorrSFS[[nbvars2select]] == max(CorrSFS[[nbvars2select]],na.rm = T))
    AssessSFS[[nbvars2select]] <- c(unlist((lapply(subSFS,'[[','AssessedVal'))[SelVar[1]]))
    # which criterion to maximize with SFS?
    WhichVar <- rownames(CorrSFS[[nbvars2select]])[SelVar[1]]
    # add selected component to selected vars
    SelectedVars <- c(SelectedVars,WhichVar)
    # delete selected component from Vars_SFS
    Vars_SFS <- Vars_SFS[-which(Vars_SFS==WhichVar)]
    EvolCorr <- c(EvolCorr,CorrSFS[[nbvars2select]][SelVar[1],])
  }
  parallel::stopCluster(cl)
  plan(sequential)
  return(list('SelectedVars' = SelectedVars, 
              'EvolCorr' = EvolCorr))
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
#' @param mtry numeric. nb of attempts for RF algo
#' @param site character. name of site
#' @param figure_path_site character. path where to save scatterplot
#' @return list 
#' @importFrom terra rast values countNA
#' @export

SFS_Regression_step <- function(algo = 'RF', NumVar_list, SelectedVars, Vars_SFS, 
                                TrainingData, TestData, ntree = NULL, mtry = NULL, 
                                site = NULL, figure_path_site = NULL) {
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
                                             mtry = mtry,
                                             importance = TRUE,
                                             do.trace = FALSE)
      
      # apply random forest on test dataset
      predictions <- predict(rf_model, newdata = Testsubset)
    }
    
    if (algo == 'liquidSVM'){
      svmModel <- liquidSVM::svmRegression(as.matrix(TrainingData[SelFeat_tmp]),
                                           as.matrix(TrainingData['target']))
      # apply SVR on test dataset
      predictions <- predict(svmModel, TestData[SelFeat_tmp])
    }
    # produce density scatterplot
    if (!is.na(figure_path_site)){
      filename <- paste0('LiDAR_LAI_model_', algo, '_vars_',paste(SelFeat_tmp, collapse = '_'))
      plot_density_scatterplot(c(predictions), Testsubset$target,
                               "Estimated LAI", "Measured LiDAR LAI",
                               title = paste("Estimated LAI vs LiDAR LAI:", site),
                               xlimits = c(0,15), ylimits = c(0,15), 
                               dirname = figure_path_site, filename = filename)
    }
    CorrVal <- cor.test(c(predictions),Testsubset$target)$estimate
    return(list('coeffcorr' = CorrVal, 
                'AssessedVal' = c(predictions)))
  }
}
