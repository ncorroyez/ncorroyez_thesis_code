# ---
# title: "function_plots.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-11-01"
# ---

library(terra)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(hrbrthemes)
library(raster)
library(rasterVis)
library(latticeExtra)
library(RColorBrewer)

#' Save a raster file in ENVI format
#' 
#' This function saves a raster object in ENVI format.
#'
#' @param raster_to_save_in_envi Raster object to be saved.
#' @param output_name Name of the output file without extension.
#' @param results_path Directory where the file will be saved.
save_envi_file <- function(raster_to_save_in_envi, output_name, results_path) {
  if (!dir.exists(results_path)) {
    dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
  }
  filename <- file.path(results_path, paste0(output_name, ".envi"))
  writeRaster(raster_to_save_in_envi, filename = filename, overwrite = TRUE)
  cat("ENVI file has been successfully created and saved at:", filename, "\n")
}

#' Save a raster file in TIF format
#' 
#' This function saves a raster object in TIF format.
#'
#' @param raster_to_save_in_tif Raster object to be saved.
#' @param output_name Name of the output file without extension.
#' @param results_path Directory where the file will be saved.
save_tif_file <- function(raster_to_save_in_tif, 
                          output_name,
                          results_path){
  if (!dir.exists(results_path)) {
    dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
  }
  filename <- file.path(results_path, paste0(output_name, "tif"))
  writeRaster(raster_to_save_in_tif, 
              filename = filename, 
              overwrite=TRUE)
  cat("TIF file has been successfully created and saved at :", 
      filename, "\n")
}

#' Open a raster file
#' 
#' This function opens a raster file with supported extensions
#' ("tif", "tiff", "envi").
#'
#' @param raster_path Path to the raster file.
#' @return Raster object.
open_raster_file <- function(raster_path){
  file_extension <- tools::file_ext(raster_path)
  valid_extensions <- c("tif", "tiff", "envi")
  if (length(file_extension) == 1 && file_extension %in% valid_extensions) {
    return(terra::rast(raster_path))
  }
  else {
    stop("Unsupported file format. 
         Please provide a file with one of these extensions: ",  
         valid_extensions)
  }
}

#' Keep positive values in a raster
#' 
#' This function replaces negative values in a raster with zero.
#'
#' @param raster Raster object.
#' @return Raster object with non-negative values.
keep_positive_values <- function(raster) {
  raster_positive_values <- values(raster)
  raster_positive_values[raster_positive_values < 0] <- 0
  return(raster_positive_values)
}

#' Open raster file and return values
#' 
#' This function opens a raster file and returns its values.
#'
#' @param raster_path Path to the raster file.
#' @return Values of the raster.
open_raster_file_as_values <- function(raster_path) {
  raster_mask <- open_raster_file(raster_path)
  return(values(raster_mask))
}

#' Open raster file and return non-negative values
#' 
#' This function opens a raster file and returns its non-negative values.
#'
#' @param raster_path Path to the raster file.
#' @return Non-negative values of the raster.
open_raster_file_as_positive_values <- function(raster_path) {
  raster_mask <- open_raster_file(raster_path)
  return(keep_positive_values(raster_mask))
}

#' Mask low vegetation in a raster
#' 
#' This function masks low vegetation areas in a raster based on a mask.
#'
#' @param var_matrix Matrix of raster values.
#' @param var_raster Raster object to be masked.
#' @param low_vegetation Logical mask for low vegetation.
#' @return List containing masked raster and masked values.
mask_low_vegetation <- function(var_matrix, var_raster, low_vegetation){
  # Mask
  index <- low_vegetation == 1 # Index where vegetation is high
  masked_values <- var_matrix
  masked_values[!index] <- NA # Mask the low vegetation (index=0)
  
  # Take the reference raster and replace its values by the masked raster 
  masked_raster <- var_raster
  values(masked_raster) <- masked_values
  setMinMax(masked_raster)
  
  return(list(masked_raster = masked_raster,
              masked_values = masked_values)
  )
}

#' Save a basic plot
#' 
#' This function saves a basic plot to a PNG file.
#'
#' @param plot_to_save Plot object to be saved.
#' @param dirname Directory where the file will be saved.
#' @param filename Name of the output file.
#' @param title Title of the plot.
#' @param width_pixels Width of the plot in pixels.
#' @param height_pixels Height of the plot in pixels.
#' @param res Resolution of the plot.
save_basic_plot <- function(plot_to_save,
                            dirname,
                            filename,
                            title = NULL,
                            width_pixels = 1920,
                            height_pixels = 1080,
                            res = 200) {
  png(file.path(dirname, filename), 
      width = width_pixels, 
      height = height_pixels, 
      units = "px", res = res)
  plot(plot_to_save, main = title)
  dev.off()
}

#' Save an x-y plot
#' 
#' This function saves an x-y plot to a PNG file.
#'
#' @param xvar X values for the plot.
#' @param yvar Y values for the plot.
#' @param dirname Directory where the file will be saved.
#' @param filename Name of the output file.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param title Title of the plot.
#' @param width_pixels Width of the plot in pixels.
#' @param height_pixels Height of the plot in pixels.
#' @param res Resolution of the plot.
save_x_y_plot <- function(xvar,
                          yvar,
                          dirname,
                          filename,
                          xlab = NULL,
                          ylab = NULL,
                          title = NULL,
                          width_pixels = 1920,
                          height_pixels = 1080,
                          res = 200) {
  png(file.path(dirname, filename), 
      width = width_pixels, 
      height = height_pixels, 
      units = "px", res = res)
  plot(xvar, yvar, main = title, xlab=xlab, ylab=ylab)
  dev.off()
}

#' Plot density scatterplot
#' 
#' This function plots a density scatterplot and optionally saves it.
#'
#' @param var_x X values for the plot.
#' @param var_y Y values for the plot.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param title Title of the plot.
#' @param dirname Directory where the file will be saved.
#' @param filename Name of the output file.
#' @param xlimits Limits for the x-axis.
#' @param ylimits Limits for the y-axis.
plot_density_scatterplot <- function(var_x,
                                     var_y,
                                     xlab, 
                                     ylab,
                                     title = NULL,
                                     dirname = NULL, 
                                     filename = NULL,
                                     xlimits = NULL,
                                     ylimits = NULL) {
  
  if (!is.vector(var_x)) {
    var_x <- as.vector(var_x)
  }
  if (!is.vector(var_y)) {
    var_y <- as.vector(var_y)
  }
  
  data <- data.frame(
    my_x = var_x,
    my_y = var_y
  )
  
  if (is.null(xlimits)){
    xlimits <- c(min(data$my_x - 0.5), max(data$my_x + 0.5))
  }
  if (is.null(ylimits)){
    ylimits <- c(min(data$my_y - 0.5), max(data$my_y + 0.5))
  }
  
  lm_model <- lm(data$my_y ~ data$my_x, data = data)
  # print(coef(lm_model))
  
  # equation <- sprintf("italic(y) == %.3f * italic(x) + %.3f ",
  #                     coef(lm_model)[2],
  #                     coef(lm_model)[1])
  equation <- sprintf("y2 = %.3f x + %.3f",
                      round(coef(lm_model)[2],2),
                      round(coef(lm_model)[1],2))
  
  statsReg1 <- cor.test(data$my_x, data$my_y)$estimate
  statsReg21 <- summary(lm_model)$r.squared
  
  ssres <- sum((data$my_y - data$my_x)^2)
  sstot <- sum((data$my_y - mean(data$my_y))^2)
  # print(ssres)
  # print(sstot)
  r_squared <- 1 - (ssres / sstot)
  
  statsReg2 <- r_squared
  statsReg3 <- Metrics::rmse(na.omit(data$my_y), na.omit(data$my_x))
  
  # Set the text size parameters
  text_size <- 12 
  
  p3 <- ggplot(data, aes(x=var_x, y=var_y)) +
    # ggtitle(paste(sprintf("R: %.3f,", statsReg1), 
    #               sprintf("R² (model): %.3f,", statsReg21),
    #               sprintf("R² (formula): %.3f,", statsReg2),
    #               sprintf("RMSE: %.3f,", statsReg3),
    #               sprintf("Equation: %s", equation))) +
    ggtitle(title) +
    xlab(xlab) + ylab(ylab) +
    scale_fill_continuous(type = "viridis") +
    scale_x_continuous(limits = xlimits) +
    scale_y_continuous(limits = ylimits) +
    # geom_bin2d(bins = 400) +
    # geom_smooth(method="lm" , color="blue", se=TRUE, lwd=1.5) +
    geom_bin2d(bins = 400) +  # Fill aesthetic mapped to bin count
    geom_smooth(method="lm" , color="blue", se=TRUE, lwd=1.5) +
    geom_abline(slope = 1, intercept = 0, color="red", lwd=1.5) + 
    theme(
      # Adjust main title size
      plot.title = element_text(size = text_size * 1.2, face = "bold", hjust = 0.5),
      
      # Adjust axis label size
      axis.title.x = element_text(size = text_size * 1.2, color = "black"),
      axis.title.y = element_text(size = text_size * 1.2, color = "black"),
      axis.text = element_text(size = text_size * 1.2, color = "black"),
      legend.text = element_text(size = text_size * 1.2, color = "black")
      
      # Adjust text annotations size
      # Uncomment the lines below if you want to add text annotations
      # axis.text = element_text(size = text_size * 0.8), 
      # axis.text.x = element_text(size = text_size * 0.8),
      # axis.text.y = element_text(size = text_size * 0.8)
    ) +
    annotate("text", x = 8.5, y = 5,
             label = sprintf("%s", equation), size = text_size * 0.5, color = "black") +
    annotate("text", x = 9, y = 8,
             label = "y1 = x", size = text_size * 0.5, color = "black")
    # annotate("text", x = 1, y = 9.75,
    #          label = sprintf("R: %.2f", statsReg1), size = text_size * 0.5, color = "black") +
    # annotate("text", x = 1.0685, y = 9,
    #          label = sprintf("R²: %.2f", statsReg21), size = text_size * 0.5, color = "black") +
    # coord_fixed(ratio = 1)
  
  plot(p3)
  
  if (!is.null(dirname) && !is.null(filename)){
    file_path <- file.path(dirname, filename)
    ggsave(filename = paste0(file_path, ".png"),
           plot = p3, device = "png", scale = 1, 
           width = 1920, height = 1080, units = "px", dpi = 100)
    cat("Plots saved at: ", paste0(file_path, ".png"), "\n")
  }
}

#' Plot Histogram
#'
#' This function creates and optionally saves a histogram for one or more
#' sets of variables.
#'
#' @param vars_list A list or vector of numeric values to be plotted.
#' @param title A character string specifying the title of the plot. 
#' Default is "LAI Distributions".
#' @param xlab A character string specifying the x-axis label. Default is "LAI".
#' @param var_labs A vector of character strings specifying the labels 
#' for the variables.
#' @param dirname A character string specifying the directory name 
#' where the plot will be saved.
#' @param filename A character string specifying the name of the file 
#' to save the plot.
#' @param limits A numeric vector of length 2 specifying the limits of 
#' the x-axis. Default is the range of the data.
#'
#' @return A histogram plot. Optionally saves the plot as a PNG file 
#' if `dirname` and `filename` are provided.
#'
#' @examples
#' vars_list <- list(rnorm(1000), rnorm(1000, mean = 3))
#' var_labs <- c("Group 1", "Group 2")
#' plot_histogram(vars_list, title = "My Histogram", var_labs = var_labs, 
#' dirname = "plots", filename = "histogram")
#'
#' @export
plot_histogram <- function(vars_list, 
                           title = NULL, 
                           xlab = NULL, 
                           var_labs, 
                           dirname = NULL, 
                           filename = NULL,
                           limits = NULL # c(0.1, 10)
) {
  
  if (length(vars_list) == 0) {
    cat("No variables provided. Exiting.\n")
    return(NULL)
  }
  if (is.null(title)) {
    title <- "LAI Distributions"
  }
  if (is.null(xlab)) {
    xlab <- "LAI"
  }
  if (is.matrix(vars_list) 
      || is.array(vars_list)
      || is.double(vars_list)) { # A single variable (i.e. not a list) was provided
    vars_list <- as.vector(vars_list)
    fill_colors = brewer.pal(3, "Set2")
    data <- data.frame(
      type = rep(as.character(var_labs), each = length(vars_list)),
      value = vars_list
    )
    if (is.null(limits)){
      limits = c(min(data$value, na.rm = TRUE), max(data$value, na.rm = TRUE))
    }
    
    h <- ggplot(data, aes(x = value, fill = type)) +
      geom_histogram(bins = 500, alpha = 0.4, position = "identity") +
      scale_fill_manual(values = fill_colors[1:length(vars_list)]) +
      labs(x = xlab, y = "Frequency", fill = "") +
      ggtitle(title) +
      scale_x_continuous(limits = limits)
  }
  if (is.list(vars_list)) { # A list is provided
    fill_colors = brewer.pal(length(vars_list), "Set2")
    data <- data.frame(
      type = rep(var_labs, each = length(vars_list[[1]])),
      value = do.call(c, vars_list) # do.call(c, vars_list) # unlist(vars_list)
    )
    if (is.null(limits)){
      limits = c(min(data$value, na.rm = TRUE), max(data$value, na.rm = TRUE))
    }
    
    h <- ggplot(data, aes(x = value, fill = type)) +
      geom_histogram(bins = 500, alpha = 0.6, position = "identity") +
      scale_fill_manual(values = fill_colors[1:length(vars_list)]) +
      labs(x = xlab, y = "Frequency", fill = "") +
      ggtitle(title) +
      scale_x_continuous(limits = limits)
  }
  print(h)
  
  if (!is.null(dirname) && !is.null(filename)) {
    file_path <- file.path(dirname, filename)
    ggsave(filename = paste0(file_path, ".png"),
           plot = h, device = "png", scale = 1, 
           width = 1920, height = 1080, units = "px", dpi = 200)
    cat("Plot is successfully saved at: ", paste0(file_path, ".png"), "\n")
  } else {
    cat("Plot is not saved because", "dirname =", dirname, 
        "and filename = ", filename, "\n")
  }
}

#' Correlation Test Function
#'
#' This function performs a correlation test between two numeric variables 
#' and prints the result.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param method A character string specifying the method for the correlation
#'  test. Default is "pearson". Other options include "kendall" and "spearman".
#'
#' @return The estimated correlation coefficient.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' correlation <- correlation_test_function(x, y)
#'
#' @export
correlation_test_function <- function(x, y, method = "pearson") {
  correlation_test <- cor.test(x, y, method = method)
  print(correlation_test)
  
  p_value <- correlation_test$p.value
  if (p_value < 0.05) {
    cat("Correlation is statistically significant (p < 0.05)\n")
  } else {
    cat("Correlation is not statistically significant (p >= 0.05)\n")
  }
  return(correlation_test$estimate)
}

#' Plot Distributions
#'
#' This function creates and optionally saves histograms 
#' for each variable in a data frame.
#'
#' @param InputPROSAIL A data frame containing the variables to be plotted.
#' @param save_plots A logical value indicating whether to save the plots. 
#' Default is \code{FALSE}.
#' @param dirname A character string specifying the directory name 
#' where the plots will be saved.
#' @param filename A character string specifying the name of the file 
#' to save the plots.
#'
#' @return A grid of histogram plots. Optionally saves the plots as a PNG file
#'  if \code{save_plots} is \code{TRUE} and 
#'  both \code{dirname} and \code{filename} are provided.
#'
#' @examples
#' InputPROSAIL <- data.frame(var1 = rnorm(1000), var2 = rnorm(1000, mean = 3))
#' plot_distributions(InputPROSAIL, save_plots = TRUE, dirname = "plots", 
#' filename = "distributions")
#'
#' @export
plot_distributions <- function(InputPROSAIL,
                               save_plots = FALSE,
                               dirname = NULL,
                               filename = NULL) {
  melted_df <- melt(InputPROSAIL)
  custom_green <- rgb(0, 200, 30, maxColorValue = 255)
  
  plot_grid <- lapply(unique(melted_df$variable), function(var) {
    ggplot(data = subset(melted_df, variable == var), aes(value)) +
      geom_histogram(bins = 30, fill = custom_green, color = "black") +
      ggtitle(paste("Histogram of", var)) +
      xlab(var) + ylab("Frequency")
  })
  
  final_plot <- plot_grid(plotlist = plot_grid, ncol = 4)
  
  if (save_plots) {
    if (!is.null(dirname) & !is.null(filename)) {
      if (!dir.exists(dirname)) {
        dir.create(dirname)
      }
      file_path <- file.path(dirname, filename)
      ggsave(filename = paste0(file_path, ".png"), final_plot, 
             width = 1920, height = 1080, units = "px", dpi = 100)
      cat("Plots saved at: ", paste0(file_path, ".png"), "\n")
    } else {
      warning("Please provide both directory name and file name for saving plots.")
    }
  } else {
    print(final_plot)
  }
}

# atbd_resample_and_mask <- function(biophysical_variable,
#                                    noise, 
#                                    band_res, 
#                                    input_base_path, 
#                                    output_directory, 
#                                    mask_raster # single_layer_forest
#                                    ) {
#   
#   input_path <- file.path(
#     input_base_path, 
#     sprintf('PRO4SAIL_INVERSION_atbd/%s/%s/L2A_T31UER_A031222_20210614T105443_Refl_%.envi', 
#             noise, band_res, biophysical_variable))
#   output_name <- sprintf('atbd_%s_%s', noise, band_res) # variable and value are modified separately
#   rast <- terra::rast(input_path)
#   
#   rast_resampled <- resample(rast, mask_raster)
#   rast_masked <- mask(rast_resampled, mask_raster)
#   save_envi_file(rast_masked, output_name, output_directory)
#   
#   rast_masked_values <- values(rast_masked)
#   rast_masked_values[rast_masked_values < 0] <- 0
#   
#   return(rast_masked_values)
# }
# 
# variable_resample_and_mask <- function(variable, 
#                                        value, 
#                                        input_base_path, 
#                                        output_directory, 
#                                        mask_raster # single_layer_forest
#                                        ) { 
#  
#   input_path <- file.path(input_base_path, sprintf('PRO4SAIL_INVERSION_atbd_%s_%s/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi', variable, value))
#   output_name <- sprintf('atbd_%s_%s', variable, value) # variable and value are modified separately
#   rast <- terra::rast(input_path)
#   
#   rast_resampled <- resample(rast, mask_raster)
#   rast_masked <- mask(rast_resampled, mask_raster)
#   save_envi_file(rast_masked, output_name, output_directory)
#   
#   rast_masked_values <- values(rast_masked)
#   rast_masked_values[rast_masked_values < 0] <- 0
#   
#   return(rast_masked_values)
# }
# 
# distribution_resample_and_mask <- function(distribution, 
#                                            input_base_path, 
#                                            output_directory, 
#                                            mask_raster # single_layer_forest
#                                            ) {
#   
#   input_path <- file.path(input_base_path, sprintf('PRO4SAIL_INVERSION_atbd_%s/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi', variable, value))
#   output_name <- sprintf('atbd_%s_%s', distribution)
#   rast <- terra::rast(input_path)
#   
#   rast_resampled <- resample(rast, mask_raster)
#   rast_masked <- mask(rast_resampled, mask_raster)
#   save_envi_file(rast_masked, output_name, output_directory)
#   
#   rast_masked_values <- values(rast_masked)
#   rast_masked_values[rast_masked_values < 0] <- 0
#   
#   return(rast_masked_values)
# }
# 
# generic_resample_and_mask <- function(input_path, 
#                                       output_name,
#                                       output_directory, 
#                                       site_mask # single_layer_forest
#                                       ) {
#   
#   rast <- open_envi_file(input_path)
#   
#   rast_resampled <- resample(rast, site_mask)
#   rast_masked <- mask(rast_resampled, site_mask)
#   save_envi_file(rast_masked, output_name, output_directory)
#   
#   rast_masked_values <- values(rast_masked)
#   rast_masked_values[rast_masked_values < 0] <- 0
#   
#   return(rast_masked_values)
# }
# 
# low_vegetation_resample_and_mask <- function(input_path, 
#                                              output_name,
#                                              output_directory, 
#                                              site_mask, # single_layer_forest
#                                              low_vegetation
# ) {
#   
#   rast <- open_envi_file(input_path)
#   
#   rast_resampled <- resample(rast, site_mask)
#   rast_masked <- mask(rast_resampled, site_mask)
#   
#   rast_masked_values <- values(rast_masked)
#   rast_masked_values[rast_masked_values < 0] <- 0
#   
#   masked <- mask_low_vegetation(rast_masked_values, rast_masked, low_vegetation)
#   masked_raster <- masked$masked_raster
#   save_envi_file(masked_raster, output_name, output_directory)
#   masked_values <- masked$masked_values
#   
#   return(masked_values)
# }
# 
# analyze_one_parameter <- function(parameter, 
#                                   output_directory,
#                                   plots_dir){
#   
#   # For a given parameter: atbd and param high/low/full in 10m/20m 
#   # -> 8 possibilities
#   
#   atbd_addmult_10m <- open_envi_file_as_values(
#     file.path(output_directory, "atbd",
#               paste0("atbd", "_", "addmult", "_", "10m", ".envi")))
#   
#   atbd_addmult_20m <- open_envi_file_as_values(
#     file.path(output_directory, "atbd",
#               paste0("atbd", "_", "addmult", "_", "20m", ".envi")))
#   
#   atbd_high_addmult_10m <- open_envi_file_as_values(
#     file.path(output_directory, paste0("atbd", "_", parameter),
#               paste0("atbd", "_", parameter, "_",
#                      "high", "_", "addmult", "_",
#                      "10m", ".envi")))
#   
#   atbd_high_addmult_20m <- open_envi_file_as_values(
#     file.path(output_directory, paste0("atbd", "_", parameter),
#               paste0("atbd", "_", parameter, "_",
#                      "high", "_", "addmult", "_",
#                      "20m", ".envi")))
#   
#   atbd_low_addmult_10m <- open_envi_file_as_values(
#     file.path(output_directory, paste0("atbd", "_", parameter),
#               paste0("atbd", "_", parameter, "_",
#                      "low", "_", "addmult", "_",
#                      "10m", ".envi")))
#   
#   atbd_low_addmult_20m <- open_envi_file_as_values(
#     file.path(output_directory, paste0("atbd", "_", parameter),
#               paste0("atbd", "_", parameter, "_",
#                      "low", "_", "addmult", "_",
#                      "20m", ".envi")))
#   atbd_full_addmult_10m <- open_envi_file_as_values(
#     file.path(output_directory, paste0("atbd", "_", parameter),
#               paste0("atbd", "_", parameter, "_",
#                      "full", "_", "addmult", "_",
#                      "10m", ".envi")))
#   
#   atbd_full_addmult_20m <- open_envi_file_as_values(
#     file.path(output_directory, paste0("atbd", "_", parameter),
#               paste0("atbd", "_", parameter, "_",
#                      "full", "_", "addmult", "_",
#                      "20m", ".envi")))
#   
#   param_10m_list <- list(atbd_addmult_10m,
#                          atbd_high_addmult_10m,
#                          atbd_low_addmult_10m,
#                          atbd_full_addmult_10m)
#   
#   param_10m_labs <- c("ATBD 10m",
#                       paste(parameter, "High", "10m", sep = " "),
#                       paste(parameter, "Low", "10m", sep = " "),
#                       paste(parameter, "Full", "10m", sep = " "))
#   
#   param_20m_list <- list(atbd_addmult_20m,
#                          atbd_high_addmult_20m,
#                          atbd_low_addmult_20m,
#                          atbd_full_addmult_20m)
#   
#   param_20m_labs <- c("ATBD 20m",
#                       paste(parameter, "High", "20m", sep = " "),
#                       paste(parameter, "Low", "20m", sep = " "),
#                       paste(parameter, "Full", "20m", sep = " "))     
#   
#   plot_histogram_4_vars(param_10m_list, "LAI Distributions", "LAI", 
#                         param_10m_labs, plots_dir, paste0(parameter, "10m", "_", "comp"))
#   
#   plot_histogram_4_vars(param_20m_list, "LAI Distributions", "LAI", 
#                         param_20m_labs, plots_dir, paste0(parameter, "20m", "_", "comp"))
# }