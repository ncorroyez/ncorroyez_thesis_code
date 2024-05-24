# ---
# title: "functions_lidar"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-01-31"
# ---

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("sf")

# 
# myPAI <- function(z, zmin, k=0.5){
#   Nout = sum(z<zmin)
#   Nin = length(z)
#   PAI = -log(Nout/Nin)/k
# }

mask_site_edges <- function(reflectance,
                            site_edges_path, # utm format
                            resolution,
                            masks_dir
){
  site_edges_mask <- terra::vect(site_edges_path)
  site_edges_mask <- mask(reflectance, site_edges_mask)
  
  viz <- site_edges_mask
  viz[is.na(viz)] <- 0
  viz[viz == 1] <- NA
  
  save_basic_plot(plot_to_save = site_edges_mask,
                  dirname = masks_dir,
                  filename = sprintf("site_edges_mask_res_%s_m.png", resolution),
                  title = sprintf("site_edges_mask_res_%s_m", resolution))
  plot(site_edges_mask, main=sprintf("site_edges_mask_res_%s_m", resolution))
  
  save_envi_file(viz, sprintf("site_edges_mask_res_%s_m_viz", resolution), masks_dir)
  save_envi_file(site_edges_mask, sprintf("site_edges_mask_res_%s_m", resolution), masks_dir)
  
  cat("Site edges mask has been successfully saved at :", 
      masks_dir, "\n")
  return(site_edges_mask)
}

mask_clouds <- function(reflectance,
                        cloud_path,
                        site_edges,
                        resolution,
                        masks_dir){
  cloud_mask <- terra::rast(cloud_path)

  
  cloud_mask <- terra::project(cloud_mask, 
                               site_edges)
  
  cloud_mask <- reflectance * cloud_mask
  
  save_basic_plot(plot_to_save = cloud_mask,
                  dirname = masks_dir,
                  filename = sprintf("cloud_mask_res_%s_m.png", resolution),
                  title = sprintf("cloud_mask_res_%s_m", resolution))
  plot(cloud_mask, main=sprintf("cloud_mask_res_%s_m", resolution))
  
  save_envi_file(cloud_mask, sprintf("cloud_mask_res_%s_m", resolution), masks_dir)
  
  cat("Clouds mask has been successfully saved at :", 
      masks_dir, "\n")
  return(cloud_mask)
}

create_mnc <- function(mnc,
                       masks_dir){
  
  # mnc <- dts_charged - dtm_charged
  save_basic_plot(plot_to_save = mnc,
                  dirname = masks_dir,
                  filename = "mnc.png",
                  title = "mnc")
  save_envi_file(mnc, "mnc", masks_dir)
  cat("CHM has been successfully saved at :", 
      masks_dir, "\n")
  return(mnc)
}

mask_mnc_threshold <- function(mnc,
                               threshold,
                               resolution,
                               masks_dir){
  mnc_thresholded <- mnc > threshold
  
  save_basic_plot(plot_to_save = mnc_thresholded,
                  dirname = masks_dir,
                  filename = sprintf("mnc_thresholded_%d_m_res_%s_m.png", 
                                     threshold, 
                                     resolution),
                  title = sprintf("mnc_thresholded_%d_m_res_%s_m", 
                                  threshold, 
                                  resolution))
  plot(mnc_thresholded, 
       main=sprintf("mnc_thresholded_%d_m_res_%s_m", 
                    threshold, 
                    resolution))
  
  save_envi_file(mnc_thresholded, 
                 sprintf("mnc_thresholded_%d_m_res_%s_m", 
                         threshold, 
                         resolution), 
                 masks_dir)
  cat(sprintf("MNC thresholded %d m has been successfully saved at :", 
              threshold), 
      masks_dir, "\n")
  return(mnc_thresholded)
}

project_mnc_threshold <- function(mnc_thresholded,
                                  raster_with_right_coordinates,
                                  resolution,
                                  masks_dir){
  average_mnc_thresholded <- terra::project(mnc_thresholded, 
                                            raster_with_right_coordinates, 
                                            method = 'average')
  save_basic_plot(plot_to_save = average_mnc_thresholded,
                  dirname = masks_dir,
                  filename = sprintf("mnc_1_m_thresholded_average_to_res_%s_m.png", resolution),
                  title = sprintf("mnc_1_m_thresholded_average_to_res_%s_m", resolution))
  plot(average_mnc_thresholded, 
       main=sprintf("mnc_1_m_thresholded_average_to_res_%s_m", resolution))
  
  save_envi_file(mnc_thresholded, sprintf("mnc_1_m_thresholded_average_to_res_%s_m", resolution), masks_dir)
  cat("MNC thresholded projected has been successfully saved at :", 
      masks_dir, "\n")
  return(average_mnc_thresholded)
}

mask_majority_project_mnc_threshold <- function(mnc,
                                                average_mnc_thresholded,
                                                percentage,
                                                resolution,
                                                masks_dir){
  
  # Majority
  average_mnc_thresholded_percentage <- average_mnc_thresholded > percentage
  save_basic_plot(plot_to_save = average_mnc_thresholded_percentage,
                  dirname = masks_dir,
                  filename = sprintf("average_mnc_thresholded_%s_p_kept_res_%s_m.png", 
                                     percentage*100, resolution),
                  title = sprintf("average_mnc_thresholded_%s_p_kept_res_%s_m", 
                                  percentage*100, resolution))
  plot(average_mnc_thresholded_percentage, 
       main = sprintf("average_mnc_thresholded_%s_p_kept_res_%s_m", 
                      percentage*100, resolution))
  
  save_envi_file(average_mnc_thresholded_percentage, 
                 sprintf("average_mnc_thresholded_%s_p_kept_res_%s_m", 
                         percentage*100, resolution), 
                 masks_dir)
  
  
  # Project MNC
  average_mnc <- terra::project(mnc, 
                                average_mnc_thresholded, 
                                method = 'average')
  mnc_val <- values(average_mnc)
  mnc_val[mnc_val < 0] <- 0
  
  index <- average_mnc_thresholded_percentage == 1
  mnc_final_mask <- average_mnc
  mnc_final_mask[!index] <- NA
  
  save_basic_plot(plot_to_save = mnc_final_mask,
                  dirname = masks_dir,
                  filename = sprintf("mnc_masked_with_average_mnc_thresholded_%s_p_kept_res_%s_m.png", 
                                     percentage*100, resolution),
                  title = sprintf("mnc_masked_with_average_mnc_thresholded_%s_p_kept_res_%s_m", 
                                  percentage*100, resolution))
  plot(mnc_final_mask,
       main = sprintf("mnc_masked_with_average_mnc_thresholded_%s_p_kept_res_%s_m", 
                      percentage*100, resolution))
  
  save_envi_file(mnc_final_mask, 
                 sprintf("mnc_masked_with_average_mnc_thresholded_%s_p_kept_res_%s_m", 
                         percentage*100, resolution), 
                 masks_dir)
  
  # Comparison to assess the mask
  mnc_final_mask_val <- values(mnc_final_mask)
  mnc_final_mask_val[mnc_final_mask_val < 0] <- 0
  
  mnc_list <- list(mnc_val, mnc_final_mask_val)
  labs <- c("mnc", "mnc_final_mask")
  plot_histogram(vars_list = mnc_list, 
                 title = sprintf("MNC vs MNC masked with threshold majorated %s p res %s m", percentage*100, resolution),
                 xlab = "heights",
                 var_labs = labs,
                 dirname = masks_dir,
                 filename = sprintf("hist_mnc_vs_mnc_masked_threshold_majorated_%s_p_kept_res_%s_m", 
                                    percentage*100, resolution))
  return(average_mnc_thresholded_percentage)
}

calculate_lad_profiles <- function(las_i, z0, zmax){
  
  las_i <- filter_poi(las_i, Z >= z0)
  Hmax_i <- max(las_i$Z)
  if (Hmax_i < z0) {
    # Set all values to 0 or NA or NULL
    rumple_i <- NA
    VCI_i <- NA
    CV_PAD_i <- NA
    # Hmax_i <- NA
    PAI_i <- NA
    max_PAD_i  <- NA
    H_max_PAD_i <- NA
    pad_i_translated <- NULL
    pad_i_translated$z <- 0
    pad_i_translated$lad <- NA
  } else {
    
    # CHM
    # chm_i <- grid_canopy(las_i, res = 1,
    #                      pitfree(thresholds = c(0, 2, 5, 10, 15),
    #                              max_edge = c(0, 1),
    #                              subcircle = 0.5))
    
    # DTM
    # dtm_i <- rasterize_terrain(las_i, res = 1, algorithm = tin())
    
    # DTS
    # dts_i <- rasterize_canopy(las_i, res = 1, p2r())
    
    # CHM
    # chm_i <- dts_i # - dtm_i
    
    # Calculations before translation
    # rumple_i <- rumple_index(chm_i)
    VCI_i <- VCI(las_i$Z, zmax = Hmax_i)
    
    # Translation
    las_i$Z <- las_i$Z + zmax - max(las_i$Z)
    pad_i_translated <- LAD(las_i@data$Z, z0=0)
    
    # Calculations after translation
    CV_PAD_i <- 100 * sd(pad_i_translated$lad, 
                         na.rm = TRUE) / mean(pad_i_translated$lad, 
                                              na.rm = TRUE) # 100*: Expressed in %
    
    PAI_i <- sum(pad_i_translated$lad, na.rm = TRUE)
    max_PAD_i <- max(pad_i_translated$lad, na.rm = TRUE)
    H_max_PAD_i <- pad_i_translated$z[which(pad_i_translated$lad==max_PAD_i)]
  }
  
  # Output
  metrics = list(
    # rumple_i,
    Hmax_i,
    VCI_i,
    CV_PAD_i, 
    PAI_i, 
    max_PAD_i, 
    H_max_PAD_i, 
    pad_i_translated)
  
  names(metrics)=c(
    # "Rumple_Index",
    "Hmax",
    "VCI",
    "CV_PAD",
    "PAI",
    "Max_PAD",
    "H_maxPAD",
    "PAD_Profile")
  
  return(metrics)
}

flatten_and_merge <- function(metrics_list, nb_z) {
  flattened_list <- lapply(seq_along(metrics_list), function(i) {
    las <- metrics_list[[i]]
    pad_profile <- as.data.frame(las$PAD_Profile)
    
    # Flatten the PAD_Profile data frame before combining
    pad_values <- unlist(pad_profile$lad)
    
    # Ensure that the number of columns in pad_profile matches nb_z
    if (length(pad_values) < nb_z) {
      pad_values <- c(pad_values, rep(0, nb_z - length(pad_values)))
    } else if (length(pad_values) > nb_z) {
      pad_values <- pad_values[1:nb_z]  # Trim excess values
    }
    
    las_num <- paste0("las_", i)
    
    return(c(las_num, 
             as.numeric(unlist(las[c(
               # "Rumple_Index", 
               "Hmax", 
               "VCI", 
               "CV_PAD", 
               "PAI", 
               "Max_PAD", 
               "H_maxPAD")], 
               use.names = FALSE)), 
             pad_values))
  })
  
  return(do.call(rbind, flattened_list))
}

calculate_lidar_metrics <- function(las_i, z0, zmax){
  
  # Rumple
  chm_i <- grid_canopy(las_i, res = 1, 
                       pitfree(thresholds = c(0, 2, 5, 10, 15), 
                               max_edge = c(0, 1), 
                               subcircle = 0.5))
  
  rumple_i <- rumple_index(chm_i)
  
  # Z
  las_i$Z <- las_i$Z + zmax - max(las_i$Z)
  
  lad_profiles <- list()
  for (z_start in seq(zmax - 1, by = -1)) {
    z_end <- zmax
    if (z_start >= z0) {
      # Calculate LAD profile for the current height stratum
      pad_i_norm <- LAD(las_i$Z)
      
      # Calculate the total LAD within the interval
      total_lad <- sum(pad_i_norm$lad[pad_i_norm$z >= z_end & pad_i_norm$z <= zmax])
      
      # pad_i_norm_sum <- data.frame(z = paste0(z_start, "-", z_end), 
      #                          lad = sum(pad_i_norm$lad, na.rm = TRUE))
      # pad_i_norm_sum <- sum(pad_i_norm$lad, na.rm = TRUE)
      
      # Store the LAD profile
      lad_profiles[[paste0(z_start, "_", z_end)]] <- pad_i_norm$lad
    }
  }
  
  las_subset <- filter_poi(las_i, Z >= z0)
  if (length(las_subset$Z) == 0) {   
    Hmax_i <- max(las_i$Z)
    VCI_i <-NA
    
    pad_i <-NULL
    pad_i$z <-0
    pad_i$lad <- NA 
    
    pad_i_norm <- data.frame(z = seq(z_start, z_end), 
                             lad = rep(0, z_end - z_start + 1)) 
    
    CV_PAD_i <- NA
    PAI_i <- NA
    max_PAD_i <- NA
    H_max_PAD_i <- NA
  } 
  else {
    
    Hmax_i <- max(las_i$Z)
    VCI_i <- VCI(las_subset$Z, by = 1, zmax= Hmax_i)
    
    step <- 1
    pad_i <-LAD(las_i$Z, dz=step, k=0.5, z0=z0)
    
    if (dim(pad_i)[1] == 0) {
      
      pad_i <- NULL
      nb_ind_z <- ceiling(Hmax_i - (z0 + 0.5*step )) 
      
      pad_i <- data.frame(matrix(NA,ncol= 2, nrow=nb_ind_z)) 
      colnames(pad_i) = c("z", "lad")
      pad_i$z <- seq(from = (z0+0.5*step),
                     to = (z0+ (nb_indice_z - 0.5)*step), by =step) 
      pad_i$lad <- rep(NA, nb_ind_z)
      
      CV_PAD_i <- NA
      PAI_i <- NA
      max_PAD_i <- NA
      H_max_PAD_i <- NA
      
    } else {
      
      CV_PAD_i <- 100 * sd(pad_i$lad, 
                           na.rm = TRUE) / mean(pad_i$lad, 
                                                na.rm = TRUE) 
      PAI_i <- sum(pad_i$lad, na.rm = TRUE)  
      
      max_PAD_i <- max(pad_i$lad, na.rm = TRUE)
      H_max_PAD_i <- pad_i$z[which(pad_i$lad==max_PAD_i)]
      
      # Norm
      CV_PAD_norm_i <- 100 * sd(pad_i$lad, 
                                na.rm = TRUE) / mean(pad_i$lad, 
                                                     na.rm = TRUE) 
      PAI_norm_i <- sum(pad_i$lad, na.rm = TRUE)  
      
      max_PAD_norm_i <- max(pad_i$lad, na.rm = TRUE)
      H_max_PAD_norm_i <- pad_i$z[which(pad_i$lad==max_PAD_i)]
      
    }
  }
  
  # # gap fraction raster ave seuil de hauteur pour definir les trouees
  # # raster moins lisse que pour rumple, restitue mieux les trouees 
  # chm_tr_i <- grid_canopy(las_i_h, res=0.5, p2r() )    
  # 
  # gapf_raster_i <- 100 * (length(
  #   chm_tr_i@data@values[!is.na(chm_tr_i@data@values) 
  #                        & chm_tr_i@data@values < seuil_gap ]))
  # / ( length( chm_tr_i@data@values[!is.na(chm_tr_i@data@values)] ) )
  # 
  # # gap fraction issu du ratio des premiers retours parvenant sous le seuil 
  # # des hauteurs par rapport au nbre total de premier retours
  # 
  # gapf_pts_i <- 100 * dim(filter_poi(
  #   las_i_h, ReturnNumber ==1L & Z <seuil_gap)@data)[1] /  
  #   dim(filter_first(las_i_h)@data)[1]
  
  
  metrics = list(rumple_i,
                 VCI_i, 
                 Hmax_i,
                 CV_PAD_i, 
                 PAI_i, 
                 max_PAD_i, 
                 H_max_PAD_i, 
                 # gapf_raster_i, 
                 # gapf_pts_i,
                 pad_i,
                 pad_i_norm_sum
  )
  names(metrics)=c("Rumple_index",
                   "VCI",
                   "Hmax", 
                   "CV_PAD",
                   "PAI",
                   "Max_PAD",
                   "H_maxPAD",
                   # "Gap_fraction_raster",
                   # "Gap_fraction_pts",
                   "Profil_pad",
                   "Profil_pad_norm_from_top")
  
  return(metrics)
  
}

cv <- function(x) {
  sd(x) / mean(x)
}
