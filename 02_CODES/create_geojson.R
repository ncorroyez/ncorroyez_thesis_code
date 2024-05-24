# ---
# title: "create_geojson.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-12-15"
# ---

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

library(sf)
library(tidyr)

# --------------------------------------------- data_plot.csv --------------------------------------------------------

data <- read.csv("../../MaCCMic_Imprint/Data_placettes_Aigoual/data_plots_202402281200.csv", sep = ",")

data$coord_x_l93 <- as.numeric(gsub(",", ".", gsub("�", "", data$coord_x_l93)))
data$coord_y_l93 <- as.numeric(gsub(",", ".", gsub("�", "", data$coord_y_l93)))
data$long_wgs84 <- as.numeric(gsub(",", ".", data$long_wgs84))
data$lat_wgs84 <- as.numeric(gsub(",", ".", data$lat_wgs84))
data$elevation <- as.numeric(gsub(",", ".", data$elevation))

data_blois <- data[data$forest == "Blois", ]
data_mormal <- data[data$forest == "Mormal", ]
data_aigoual <- data[data$forest == "Aigoual", ]

coordinates_blois <- st_as_sf(data_blois, 
                              coords = c("long_wgs84", "lat_wgs84"), 
                              crs = 4326)
coordinates_mormal <- st_as_sf(data_mormal, 
                               coords = c("long_wgs84", "lat_wgs84"), 
                               crs = 4326) 
coordinates_aigoual <- st_as_sf(data_aigoual, 
                                coords = c("long_wgs84", "lat_wgs84"), 
                                crs = 4326) 

st_write(coordinates_blois, "data_blois.geojson", driver = "GeoJSON")
st_write(coordinates_mormal, "data_mormal.geojson", driver = "GeoJSON")
st_write(coordinates_aigoual, "data_aigoual.geojson", driver = "GeoJSON")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------------- data_lidar_pad.csv -----------------------------------------------------

# data <- read.csv("../../MaCCMic_Imprint/MaccMIC_dec2023/data_lidar_pad.csv", sep = "\t")
# 
# data$coord_x_l93 <- as.numeric(gsub(",", ".", gsub("�", "", data$coord_x_l93)))
# data$coord_y_l93 <- as.numeric(gsub(",", ".", gsub("�", "", data$coord_y_l93)))
# data$long_wgs84 <- as.numeric(gsub(",", ".", data$long_wgs84))
# data$lat_wgs84 <- as.numeric(gsub(",", ".", data$lat_wgs84))
# data$elevation <- as.numeric(gsub(",", ".", data$elevation))
# 
# data_blois <- data[data$forest == "Blois", ]
# data_mormal <- data[data$forest == "Mormal", ]
# 
# coordinates_blois <- st_as_sf(data_blois, 
#                               coords = c("long_wgs84", "lat_wgs84"), 
#                               crs = 4326)
# coordinates_mormal <- st_as_sf(data_mormal, 
#                                coords = c("long_wgs84", "lat_wgs84"), 
#                                crs = 4326) 
# 
# st_write(coordinates_blois, "data_blois.geojson", driver = "GeoJSON")
# st_write(coordinates_mormal, "data_mormal.geojson", driver = "GeoJSON")
