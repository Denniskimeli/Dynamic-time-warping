library(raster)
library(dtwSat) 
library(rgdal)
library(caret)
library(sp)
library(sf)
library(zoo)
library(rgeos)
library(coin)
library(rlang)
library(terra)
library(devtools)

install.packages('dtwSat')
install.packages('snow')
install.packages('htmltools')

# load AOI
AOI <- readOGR("C:/Msc/research progress/Data collection/ADC data/AOI/min_box.shp")
AOI
AOI <-spTransform(AOI,proj_str)
AOI

# function to calculate mean 
# calculate_mean_of_consecutive_bands <- function(x) {
#   num_layers <- nlayers(x)
#   aggregated_stack <- stack()
#   for (i in seq(1, num_layers - 2, by = 3)) {
#     aggregated_mean <- (x[[i]] + x[[i + 1]] + x[[i + 2]]) / 3
#     aggregated_stack <- addLayer(aggregated_stack, aggregated_mean)
#   }
#   names(aggregated_stack) <- names(x)[seq(3, num_layers, by = 3)]
#   return(aggregated_stack)
# }

calculate_mean_of_consecutive_bands <- function(x) {
  num_layers <- nlayers(x)
  aggregated_stack <- stack()
  for (i in seq(1, num_layers - 1, by = 2)) {
    aggregated_mean <- (x[[i]] + x[[i + 1]]) / 2
    aggregated_stack <- addLayer(aggregated_stack, aggregated_mean)
  }
  names(aggregated_stack) <- names(x)[seq(2, num_layers, by = 2)]
  return(aggregated_stack)
}

timeline <- scan(file="D:/MSC_data/CROP MAPPING/shapefiles/dates_VV_VH.txt",what='date')


# function to remove NA
remove_na_bands <- function(raster_stack) {
  num_layers <- nlayers(raster_stack)
  # Calculate the fraction of NA values in each band
  fraction_na <- cellStats(is.na(raster_stack), sum) / ncell(raster_stack)
  # Identify the bands with less than 50% NA values
  valid_bands <- which(fraction_na < 0.5)
  # Extract and return the valid bands in a new raster stack
  return(raster_stack[[valid_bands]])
}

proj_str = scan(file="D:/Compressed/Landasat TZ/timeline/samples_projection.txt", what = "character") 
# proj_str<- "+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs"
beginCluster()
#  load evi 2022
evi_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/EVI/0/stack.vrt")
# ev_test <- evi_2022[[1]]
# plot(ev_test)
evi_no_na <- remove_na_bands(evi_2022)
evi_no_na
# names_with_na <- names(evi_2022)
names_withouth_na <- names(evi_no_na)
names_withouth_na
aggregated_evi <- calculate_mean_of_consecutive_bands(evi_no_na)
aggregated_evi
# evi_2022_norm <- aggregated_evi/10000
evi_brick <- brick(aggregated_evi)
evi_proj <- projectRaster(evi_brick, crs = proj_str)
# evi_proj <- reclassify(evi_proj, cbind(NA, 0))
evi_proj <- evi_proj[[1:17]]
plot(evi_proj)
# load swir1
swir_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/SWIR1/0/stack.vrt")
swir_sel <- swir_2022[[names_withouth_na]]
swir_sel
aggregated_swir <- calculate_mean_of_consecutive_bands(swir_sel)
# swir_2022_norm <- aggregated_swir/10000
swir_brick <- brick(aggregated_swir)
swir_proj <- projectRaster(swir_brick, crs = proj_str)
swir_proj <- swir_proj[[1:17]]
# swir_proj <- reclassify(swir_proj, cbind(NA, 0))
plot(swir_proj)
# loadd NDVI 
ndvi_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/NDVI/0/stack.vrt")
ndvi_sel <- ndvi_2022[[names_withouth_na]]
aggregated_ndvi <- calculate_mean_of_consecutive_bands(ndvi_sel)
# ndvi_2022_norm <- aggregated_ndvi/10000
ndvi_brick <- brick(aggregated_ndvi)
ndvi_proj <- projectRaster(ndvi_brick, crs = proj_str)
ndvi_proj <- ndvi_proj[[1:17]]
# ndvi_proj <- reclassify(ndvi_proj, cbind(NA, 0))
# load swir2
swir2_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/SWIR2/0/stack.vrt")
swir2_sel <- swir2_2022[[names_withouth_na]]
aggregated_swir2 <- calculate_mean_of_consecutive_bands(swir2_sel)
# swir2_2022_norm <- aggregated_swir2/10000
swir2_brick <- brick(aggregated_swir2)
swir2_proj <- projectRaster(swir2_brick, crs = proj_str)
swir2_proj <- swir2_proj[[1:17]]
# swir2_proj <- reclassify(swir2_proj, cbind(NA, 0))

# loadd blue 
blue_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/BLUE/0/stack.vrt")
blue_sel <- blue_2022[[names_withouth_na]]
aggregated_blue <- calculate_mean_of_consecutive_bands(blue_sel)
# blue_2022_norm <- aggregated_blue/10000
blue_brick <- brick(aggregated_blue)
blue_proj <- projectRaster(blue_brick, crs = proj_str)
blue_proj <- blue_proj[[1:17]]
# blue_proj <- reclassify(blue_proj, cbind(NA, 0))
# load green
green_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/GREEN/0/stack.vrt")
green_sel <- green_2022[[names_withouth_na]]
aggregated_green <- calculate_mean_of_consecutive_bands(green_sel)
# green_2022_norm <- aggregated_green/10000
green_brick <- brick(aggregated_green)
green_proj <- projectRaster(green_brick, crs = proj_str)
green_proj <- green_proj[[1:17]]
# green_proj <- reclassify(green_proj, cbind(NA, 0))
# load red
red_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/RED/0/stack.vrt")
red_sel <- red_2022[[names_withouth_na]]
aggregated_red <- calculate_mean_of_consecutive_bands(red_sel)
# red_2022_norm <- aggregated_red/10000
red_brick <- brick(aggregated_red)
red_proj <- projectRaster(red_brick, crs = proj_str)
red_proj <- red_proj[[1:17]]
# red_proj <- reclassify(red_proj, cbind(NA, 0))
# load NIR
nir_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/NIR/0/stack.vrt")
nir_sel <- nir_2022[[names_withouth_na]]
aggregated_nir <- calculate_mean_of_consecutive_bands(nir_sel)
# nir_2022_norm <- aggregated_nir/10000
nir_brick <- brick(aggregated_nir)
nir_proj <- projectRaster(nir_brick, crs = proj_str)
nir_proj <- nir_proj[[1:17]]
# nir_proj <- reclassify(nir_proj, cbind(NA, 0))

# load NDMI
ndmi_2022 <- raster::stack("D:/MSC_data/CROP MAPPING/Sepal/2022/NDMI/0/stack.vrt")
ndmi_sel <- ndmi_2022[[names_withouth_na]]
aggregated_ndmi <- calculate_mean_of_consecutive_bands(ndmi_sel)
# ndmi_2022_norm <- aggregated_ndmi/10000
ndmi_brick <- brick(aggregated_ndmi)
ndmi_proj <- projectRaster(ndmi_brick, crs = proj_str)
ndmi_proj <- ndmi_proj[[1:17]]
# ndmi_proj <-reclassify(ndmi_proj, cbind(NA, 0))
plot(ndmi_proj)
ndmi_proj[[1:17]]
# Replace NA values with zeros in the raster stack
# aggregated_ndmi_no_na <- reclassify(ndmi_proj, cbind(NA, 0))

tbs <- gsub("^X", "", names(red_proj)) # Remove the leading "X"
timeline <- as.Date(tbs, format = "%Y.%m.%d")
timeline
samples_22 <- read.csv("D:/MSC_data/CROP MAPPING/DTW/Samples_22.csv")
samples_22 <- subset(samples_22, label != 'Barley')

# sentinel one
files_s1_22 <- list.files("D:/MSC_data/CROP MAPPING/Sentinel1/Via_Drive/staceked", pattern = "\\.tif$", full.names = T)

vvvh_2022 <- lapply(files_s1_22,raster::stack)
vvvh_2022

# clipped_stacks <- lapply(vvvh_2022, function(stack) {
#   clipped_layers <- lapply(1:nlayers(stack), function(i) {
#     mask(stack[[i]], AOI)
#   })
#   stack(clipped_layers)
# })
# 
# clipped_stacks

vvvh_resampled_2022 <- lapply(vvvh_2022, function(x) resample(x, vvvh_2022[[11]], method = "bilinear"))
vvvh_resampled_2022

extract_vv_layer <- function(stack) {
  return(stack$VV)  # Access the VV layer using $VV
}

extract_vh_layer <- function(stack) {
  return(stack$VH)  # Access the VV layer using $VV
}


vv_unsampled <- lapply(vvvh_resampled_2022, extract_vv_layer)
vv_unsampled
vv_unsampled <- crop(vv_unsampled,AOI)

vv_layers <- lapply(vvvh_resampled_2022, extract_vv_layer)
vv_2022 <- crop(stack(vv_layers),AOI)

# crop_to_aoi <- function(layer) {
#   cropped_layer <- mask(layer, AOI)
#   return(cropped_layer)
# }

# Apply the cropping function to each extracted layer in vv_layers
# cropped_vv_layers <- lapply(vv_layers, crop_to_aoi)

vh_layers <- lapply(vvvh_resampled_2022, extract_vh_layer)
vh_2022 <- crop(stack(vh_layers),AOI)
vh_2022

vv_2022_df <- as.data.frame(vv_2022) 
print(vh_2022_df[1:10, ])
na_matrix <- is.na(vv_2022_df)
na_count <- colSums(na_matrix)
# Print the count of NA values in each column
print(na_count)


evi_proj_sample <- resample(evi_proj, vv_2022, method = "bilinear")
extent(evi_proj_sample)

plot(vh_2022)
na_count <- cellStats(is.na(vh_2022_df), sum)
na_count
rts2022_optical <- twdtwRaster(ndmi_proj,swir_proj,swir2_proj,blue_proj,green_proj,red_proj,ndvi_proj,evi_proj,nir_proj,vh_2022,vv_2022, timeline= timeline)
rts2022_radar <- twdtwRaster(vv_2022,vh_2022,evi_proj_sample, timeline = timeline)
# rts2022_optical@layers

training_ts_radar <- getTimeSeries(rts2022_radar,y=samples_22,proj4string = proj_str)

training_ts_radar@timeseries
training_ts_optical <- getTimeSeries(rts2022_optical, y = samples_22, proj4string = proj_str)
temporal_patterns_S1 <- createPatterns(training_ts_radar, freq = 8, formula = y ~ s(x))

plot(temporal_patterns, type = "patterns")

log_funO <- logisticWeight(alpha = -0.1, beta = 50)

# Get sample time series from raster time series 

field_samples_ts
# field_samples_filtered <- subset(samples_22, label != 'Wheat')
field_samples_ts = getTimeSeries(rts2022_radar,
                                 y = field_samples_filtered, proj4string = proj_str)
# Run cross validation
set.seed(1)
# Define TWDTW weight function
log_fun = logisticWeight(alpha=-0.1, beta=50)
cross_validation = twdtwCrossValidate(field_samples_ts, times=3, p=0.1,
                                      freq = 8, formula = y ~ s(x, bs="cc"), weight.fun = log_fun)
cross_validation

summary(cross_validation)

plot(cross_validation)

r_dist_2022 <-twdtwApply(x = rts2022_radar, y = temporal_patterns, weight.fun = log_fun, progress = 'text')
crops2022 <- twdtwClassify(r_dist_2022, format="GTiff", overwrite=TRUE)
plot(x = crops2022, type = "maps")

plot(AOI)


