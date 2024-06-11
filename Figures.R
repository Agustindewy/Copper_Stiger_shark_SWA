
# Figures

library(viridis)
library(ggplot2)
library(rworldmap)
library(maps)
library(ggmap)
library(ggthemes)
library(sf)
library(rgdal)
library(dplyr)
library(rgl)
library(raster)
library(rgeos)
library(ggsn)
library(marmap)
library(SDMtune)

setwd('SET YOUR WORKING DIRECTORY')

# Coastline (spatial polygons) - taken from GSHHG coast database as of June 15, 2017
coastline <- 'Coast shape'
coast0 <- readOGR(dsn = coastline, layer = 'GSHHS_f_L1_SouthAmerica')
coast <- readOGR(dsn = coastline, layer = 'GSHHS_f_L1_World')

# Political borders
Borders <- 'Political borders'
Border <- readOGR(dsn = Borders, layer = 'WDBII_border_f_L1')
Border <- st_read(dsn = Borders, layer = 'WDBII_border_f_L1')

# Lat-lon projection
crs <- CRS('+proj=longlat +datum=WGS84')

# Functions to compute 5% and 10% Minimum Training Presence (taken from https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/)
sdm_threshold.5 <- function(sdm, occs, type = 'mtp', binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == 'mtp'){
    thresh <- min(na.omit(occPredVals))
  } else if(type == 'p05'){
    if(length(occPredVals) < 10){
      p05 <- floor(length(occPredVals) * 0.95)
    } else {
      p05 <- ceiling(length(occPredVals) * 0.95)
    }
    thresh <- rev(sort(occPredVals))[p05]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
} # MTP5
sdm_threshold.10 <- function(sdm, occs, type = 'mtp', binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == 'mtp'){
    thresh <- min(na.omit(occPredVals))
  } else if(type == 'p10'){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
} # MTP10

# Functions for response curves
.get_presence <- function(swd) {
  return(swd@data[swd@pa == 1, , drop = FALSE])
}
.get_absence <- function(swd) {
  return(swd@data[swd@pa == 0, , drop = FALSE])
}


#------------------------------------ Figure 1a ---------------------------------------

# Continuous habitat suitability globally at annual scale

pred <- raster('Carcharhinus brachyurus/Predictions.tif')

# Strict extrapolation areas
global_mop <- raster('Carcharhinus brachyurus/Mop.tif') 
global_mop[global_mop > 0] <- NA
global_mop[global_mop == 0] <- 2
pred <- crop(pred, global_mop)

# Remove MOP areas from prediction
pred_mop <- mosaic(pred, global_mop, fun = max) 
pred_mop[pred_mop > 1] <- NA

df <- data.frame(coordinates(pred_mop), as.data.frame(pred_mop))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = layer)) + 
  geom_sf(data = coast, col = 'grey10', fill = 'grey20', linewidth = 0.1) +
  scale_fill_viridis(option = 'G', na.value = '#0B0405FF') + coord_sf(expand = F) +
  scale_y_continuous(name = NULL, breaks = c(-50, 0, 50), labels = c('50ºS', '0º', '50ºN')) + 
  scale_x_continuous(name = NULL, breaks = c(-100, 0, 100), labels = c('100ºW', '0º', '100ºE')) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(),
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5))
ggsave('Figure 1a.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm', device = grDevices::tiff)


#------------------------------------ Figure 1b ---------------------------------------

# Calibration area + calibration points
library(stars)
occ_cal <- read.csv('Carcharhinus brachyurus/Calibration_points.csv') 
env.M <- read_stars('Carcharhinus brachyurus/Calibration_areas.tif')
env_buff <- env.M > -Inf
env_buff <- st_as_sf(env_buff, as_points = F, merge = T)

ggplot() +
  geom_sf(data = env_buff, col = 'black', fill = '#357BA2FF', linewidth = 0.5) +
  geom_sf(data = coast, col = 'grey30', fill = 'grey50', linewidth = 0.1) +
  geom_point(data = occ_cal, aes(x = Longitude, y = Latitude), shape = 21, size = 0.5, color = '#A0DFB9FF', fill = '#0B0405FF') +
  scale_fill_manual(values = '#357BA2FF', na.value = 'grey95') + coord_equal(expand = 0) +
  scale_y_continuous(name = NULL, breaks = c(-50, 0, 50), labels = c('50ºS', '0º', '50ºN')) + 
  scale_x_continuous(name = NULL, breaks = c(-100, 0, 100), labels = c('100ºW', '0º', '100ºE')) +
  coord_sf(xlim = c(-180, 180), ylim = c(-68.9, 83.6), expand = F) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(),
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5))
ggsave('Figure 1b.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm', device = grDevices::tiff)


#------------------------------------ Figure 1c & S1.5 ---------------------------------------

# Annual model response plots

# Feature class and regularization multiplier selected
fc <- 'lq' 
rm <- 0.1

occ_cal <- read.csv('Carcharhinus brachyurus/Calibration_points.csv') 
occ_cal <- occ_cal[, c('Longitude', 'Latitude')]
env.M <- stack('Carcharhinus brachyurus/Calibration_areas.tif') 
var_set <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity', 'Distance_to_coast', 'Thermal_fronts')
names(env.M) <- var_set

p_df = data.frame()
a_df = data.frame()
plot_data_df = data.frame()

source('prepareSWD_adw.R') # Function from the 'SDMtune' package with a bug fixed

for(j in names(env.M)) {
  
  # Prepare data
  env <- env.M[[j]]
  env <- stack(env)
  Vars <- j
  
  # Background points
  set.seed(111)
  notna <- which(complete.cases(values(env)))
  samp <- sample(notna, 10000, replace = F)
  samplocs <- as.data.frame(xyFromCell(env, samp))
  
  # SWD object
  data <- prepareSWD_adw(species = 'Carcharhinus brachyurus', p = occ_cal, a = samplocs, env = env)
  
  # Run maxent replicates
  folds = randomFolds(data, k = 10, only_presence = T, seed = 111)
  default_model <- train(method = 'Maxent', data = data, fc = fc, reg = rm, iter = 1000, folds = folds)
  
  # Presences and pseudo-absences with environmental data
  p <- .get_presence(default_model@data)
  p$var <- j
  names(p)[names(p) == j] <- 'values'
  a <- .get_absence(default_model@data)
  a$var <- j
  names(a)[names(a) == j] <- 'values'
  
  pred <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 10))
  for(i in 1:10){
    pred[, i] <- predict(default_model@models[[i]], data = data@data, type = 'logistic')
  }
  
  # Plot data
  plot_data <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 4))
  names(plot_data) <- c('mean', 'sd', 'max', 'min')
  plot_data$mean <- rowMeans(pred)
  plot_data$sd <- apply(pred, 1, sd)
  plot_data$max <- plot_data$mean + plot_data$sd
  plot_data$min <- plot_data$mean - plot_data$sd
  plot_data$var <- j
  plot_data$values <- data@data[, j]
  
  # Data frames
  p_df <- rbind(p_df, p)
  a_df <- rbind(a_df, a)
  plot_data_df <- rbind(plot_data_df, plot_data)
  
}

# Re-order
order_i_want1 <- c('Distance_to_coast', 'Temperature', 'Surface_temperature', 'Salinity', 'Kd490', 'Primary_productivity', 'Thermal_fronts')
p_df <- p_df[!(p_df$var == 'Salinity' & p_df$values < 20), ]
p_df <- p_df[!(p_df$var == 'Distance_to_coast' & p_df$values > 750), ]
p_df <- p_df[!(p_df$var == 'Kd490' & p_df$values > 0.5), ]
p_df <- p_df[!(p_df$var == 'Primary_productivity' & p_df$values > 0.07), ]
p_df <- p_df[!(p_df$var == 'Thermal_fronts' & p_df$values > 0.6), ]
p_df <- transform(p_df, var = factor(var, levels = order_i_want1))

a_df <- a_df[!(a_df$var == 'Salinity' & a_df$values < 20), ]
a_df <- a_df[!(a_df$var == 'Distance_to_coast' & a_df$values > 750), ]
a_df <- a_df[!(a_df$var == 'Kd490' & a_df$values > 0.5), ]
a_df <- a_df[!(a_df$var == 'Primary_productivity' & a_df$values > 0.07), ]
a_df <- a_df[!(a_df$var == 'Thermal_fronts' & a_df$values > 0.6), ]
a_df <- transform(a_df, var = factor(var, levels = order_i_want1))

plot_data_df <- plot_data_df[!(plot_data_df$var == 'Salinity' & plot_data_df$values < 20), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Distance_to_coast' & plot_data_df$values > 750), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Kd490' & plot_data_df$values > 0.5), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Primary_productivity' & plot_data_df$values > 0.07), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Thermal_fronts' & plot_data_df$values > 0.6), ]
plot_data_df <- transform(plot_data_df, var = factor(var, levels = order_i_want1))

# Labeller
VAR_names = as_labeller(c(Temperature = 'Bottom~temperature~(ºC)', Distance_to_coast = 'Distance~to~coast~(km)',
                          Surface_temperature = 'Surface~temperature~(ºC)', Kd490 = 'Kd490~coefficient~(m^-1)', 
                          Salinity = 'Bottom~salinity~(psu)', Primary_productivity = 'Primary~prod.~(gCm^-2~día^-1)',
                          Thermal_fronts = 'Thermal~fronts~(ºC)'), default = label_parsed)

# Plot
ggplot(data = plot_data_df, aes(x = values, y = mean, ymin = min, ymax = max)) + 
  geom_line(colour = '#0B0405FF') + 
  geom_ribbon(fill = '#0B0405FF', alpha = 0.2) +
  geom_rug(data = p_df, inherit.aes = F, aes(values), sides = 't', color = '#A0DFB9FF', linewidth = 0.3) + 
  geom_rug(data = a_df, inherit.aes = F, aes(values), sides = 'b', color = '#357BA2FF', linewidth = 0.3) + 
  labs(x = NULL, y = 'Logistic output') + ylim(0, 1) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.2),
        panel.grid.major = element_line(linewidth = 0.2, colour = 'grey90'),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(vjust = -0.5, size = 9),
        plot.margin = unit(c(0, 0, 0.2, 0.2), 'cm'),
        axis.title = element_text(size = 10)) +
  facet_wrap(~var, scales = 'free_x', labeller = VAR_names, nrow = 3) 
ggsave('Figure 1c.pdf', width = 15, height = 15, units = 'cm')


#------------------------------------ Figure 2 & Figure 5a ---------------------------------------

# Annual continuous habitat suitability at a population scale in the Southwest Atlantic

pred <- raster('Carcharhinus brachyurus/Predictions.tif')
swa_ext <- extent(-68, -39, -48, -19)
pred <- crop(pred, swa_ext)

occ <- read.csv('Free_occurrence_records.csv') 
occ <- subset(occ, Species == 'Carcharhinus brachyurus')
occ <- subset(occ, Region == 'Southwest Atlantic')

df <- data.frame(coordinates(pred), as.data.frame(pred))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = layer)) + 
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), color = 'grey10', fill = 'grey50', linewidth = 0.15) +
  geom_point(data = occ, aes(x = Longitude, y = Latitude), shape = 21, size = 0.9, stroke = 0.2, col = 'darkred', fill = 'darkorange2') + 
  scale_fill_viridis(option = 'G', na.value = '#0B0405FF') + 
  scale_y_continuous(name = NULL, breaks = c(-44, -34, -24), labels = c('44º', '34º', '24º')) +
  scale_x_continuous(name = NULL, breaks = c(-66, -58, -50, -42), labels = c('66º', '58º', '50º', '42º')) +
  coord_equal(xlim = c(-68, -39), ylim = c(-48, -19), expand = F) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(),
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5))
ggsave('Figure 2a.tiff', dpi = 900, width = 10, height = 10, units = 'cm', device = grDevices::tiff)


# Annual binary habitat suitability at a population scale in the Southwest Atlantic

pred <- raster('Carcharhinus brachyurus/Predictions.tif')
bg_binary <- pred
bg_binary[bg_binary >= 0] <- 1

occ_cal <- read.csv('Carcharhinus brachyurus/Calibration_points.csv') 

# Threshold (5% minimum training presence)
mtp5 <- sdm_threshold.5(pred, occ_cal[, c('Longitude', 'Latitude')], 'p05', binary = F)
bin_5 <- pred
bin_5[bin_5 < minValue(mtp5)] <- NA
bin_5[bin_5 >= minValue(mtp5)] <- 2

# Threshold (10% minimum training presence)
mtp10 <- sdm_threshold.10(pred, occ_cal[, c('Longitude', 'Latitude')], 'p10', binary = F)
bin_10 <- pred
bin_10[bin_10 < minValue(mtp10)] <- NA
bin_10[bin_10 >= minValue(mtp10)] <- 3

# Final mosaic
mod_mosaic <- mosaic(bg_binary, bin_10, bin_5, fun = sum)
mod_mosaic[mod_mosaic == 1] <- 1 # Unsuitable areas
mod_mosaic[mod_mosaic == 3] <- 2 # MTP5
mod_mosaic[mod_mosaic == 6] <- 3 # MTP overlap
mod_mosaic <- crop(mod_mosaic, swa_ext)

Col <- c('grey95', '#0B0405FF', '#357BA2FF', '#A0DFB9FF', 'grey90')

df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), color = 'grey10', fill = 'white', linewidth = 0.15) +
  geom_path(data = Border, aes(x = long, y = lat, group = group), color = 'grey20', linewidth = 0.2) +
  scale_y_continuous(name = NULL, breaks = c(-44, -34, -24), labels = c('44º', '34º', '24º')) +
  scale_x_continuous(name = NULL, breaks = c(-66, -58, -50, -42), labels = c('66º', '58º', '50º', '42º')) +
  coord_equal(xlim = c(-68, -39), ylim = c(-48, -19), expand = F) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5)) 
ggsave('Figure 2b.tiff', dpi = 900, width = 10, height = 10, units = 'cm', device = grDevices::tiff)


# Southern limit

# Provincial borders
provincias <- readOGR(dsn = Borders, layer = 'provincias_argentina')

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), color = 'grey10', fill = 'white', linewidth = 0.15) +
  geom_path(data = provincias, aes(x = long, y = lat, group = group), color = 'grey70', linewidth = 0.15, linetype = "dashed") +
  scale_y_continuous(name = NULL, breaks = c(-46, -42, -38), labels = c('46º', '42º', '38º')) +
  scale_x_continuous(name = NULL, breaks = c(-66, -62, -58), labels = c('66º', '62º', '58º')) +
  coord_equal(xlim = c(-68, -56), ylim = c(-48, -36), expand = F) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5)) 
ggsave('Figure 5a.tiff', dpi = 900, width = 10, height = 10, units = 'cm', device = grDevices::tiff)


#------------------------------------ Figure 3a ---------------------------------------

# Annual continuous habitat suitability at a global scale

pred <- raster('Carcharias taurus/Predictions.tif')

# Strict extrapolation areas 
global_mop <- raster('Carcharias taurus/Mop.tif') 
global_mop[global_mop > 0] <- NA
global_mop[global_mop == 0] <- 2
pred <- crop(pred, global_mop)

# Remove MOP areas from the prediction
pred_mop <- mosaic(pred, global_mop, fun = max) 
pred_mop[pred_mop > 1] <- NA

df <- data.frame(coordinates(pred_mop), as.data.frame(pred_mop))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = layer)) + 
  geom_sf(data = coast, col = 'grey10', fill = 'grey20', linewidth = 0.1) +
  scale_fill_gradientn(colors = pals::parula(1000), na.value = '#352A87') + coord_sf(expand = F) +
  scale_y_continuous(name = NULL, breaks = c(-50, 0, 50), labels = c('50ºS', '0º', '50ºN')) + 
  scale_x_continuous(name = NULL, breaks = c(-100, 0, 100), labels = c('100ºW', '0º', '100ºE')) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(),
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5))
ggsave('Figure 3a.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm', device = grDevices::tiff)


#------------------------------------ Figure 3b ---------------------------------------

# Calibration area + calibration points

library(stars)
occ_cal <- read.csv('Carcharias taurus/Calibration_points.csv') 
env.M <- read_stars('Carcharias taurus/Calibration_areas.tif')
env_buff <- env.M > -Inf
env_buff <- st_as_sf(env_buff, as_points = F, merge = T)

ggplot() +
  geom_sf(data = env_buff, col = 'black', fill = '#33B7A0', linewidth = 0.5) +
  geom_sf(data = coast, col = 'grey30', fill = 'grey50', linewidth = 0.1) +
  geom_point(data = occ_cal, aes(x = Longitude, y = Latitude), shape = 21, size = 0.5, color = '#F9FB0E', fill = '#352A87') +
  scale_fill_manual(values = '#33B7A0', na.value = 'grey95') + coord_equal(expand = 0) +
  scale_y_continuous(name = NULL, breaks = c(-50, 0, 50), labels = c('50ºS', '0º', '50ºN')) + 
  scale_x_continuous(name = NULL, breaks = c(-100, 0, 100), labels = c('100ºW', '0º', '100ºE')) +
  coord_sf(xlim = c(-180, 180), ylim = c(-68.9, 83.6), expand = F) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(),
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5))
ggsave('Figure 3b.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm', device = grDevices::tiff)


#------------------------------------ Figure 3c & S1.7 ---------------------------------------

# Annual model response plots

# Selected feature class and regularization multiplier
fc <- 'lq' 
rm <- 0.1

occ_cal <- read.csv('Carcharias taurus/Calibration_points.csv') 
occ_cal <- occ_cal[, c('Longitude', 'Latitude')]
env.M <- stack('Carcharias taurus/Calibration_areas.tif') 
var_set <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity', 'Distance_to_coast', 'Thermal_fronts')
names(env.M) <- var_set

p_df = data.frame()
a_df = data.frame()
plot_data_df = data.frame()

source('prepareSWD_adw.R') # Function from the 'SDMtune' package with a corrected bug

for(j in names(env.M)) {
  
  # Prepare data
  env <- env.M[[j]]
  env <- stack(env)
  Vars <- j
  
  # Background points
  set.seed(111)
  notna <- which(complete.cases(values(env)))
  samp <- sample(notna, 10000, replace = F)
  samplocs <- as.data.frame(xyFromCell(env, samp))
  
  # SWD object
  data <- prepareSWD_adw(species = 'Carcharias taurus', p = occ_cal, a = samplocs, env = env)
  
  # Run maxent replicas
  folds = randomFolds(data, k = 10, only_presence = T, seed = 111)
  default_model <- train(method = 'Maxent', data = data, fc = fc, reg = rm, iter = 1000, folds = folds)
  
  # Presences and pseudo-absences with environmental data
  p <- .get_presence(default_model@data)
  p$var <- j
  names(p)[names(p) == j] <- 'values'
  a <- .get_absence(default_model@data)
  a$var <- j
  names(a)[names(a) == j] <- 'values'
  
  pred <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 10))
  for(i in 1:10){
    pred[, i] <- predict(default_model@models[[i]], data = data@data, type = 'logistic')
  }
  
  # Plot data
  plot_data <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 4))
  names(plot_data) <- c('mean', 'sd', 'max', 'min')
  plot_data$mean <- rowMeans(pred)
  plot_data$sd <- apply(pred, 1, sd)
  plot_data$max <- plot_data$mean + plot_data$sd
  plot_data$min <- plot_data$mean - plot_data$sd
  plot_data$var <- j
  plot_data$values <- data@data[, j]
  
  # Data frames
  p_df <- rbind(p_df, p)
  a_df <- rbind(a_df, a)
  plot_data_df <- rbind(plot_data_df, plot_data)
  
}

# Reorder
order_i_want1 <- c('Distance_to_coast', 'Temperature', 'Surface_temperature', 'Salinity', 'Kd490', 'Primary_productivity', 'Thermal_fronts')
p_df <- p_df[!(p_df$var == 'Distance_to_coast' & p_df$values > 600), ]
p_df <- p_df[!(p_df$var == 'Primary_productivity' & p_df$values > 0.07), ]
p_df <- p_df[!(p_df$var == 'Thermal_fronts' & p_df$values > 0.7), ]
p_df <- transform(p_df, var = factor(var, levels = order_i_want1))

a_df <- a_df[!(a_df$var == 'Distance_to_coast' & a_df$values > 600), ]
a_df <- a_df[!(a_df$var == 'Primary_productivity' & a_df$values > 0.07), ]
a_df <- a_df[!(a_df$var == 'Thermal_fronts' & a_df$values > 0.7), ]
a_df <- transform(a_df, var = factor(var, levels = order_i_want1))

plot_data_df <- plot_data_df[!(plot_data_df$var == 'Distance_to_coast' & plot_data_df$values > 600), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Primary_productivity' & plot_data_df$values > 0.07), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Thermal_fronts' & plot_data_df$values > 0.7), ]
plot_data_df <- transform(plot_data_df, var = factor(var, levels = order_i_want1))

# Labeller
VAR_names = as_labeller(c(Temperature = 'Bottom~temperature~(ºC)', Distance_to_coast = 'Distance~to~coast~(km)',
                          Surface_temperature = 'Surface~temperature~(ºC)', Kd490 = 'Kd490~coefficient~(m^-1)', 
                          Salinity = 'Bottom~salinity~(psu)', Primary_productivity = 'Primary~prod.~(gCm^-2~día^-1)',
                          Thermal_fronts = 'Thermal~fronts~(ºC)'), default = label_parsed)

# Plot
ggplot(data = plot_data_df, aes(x = values, y = mean, ymin = min, ymax = max)) + 
  geom_line(colour = '#352A87') + 
  geom_ribbon(fill = '#352A87', alpha = 0.2) +
  scale_color_manual(values = rev(pals::parula(7))) +
  geom_rug(data = p_df, inherit.aes = F, aes(values), sides = 't', color = '#F9FB0E', linewidth = 0.3) + 
  geom_rug(data = a_df, inherit.aes = F, aes(values), sides = 'b', color = '#33B7A0', linewidth = 0.3) + 
  labs(x = NULL, y = 'Logistic output') + ylim(0, 1) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.2),
        panel.grid.major = element_line(linewidth = 0.2, colour = 'grey90'),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(vjust = -0.5, size = 9),
        plot.margin = unit(c(0, 0, 0.2, 0.2), 'cm'),
        axis.title = element_text(size = 10)) +
  facet_wrap(~var, scales = 'free_x', labeller = VAR_names, nrow = 3) 
ggsave('Figure 3c.pdf', width = 15, height = 15, units = 'cm')


#------------------------------------ Figure 4 & Figure 5b ---------------------------------------

# Annual population-scale habitat suitability in the Southwest Atlantic

pred <- raster('Carcharias taurus/Predictions.tif')
swa_ext <- extent(-68, -37, -46, -12)
pred <- crop(pred, swa_ext)

occ <- read.csv('Free_occurrence_records.csv') 
occ <- subset(occ, Species == 'Carcharias taurus')
occ <- subset(occ, Region == 'Southwest Atlantic')

df <- data.frame(coordinates(pred), as.data.frame(pred))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = layer)) + 
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), color = 'grey10', fill = 'grey50', linewidth = 0.15) +
  geom_point(data = occ, aes(x = Longitude, y = Latitude), shape = 21, size = 1, stroke = 0.2, col = 'darkred', fill = 'darkorange2') + 
  scale_fill_gradientn(colors = pals::parula(1000), na.value = '#352A87') + 
  scale_y_continuous(name = NULL, breaks = c(-44, -34, -24, -14), labels = c('44º', '34º', '24º', '14º')) +
  scale_x_continuous(name = NULL, breaks = c(-66, -58, -50, -42), labels = c('66º', '58º', '50º', '42º')) +
  coord_equal(xlim = c(-68, -37), ylim = c(-46, -12), expand = F) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(),
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5))
ggsave('Figure 4a.tiff', dpi = 900, width = 10, height = 10, units = 'cm', device = grDevices::tiff)


# Annual population-scale binary habitat suitability in the Southwest Atlantic

pred <- raster('Carcharias taurus/Predictions.tif')
bg_binary <- pred
bg_binary[bg_binary >= 0] <- 1

occ_cal <- read.csv('Carcharias taurus/Calibration_points.csv') 

# Threshold (%5 minimum training presence)
mtp5 <- sdm_threshold.5(pred, occ_cal[, c('Longitude', 'Latitude')], 'p05', binary = F)
bin_5 <- pred
bin_5[bin_5 < minValue(mtp5)] <- NA
bin_5[bin_5 >= minValue(mtp5)] <- 2

# Threshold (%10 minimum training presence)
mtp10 <- sdm_threshold.10(pred, occ_cal[, c('Longitude', 'Latitude')], 'p10', binary = F)
bin_10 <- pred
bin_10[bin_10 < minValue(mtp10)] <- NA
bin_10[bin_10 >= minValue(mtp10)] <- 3

# Final mosaic
mod_mosaic <- mosaic(bg_binary, bin_10, bin_5, fun = sum)
mod_mosaic[mod_mosaic == 1] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 3] <- 2 # MTP5
mod_mosaic[mod_mosaic == 6] <- 3 # MTP overlap
mod_mosaic <- crop(mod_mosaic, swa_ext)

Col <- c('grey95', '#352A87', '#33B7A0', '#F9FB0E', 'grey90')

df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), color = 'grey10', fill = 'white', linewidth = 0.15) +
  geom_path(data = Border, aes(x = long, y = lat, group = group), color = 'grey20', linewidth = 0.2) +
  scale_y_continuous(name = NULL, breaks = c(-44, -34, -24, -14), labels = c('44º', '34º', '24º', '14º')) +
  scale_x_continuous(name = NULL, breaks = c(-66, -58, -50, -42), labels = c('66º', '58º', '50º', '42º')) +
  coord_equal(xlim = c(-68, -37), ylim = c(-46, -12), expand = F) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5)) 
ggsave('Figure 4b.tiff', dpi = 900, width = 10, height = 10, units = 'cm', device = grDevices::tiff)


# Southern limit

# Provincial borders
provincias <- readOGR(dsn = Borders, layer = 'provincias_argentina')

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), color = 'grey10', fill = 'white', linewidth = 0.15) +
  geom_path(data = provincias, aes(x = long, y = lat, group = group), color = 'grey70', linewidth = 0.15, linetype = "dashed") +
  scale_y_continuous(name = NULL, breaks = c(-46, -42, -38), labels = c('46º', '42º', '38º')) +
  scale_x_continuous(name = NULL, breaks = c(-66, -62, -58), labels = c('66º', '62º', '58º')) +
  coord_equal(xlim = c(-68, -56), ylim = c(-48, -36), expand = F) +
  theme(panel.background = element_rect(fill = NULL), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5)) 
ggsave('Figure 5b.tiff', dpi = 900, width = 10, height = 10, units = 'cm', device = grDevices::tiff)


#------------------------------------ Figure S1.1 ---------------------------------------

# Global map of freely accessible occurrence distribution of Carcharhinus brachyurus

dat <- read.csv('Free_occurrence_records.csv')
dat <- subset(dat, Species == 'Carcharhinus brachyurus')
dat <- subset(dat, Source1 %in% c('GBIF', 'OBIS', 'Published literature', 'Grey literature', 'Social media'))
dat$Region <- as.factor(dat$Region)
dat <- transform(dat, Region = factor(Region, levels = c('Northeast Pacific', 'Southeast Pacific', 'Southwest Atlantic', 'Northeast Atlantic',
                                                         'Southeast Atlantic', 'Northwest Pacific', 'Australia', 'New Zealand')))
levels(dat$Region) <- c('Northeast Pacific', 'Southeast Pacific', 'Southwest Atlantic', 'Northeast Atlantic',
                        'Southeast Atlantic', 'Northwest Pacific', 'Australia', 'New Zealand')

ggplot() +
  geom_point(data = dat, aes(x = Longitude, y = Latitude, color = Region), size = 5) + 
  scale_color_viridis_d(option = 'G', direction = -1) +
  geom_point(data = dat, aes(x = Longitude, y = Latitude), shape = 21, size = 1.25, stroke = 0.001, col = 'darkred', fill = 'darkorange') + 
  borders('world', colour = 'black', fill = 'black', size = 0.2) +
  borders('world', colour = 'white', fill = 'white', size = 0.01) +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85), expand = 0) + 
  scale_y_continuous(name = 'Latitude', breaks = c(-50, 0, 50), labels = c('50ºS', '0º', '50ºN')) + 
  scale_x_continuous(name = 'Longitude', breaks = c(-100, 0, 100), labels = c('100ºW', '0º', '100ºE')) +
  theme(panel.background = element_rect(fill = 'grey95'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5))
ggsave('Figure S1.1.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm', device = grDevices::tiff)


# Global map of freely accessible occurrence distribution of Carcharias taurus

dat <- read.csv('Free_occurrence_records.csv') 
dat <- subset(dat, Species == 'Carcharias taurus')
dat <- subset(dat, Source1 %in% c('GBIF', 'OBIS', 'Published literature', 'Grey literature', 'Social media'))
dat$Region <- as.factor(dat$Region)
dat <- transform(dat, Region = factor(Region, levels = c('Northwest Atlantic', 'Southwest Atlantic', 'Northeast Atlantic', 'Southeast Atlantic',
                                                         'Northwest Indian', 'Northwest Pacific', 'Australia')))
levels(dat$Region) <- c('Northwest Atlantic', 'Southwest Atlantic', 'Northeast Atlantic', 'Southeast Atlantic',
                        'Northwest Indian', 'Northwest Pacific', 'Australia')

ggplot() +
  geom_point(data = dat, aes(x = Longitude, y = Latitude, color = Region), size = 5) + 
  scale_color_manual(values = rev(pals::parula(7))) +
  geom_point(data = dat, aes(x = Longitude, y = Latitude), shape = 21, size = 1.25, stroke = 0.001, col = 'darkred', fill = 'darkorange') + 
  borders('world', colour = 'black', fill = 'black', size = 0.2) +
  borders('world', colour = 'white', fill = 'white', size = 0.01) +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85), expand = 0) + 
  scale_y_continuous(name = 'Latitude', breaks = c(-50, 0, 50), labels = c('50ºS', '0º', '50ºN')) + 
  scale_x_continuous(name = 'Longitude', breaks = c(-100, 0, 100), labels = c('100ºW', '0º', '100ºE')) +
  theme(panel.background = element_rect(fill = 'grey95'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5))
ggsave('Figure S1.2.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm', device = grDevices::tiff)


#------------------------------------ Figure S1.2 ---------------------------------------

# Total occurrences compiled per species

# Occurrences
dat <- read.csv('Free_occurrence_records.csv')
dat <- subset(dat, Region == 'Southwest Atlantic')
dat <- transform(dat, Species = factor(Species, levels = c('Galeorhinus galeus', 'Notorynchus cepedianus', 'Carcharhinus brachyurus', 'Carcharias taurus')))
coast <- crop(coast0, extent(c(-70, -38, -53.2, -18)))

ggplot() +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'grey50', fill = 'white', linewidth = 0.01) +
  geom_point(data = dat, aes(x = Longitude, y = Latitude), shape = 21, size = 1.25, stroke = 0.01, color = 'black', fill = 'darkblue') +
  scale_y_continuous(name = 'Latitude (S)', breaks = c(-48, -36, -24), labels = c('48º', '36º', '24º')) + 
  scale_x_continuous(name = 'Longitude (W)', breaks = c(-68, -56, -44), labels = c('68º', '56º', '44º')) +
  facet_wrap(~ Species, ncol = 2) +
  coord_equal(expand = 0) + 
  theme(panel.background = element_rect(fill = 'grey95'),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.25))
ggsave('Figure S1.2.pdf', width = 15, height = 15, units = 'cm')


#------------------------------------ Figure S1.3 ---------------------------------------

# Discrepancies in depth variable

# Read predictors and name them
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope', 'Thermal_fronts')
names(env) <- var_names

# Read occurrences
dat <- read.csv('Free_occurrence_records.csv') # not all are here as some are not freely accessible
dat <- subset(dat, Species == 'Carcharhinus brachyurus')
dat <- subset(dat, Source1 %in% c('GBIF', 'OBIS', 'Government data', 'Published literature', 'Grey literature', 'Social media'))

# Remove duplicates
dups <- duplicated(dat[c('Longitude', 'Latitude')])
occ_cal <- dat[!dups, ]

# Extract values
values <- data.frame(Depth = extract(env[['Bathymetry']], SpatialPoints(occ_cal[, c('Longitude', 'Latitude')], crs)))

ggplot(data = values, aes(sample(seq_along(Depth)), Depth)) +
  geom_point(col = "#0B0405FF", fill = "#A0DFB9FF", shape = 21, stroke = 0.75) +
  geom_hline(yintercept = min(-dat$Depth, na.rm = TRUE)) +
  ylab('Depth (m)') + xlab('Number of records') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_line(color = 'grey95'),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure S1.3.pdf', width = 15, height = 10, units = 'cm')


#------------------------------------ Figure S1.4 ---------------------------------------
  
# Discrepancies in the depth variable

# Read predictors and name them
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope', 'Thermal_fronts')
names(env) <- var_names

# Read occurrences
dat <- read.csv('Free_occurrence_records.csv') # not all are here as some are not freely accessible
dat <- subset(dat, Species == 'Carcharias taurus')
dat <- subset(dat, Source1 %in% c('GBIF', 'OBIS', 'Government data', 'Published literature', 'Grey literature', 'Social media'))

# Remove duplicates
dups <- duplicated(dat[c('Longitude', 'Latitude')])
occ_cal <- dat[!dups, ]

# Extract values
values <- data.frame(Depth = extract(env[['Bathymetry']], SpatialPoints(occ_cal[, c('Longitude', 'Latitude')], crs)))

ggplot(data = values, aes(sample(seq_along(Depth)), Depth)) +
  geom_point(col = "#352A87", fill = "#F9FB0E", shape = 21, stroke = 0.75) +
  geom_hline(yintercept = min(-dat$Depth, na.rm = TRUE)) +
  ylab('Depth (m)') + xlab('Number of records') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_line(color = 'grey95'),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure S1.4.pdf', width = 15, height = 10, units = 'cm')


#------------------------------------ Figure S1.6 ---------------------------------------

# Occurrences in unsuitable areas
occ <- read.csv('Free_occurrence_records.csv') 
occ <- subset(occ, Species == 'Carcharhinus brachyurus')
occ <- subset(occ, Region == 'Southwest Atlantic')
occ <- subset(occ, Source2 %in% c('Mas Bervejillo (2012)', 'Soto & Mincarone (2004)*'))

# Calibration points
occ_cal <- read.csv('Carcharhinus brachyurus/Calibration_points.csv') 

# Environmental space
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope', 'Thermal_fronts')
names(env) <- var_names

# Extract values
values1 <- data.frame(Distance = extract(env[['Distance_to_coast']], SpatialPoints(occ_cal[, c('Longitude', 'Latitude')], crs)))
values2 <- data.frame(Distance = extract(env[['Distance_to_coast']], SpatialPoints(occ[, c('Longitude', 'Latitude')], crs)))

ggplot(data = values1, aes(sample(seq_along(Distance)), Distance)) + ylim(0, 495) +
  geom_point(col = "#0B0405FF", fill = "#357BA2FF", shape = 21, stroke = 0.75) +
  geom_point(data = values2, aes(sample(seq_along(Distance)), Distance), 
             col = "#0B0405FF", fill = "#A0DFB9FF", shape = 21, stroke = 0.75) +
  ylab('Distance to coast (km)') + xlab('Number of records') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_line(color = 'grey95'),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5)) 
ggsave('Figure S1.6.pdf', width = 15, height = 10, units = 'cm')


#------------------------------------ END ---------------------------------------




