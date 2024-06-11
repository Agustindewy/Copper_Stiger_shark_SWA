
# De Wysiecki et al. - Global ENM projection for the copper shark (Carcharhinus brachyurus)

# Preparation of data and global niche modeling

library(sp)
library(raster)
library(sf)
library(rgdal)
library(sqldf)
library(maps)
library(ggplot2)
library(rgeos)
library(dismo)
library(blockCV)
library(kuenm)
library(ntbox)# available on GitHub (luismurao/ntbox)
library(rSDM) # available on GitHub (Pakillo/rSDM)
library(ENMeval) # available on GitHub (jamiemkass/ENMeval)
library(SDMtune) # available on GitHub (ConsBiol-unibern/SDMtune)

setwd('SET YOUR WORKING DIRECTORY')

# Projection
crs <- CRS('+proj=longlat +datum=WGS84')

# Specie
species <- 'Carcharhinus brachyurus'

# Seed
set.seed(111)

# Create folder for the species
dir.create(paste(species))

#-------------------------------- G space --------------------------------------------

# Read predictors and name them
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope', 'Thermal_fronts')
names(env) <- var_names

# Discard bathymetry and slope due to discrepancies on the continental shelf margin (De Wysiecki et al. 2022) 
var_set <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity', 
             'Distance_to_coast', 'Thermal_fronts')

# G space
env.G <- env[[var_set]]

#-------------------------------- Preparation ----------------------------------

# Read occurrences
dat <- read.csv('Free_occurrence_records.csv') # not all are included as some are not publicly accessible
dat <- subset(dat, Species == species)
dat <- subset(dat, Source1 %in% c('GBIF', 'OBIS', 'Government data', 'Published literature', 'Grey literature', 'Social media'))

# Remove duplicates
dups <- duplicated(dat[c('Longitude', 'Latitude')])
occ_cal <- dat[!dups, c(1, 4, 5)]

# Spatial filtering
train_thin <- thin(occ_cal, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                   thin.par = 50, reps = 10, locs.thinned.list.return = T, write.files = F, write.log.file = F)
occ_cal <- data.frame(Species = species, Longitude = round(train_thin[[1]]$Longitude, 2),
                      Latitude = round(train_thin[[1]]$Latitude, 2))

# Initial M space for environmental filtering
# some data points fall in New Zealand, requiring buffers that exceed the 180° meridian
# therefore, recenter the points to 110°
occ_cal$Longitude2 <- occ_cal$Longitude + 110
occ_cal$Longitude2 <- ifelse(occ_cal$Longitude2 > 180, (occ_cal$Longitude2 - 180) - 180, occ_cal$Longitude2)
coord <- data.frame(Longitude = occ_cal$Longitude2, Latitude = occ_cal$Latitude)
occ.tot = st_as_sf(coord, coords = c('Longitude', 'Latitude'), crs = 4326)
occ.buff <- st_buffer(occ.tot, 1000000) # 1000 km buffer
# recenter the raster also to 110°
env1 <- crop(env.G, extent(-180, 70, -90, 90))
env2 <- crop(env.G, extent(70, 180, -90, 90))   
extent(env1) <- c(-70, 180, -90, 90)
extent(env2) <- c(-180, -70, -90, 90)
env0 <- merge(env1, env2)
names(env0) <- var_set
env.M <- crop(env0, occ.buff) # crop the raster to the buffer
env.M <- mask(env.M, occ.buff) # mask the raster to the buffer
env.M <- stack(env.M)

# Rescue records that fall outside any raster layer and move them to the nearest cell with data
for(i in 1:length(env.M@layers)){
  spp <- SpatialPoints(occ_cal[, c('Longitude2', 'Latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M, layer = i, move = T)
  occ_cal$Longitude2 <- round(spp_corrected@coords[, 1], 3)  # replace coordinates, including new ones if any
  occ_cal$Latitude <- round(spp_corrected@coords[, 2], 3)
} # plots appear if any corrections are applied

# Remove duplicates again in case any have emerged
dups <- duplicated(occ_cal[c('Longitude2', 'Latitude')])
occ_cal <- occ_cal[!dups, ]

# Environmental filtering
source('envSample.R')
coords <- SpatialPoints(occ_cal[, c('Longitude2', 'Latitude')], crs)
env.data <- extract(env.M, coords)
env.data <- as.data.frame(env.data)
coords <- envSample(coords, filters = list(env.data$Surface_temperature, env.data$Distance_to_coast),
                    res = list(0.5, 1), do.plot = T)
occ_cal <- merge(occ_cal, coords, by.x = c('Longitude2', 'Latitude'), by.y = c('lon', 'lat'))

# Outliers in environmental space (Cobos et al., 2018)
variables_values <- na.omit(values(env.M))
occ_variables <- na.omit(cbind(as.numeric(row.names(occ_cal)), extract(env.M, SpatialPoints(occ_cal[, c('Longitude2', 'Latitude')], crs))))
occ_variables <- as.data.frame(occ_variables)
if(dim(variables_values)[1] > 10000) {
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ]
}

# Plot all combinations and look for outliers, writing the number below
for(i in 1:dim(env.M)[3]){
  for(j in 1:dim(env.M)[3]){
    par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
    plot(variables_values[, i], variables_values[, j], col = 'grey65', pch = 1,
         xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    points(occ_variables[, i + 1], occ_variables[, j + 1], col = 'blue', pch = 19, cex = 1.5)
    legend('bottomright', legend = c('Region of interest', 'Occurrences'),
           pch = c(1, 19), col = c('grey65', 'black'), bty = 'n')
    plot(variables_values[, i], variables_values[, j], col = 'grey65',
         pch = 1, xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    text(occ_variables[, i + 1], occ_variables[, j + 1], occ_variables[, 1], cex = 1, col = 'blue')
    legend('bottomright', legend = 'Occurrence ID', bty = 'n')
  }
}

# Remove outliers
occ_variables <- subset(occ_variables, !V1 %in% c(99999)) # enter the numbers here
occ_cal <- occ_cal[which(row.names(occ_cal) %in% as.character(occ_variables$V1)), ]

# Final calibration points
occ_cal <- rbind(occ_cal[, c(3, 4, 2)])
write.csv(occ_cal, paste(species, '/', 'Calibration_points.csv', sep = ''), row.names = F)


#-------------------------------- Calibration area ----------------------------------------------

occ_cal <- read.csv(paste(species, '/', 'Calibration_points.csv', sep = ''))

# Final M space (calibration area)
# Recenter due to data in New Zealand
occ_cal$Longitude2 <- occ_cal$Longitude + 110
occ_cal$Longitude2 <- ifelse(occ_cal$Longitude2 > 180, (occ_cal$Longitude2 - 180) - 180, occ_cal$Longitude2)
coord <- data.frame(Longitude = occ_cal$Longitude2, Latitude = occ_cal$Latitude)
occ.tot = st_as_sf(coord, coords = c('Longitude', 'Latitude'), crs = 4326)
occ.buff <- st_buffer(occ.tot, 1000000) # 1000 km buffer
# Recenter the raster
env1 <- crop(env.G, extent(-180, 70, -90, 90))
env2 <- crop(env.G, extent(70, 180, -90, 90))   
extent(env1) <- c(-70, 180, -90, 90)
extent(env2) <- c(-180, -70, -90, 90)
env0 <- merge(env1, env2)
names(env0) <- var_set
env.M <- crop(env0, occ.buff) # crop the raster to the buffer
env.M <- mask(env.M, occ.buff) # mask the raster to the buffer

# Once cropped, recenter to 0°
env.M <- extend(env.M, extent(-180, 180, -52.66666, 54.08333)) # back to -180°, 180°
env.M1 <- crop(env.M, extent(-180, -70, -52.66666, 54.08333))
env.M2 <- crop(env.M, extent(-70, 180, -52.66666, 54.08333))   
extent(env.M1) <- c(70, 180, -52.66666, 54.08333)
extent(env.M2) <- c(-180, 70, -52.66666, 54.08333)
env.M <- merge(env.M2, env.M1)
names(env.M) <- var_set
env.M <- stack(env.M)

# Remove some incoherent areas within the buffers that cross land, e.g., Dead Sea
mar_muerto <- matrix(c(23.55, 31.32, 47.24, 37.60, 23.55, 43.28, 51.13, 43.70, 35.53, 43.28), 5, 2)
mar_muerto <- SpatialPolygons(list(Polygons(list(Polygon(mar_muerto)), ID = 1)), proj4string = crs)
golfo_mx <- matrix(c(-102.82, -100.20, -82.00, -86.36, -102.82, 19.54, 30.26, 24.15, 14.65, 19.54), 5, 2)
golfo_mx <- SpatialPolygons(list(Polygons(list(Polygon(golfo_mx)), ID = 1)), proj4string = crs)
chile <- matrix(c(-85.19, -67.62, -72.24, -94.47, -85.19, -30.61, -31.32, -56.03, -50.62, -30.61), 5, 2)
chile <- SpatialPolygons(list(Polygons(list(Polygon(chile)), ID = 1)), proj4string = crs)
costa_rica <- matrix(c(-81.21, -86.61, -81.42, -79.12, -73.80, -67.91, -81.21, 18.82, 12.48, 8.20, 9.36, 5.56, 11.92, 18.82), 7, 2)
costa_rica <- SpatialPolygons(list(Polygons(list(Polygon(costa_rica)), ID = 1)), proj4string = crs)
env.M <- mask(env.M, mar_muerto, inverse = T)
env.M <- mask(env.M, golfo_mx, inverse = T)
env.M <- mask(env.M, chile, inverse = T)
env.M <- mask(env.M, costa_rica, inverse = T)

# Check correlation in the M space
layerStats(env.M, 'pearson', na.rm = T) # no significant correlation

# Save to directory
writeRaster(env.M, filename = paste(species, '/', 'Calibration_areas.tif', sep = ''))


#-------------------------------- Calibration, evaluation, and model selection ---------------------

# Calibration points
occ_cal <- read.csv(paste(species, '/', 'Calibration_points.csv', sep = ''))
occ_cal <- occ_cal[, c('Longitude', 'Latitude')]
colnames(occ_cal) <- c('longitude', 'latitude')

# Calibration areas
env.M <- stack(paste(species, '/','Calibration_areas.tif', sep = ''))
names(env.M) <- var_set

# Background points
bg_points <- as.data.frame(randomPoints(env.M[[1]], n = 10000))
colnames(bg_points) <- c('longitude', 'latitude')

# Calibration points + background
all_pts <- rbind(occ_cal, bg_points)
all_pts <- SpatialPoints(all_pts, proj4string = crs)

# Check effective range of spatial autocorrelation
eff.range <- spatialAutoRange(rasterLayer = env.M, sampleNumber = 10000, doParallel = T, showPlots = T)
median_range <- 99999 # enter the optimal range here 

# Spatial block design based on 'blockCV' by Valavi et al. (2019)
sp_block <- spatialBlock(speciesData = all_pts, rasterLayer = env.M[[1]], theRange = median_range, k = 5, selection = 'random')

# Get block information for modeling
user.grp <- list(occs.grp = sp_block$foldID[1:nrow(occ_cal)],
                 bg.grp = sp_block$foldID[(nrow(occ_cal) + 1):length(sp_block$foldID)])

# Function to include partial ROC statistic in ENM evaluation
proc <- function(vars) {
  proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(proc_auc_ratio = proc$pROC_summary[1], 
                    proc_pval = proc$pROC_summary[2], row.names = NULL)
  return(out)
} 

# Run and evaluate candidate models
mod_global <- ENMevaluate(occs = occ_cal, envs = env.G, bg = bg_points, user.grp = user.grp, 
                          algorithm = 'maxent.jar', partitions = 'user', user.eval = proc, doClamp = F, 
                          tune.args = list(fc = c('LQ', 'LP', 'LQP'), rm = c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10)))

# Model selection following the 3-step criterion by Cobos et al. (2019)
res <- eval.results(mod_global)
opt.seq <- res %>% 
  filter(delta.AICc == min(delta.AICc)) %>%
  filter(proc_pval.avg == min(proc_pval.avg)) %>%
  filter(or.10p.avg == min(or.10p.avg))

# Optimal model
mod.seq <- eval.models(mod_global)[[opt.seq$tune.args]]

# Predictor contribution to the model
mod_global@variable.importance

# Model prediction
pred.seq <- eval.predictions(mod_global)[[opt.seq$tune.args]]
writeRaster(pred.seq, filename = paste(species, '/','Predictions.tif', sep = ''))


#-------------------------------- Extrapolation risk ---------------------------------


# Based on the Mobility Oriented Parity (MOP) analysis by Owens et al. (2013)

# G space 
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope', 'Thermal_fronts')
names(env) <- var_names
var_set <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity', 
             'Distance_to_coast', 'Thermal_fronts')
env.G <- env[[var_set]]

# Coastline (spatial polygons) - taken from the GSHHG coast database as of June 15, 2017
coast <- st_read(dsn = '', layer = 'GSHHS_f_L1_World')
env.G <- crop(env.G, coast) # crop regions of the Arctic and Antarctica

# M space
env.M <- stack(paste(species, '/', 'Calibration_areas.tif', sep = ''))
names(env.M) <- var_set
mop_analysis <- mop(M_stack = env.M, G_stack = env.G, percent = 10, comp_each = 2000, parallel = T)
writeRaster(mop_analysis, filename = paste(species, '/', 'Mop.tif', sep = ''))


#-------------------------------- END -------------------------------------
