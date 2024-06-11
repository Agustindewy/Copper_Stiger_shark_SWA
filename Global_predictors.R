
# Global scale scenopoetic predictors for C. brachyurus and C. taurus: download and preparation

library(sdmpredictors)
library(raster)
library(grec)

setwd('WORKING DIRECTORY')

#----------------------------------------- Predictors ------------------------------------------------

# Bio-Oracle
list_layers(c('Bio-ORACLE'))[, 2:4] # explore layers

# Check correlation at global scale
layers_correlation(c('BO2_tempmean_bdmean', 'BO_sstmean', 'BO2_ppmean_bdmean', 'BO_damean', 'BO_chlomean',
                     'BO2_salinitymean_bdmean', 'BO_bathymean'))
# BO_chlomean vs BO_damean corr > 0.8, we discard BO_chlomean

bio <- load_layers(c('BO2_tempmean_bdmean', 'BO_sstmean', 'BO2_ppmean_bdmean', 'BO_damean',
                     'BO2_salinitymean_bdmean', 'BO_bathymean'))

# Distance to coastline from Global Self-consistent, Hierarchical, High-resolution Geography Database
# Downloaded from https://www.soest.hawaii.edu/pwessel/gshhg/
dis <- raster('DOWNLOAD AND READ THE RASTER.tif') # read the tif
dis <- aggregate(dis, fact = 2, fun = mean) # reduce resolution to match Bio-ORACLE predictors

# Slope (created from bathymetry raster)
slo <- terrain(bio[['BO_bathymean']], opt = 'slope', unit = 'degrees', neighbors = 4)
slo <- projectRaster(slo, bio)

# Thermal fronts (created from sea surface temperature raster)
sst <- bio[['BO_sstmean']]
tfr <- detectFronts(sst, method = 'median_filter')
m <- matrix(1, ncol = 5, nrow = 5)
tfr <- focal(tfr, m, fun = mean, NAonly = T, pad = T, na.rm = T) 
tfr <- mask(tfr, sst)

# Final stack of predictors
env <- stack(bio, dis, slo, tfr)

# Check correlation again with new predictors
layerStats(env, 'pearson', na.rm = T) # no correlation

# Check if environmental values make sense in each case
env # bathymetry has some positive values, convert them to NA
env[[6]][env[[6]] > -1] <- NA

writeRaster(env, filename = 'predictors.tif') # save as tif


#-------------------------------------- END ------------------------------------------
