##############################################################################
# Predict various models to full study area and write to raster
# Requires functions in the *HOTR_model_functions.r* file
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 03/12/2020
#
# Code licensed under the GNU General Public License version 3.
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/
#############################################################################

library(data.table)
library(stringr)
library(raster)
library(rgdal)
library(gbm)
library(randomForest)
library(lme4)
library(doSNOW)

pth <- "H:/HOTR_models/"
setwd(pth)
source("HOTR_model_functions.r")

rasterOptions(tmpdir = "E:/MMF/temp",
              chunksize = 5e+08)
ncores <- 16  # Set to number of physical cores, not logical processors
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=4", "BIGTIFF=YES")

### Predict each model - have to do one at a time ###
modfile <- "RF_Oct22_full.rds"  #"BRT_Oct28_full.rds" "GLMER_Oct22_full.rds"

# Load all raster names to start
some_inputs <- fread("env_inputs_11252019.csv", header = TRUE, sep = ",")
add_inputs <- fread("onehot_rasters.csv", header = TRUE, sep = ",")
env_inputs <- rbind(some_inputs[, .(raster, label)], add_inputs[, .(raster, label)])

## Pair down list to only those inputs used in each model. ##
# How to do this changes with model type. Uncomment the correct one.
# - Random Forest
mod <- readRDS(paste0(pth,modfile))
subdf <- env_inputs[label %in% names(mod$forest$ncat)] #RF

# - Boosted Regression Tree
#mod <- readRDS(paste0(pth,modfile))
#subdf <- env_inputs[label %in% mod$var.names] #BRT

# - Generalized Linear Mixed Model
#mod <- readRDS(file.path(pth,modfile))
#subdf <- env_inputs[label %in% names(mod@frame)] #GLMM

### Only need to run tilebuilder once ###
tiles_index <- tilebuilder(raster = raster(env_inputs$raster[1]), out = 'data.frame')
# Tiles are in a rectangular matrix, but study area does not completely fill it.
# Save tiles_index to a file, edit out those tiles that are all NAs
fwrite(tiles_index, file="tiles_index.csv")

### Once the above has been done once, start from Here ###
# Reload the edited file (went from 322 to 213 tiles)
tiles_index <- fread("tiles_index.csv", header = TRUE, sep = ",")

# Make raster stack
layerStk <- raster::stack(subdf$raster)
names(layerStk) <- subdf$label

# Start parallel cluster
cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

## Create the raster tiles ##
foreach(i_tile = 1:nrow(tiles_index), .packages = "raster") %dopar% {
  ext_i <- extent(tiles_index$xmin[i_tile], tiles_index$xmax[i_tile],
                  tiles_index$ymin[i_tile], tiles_index$ymax[i_tile])
  fname <- file.path(pth, paste0("tile_", as.character(tiles_index$tile[i_tile]), ".tif"))
  if(!file.exists(fname)){
    tile <- crop(layerStk, ext_i, filename=fname, format="GTiff",
                 overwrite=TRUE, options = gopts)
  }
}

## Use the tiles to create a prediction raster ##
# Make sure you're using the correct tile_prefix and package
tile_prefix <- "GLMER_Oct22_" # "RF_Oct22_" "BRT_Oct28_"
pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "randomForest")) %dopar% {  #"randomForest" "gbm" "lme4"
                        raspred(tiles_index$tile[i_tile], model = mod,
                                lyrnames = subdf$label, oname = tile_prefix, modtype = "RF")
                      }

# Release cluster to free up memory
snow::stopCluster(cl)

## Merge the prediction tiles into single raster ##
# Now that cluster is released, can use more threads
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=14", "BIGTIFF=YES")
pred_tiles$filename <- str_replace(modfile, "rds", "tif")
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
