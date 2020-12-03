##############################################################################
# Predict various models to full study area and write to raster
# Requires functions in the *HOTR_model_functions.r* file
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 05/21/2020
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
rasterOptions(tmpdir = "E:/MMF/temp")

source("HOTR_model_functions.r")
ncores <- 14  # physical cores = 16, setting to 14 prevents corruption of tile rasters
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=4", "BIGTIFF=YES")

# The model names
BRTm <- "BRT_gbmholdout_50ktMay13.rds"
GLMm <- "GLMER_May14_70pct.rds"
RFm <- "RF_4800trees_May12.rds"

# Load raster names
env_inputs <- fread("env_inputs_04022020.csv", header = TRUE, sep = ",")

### Only need to run tilebuilder once ###
tiles_index <- tilebuilder(raster = raster(env_inputs$raster[1]), out = 'data.frame')
# Tiles are in a rectangular matrix, but study area does not completely fill it.
# Save tiles_index to a file, edit out those tiles that are all NAs
fwrite(tiles_index, file="tiles_index.csv")

### Once the above has been done once, start from Here ###
# Reload the edited file (went from 322 to 213 tiles)
tiles_index <- fread("tiles_index.csv", header = TRUE, sep = ",")
######

### Make the raster tiles for each set of inputs ###
# RF and BRT need these
modR <- readRDS(file.path(pth,RFm))
subdfR <- env_inputs[label %in% names(modR$forest$ncat)] #RF
layerStk <- raster::stack(subdfR$raster)
names(layerStk) <- subdfR$label

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

foreach(i_tile = 1:nrow(tiles_index), .packages = "raster") %dopar% {
  maketiles(tiles_index, layerStk, gopts, i_tile)
}

snow::stopCluster(cl)

# GLM needs these.
# Note: move the previously created tiles to another directory first
modG <- readRDS(file.path(pth,GLMm))
subdfG <- env_inputs[label %in% names(modG@frame)] #GLMM
layerStk <- raster::stack(subdfG$raster)
names(layerStk) <- subdfG$label

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

foreach(i_tile = 1:nrow(tiles_index), .packages = "raster") %dopar% {
  maketiles(tiles_index, layerStk, gopts, i_tile)
}

snow::stopCluster(cl)
######

### Predict each model###
## RF
# Here's where the raster stack Tiles have been moved to
tile_pth <- "H:/HOTR_models/RF_tiles"
tile_prefix <- "RF_70pct_"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)
# this bit takes anywhere from ~14.5 hours to 2.5 days
pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "randomForest")) %dopar% {
                        raspred(tiles_index$tile[i_tile], tiles_pth = tile_pth, model = modR,
                                lyrnames = subdfR$label, oname = tile_prefix, modtype = "RF")
                      }

snow::stopCluster(cl)

# Merge the prediction tiles into single raster
# Now that cluster is released, can use more threads
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=14", "BIGTIFF=YES")
pred_tiles$filename <- str_replace(RFm, "rds", "tif")
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)

## BRT
modB <- readRDS(file.path(pth,BRTm))
# not really necessary (same as RF), but good to remember how to read the BRT model
subdfB <- env_inputs[label %in% modB$var.names] #BRT
tile_prefix <- "BRT_70pct_"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "gbm")) %dopar% {
                        raspred(tiles_index$tile[i_tile], tiles_pth = tile_pth, model = modB,
                                lyrnames = subdfB$label, oname = tile_prefix, modtype = "BRT")
                      }

snow::stopCluster(cl)

pred_tiles$filename <- str_replace(BRTm, "rds", "tif")
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)

## GLMM
tile_pth <- "H:/HOTR_models/GLM_tiles"
tile_prefix <- "GLMER_70pct_"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "lme4")) %dopar% {
                        raspred(tiles_index$tile[i_tile], tiles_pth = tile_pth, model = modG,
                                lyrnames = subdfG$label, oname = tile_prefix, modtype = "GLM")
                      }

snow::stopCluster(cl)

pred_tiles$filename <- str_replace(GLMm, "rds", "tif")
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
#####