##### Raster Creation #####
# Work In Progress #
library(data.table)
library(stringr)
library(raster)
library(rgdal)
#library(gbm)
library(randomForest)
library(doSNOW)

pth <- "H:/HOTR_models/"
setwd(pth)
source("HOTR_model_functions.r")

rasterOptions(tmpdir = "E:/MMF/temp",
              chunksize = 5e+08)
ncores <- 16  # Set to number of physical cores, not logical processors
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=4", "BIGTIFF=YES")

### Predict each model - have to do one at a time ###
#BRT1 <- "BRT_Oct28_full.rds"
RF2 <- "RF_Oct22_full.rds"

# Load all raster names to start
some_inputs <- fread("env_inputs_11252019.csv", header = TRUE, sep = ",")
add_inputs <- fread("onehot_rasters.csv", header = TRUE, sep = ",")
env_inputs <- rbind(some_inputs[, .(raster, label)], add_inputs[, .(raster, label)])

# Pair down list to only those inputs used in each model.
# How to do this changes with model type
# - Random Forest
mod <- readRDS(paste0(pth,RF2))
subdf <- env_inputs[label %in% names(mod$forest$ncat)] #RF

# - Boosted Regression Tree
#mod <- readRDS(paste0(pth,BRT1))
#subdf <- env_inputs[label %in% mod$var.names] #BRT

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

# Create the raster tiles
foreach(i_tile = 1:nrow(tiles_index), .packages = "raster") %dopar% {
  ext_i <- extent(tiles_index$xmin[i_tile], tiles_index$xmax[i_tile],
                  tiles_index$ymin[i_tile], tiles_index$ymax[i_tile])
  fname <- file.path(pth, paste0("tile_", as.character(tiles_index$tile[i_tile]), ".tif"))
  tile <- crop(layerStk, ext_i, filename=fname, format="GTiff",
               overwrite=TRUE, options = gopts)
}

# Use the tiles to create a prediction raster
pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "randomForest")) %dopar% {
  raspred(tiles_index$tile[i_tile], model = mod, lyrnames = subdf$label)
}

# Release cluster to free up memory
snow::stopCluster(cl)

# Merge the prediction tiles into single raster
pred_tiles$filename <- str_replace(RF2, "rds", "tif")
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
