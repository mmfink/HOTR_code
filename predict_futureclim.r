##############################################################################
# Predict future climate models to full study area and write to raster
# Models created in Jan2021_models.r, uses Performance_functions.r
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 03/03/2021
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
.libPaths("E:/MMF/R/win-library/3.6")
Sys.setenv(R_USER = "E:/MMF")
library(stringr)
library(data.table)
library(dplyr)
library(doSNOW)
library(randomForest)
library(lme4)
library(gbm)
library(dismo)
library(ROCR)
library(ggplot2)
library(raster)
library(rgdal)

pth <- "H:/HOTR_models"
setwd(pth)
source("Performance_functions.r")
set.seed(855)

##Prediction Function##
raspred <- function(tile_i, tiles_pth=pth, model, lyrnames,
                    oname, modtype){
  fname <- file.path(tiles_pth,
                     paste0("tile_", as.character(tile_i), ".tif"))
  oname <- paste0(oname, as.character(tile_i), ".tif")
  if(file.exists(oname)){
    outras <- raster(oname)
  } else {
    ras <- brick(fname)
    names(ras) <- lyrnames
    opt <- c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES")
    if(modtype == "RF"){
      #Random Forest:
      outras <- raster::predict(ras, model=model, type = "prob", index = 2,
                                filename = oname, format = "GTiff", overwrite = TRUE,
                                options = opt)
    }
    if(modtype == "BRT"){
      #Boosted Regression Tree:
      outras <- raster::predict(ras, model=model, index=1, filename=oname, na.rm=FALSE,
                                format="GTiff", overwrite=TRUE, options=opt,
                                n.trees=26270, type="response", single.tree=FALSE)
    }
    if(modtype == "GLM"){
      #Generalized Linear Mixed Model
      outras <- raster::predict(ras, model=model, index=1, filename=oname,
                                format="GTiff", overwrite=TRUE, options=opt,
                                type="response", re.form=NA)
    }
  }
  return(outras)
}
########

BRT70 <- "BRT_gbmholdout_50ktJan11_2021.rds"
GLM70 <- "GLMER_Jan07_2021.rds"
RF70 <- "RF_4200trees_Jan12_2021.rds"

ncores <- 10
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=6", "BIGTIFF=YES")

ww_inputs <- fread("WarmWet_inputs_02042021.csv", header = TRUE, sep = ",")
hd_inputs <- fread("HotDry_inputs_02042021.csv", header = TRUE, sep = ",")
tiles_index <- fread("tiles_index.csv", header = TRUE, sep = ",")

####### ******* WarmWet ******* #######
### pair down list to only those inputs used in each model. ###
# How to do this changes with model type
mod <- readRDS(file.path(pth,RF70))
subdf <- ww_inputs[label %in% names(mod$forest$ncat)] #RF

# make raster stack #
layerStk <- raster::stack(subdf$raster)
names(layerStk) <- subdf$label

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

# Create the raster tiles #
foreach(i_tile = 1:nrow(tiles_index), .packages = "raster") %dopar% {
  ext_i <- extent(tiles_index$xmin[i_tile], tiles_index$xmax[i_tile],
                  tiles_index$ymin[i_tile], tiles_index$ymax[i_tile])
  fname <- file.path(pth, paste0("tile_", as.character(tiles_index$tile[i_tile]), ".tif"))
  if(!file.exists(fname)){
    tile <- crop(layerStk, ext_i, filename=fname, format="GTiff",
                 overwrite=TRUE, options = gopts)
  }
}

snow::stopCluster(cl)

### Use the tiles to create a prediction raster ###
### RF
tile_prefix <- "RF_WW_"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "randomForest")) %dopar% {
  raspred(tiles_index$tile[i_tile], model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "RF")
}

snow::stopCluster(cl)

pred_tiles$filename <- "RF_WarmWet.tif"
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
#####

### BRT
mod <- readRDS(file.path(pth,BRT70))
subdf <- ww_inputs[label %in% mod$var.names] #BRT
tile_prefix <- "BRT_WW_"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "gbm")) %dopar% {
  raspred(tiles_index$tile[i_tile], model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "BRT")
}

snow::stopCluster(cl)

pred_tiles$filename <- "BRT_WarmWet.tif"
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
########

### GLM
mod <- readRDS(file.path(pth,GLM70))
subdf <- ww_inputs[label %in% names(mod@frame)] #GLMM
# make raster stack #
layerStk <- raster::stack(subdf$raster)
names(layerStk) <- subdf$label

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

# Create the raster tiles
foreach(i_tile = 1:nrow(tiles_index), .packages = "raster") %dopar% {
  ext_i <- extent(tiles_index$xmin[i_tile], tiles_index$xmax[i_tile],
                  tiles_index$ymin[i_tile], tiles_index$ymax[i_tile])
  fname <- file.path(pth, paste0("tile_", as.character(tiles_index$tile[i_tile]), ".tif"))
  if(!file.exists(fname)){
    tile <- crop(layerStk, ext_i, filename=fname, format="GTiff",
                 overwrite=TRUE, options = gopts)
  }
}

snow::stopCluster(cl)

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

### Use the tiles to create a prediction raster
tile_prefix <- "GLMM_WW_"
pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "lme4")) %dopar% {
  raspred(tiles_index$tile[i_tile], model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "GLM")
}

snow::stopCluster(cl)

gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=14", "BIGTIFF=YES")
pred_tiles$filename <- "GLMM_WarmWet.tif"
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
#####

####### ******* HotDry ******* #######
### pair down list to only those inputs used in each model. ###
# How to do this changes with model type
mod <- readRDS(file.path(pth,RF70))
subdf <- hd_inputs[label %in% names(mod$forest$ncat)] #RF

# make raster stack #
layerStk <- raster::stack(subdf$raster)
names(layerStk) <- subdf$label

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

# Create the raster tiles
foreach(i_tile = 1:nrow(tiles_index), .packages = "raster") %dopar% {
  ext_i <- extent(tiles_index$xmin[i_tile], tiles_index$xmax[i_tile],
                  tiles_index$ymin[i_tile], tiles_index$ymax[i_tile])
  fname <- file.path(pth, paste0("tile_", as.character(tiles_index$tile[i_tile]), ".tif"))
  if(!file.exists(fname)){
    tile <- crop(layerStk, ext_i, filename=fname, format="GTiff",
                 overwrite=TRUE, options = gopts)
  }
}

snow::stopCluster(cl)

### Use the tiles to create a prediction raster ###
### RF
tile_prefix <- "RF_HD_"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "randomForest")) %dopar% {
  raspred(tiles_index$tile[i_tile], model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "RF")
}

snow::stopCluster(cl)

pred_tiles$filename <- "RF_HotDry.tif"
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
#####

### BRT
mod <- readRDS(file.path(pth,BRT70))
subdf <- hd_inputs[label %in% mod$var.names] #BRT
tile_prefix <- "BRT_HD_"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "gbm")) %dopar% {
  raspred(tiles_index$tile[i_tile], model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "BRT")
}

snow::stopCluster(cl)

pred_tiles$filename <- "BRT_HotDry.tif"
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
########

### GLM
mod <- readRDS(file.path(pth,GLM70))
subdf <- hd_inputs[label %in% names(mod@frame)] #GLMM
# make raster stack #
layerStk <- raster::stack(subdf$raster)
names(layerStk) <- subdf$label

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

# Create the raster tiles
foreach(i_tile = 1:nrow(tiles_index), .packages = "raster") %dopar% {
  ext_i <- extent(tiles_index$xmin[i_tile], tiles_index$xmax[i_tile],
                  tiles_index$ymin[i_tile], tiles_index$ymax[i_tile])
  fname <- file.path(pth, paste0("tile_", as.character(tiles_index$tile[i_tile]), ".tif"))
  if(!file.exists(fname)){
    tile <- crop(layerStk, ext_i, filename=fname, format="GTiff",
                 overwrite=TRUE, options = gopts)
  }
}

snow::stopCluster(cl)

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

### Use the tiles to create a prediction raster
tile_prefix <- "GLMM_HD_"
pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "lme4")) %dopar% {
  raspred(tiles_index$tile[i_tile], model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "GLM")
}

snow::stopCluster(cl)

gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=14", "BIGTIFF=YES")
pred_tiles$filename <- "GLMM_HotDry.tif"
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
#####


## Manual merging if necessary - do not run this section otherwise ##
# Identify corrupt tiles
intifs <- list.files(pattern = "BRT_HD_\\d+\\.tif")
for(t in intifs){plot(raster(t))}

# Fix
i_tile <- 131
raspred(tiles_index$tile[i_tile], model = mod,
        lyrnames = subdf$label, oname = tile_prefix, modtype = "RF")

gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=14", "BIGTIFF=YES")
outname <- "BRT_HotDry.tif"

raslist <- lapply(intifs, raster)
raslist$filename <- outname
raslist$format <- "GTiff"
raslist$options <- gopts
out <- do.call(raster::merge, raslist)
#####

### Now evaluate models with Test data and compare ###
library(dismo)
library(ROCR)
library(ggplot2)
source("Performance_functions.r")

RF_mod <- readRDS(RF70)
GLM_mod <- readRDS(GLM70)
BRT_mod <- readRDS(BRT70)

test_data <- fread("testing_data01072021.csv", header=TRUE, sep=",")
bestvars_T <- fread("vars_to_keep_Mar02.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]
bestvars_G <- fread("vars_to_keep_Feb25.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]

## **RF** ##
allcols <- c("response", varT)
testdat <- test_data[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]
testRF <- predict(RF_mod, newdata=testdat, type="prob", norm.votes=TRUE)
evRF <- data.table(prob=testRF[,2], resp=testdat[, response])
outRF <- testEval(evRF, "RF_test15pct")

## **BRT** ##
testdat <- test_data[, ..allcols][, nlcd := as.factor(nlcd)]
testBRT <- predict(BRT_mod, newdata=testdat, n.trees=26270, type="response")
evBRT <- data.table(prob=testBRT, resp=testdat[, response])
outBRT <- testEval(evBRT, "BRT_test15pct")

## **GLMM** ##
allcols <- c("response", "Grid_ID", varG)
testdat <- test_data[, ..allcols][, Grid_ID := as.factor(Grid_ID)]
testGLM <- predict(GLM_mod, newdata=testdat, type="response", allow.new.levels=TRUE)
evGLM <- data.table(prob=testGLM, resp=testdat[, response])
outGLM <- testEval(evGLM, "GLM_test15pct")

testOut <- tibble::tribble(~Name, ~AUC, ~TSS, ~err_rate, ~kappa, ~PCC, ~Sensitivity, ~Specificity, ~Threshold,
outGLM[[1]],outGLM[[2]],outGLM[[3]],outGLM[[4]],outGLM[[5]],outGLM[[6]],outGLM[[7]],outGLM[[8]],outGLM[[9]],
outRF[[1]],outRF[[2]],outRF[[3]],outRF[[4]],outRF[[5]],outRF[[6]],outRF[[7]],outRF[[8]],outRF[[9]],
outBRT[[1]],outBRT[[2]],outBRT[[3]],outBRT[[4]],outBRT[[5]],outBRT[[6]],outBRT[[7]],outBRT[[8]],outBRT[[9]])
fwrite(testOut, file = "Testing15pct_Jan192021.csv")

fullcv <- fread("CrossValidation_metrics_Sens95_01122021.csv", header=TRUE, sep=",")
tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
toadd <- melt(testOut, id.vars = "Name", measure.vars = c(2:9), variable.name = "metric")

ggplot(data=tograph, aes(x=factor(metric), y=value)) +
  geom_boxplot(aes(fill=factor(model))) +
  ylim(c(0,1)) +
  labs(x="Performance Metric", y="",
       title="10-fold Cross-Validation Plus Testing Data, Sensitivity=0.95") +
  geom_point(data=toadd, aes(x=factor(metric),y=value,
                             color=factor(Name)),
             size=4, pch=18) +
  guides(fill=guide_legend(title = "10fold Training"),
         color=guide_legend(title = "Testing data"))

## again, with Sensitiviy = Specificity
## **RF** ##
outRF <- testEval(evRF, "RF_test15pct", cutype="senspec")

## **BRT** ##
outBRT <- testEval(evBRT, "BRT_test15pct", cutype="senspec")

## **GLMM** ##
outGLM <- testEval(evGLM, "GLM_test15pct", cutype="senspec")

testOut <- tibble::tribble(~Name, ~AUC, ~TSS, ~err_rate, ~kappa, ~PCC, ~Sensitivity, ~Specificity, ~Threshold,
                           outGLM[[1]],outGLM[[2]],outGLM[[3]],outGLM[[4]],outGLM[[5]],outGLM[[6]],outGLM[[7]],outGLM[[8]],outGLM[[9]],
                           outRF[[1]],outRF[[2]],outRF[[3]],outRF[[4]],outRF[[5]],outRF[[6]],outRF[[7]],outRF[[8]],outRF[[9]],
                           outBRT[[1]],outBRT[[2]],outBRT[[3]],outBRT[[4]],outBRT[[5]],outBRT[[6]],outBRT[[7]],outBRT[[8]],outBRT[[9]])

fwrite(testOut, file = "Testing15pct_senspec_Jan192021.csv")
fullcv <- fread("CrossValidation_metrics_SensSpec_01122021.csv", header=TRUE, sep=",")
tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
toadd <- melt(testOut, id.vars = "Name", measure.vars = c(2:9), variable.name = "metric")

ggplot(data=tograph, aes(x=factor(metric), y=value)) +
  geom_boxplot(aes(fill=factor(model))) +
  ylim(c(0,1)) +
  labs(x="Performance Metric", y="",
       title="10-fold Cross-Validation Plus Testing Data, Sensitivity=Specificity") +
  geom_point(data=toadd, aes(x=factor(metric),y=value,
                             color=factor(Name)),
             size=4, pch=18) +
  guides(fill=guide_legend(title = "10fold Training"),
         color=guide_legend(title = "Testing data"))

### NOTE: the final future climate ensembles were created in ArcGIS using
# RasterCalculator "("GLMM_HotDry.tif" * 0.3042364) + ("BRT_HotDry.tif" * 0.3393725) +
# ("RF_HotDry.tif" * 0.3563911)" H:\HOTR_models\Ensemble_HotDry.tif
# and
# RasterCalculator "("GLMM_WarmWet.tif" * 0.3042364) + ("BRT_WarmWet.tif" * 0.3393725) +
# ("RF_WarmWet.tif" * 0.3563911)" H:\HOTR_models\Ensemble_WarmWet.tif
#
# Weights used were the same as those for the current climate ensemble.
