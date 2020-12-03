##############################################################################
# Predict models to full study area and write to raster
# Models created in Nov2020_models.r
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 11/30/2020
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

pth <- "H:/HOTR_models"
setwd(pth)

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
                                n.trees=20490, type="response", single.tree=FALSE)
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
####

BRT70 <- "BRT_gbmholdout_50ktNov24_2020.rds"
GLM70 <- "GLMER_Nov24_2020.rds"
RF70 <- "RF_4800trees_Nov24_2020.rds"

ncores <- 16
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=4", "BIGTIFF=YES")

env_inputs <- fread("env_inputs_04022020.csv", header = TRUE, sep = ",")
tiles_index <- fread("tiles_index.csv", header = TRUE, sep = ",")

### pair down list to only those inputs used in each model. ###
# How to do this changes with model type
mod <- readRDS(file.path(pth,RF70))
subdf <- env_inputs[label %in% names(mod$forest$ncat)] #RF

### raster stack Tiles have already been created (see Revised_TWI_models.R)
tile_pth <- "H:/HOTR_models/RF_tiles"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

### Use the tiles to create a prediction raster
tile_prefix <- "RF_70pct_"
pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "randomForest")) %dopar% {
  raspred(tiles_index$tile[i_tile], tiles_pth = tile_pth, model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "RF")
}

snow::stopCluster(cl)

gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=14", "BIGTIFF=YES")
pred_tiles$filename <- str_replace(RF70, "rds", "tif")
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
#####

### BRT
mod <- readRDS(file.path(pth,BRT70))
subdf <- env_inputs[label %in% mod$var.names] #BRT
tile_prefix <- "BRT_70pct_"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "gbm")) %dopar% { #dismo??
  raspred(tiles_index$tile[i_tile], tiles_pth = tile_pth, model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "BRT")
}

snow::stopCluster(cl)

gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=14", "BIGTIFF=YES")
pred_tiles$filename <- str_replace(BRT70, "rds", "tif")
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
#####

### GLM
mod <- readRDS(file.path(pth,GLM70))
subdf <- env_inputs[label %in% names(mod@frame)] #GLMM
tile_pth <- "H:/HOTR_models/GLM_tiles"

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

### Use the tiles to create a prediction raster
tile_prefix <- "GLMER_70pct_"
pred_tiles <- foreach(i_tile = 1:nrow(tiles_index),
                      .packages = c("raster", "lme4")) %dopar% {
  raspred(tiles_index$tile[i_tile], tiles_pth = tile_pth, model = mod,
          lyrnames = subdf$label, oname = tile_prefix, modtype = "GLM")
}

snow::stopCluster(cl)

gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=14", "BIGTIFF=YES")
pred_tiles$filename <- str_replace(GLM70, "rds", "tif")
pred_tiles$format <- "GTiff"
pred_tiles$options <- gopts

out <- do.call(raster::merge, pred_tiles)
#####

### Now evaluate models with Test data and compare ###
library(dismo)
library(ROCR)
library(ggplot2)
source("Performance_functions.r")

RF_mod <- readRDS(RF70)
GLM_mod <- readRDS(GLM70)
BRT_mod <- readRDS(BRT70)

test_data <- fread("testing_data11192020.csv", header=TRUE, sep=",")
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
testBRT <- predict(BRT_mod, newdata=testdat, n.trees=25083, type="response")
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

fwrite(testOut, file = "Testing15pct_Dec012020.csv")
fullcv <- fread("CrossValidation_metricsNov23.csv", header=TRUE, sep=",")
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

fwrite(testOut, file = "Testing15pct_senspec_Dec012020.csv")
fullcv <- fread("CrossValidation_metrics_SenSpec_Nov23.csv", header=TRUE, sep=",")
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
