##############################################################################
# Creating and evaluating a weighted-average ensemble model
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 06/22/2020
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
library(raster)
library(ggplot2)
library(ROCR)
library(randomForest)
library(dismo)
library(lme4)
library(stringr)

source("Performance_functions.r")
pth <- "H:/HOTR_models/"
setwd(pth)

### Calculate weights ###
# See kfold.r for creation of cross-validation metrics
fullcv <- fread(paste0(pth, "CrossValidation_metricsMay14.csv"),
                header=TRUE, sep=",")
lstmodl <- c("GLM_May14", "RF_May12", "BRT_May13")

cvmean <- fullcv[, lapply(.SD, mean),
                 .SDcols=c("AUC","TSS","kappa","PCC","Sensitivity","Specificity","Threshold"),
                 by=.(model)]
tpose <- transpose(cvmean, make.names = "model")
# weights are based on mean of the first 6 metrics
wght <- tpose[1:6,][,lapply(.SD, mean),.SDcols=lstmodl]
stwght <- wght[,lapply(.SD, function(x){x/sum(wght)}),.SDcols=lstmodl]

# stwght:
# GLM_May14  RF_May12 BRT_May13
# 1: 0.3019841 0.3537902 0.3442257
######

### Create ensemble geoTIFF ###
BRTras <- raster(paste0(pth,"BRT_gbmholdout_50ktMay13.tif"))
GLMras <- raster(paste0(pth,"GLMER_May14_70pct.tif"))
RFras <- raster(paste0(pth,"RF_4800trees_May12.tif"))

Allras <- brick(GLMras, RFras, BRTras)
outname <- "ensemble_wave_May26_70pct.tif"
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=4", "BIGTIFF=YES")

beginCluster(6, type="SOCK")
ENras <- clusterR(Allras, fun=function(x){
  (x[[1]] * stwght[[1]]) + (x[[2]] * stwght[[2]]) +
    (x[[3]] * stwght[[3]])},
  export="stwght", filename=outname,
  format="GTiff", options=gopts)
endCluster()

######

### Evaluate ensemble with 15% Testing data ###
test_data <- fread("test_data_15correct.csv", header=TRUE, sep=",")
xy15 <- setnames(test_data[, 37:38], c("coords.x1","coords.x2"),
                 c("X","Y"))
EN15vals <- extract(ENras, xy15)
evEN <- data.table(prob=EN15vals, resp=test_data[, response])
outEN <- testEval(evEN, "EN_test15pct")

# Combine and graph all metrics
testOut <- tibble::tribble(~Name, ~AUC, ~TSS, ~err_rate, ~kappa, ~PCC, ~Sensitivity, ~Specificity, ~Threshold,
                           outEN[[1]],outEN[[2]],outEN[[3]],outEN[[4]],outEN[[5]],outEN[[6]],outEN[[7]],outEN[[8]],outEN[[9]])

# See kfold.r for creation of testing metrics
testIn <- fread("Testing15pct_3models.csv", header=TRUE, sep=",")

fulltest <- rbind(testIn, testOut)
fwrite(fulltest, file = paste0(pth, "Testing15pct_evalMetrics.csv"))

tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
toadd <- melt(fulltest, id.vars = "Name", measure.vars = c(2:9), variable.name = "metric")

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
######

### Create and Save the comprehensive Testing Prediction object ###
# This is so I can plot ROC curves and similar. See Ensemble_Review.Rmd
RF_mod <- readRDS(paste0(pth,"RF_4800trees_May12.rds"))
GLM_mod <- readRDS(paste0(pth,"GLMER_May14_70pct.rds"))
BRT_mod <- readRDS(paste0(pth,"BRT_gbmholdout_50ktMay13.rds"))

bestvars_T <- fread(paste0(pth, "vars_to_keep_Mar02.csv"), header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]

# RF
allcols <- c("response", varT)
testdat <- test_data[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]
testRF <- predict(RF_mod, newdata=testdat, type="prob", norm.votes=TRUE)
evRF <- data.table(prob=testRF[,2], resp=testdat[, response])

# BRT
testdat <- test_data[, ..allcols][, nlcd := as.factor(nlcd)]
testBRT <- predict(BRT_mod, newdata=testdat, n.trees=20490, type="response")
evBRT <- data.table(prob=testBRT, resp=testdat[, response])

# GLMM
bestvars_G <- fread(paste0(pth,"vars_to_keep_Feb25.csv"), header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
testdat <- test_data[, ..allcols][, Grid_ID := as.factor(Grid_ID)]
testGLM <- predict(GLM_mod, newdata=testdat, type="response", allow.new.levels=TRUE)
evGLM <- data.table(prob=testGLM, resp=testdat[, response])

predAll <- prediction(data.frame("GLMprob"=evGLM[,1],
                                 "RFprob"=evRF[,1],
                                 "BRTprob"=evBRT[,1],
                                 "ENprob"=evEN[,1]),
                      data.frame("GLMresp"=evRF[,2],
                                 "RFresp"=evRF[,2],
                                 "BRTresp"=evRF[,2],
                                 "ENresp"=evRF[,2]))

saveRDS(predAll, file = paste0(pth, "predictions15pct_allMod.rds"))
######