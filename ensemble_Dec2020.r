##############################################################################
# Create an Ensemble weighted average model from the three models.
# Then evaluate Ensemble with the withheld Testing dataset.
# Models created in Nov2020_models.r and predicted to raster in predict_70pct.r
# Requires functions in Performance_functions.r
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 12/02/2020
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
library(PresenceAbsence) # for Kappa and PCC
library(stringr)
library(doSNOW)

pth <- "H:/HOTR_models/"
setwd(pth)
source("Performance_functions.r")

BRT70 <- "BRT_gbmholdout_50ktNov24_2020.rds"
GLM70 <- "GLMER_Nov24_2020.rds"
RF70 <- "RF_4800trees_Nov24_2020.rds"

ncores <- 12
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=4", "BIGTIFF=YES")

RF_mod <- readRDS(RF70)
GLM_mod <- readRDS(GLM70)
BRT_mod <- readRDS(BRT70)

test_data <- fread("testing_data11192020.csv", header=TRUE, sep=",")
fullcv <- fread("CrossValidation_metricsNov23.csv", header=TRUE, sep=",")
lstmodl <- c("GLM_Nov19", "RF_Nov23", "BRT_Nov23")

cvmean <- fullcv[, lapply(.SD, mean),
                    .SDcols=c("AUC","TSS","kappa","PCC","Sensitivity","Specificity","Threshold"),
                    by=.(model)]
tpose <- transpose(cvmean, make.names = "model")
wght <- tpose[1:6,][,lapply(.SD, mean),.SDcols=lstmodl]
stwght <- wght[,lapply(.SD, function(x){x/sum(wght)}),.SDcols=lstmodl]

### Create ensemble geoTIFF ###
BRTras <- raster("BRT_gbmholdout_50ktNov24_2020.tif")
GLMras <- raster("GLMER_Nov24_2020.tif")
RFras <- raster("RF_4800trees_Nov24_2020.tif")

Allras <- brick(GLMras, RFras, BRTras)
outname <- "ensemble_wave_Dec012020.tif"

beginCluster(ncores, type="SOCK")
ENras <- clusterR(Allras, fun=function(x){
  (x[[1]] * stwght[[1]]) + (x[[2]] * stwght[[2]]) +
    (x[[3]] * stwght[[3]])},
  export="stwght", filename=outname, overwrite=TRUE,
  format="GTiff", options=gopts)
endCluster()

####

# Evaluate ensemble with Testing data
xy15 <- setnames(test_data[, 37:38], c("coords.x1","coords.x2"),
                 c("X","Y"))
EN15vals <- extract(ENras, xy15)
evEN <- data.table(prob=EN15vals, resp=test_data[, response])
outEN <- testEval(evEN, "EN_test15pct")

## Combine and graph
testEn <- tibble::tribble(~Name, ~AUC, ~TSS, ~err_rate, ~kappa, ~PCC, ~Sensitivity, ~Specificity, ~Threshold,
                           outEN[[1]],outEN[[2]],outEN[[3]],outEN[[4]],outEN[[5]],outEN[[6]],outEN[[7]],outEN[[8]],outEN[[9]])
modTest <- fread("Testing15pct_Dec012020.csv", header=TRUE, sep=",")
testOut <- rbind(modTest, testEn)

fwrite(testOut, file = "Testing15pct_allModels_Dec2020.csv")

# TODO - the ev... files were already in memory when I ran this. Should
# add the code to create them here for completeness' sake.
predAll <- prediction(data.frame("GLMprob"=evGLM[,1],
                                 "RFprob"=evRF[,1],
                                 "BRTprob"=evBRT[,1],
                                 "ENprob"=evEN[,1]),
                      data.frame("GLMresp"=evRF[,2],
                                 "RFresp"=evRF[,2],
                                 "BRTresp"=evRF[,2],
                                 "ENresp"=evRF[,2]))

saveRDS(predAll, file = "predictions15pct_allMod_Dec2020.rds")

tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
toadd <- melt(testOut, id.vars = "Name", measure.vars = c(2:9), variable.name = "metric")

m4Palette <- c("#cf6329", "#a921ff", "#00ba38", "#619cff")

ggplot(data=tograph, aes(x=factor(metric), y=value)) +
  geom_boxplot(aes(fill=factor(model))) +
  ylim(c(0,1)) +
  labs(x="Performance Metric", y="",
       title="10-fold Cross-Validation Plus Testing Data, Sensitivity=0.95") +
  geom_jitter(data=toadd, aes(x=factor(metric),y=value,
                             color=factor(Name)),
             size=3, pch=10) +
  scale_colour_manual(values=m4Palette) +
  guides(fill=guide_legend(title = "10fold Training"),
         color=guide_legend(title = "Testing data"))
