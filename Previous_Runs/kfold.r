##############################################################################
# k-fold crossvalidation of models, stratifying by Grid_ID
# Requires functions in 'Performance_functions.r'
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 11/23/2020
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

source("Performance_functions.r")

pth <- "H:/HOTR_models"
setwd(pth)
set.seed(5498)

## Ran this first with Sensitivity = 0.95, and numbers pretty low,
## so running a second time with Sensitivity=Specificity just to see.

train_data <- fread("training_kfold11192020.csv", header=TRUE, sep=",")
cvdat <- train_data[, Grid_ID := as.factor(Grid_ID)]

## **GLM** ##
ncores <- 10
bestvars_G <- fread("vars_to_keep_Feb25.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(cvdat[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "I(GDD5^2)", "TWI:ppt_ws", "TWI:ppt_sf", "GDD5:clay",
             "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

allcols <- c("KFsplit", "response", "Grid_ID", varG)
dat <- cvdat[, ..allcols]

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

cv_metrics <- foreach(i = 1:10,
                      .combine = "rbind",
                      .packages = c("data.table", "lme4", "PresenceAbsence", "ROCR")) %dopar% {
                        datrain <- dat[KFsplit != i]
                        datest <- dat[KFsplit == i]
                        fit <- glmer(f_glm, data = datrain,
                                     family = "binomial",
                                     control = glmerControl(optimizer = "bobyqa"),
                                     nAGQ = 0,
                                     subset = NULL,
                                     weights = NULL,
                                     na.action = na.omit,
                                     offset = NULL,
                                     mustart = NULL,
                                     etastart = NULL)
                        evfit <- predict(fit, newdata=datest, type="response",
                                         allow.new.levels=TRUE)
                        evdt <- data.table(prob=evfit, resp=datest[, response])
                        testEval(evdt, paste("fold", as.character(i)), cutype="senspec")
                      }

snow::stopCluster(cl)

cvout <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                    TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                    kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                    Sensitivity=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                    Threshold=unlist(cv_metrics[,9]))

cvout
cvout$model <- "GLM_Nov19"

## **RF** ##
# Tune first
bestvars_T <- fread("vars_to_keep_Mar02.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]
allcols <- c("KFsplit", "response", varT)
dat <- cvdat[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]
# 1) Tune mtry
sub0 <- dat %>% dplyr::filter(response == 0) %>% sample_frac(0.3)
sub1 <- dat %>% dplyr::filter(response == 1) %>% sample_frac(0.3)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varT], y=subtune[, response],
            ntreeTry = 501,
            stepFactor = 1.5,
            improve=0.01)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
# mt <- 9 same as last run

ntrees <- 4800
mt <- 9
#runs out of memory with anything higher than 2 cores
ncores <- 2

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

cv_metrics <- foreach(i = 1:10,
                      .combine = "rbind",
                      .packages = c("data.table", "randomForest", "PresenceAbsence", "ROCR")) %dopar% {
                        datrain <- dat[KFsplit != i]
                        datest <- dat[KFsplit == i]
                        fit <- randomForest(datrain[, 3:15],
                                            y=datrain[, response],
                                            ntree=ntrees,
                                            mtry=mt,
                                            replace = TRUE,
                                            norm.votes = TRUE)
                        evfit <- predict(fit, newdata=datest, type="prob",
                                         norm.votes=TRUE)
                        evdt <- data.table(prob=evfit[,2], resp=datest[, response])
                        testEval(evdt, paste("fold", as.character(i)), cutype="senspec")
                      }

snow::stopCluster(cl)

cvout2 <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                     TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                     kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                     Sensitivity=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                     Threshold=unlist(cv_metrics[,9]))

cvout2$model <- "RF_Nov23"
addcv <- rbind(cvout, cvout2)

# **BRT**
ncores <- 8
ntrees <- 20490
dat <- cvdat[, ..allcols][, nlcd := as.factor(nlcd)]

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

cv_metrics <- foreach(i = 1:10,
                      .combine = "rbind",
                      .packages = c("data.table", "gbm", "dismo", "PresenceAbsence", "ROCR")) %dopar% {
                        datrain <- dat[KFsplit != i]
                        datest <- dat[KFsplit == i]
                        fit <- gbm.fixed(datrain,
                                         gbm.x = c(3:15),
                                         gbm.y = 2,
                                         tree.complexity = 2,
                                         verbose = FALSE,
                                         learning.rate = 0.11,
                                         n.trees = ntrees,
                                         bag.fraction = 0.5)
                        evfit <- predict(fit, newdata=datest,
                                         n.trees=ntrees,
                                         type="response")
                        evdt <- data.table(prob=evfit, resp=datest[, response])
                        testEval(evdt, paste("fold", as.character(i)), cutype="senspec")
                      }

snow::stopCluster(cl)

cvout3 <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                     TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                     kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                     Sensitivity=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                     Threshold=unlist(cv_metrics[,9]))

cvout3$model <- "BRT_Nov23"

fullcv <- rbind(addcv, cvout3)
fwrite(fullcv, file = "CrossValidation_metrics_SenSpec_Nov23.csv")
tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
ggplot(data=tograph, aes(x=factor(metric), y=value)) +
  geom_boxplot(aes(fill=factor(model))) +
  ylim(c(0,1)) +
  labs(x="Performance Metric", y="",
       title="10-fold Cross-Validation of 70% Training, Sensitivity=Specificity")
####################################stopped here##############################

### Now evaluate models with Test data and compare ###
RF_mod <- readRDS(paste0(pth,"RF_4800trees_May12.rds"))
GLM_mod <- readRDS(paste0(pth,"GLMER_May14_70pct.rds"))
BRT_mod <- readRDS(paste0(pth,"BRT_gbmholdout_50ktMay13.rds"))

test_data <- fread(paste0(pth, "test_data_15correct.csv"), header=TRUE, sep=",")
bestvars_T <- fread(paste0(pth, "vars_to_keep_Mar02.csv"), header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]

# RF
allcols <- c("response", varT)
testdat <- test_data[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]
testRF <- predict(RF_mod, newdata=testdat, type="prob", norm.votes=TRUE)
evRF <- data.table(prob=testRF[,2], resp=testdat[, response])
outRF <- testEval(evRF, "RF_test15pct")

# BRT
testdat <- test_data[, ..allcols][, nlcd := as.factor(nlcd)]
testBRT <- predict(BRT_mod, newdata=testdat, n.trees=20490, type="response")
evBRT <- data.table(prob=testBRT, resp=testdat[, response])
outBRT <- testEval(evBRT, "BRT_test15pct")

# GLMM
allcols <- c("response", "Grid_ID", varG)
testdat <- test_data[, ..allcols][, Grid_ID := as.factor(Grid_ID)]
testGLM <- predict(GLM_mod, newdata=testdat, type="response", allow.new.levels=TRUE)
evGLM <- data.table(prob=testGLM, resp=testdat[, response])
outGLM <- testEval(evGLM, "GLM_test15pct")

testOut <- tibble::tribble(~Name, ~AUC, ~TSS, ~err_rate, ~kappa, ~PCC, ~Sensitivity, ~Specificity, ~Threshold,
                           outGLM[[1]],outGLM[[2]],outGLM[[3]],outGLM[[4]],outGLM[[5]],outGLM[[6]],outGLM[[7]],outGLM[[8]],outGLM[[9]],
                           outRF[[1]],outRF[[2]],outRF[[3]],outRF[[4]],outRF[[5]],outRF[[6]],outRF[[7]],outRF[[8]],outRF[[9]],
                           outBRT[[1]],outBRT[[2]],outBRT[[3]],outBRT[[4]],outBRT[[5]],outBRT[[6]],outBRT[[7]],outBRT[[8]],outBRT[[9]])

fwrite(testOut, file = paste0(pth, "Testing15pct_3models.csv"))
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
