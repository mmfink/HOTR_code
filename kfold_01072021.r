##############################################################################
# k-fold crossvalidation of models, stratifying by Grid_ID
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 01/12/2021
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

source("Performance_functions.r")

pth <- "H:/HOTR_models"
setwd(pth)
set.seed(7981)

## Run this first with Sensitivity = 0.95
train_data <- fread("training_kfold01072021.csv", header=TRUE, sep=",")
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
                        testEval(evdt, paste("fold", as.character(i)), cutype="sens95")
                      }

snow::stopCluster(cl)

cvout <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                    TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                    kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                    Sensitivity=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                    Threshold=unlist(cv_metrics[,9]))

cvout
cvout$model <- "GLM_Jan072021"

## **RF** ##
# Tune first
bestvars_T <- fread("vars_to_keep_Mar02.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]
allcols <- c("KFsplit", "response", varT)
dat <- train_data[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]
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
# mt <- 6, lower than previous runs

ntrees <- 4200
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
                        testEval(evdt, paste("fold", as.character(i)), cutype="sens95")
                      }

snow::stopCluster(cl)

cvout2 <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                     TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                     kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                     Sensitivity=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                     Threshold=unlist(cv_metrics[,9]))

cvout2$model <- "RF_Jan122021"
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
                        testEval(evdt, paste("fold", as.character(i)), cutype="sens95")
                      }

snow::stopCluster(cl)

cvout3 <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                     TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                     kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                     Sensitivity=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                     Threshold=unlist(cv_metrics[,9]))

cvout3$model <- "BRT_Jan072021"

fullcv <- rbind(addcv, cvout3)
fwrite(fullcv, file = "CrossValidation_metrics_Sens95_01072021.csv")
tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
ggplot(data=tograph, aes(x=factor(metric), y=value)) +
  geom_boxplot(aes(fill=factor(model))) +
  ylim(c(0,1)) +
  labs(x="Performance Metric", y="",
       title="10-fold Cross-Validation of 70% Training, Sensitivity=0.95")
