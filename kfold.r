##############################################################################
# k-fold crossvalidation of models, stratifying by Grid_ID
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 05/19/2020
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
library(PresenceAbsence) # for Kappa and PCC
library(ggplot2)

opt.cut <- function(perf, pred, type="senspec"){
  # Returns [Sensitivity, Specificity, Threshold]
  # adapted from:
  #https://www.r-bloggers.com/a-small-introduction-to-the-rocr-package/
  #type="senspec" calculates sensitivity = specificity threshold
  #type="sens95" caclulates sensitivity = 0.95 threshold

  mapply(FUN=function(x, y, p){
    if(type == "senspec"){
      d = (x - 0)^2 + (y-1)^2
      ind = which(d == min(d))
      c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
        cutoff = p[[ind]])
    } else if(type == "sens95"){
      d = which(y >= 0.95)
      ind = which(y == min(y[d]))
      if(length(ind > 1)){ # if ties, take first
        idx <- ind[1]
      } else {
        idx <- ind
      }
      c(sensitivity = y[[idx]], specificity = 1-x[[idx]],
        cutoff = p[[idx]])
    }
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

testEval <- function(evobj, run_name){
  #evobj == data.table(prob, resp)
  pred.mod <- prediction(evobj$prob, evobj$resp)
  perf.roc <- performance(pred.mod, "tpr", "fpr")
  perf.err <- performance(pred.mod, "err")
  AUC <- performance(pred.mod, measure = "auc")@y.values[[1]]
  sst <- opt.cut(perf.roc, pred.mod, type = "sens95")
  TSS <- sst[1] + sst[2] - 1
  err_rate <- perf.err@y.values[[1]][which(perf.err@x.values[[1]] == sst[3])]

  #Presence.Absence insists on its own data structure
  padat <- data.frame(ID=seq(1:length(pred.mod@predictions[[1]])),
                      OBS=as.numeric(levels(pred.mod@labels[[1]]))[pred.mod@labels[[1]]],
                      PRED=pred.mod@predictions[[1]])

  pacmx <- cmx(padat, threshold = sst[3])
  kappa <- Kappa(pacmx, st.dev = FALSE)
  PCC <- pcc(pacmx, st.dev = FALSE)

  return(list(run_name, AUC, TSS, err_rate, kappa, PCC, sst[1], sst[2], sst[3]))
}

pth <- "H:/HOTR_models"
setwd(pth)
set.seed(135)
folds <- 10

train_data <- fread(paste0(pth, "train_data_70pct.csv"), header=TRUE, sep=",")
indata <- train_data[, Grid_ID := as.factor(Grid_ID)]

# Assign a K-fold Split number to each Grid_ID.
# NOTE while this technique selects whole grids, it selects them randomly
# across the region and not stratified geographically.
# Also note that the number of presence vs absence points in each fold
# will not be strictly equal.

kf <- kfold(levels(indata$Grid_ID), k=folds)
names(kf) <- levels(indata$Grid_ID)
kfdt <- data.table(Grid_ID=names(kf), KFsplit=kf, stringsAsFactors = TRUE)

# Add the Split column to the dataset
cvdat <- train_data[kfdt, on=.(Grid_ID=Grid_ID)]

# Let's run some summary stats on the folds to check for weirdness.
inroll <- rollup(cvdat, j = c(list(cnt=.N), lapply(.SD, mean)),
                 by = c("KFsplit", "response"), .SDcols = c(4:25))
subresp <- inroll[!is.na(response)]
setorder(subresp, KFsplit, response)
fwrite(subresp, file = "KFSplit_byResponse_May14.csv")

qplot(x=factor(KFsplit), data=cvdat, geom="bar", fill=factor(response),
      xlab="K-fold Split", ylab="Count",
      main="Presence/Absence per Split")

tograph <- melt(subresp, id.vars = c("KFsplit", "response"),
                measure.vars = c(4:25), variable.name = "metric")

qplot(factor(response), value, data=tograph,
      geom="boxplot", facets=.~factor(metric), color=factor(response),
      xlab="K-fold mean value spread for each metric by response")
####

# **GLM**
ncores <- 8
bestvars_G <- fread("vars_to_keep_Feb25.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("KFsplit", "response", "Grid_ID", varG)
dat <- cvdat[, ..allcols]

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

cv_metrics <- foreach(i = 1:10,
                      .combine = "rbind",
                      .packages = c("data.table", "lme4", "PresenceAbsence", "ROCR")) %dopar% {
                        datrain <- dat[KFsplit != i]
                        datest <- dat[KFsplit == i]
                        # Note: error in lme4 v1.1-23 making optional args required
                        # so glad I updated. :-/
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
                        testEval(evdt, paste("fold", as.character(i)))
                      }

snow::stopCluster(cl)

cvout <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                    TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                    kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                    Sensitivity=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                    Threshold=unlist(cv_metrics[,9]))

cvout
cvout$model <- "GLM_May14"

# **RF**
ntrees <- 4800
mt <- 9
#runs out of memory with anything higher than 2 cores
ncores <- 2
bestvars_T <- fread("vars_to_keep_Mar02.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]
allcols <- c("KFsplit", "response", varT)
dat <- cvdat[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]

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
                        testEval(evdt, paste("fold", as.character(i)))
                      }

snow::stopCluster(cl)

cvout2 <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                     TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                     kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                     Sensitivty=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                     Threshold=unlist(cv_metrics[,9]))

cvout2$model <- "RF_May12"
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
                        testEval(evdt, paste("fold", as.character(i)))
                      }

snow::stopCluster(cl)

cvout3 <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                     TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                     kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                     Sensitivty=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                     Threshold=unlist(cv_metrics[,9]))

cvout3$model <- "BRT_May13"

fullcv <- rbind(addcv, cvout3)
fwrite(fullcv, file = "CrossValidation_metricsMay14.csv")
tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
ggplot(data=tograph, aes(x=factor(metric), y=value)) +
  geom_boxplot(aes(fill=factor(model))) +
  ylim(c(0,1)) +
  labs(x="Performance Metric", y="",
       title="10-fold Cross-Validation of 70% Training, Sensitivity=0.95")
######

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
