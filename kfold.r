##############################################################################
# k-fold crossvalidation of models, stratifying by Grid_ID
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 03/03/2020
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

#.libPaths("E:/MMF/R/win-library/3.6") #specific to the modeling server
library(stringr)
library(data.table)
library(dplyr)
library(doSNOW)
library(randomForest)
library(lme4)
library(dismo)
library(ROCR)
library(PresenceAbsence) # for Kappa and PCC
library(ggplot2)

opt.cut <- function(perf, pred){
  #https://www.r-bloggers.com/a-small-introduction-to-the-rocr-package/
  #calculates sensitivity = specificity threshold
  mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

testEval <- function(evobj, run_name){
  #evobj == data.table(prob, resp)
  pred.mod <- prediction(evobj$prob, evobj$resp)
  perf.roc <- performance(pred.mod, "tpr", "fpr")
  perf.err <- performance(pred.mod, "err")
  AUC <- performance(pred.mod, measure = "auc")@y.values[[1]]
  sst <- opt.cut(perf.roc,pred.mod)
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
rnam <- "Mar03" #not currently used
folds <- 10
ntrees <- 2000
ncores <- 6 #if modeling server: 20

train_data <- fread("training_data_seed655.csv", header=TRUE, sep=",")
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
indata <- indata[kfdt, on=.(Grid_ID=Grid_ID)]

## Let's run some summary stats on the folds to check for weirdness. ##
inroll <- rollup(indata, j = c(list(cnt=.N), lapply(.SD, mean)),
                 by = c("KFsplit", "response"), .SDcols = c(4:25))
subresp <- inroll[!is.na(response)]
setorder(subresp, KFsplit, response)
subresp

qplot(x=factor(KFsplit), data=indata, geom="bar", fill=factor(response),
      xlab="K-fold Split", ylab="Count",
      main="Presence/Absence per Split")

tograph <- melt(subresp, id.vars = c("KFsplit", "response"),
                measure.vars = c(4:25), variable.name = "metric")

qplot(factor(response), value, data=tograph,
      geom="boxplot", facets=.~factor(metric), color=factor(response),
      xlab="K-fold mean value spread for each metric by response")
####

### Model 1) GLMER_alltypes_Feb25 ###
bestvars_G <- fread("vars_to_keep_Feb25.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("KFsplit", "response", "Grid_ID", varG)
dat <- indata[, ..allcols]
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(dat[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

cv_metrics <- foreach(i = 1:10,
                      .combine = "rbind",
                      .packages = c("data.table", "lme4", "PresenceAbsence", "ROCR")) %dopar% {
                        datrain <- dat[KFsplit != i]
                        datest <- dat[KFsplit == i]
                        fit <- glmer(f_glm, data = datrain, family = "binomial",
                                     control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)
                        evfit <- predict(fit, newdata=datest, type="response",
                                         allow.new.levels=TRUE)
                        evdt <- data.table(prob=evfit, resp=datest[, response])
                        testEval(evdt, paste("fold", as.character(i)))
                      }

snow::stopCluster(cl)

# TODO: I create a results table, but then don't save the table anywhere.
cvout <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                    TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                    kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                    Sensitivty=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                    Threshold=unlist(cv_metrics[,9]))

tograph <- melt(cvout, id.vars = "cv", measure.vars = c(2:9), variable.name = "metric")
qplot(factor(metric), value, data=tograph, geom="boxplot",
      xlab = "Performance Metric", main = "10-fold Cross-Validation for GLMER_alltypes_Feb25")
######

### Model 2) GLMER_alltypes_Mar02 ###
# Same data but changes to formula
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "I(GDD5^2)", "TWI:ppt_ws", "TWI:ppt_sf", "GDD5:clay",
             "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

cv_metrics <- foreach(i = 1:10,
                      .combine = "rbind",
                      .packages = c("data.table", "lme4", "PresenceAbsence", "ROCR")) %dopar% {
                        datrain <- dat[KFsplit != i]
                        datest <- dat[KFsplit == i]
                        fit <- glmer(f_glm, data = datrain, family = "binomial",
                                     control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)
                        evfit <- predict(fit, newdata=datest, type="response",
                                         allow.new.levels=TRUE)
                        evdt <- data.table(prob=evfit, resp=datest[, response])
                        testEval(evdt, paste("fold", as.character(i)))
                      }

snow::stopCluster(cl)

# See TODO above
cvout <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                    TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                    kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                    Sensitivty=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                    Threshold=unlist(cv_metrics[,9]))

cvout
tograph <- melt(cvout, id.vars = "cv", measure.vars = c(2:9), variable.name = "metric")
qplot(factor(metric), value, data=tograph, geom="boxplot",
      xlab = "Performance Metric", main = "10-fold Cross-Validation for GLMER_alltypes_Mar02")

### Model 3) RF_alltypes_Mar02 ###
# Uses factor NLCD instead of the one-hot variables for land use

train_dataRF <- fread("training_data_seed655_forRF.csv", header=TRUE, sep=",")
indata <- train_dataRF[, Grid_ID := as.factor(Grid_ID)]

# Assign a K-fold Split number to each Grid_ID with this dataset.
kf <- kfold(levels(indata$Grid_ID), k=folds)
names(kf) <- levels(indata$Grid_ID)
kfdt <- data.table(Grid_ID=names(kf), KFsplit=kf, stringsAsFactors = TRUE)

# Add the Split column to the dataset
#Note that, with the same seed & number of records, splits are the same as for above dataset
indata <- indata[kfdt, on=.(Grid_ID=Grid_ID)]

bestvars_G <- fread("vars_to_keep_Mar02.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("KFsplit", "response", varG)
dat <- indata[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

# FIXME: The below very quickly uses all 147 GB of memory on modeling server
cv_metrics <- foreach(i = 1:10,
                      .combine = "rbind",
                      .packages = c("data.table", "randomForest", "PresenceAbsence", "ROCR")) %dopar% {
                        datrain <- dat[KFsplit != i]
                        datest <- dat[KFsplit == i]
                        fit <- randomForest(datrain[, 3:15],
                                  y=datrain[, response],
                                  ntree=ntrees,
                                  mtry=4,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                        evfit <- predict(fit, newdata=datest, type="prob",
                                         norm.votes=TRUE)
                        evdt <- data.table(prob=evfit[,2], resp=datest[, response])
                        testEval(evdt, paste("fold", as.character(i)))
                      }

snow::stopCluster(cl)

# See TODO above
cvout <- data.table(cv=unlist(cv_metrics[,1]), AUC=unlist(cv_metrics[,2]),
                    TSS=unlist(cv_metrics[,3]), err_rate=unlist(cv_metrics[,4]),
                    kappa=unlist(cv_metrics[,5]), PCC=unlist(cv_metrics[,6]),
                    Sensitivty=unlist(cv_metrics[,7]), Specificity=unlist(cv_metrics[,8]),
                    Threshold=unlist(cv_metrics[,9]))

cvout
tograph <- melt(cvout, id.vars = "cv", measure.vars = c(2:9), variable.name = "metric")
qplot(factor(metric), value, data=tograph, geom="boxplot",
      xlab = "Performance Metric", main = "10-fold Cross-Validation for RF_alltypes_Mar02")