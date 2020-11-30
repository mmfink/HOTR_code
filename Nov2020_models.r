#############################################################################
# Run models with new (Nov 2020) data sampling scheme.
# See redo_sampling.r and kfold.r for what came before.
# Requires functions in 'Performance_functions.r'
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 11/24/2020
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

train_data <- fread("training_kfold11192020.csv", header=TRUE, sep=",")
indata <- train_data[, Grid_ID := as.factor(Grid_ID)]
rnam <- "Nov24_2020"

## **GLM** ##
outnam <- paste0("GLMER_", rnam)

bestvars_G <- fread("vars_to_keep_Feb25.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "I(GDD5^2)", "TWI:ppt_ws", "TWI:ppt_sf", "GDD5:clay",
             "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

allcols <- c("response", "Grid_ID", varG)
dat <- indata[, ..allcols]

glmm.fit <- glmer(f_glm, data = indata,
                  family = "binomial",
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 0,
                  subset = NULL,
                  weights = NULL,
                  na.action = na.omit,
                  offset = NULL,
                  mustart = NULL,
                  etastart = NULL)

saveRDS(glmm.fit, paste0(outnam, ".rds"))
summary(glmm.fit)

# Evaluate with Tweak data #
tweak_data <- fread("tweaking_data11192020.csv", header=TRUE, sep=",")
twkdat <- tweak_data[, ..allcols][, Grid_ID := as.factor(Grid_ID)]
twkGLM <- predict(glmm.fit, newdata=twkdat, type="response", allow.new.levels=TRUE)
evGLM <- data.table(prob=twkGLM, resp=twkdat[, response])
outGLM <- testEval(evGLM, outnam)

## **RF** ##
ncores <- 16
outnam <- paste0("RF_4800trees_", rnam)
ntrees <- 4800
mt <- 9
treeSubs <- ntrees/ncores
bestvars_T <- fread("vars_to_keep_Mar02.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]
allcols <- c("response", varT)
indata <- train_data[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit1 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, 2:14],
                                  y=indata[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }
saveRDS(rf.fit1, paste0(outnam, ".rds"))
impvals <- importance(rf.fit1)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit1)

snow::stopCluster(cl)

# Evaluate with Tweak data #
twkdat <- tweak_data[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]
twkRF <- predict(rf.fit1, newdata=twkdat, type="prob", norm.votes=TRUE)
evRF <- data.table(prob=twkRF[,2], resp=twkdat[, response])
outRF <- testEval(evRF, "RF_4800t_mtry9")

## **BRT** ##
outnam <- paste0("BRT_gbmholdout_50kt", rnam)
indata <- train_data[, ..allcols][, nlcd := as.factor(nlcd)]

brt.fit <- gbm.holdout(indata, gbm.x = c(2:14), gbm.y = 1,
                        learning.rate = 0.11, tree.complexity = 2,
                        n.trees = 20000, add.trees = 5000,
                        max.trees = 50000,
                        train.fraction = 0.6)

saveRDS(brt.fit, paste0(outnam, ".rds"))
summary(brt.fit)
gbm.interactions(brt.fit, mask.object = brt.fit)
# stopped at 25083 trees, ~5k more than last run.

# Evaluate with Tweak data #
twkdat <- tweak_data[, ..allcols][, nlcd := as.factor(nlcd)]
twkBRT <- predict(brt.fit, newdata=twkdat, n.trees=25083, type="response")
evBRT <- data.table(prob=twkBRT, resp=twkdat[, response])
outBRT <- testEval(evBRT, "BRT_25083trees")

# Combine and save evaluation metrics
twkOut <- tibble::tribble(~Name, ~AUC, ~TSS, ~err_rate, ~kappa, ~PCC, ~Sensitivity, ~Specificity, ~Threshold,
                           outGLM[[1]],outGLM[[2]],outGLM[[3]],outGLM[[4]],outGLM[[5]],outGLM[[6]],outGLM[[7]],outGLM[[8]],outGLM[[9]],
                           outRF[[1]],outRF[[2]],outRF[[3]],outRF[[4]],outRF[[5]],outRF[[6]],outRF[[7]],outRF[[8]],outRF[[9]],
                           outBRT[[1]],outBRT[[2]],outBRT[[3]],outBRT[[4]],outBRT[[5]],outBRT[[6]],outBRT[[7]],outBRT[[8]],outBRT[[9]])

fwrite(twkOut, file = "Tweak15pct_Nov242020.csv")

# Graph a comparison of Tweak vs. Cross-Validation metrics
fullcv <- fread("CrossValidation_metricsNov23.csv", header=TRUE, sep=",")
tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
toadd <- melt(twkOut, id.vars = "Name", measure.vars = c(2:9), variable.name = "metric")

ggplot(data=tograph, aes(x=factor(metric), y=value)) +
  geom_boxplot(aes(fill=factor(model))) +
  ylim(c(0,1)) +
  labs(x="Performance Metric", y="",
       title="10-fold Cross-Validation Plus Tweaking Data, Sensitivity=0.95") +
  geom_jitter(data=toadd, aes(x=factor(metric),y=value,
                             color=factor(Name)),
             size=3, pch=10) +
  guides(fill=guide_legend(title = "10fold Training"),
         color=guide_legend(title = "Tweaking data"))

# Re-do with Sensitivity=Specificty metrics, just to see
outGLM <- testEval(evGLM, "GLMER_Nov24_2020", cutype = "senspec")
outRF <- testEval(evRF, "RF_4800t_mtry9", cutype = "senspec")
outBRT <- testEval(evBRT, "BRT_25083trees", cutype = "senspec")

twkOut <- tibble::tribble(~Name, ~AUC, ~TSS, ~err_rate, ~kappa, ~PCC, ~Sensitivity, ~Specificity, ~Threshold,
                           outGLM[[1]],outGLM[[2]],outGLM[[3]],outGLM[[4]],outGLM[[5]],outGLM[[6]],outGLM[[7]],outGLM[[8]],outGLM[[9]],
                           outRF[[1]],outRF[[2]],outRF[[3]],outRF[[4]],outRF[[5]],outRF[[6]],outRF[[7]],outRF[[8]],outRF[[9]],
                           outBRT[[1]],outBRT[[2]],outBRT[[3]],outBRT[[4]],outBRT[[5]],outBRT[[6]],outBRT[[7]],outBRT[[8]],outBRT[[9]])

fwrite(twkOut, file = "Tweak15pct_Nov242020_senspec.csv")
fullcv <- fread("CrossValidation_metrics_SenSpec_Nov23.csv", header=TRUE, sep=",")
tograph <- melt(fullcv, id.vars = c("cv", "model"), measure.vars = c(2:9), variable.name = "metric")
toadd <- melt(twkOut, id.vars = "Name", measure.vars = c(2:9), variable.name = "metric")

ggplot(data=tograph, aes(x=factor(metric), y=value)) +
  geom_boxplot(aes(fill=factor(model))) +
  ylim(c(0,1)) +
  labs(x="Performance Metric", y="",
       title="10-fold Cross-Validation Plus Tweaking Data, Sensitivity=Specificity") +
  geom_jitter(data=toadd, aes(x=factor(metric),y=value,
                             color=factor(Name)),
             size=3, pch=10) +
  guides(fill=guide_legend(title = "10fold Training"),
         color=guide_legend(title = "Tweaking data"))
