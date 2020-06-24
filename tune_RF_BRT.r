##############################################################################
# Tuning Random Forest and Boosted Regression Trees call parameters
# Purpose of this code: Models trained on 75% of the data
# showed that RF and BRT 10-fold cross-validation performance metrics
# were higher than those of the GLMM, but the pattern was strongly reversed
# when the performance metrics were calculated on the 25% withheld Test data.
#
# So, this code uses 70% of the full data set as Training data. The calling
# parameters for RF and BRT were then tuned using a 'Tweaking' data set of 15%
# of the withheld data to try to minimize over-fitting. The final 15% of data
# was used to Test all 3 models. See kfold.r for that.
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

library(stringr)
library(data.table)
library(plyr)
library(doSNOW)
library(randomForest)
library(lme4)
library(dismo)
library(gbm3)
library(ROCR)
library(ggplot2)

source("Performance_functions.r")
pth <- "H:/HOTR_models/"
setwd(pth)
set.seed(6785) # for the 70/15/15 run

test_pct <- 0.30 # proportion of response data to withhold from Training
full_data <- fread(paste0(pth, "training_data_s655Apr07.csv"), header=TRUE, sep=",")

# split out test and train records
dat_abs <- full_data[response == 0,]
dat_pres <- full_data[response == 1,]
r <- length(dat_abs[, response])
test_idx <- sample(seq_len(r), size = floor(test_pct * r))
test_abs <- dat_abs[test_idx, ]
train_abs <- dat_abs[-test_idx, ]
test_pres <- dat_pres[test_idx, ]
train_pres <- dat_pres[-test_idx, ]
train_data <- bind_rows(train_pres, train_abs)

# further split testing in half
rows <- nrow(test_pres)
tweak_data <- bind_rows(test_pres[1:(rows/2),], test_abs[1:(rows/2),])

#FIXME: this didn't work as intended - included all rows, but half were NA
# wound up manually fixing and checked, no duplication between Tweak and Test. Yay!
test_data <- bind_rows(test_pres[(rows/2)+1:rows,], test_abs[(rows/2)+1:rows,])

fwrite(train_data, file = "train_data_70pct.csv")
fwrite(test_data, file = "test_data_15pct.csv")
fwrite(tweak_data, file = "tweak_data_15pct.csv")

### Random Forest tuning ###
bestvars_T <- fread("vars_to_keep_Mar02.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]
allcols <- c("KFsplit", "response", varT)

indata <- train_data[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]

# 1) re-tune mtry with larger sub-sample and more trees
sub0 <- indata %>% dplyr::filter(response == 0) %>% sample_frac(0.3)
sub1 <- indata %>% dplyr::filter(response == 1) %>% sample_frac(0.3)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varT], y=subtune[, response],
            ntreeTry = 501,
            stepFactor = 1.5,
            improve=0.01)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
# mt <- 9 # original was 4

# 2) increase trees
rnam <- "May12"
outnam <- paste0("RF_3200trees_", rnam)

ntrees <- 3200 # original was 2000
ncores <- 16 # needs to evenly divide the trees
treeSubs <- ntrees/ncores

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit1 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, 3:15],
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
####

### Evaluate with Tweak data ###
allcols <- c("response", varT)
twkdat <- tweak_data[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]
twkRF <- predict(rf.fit1, newdata=twkdat, type="prob", norm.votes=TRUE)
evRF <- data.table(prob=twkRF[,2], resp=twkdat[, response])
outRF <- testEval(evRF, "RF_3200t_mtry9") # metrics improved over original

## keep increasing trees
outnam <- paste0("RF_4800trees_", rnam)

ntrees <- 4800
treeSubs <- ntrees/ncores

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit2 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, 3:15],
                                  y=indata[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }

saveRDS(rf.fit2, paste0(outnam, ".rds"))
impvals <- importance(rf.fit2)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit2)

snow::stopCluster(cl)
###

### Evaluate with Tweak data ###
twkRF2 <- predict(rf.fit2, newdata=twkdat, type="prob", norm.votes=TRUE)
evRF2 <- data.table(prob=twkRF2[,2], resp=twkdat[, response])
outRF2 <- testEval(evRF2, "RF_4800t_mtry9") # *very* minor improvements over 3200 trees
###

# 3) so now let's add nodesize to the mix, keeping trees at 4800 and mtry at 9
outnam <- paste0("RF_4800t_nodes3", rnam)

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit3 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, 3:15],
                                  y=indata[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  nodesize = 3, #default for classification is 1
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }

saveRDS(rf.fit3, paste0(outnam, ".rds"))
impvals <- importance(rf.fit3)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit3)

snow::stopCluster(cl)
###

### Evaluate with Tweak data ###
twkRF3 <- predict(rf.fit3, newdata=twkdat, type="prob", norm.votes=TRUE)
evRF3 <- data.table(prob=twkRF3[,2], resp=twkdat[, response])
outRF3 <- testEval(evRF3, "RF_4800t_m9_n3")
outRF3 # metrics are worse. RF_4800t_mtry9 wins.
######

### Boosted Regression Trees tuning ###
# original model used gbm library:
# brt.fit <- gbm(response~., data = indata, n.trees = 5000,
#                 interaction.depth = 2, shrinkage = 0.11,
#                 bag.fraction = 0.75, cv.folds = 3, verbose = TRUE,
#                 class.stratify.cv = TRUE, n.cores = 6)

## 1) use dismo::gbm.holdout to find optimal n.trees
outnam <- paste0("BRT_holdout", rnam)
indata <- train_data[, ..allcols][, nlcd := as.factor(nlcd)]

brt.fit <- gbm.holdout(indata, gbm.x = c(2:14), gbm.y = 1,
                       learning.rate = 0.11, tree.complexity = 2,
                       n.trees = 5000, add.trees = 2000,
                       train.fraction = 0.6)

saveRDS(brt.fit, paste0(outnam, ".rds"))
summary(brt.fit)
# 20151 trees (hit default max of 20000).
# Received warning that either learning.rate too low or max trees too low

### Evaluate with Tweak data ###
twkdat <- tweak_data[, ..allcols][, nlcd := as.factor(nlcd)]
twkBRT1 <- predict(brt.fit, newdata=twkdat, n.trees=20151, type="response")
evBRT1 <- data.table(prob=twkBRT1, resp=twkdat[, response])
outBRT1 <- testEval(evBRT1, "BRT_20151trees")
outBRT1 #performed better than original

## 2) dismo does not expose parameters like minobsinnode & cv.fold,
# so next tried package gbm3
brt.fit2 <- gbm3::gbm(response~., data = indata,
                      distribution = "bernoulli",
                      n.trees = 20000,
                      weights = rep(1, nrow(indata)),
                      n.minobsinnode = 6, #default is 10
                      interaction.depth = 2,
                      shrinkage = 0.11,
                      bag.fraction = 0.5,
                      train.fraction = 0.6,
                      cv.folds = 5,
                      verbose = FALSE, #not desirable here
                      class.stratify.cv = TRUE,
                      par.details=gbmParallel(num_threads=ncores))
#received warnings about shrinkage being to high (sigh)
saveRDS(brt.fit2, paste0(outnam, ".rds"))
gbm.perf(brt.fit2, method = "cv") #Suggests 10628 trees

### Evaluate with Tweak data ###
twkBRT2 <- predict(brt.fit2, newdata=twkdat, n.trees=20000, type="response")
evBRT2 <- data.table(prob=twkBRT2, resp=twkdat[, response])
outBRT2 <- testEval(evBRT2, "BRT_gbm3_20kt")
outBRT2 #performed worse than original

## 3) try lowering shrinkage
outnam <- paste0("BRT_gbm3_08shrink", rnam)
#using gbm3 and lowering shrinkage
brt.fit3 <- gbm3::gbm(response~., data = indata,
                      distribution = "bernoulli",
                      n.trees = 25000,
                      weights = rep(1, nrow(indata)),
                      n.minobsinnode = 6, #default is 10
                      interaction.depth = 2,
                      shrinkage = 0.08,
                      bag.fraction = 0.5,
                      train.fraction = 0.6,
                      cv.folds = 5,
                      verbose = FALSE, #not desirable here
                      class.stratify.cv = TRUE,
                      par.details=gbmParallel(num_threads=ncores))
#same warnings about shrinkage being to high
saveRDS(brt.fit3, paste0(outnam, ".rds"))
gbm.perf(brt.fit3, method = "cv") #Suggests 14130 trees

### Evaluate with Tweak data ###
twkBRT3 <- predict(brt.fit3, newdata=twkdat, n.trees=25000, type="response")
evBRT3 <- data.table(prob=twkBRT3, resp=twkdat[, response])
outBRT3 <- testEval(evBRT3, "BRT_gbm3_25kt_08sh")
outBRT3 #even worse

## 4) Never mind, going back to dismo
outnam <- paste0("BRT_gbmholdout_50kt", rnam)
brt.fit4 <- gbm.holdout(indata, gbm.x = c(2:14), gbm.y = 1,
                        learning.rate = 0.11, tree.complexity = 2,
                        n.trees = 20000, add.trees = 5000,
                        max.trees = 50000,
                        train.fraction = 0.6)

saveRDS(brt.fit4, paste0(outnam, ".rds"))
summary(brt.fit4)
gbm.interactions(brt.fit4, mask.object = brt.fit4)
# stopped at 20490 trees, looks like finally hit optimal number

### Evaluate with Tweak data ###
twkBRT4 <- predict(brt.fit4, newdata=twkdat, n.trees=20490, type="response")
evBRT4 <- data.table(prob=twkBRT4, resp=twkdat[, response])
outBRT4 <- testEval(evBRT4, "BRT_20490trees")
outBRT4 #ever so slightly better than BRT_20151trees. stopping here.
######

### Need a GLM using same 70% dataset as RF and BRT ###
train_data <- fread(paste0(pth, "train_data_70pct.csv"), header=TRUE, sep=",")
rnam <- "May14_70pct"
outnam <- paste0("GLMER_", rnam)
bestvars_G <- fread("vars_to_keep_Feb25.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect

f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "I(GDD5^2)", "TWI:ppt_ws", "TWI:ppt_sf", "GDD5:clay",
             "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

# Note: error in most recent lme4 making optional args required
# so glad I updated. :-/
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
######
