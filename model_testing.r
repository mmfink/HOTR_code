##############################################################################
# Run various testing models on the now normalized data
# (see Dataset_normalization.r)
# NOTE that this is a compilation of code run over several weeks and is not
# intended to be run from start to finish as a stand alone script.
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 03/02/2020
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
library(lme4)
library(randomForest)
library(doSNOW)

pth <- "H:/HOTR_models"
setwd(pth)


### Run Exploratory GLMER models ###

set.seed(1157)
rnam <- "Jan21_norm"
train_data <- fread(file = "training_data_seed655.csv", header = TRUE, sep=",")

# 1) Model with all the variables used in GLMER_Oct22, but normalized
outnam <- paste0("GLMER_all_", rnam)
bestvars_G <- fread("vars_to_keep_GLM_GAM.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit1 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit1, paste0(outnam, ".rds"))
summary(glmm.fit1)
###

# 2) Same as 1 but without DEM
outnam <- paste0("GLMER_noDEM_", rnam)
bestvars_G <- fread("vars_to_keep_GLMER_notDEM.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit2 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit2, paste0(outnam, ".rds"))
summary(glmm.fit2)
###

# 3) Only those inputs that were significant in GLMER_Oct22 (including DEM)
outnam <- paste0("GLMER_sig_", rnam)
bestvars_G <- fread("vars_to_keep_GLMER_sig.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit3 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit3, paste0(outnam, ".rds"))
summary(glmm.fit3)
###

# 4) Only those inputs that were significant (excluding DEM)
outnam <- paste0("GLMER_sigNoDEM_", rnam)
bestvars_G <- fread("vars_to_keep_GLMER_sigNoDEM.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit4 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit4, paste0(outnam, ".rds"))
summary(glmm.fit4)
########

### Run Comparison GLMER / RF models ###
set.seed(655)
rnam <- "Feb24_norm"

# 1) Climate Variables #
##GLM
outnam <- paste0("GLMER_climate_", rnam)
train_data <- fread("training_data_seed655.csv", header=TRUE, sep=",")
bestvars_G <- fread("vars_to_keep_CLIMATE.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit1 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit1, paste0(outnam, ".rds"))
summary(glmm.fit1)

##randomForest
outnam <- paste0("RF_climate_", rnam)
ntrees <- 2000
ncores <- 20 # needs to evenly divide the trees
allcols <- c("response", varG)
indata <- train_data[, ..allcols][, response := as.factor(response)]

# tuning mtry #
sub0 <- indata %>% dplyr::filter(response == 0) %>% sample_frac(0.1)
sub1 <- indata %>% dplyr::filter(response == 1) %>% sample_frac(0.1)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varG], y=subtune[, response],
            ntreeTry = 200,
            stepFactor = 1.5)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
# mt <- 2 #

treeSubs <- ntrees/ncores

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit1 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, -1],
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

# 2) Topography and Soils #
##GLM
outnam <- paste0("GLMER_toposoils_", rnam)
train_data <- fread("training_data_seed655.csv", header=TRUE, sep=",")
bestvars_G <- fread("vars_to_keep_TopoSoils.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit2 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit2, paste0(outnam, ".rds"))
summary(glmm.fit2)

##randomForest
outnam <- paste0("RF_toposoils_", rnam)
allcols <- c("response", varG)
indata <- train_data[, ..allcols][, response := as.factor(response)]

# tuning mtry #
sub0 <- indata %>% dplyr::filter(response == 0) %>% sample_frac(0.1)
sub1 <- indata %>% dplyr::filter(response == 1) %>% sample_frac(0.1)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varG], y=subtune[, response],
            ntreeTry = 200,
            stepFactor = 1.5)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
# mt <- 4 #

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit2 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, -1],
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
####

# 3) NLCD #
##GLM
outnam <- paste0("GLMER_nlcd_", rnam)
train_data <- fread("training_data_seed655.csv", header=TRUE, sep=",") #re-reading in case of interim changes while I futz
bestvars_G <- fread("vars_to_keep_nlcd.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit3 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit3, paste0(outnam, ".rds"))
summary(glmm.fit3)

##randomForest
outnam <- paste0("RF_NLCD_onehot_", rnam)
allcols <- c("response", "nlcd.Shrubland", varG) # Adding Shrubland back in
indata <- train_data[, ..allcols][, response := as.factor(response)]

# tuning mtry #
sub0 <- indata %>% dplyr::filter(response == 0) %>% sample_frac(0.3) #smaller fraction results in ties
sub1 <- indata %>% dplyr::filter(response == 1) %>% sample_frac(0.3)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varG], y=subtune[, response],
            ntreeTry = 200,
            stepFactor = 1.5)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
# mt <- 3 #

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit3 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, -1],
                                  y=indata[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }

saveRDS(rf.fit3, paste0(outnam, ".rds"))
impvals <- importance(rf.fit3)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit3)

snow::stopCluster(cl)

# Another RF model Using just NLCD as factor (? does this even make sense ?)
# ***NOPE*** Moving on...
####

# 4) Climate + Dirt + Interactions #
rnam <- "Feb25_norm"

##GLM
outnam <- paste0("GLMER_climdirt_", rnam)
train_data <- fread("training_data_seed655.csv", header=TRUE, sep=",")
bestvars_G <- fread("vars_to_keep_climdirt.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit4 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit4, paste0(outnam, ".rds"))
summary(glmm.fit4)

##randomForest
outnam <- paste0("RF_climdirt_", rnam)
allcols <- c("response", varG)
indata <- train_data[, ..allcols][, response := as.factor(response)]

# tuning mtry #
sub0 <- indata %>% dplyr::filter(response == 0) %>% sample_frac(0.1)
sub1 <- indata %>% dplyr::filter(response == 1) %>% sample_frac(0.1)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varG], y=subtune[, response],
            ntreeTry = 200,
            stepFactor = 1.5)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
# mt <- 4 #

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit4 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, -1],
                                  y=indata[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }

saveRDS(rf.fit4, paste0(outnam, ".rds"))
impvals <- importance(rf.fit4)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit4)

snow::stopCluster(cl)
####

# 5) The Big Kahuna #

##GLM
outnam <- paste0("GLMER_alltypes_", rnam)
bestvars_G <- fread("vars_to_keep_Feb25.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_data[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit5 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit5, paste0(outnam, ".rds"))
summary(glmm.fit5)

##randomForest
outnam <- paste0("RF_alltypes_onehot_", rnam)
allcols <- c("response", varG)
indata <- train_data[, ..allcols][, response := as.factor(response)]

# tuning mtry #
sub0 <- indata %>% dplyr::filter(response == 0) %>% sample_frac(0.1)
sub1 <- indata %>% dplyr::filter(response == 1) %>% sample_frac(0.1)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varG], y=subtune[, response],
            ntreeTry = 200,
            stepFactor = 1.5)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
# mt <- 6 #

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit5 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata[, -1],
                                  y=indata[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }

saveRDS(rf.fit5, paste0(outnam, ".rds"))
impvals <- importance(rf.fit5)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit5)

snow::stopCluster(cl)
####

## 5a) Random Forest with NLCD as factor (no one-hot) ##
rnam <- "Mar02_norm"
outnam <- paste0("RF_alltypes_", rnam)

# The NLCD categorical variable did not get included in the training data before
# so have to add it in now.
respdata <- fread(file.path(pth, "response.csv"), header=TRUE, sep=",")
newtrain <- train_data[, c(1:25, 36:37)][respdata[, .(ID, nlcd)], on=.(ID=ID)]
newtrain <- newtrain[!is.na(response),]
fwrite(newtrain, file = "training_data_seed655_forRF.csv", row.names = TRUE)

bestvars_G <- fread("vars_to_keep_Mar02.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", varG)

indata <- newtrain[, ..allcols][, response := as.factor(response)][, nlcd := as.factor(nlcd)]

# tuning mtry #
sub0 <- indata %>% dplyr::filter(response == 0) %>% sample_frac(0.1)
sub1 <- indata %>% dplyr::filter(response == 1) %>% sample_frac(0.1)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varG], y=subtune[, response],
            ntreeTry = 200,
            stepFactor = 1.5)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
# mt <- 4 #

ntrees <- 2000
ncores <- 20 # needs to evenly divide the trees
treeSubs <- ntrees/ncores

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

rf.fit5a <- foreach(tree = rep(treeSubs,ncores),
                    .combine = randomForest::combine,
                    .packages = c("randomForest", "data.table"),
                    .multicombine = TRUE) %dopar% {
                      randomForest(indata[, -1],
                                   y=indata[, response],
                                   importance=TRUE,
                                   ntree=tree,
                                   mtry=mt,
                                   replace = TRUE,
                                   norm.votes = TRUE)
                    }

saveRDS(rf.fit5a, paste0(outnam, ".rds"))
impvals <- importance(rf.fit5a)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit5a)

snow::stopCluster(cl)
####

## 6) Tweaks to GLM formula March 02 ##
outnam <- paste0("GLMER_alltypes_", rnam)
train_data <- fread("training_data_seed655.csv", header=TRUE, sep=",")
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

glmm.fit6 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit6, paste0(outnam, ".rds"))
summary(glmm.fit6)
########
