##############################################################################
# Fit Multiple Species Distribution Models
# This uses the data tables already created in Variable_selection.r
# Separate code chunks create model objects for:
#    Generalized Linear Mixed Model (GLMM)
#    Random Forest (RF)
#    Boosted Regression Trees (BRT)
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 11/06/2019
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

pth <- "H:/HOTR_models/"
setwd(pth)
set.seed(9832)
rnam <- "Nov04_full" # Use to uniquely name each model run

# Load response table
train_dat <- fread("training_data_seed717.csv", header=TRUE, sep=",")


##### Random Forest - No random effects #####
library(randomForest)
library(doSNOW)

ntrees <- 2000
ncores <- 20 # needs to evenly divide the trees
outnam <- paste0("RF_", rnam)

# Filter response table by variables chosen in HOTR_model_functions.r
bestvars_T <- fread("vars_to_keep_RF_Oct22.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]

### only need to run this section once ###
sub0 <- train_dat %>% dplyr::filter(response == 0) %>% sample_frac(0.1)
sub1 <- train_dat %>% dplyr::filter(response == 1) %>% sample_frac(0.1)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varT], y=subtune[, response],
            ntreeTry = 200,
            stepFactor = 1.5)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
### mt <- 4 ### Result of tuning, just use this

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

treeSubs <- ntrees/ncores
allcols <- c("response", varT)
indata <- train_dat[, ..allcols][, nlcd := as.factor(nlcd)][, nearType := as.factor(nearType)][, response := as.factor(response)]

# Create model object
rf.fit1 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata,
                                  y=indata[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }

# Save model object
saveRDS(rf.fit1, paste0(outnam, ".rds"))

# Output some basic model evaluations while I have it loaded
# More can be found in Model_review_Nov.Rmd
impvals <- importance(rf.fit1)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit1)

##### Random Forest - try with random effects #####
### NOT YET RUN ###
# outnam <- paste0("RFM_", rnam)
# allcols <- c("response", "Grid_ID", varT)
# indata <- train_dat[, ..allcols][, nlcd := as.factor(nlcd)][, nearType := as.factor(nearType)]
# indata <- indata[, response := as.factor(response)][, Grid_ID := as.factor(Grid_ID)]
#
# rf.fit2 <- foreach(tree = rep(treeSubs,ncores),
#                    .combine = randomForest::combine,
#                    .packages = c("randomForest", "data.table"),
#                    .multicombine = TRUE) %dopar% {
#                        randomForest(indata,
#                                     y=indata[, response],
#                                     importance=TRUE,
#                                     ntree=tree,
#                                     mtry=mt,
#                                     strata = train_dat[, Grid_ID],
#                                     #sampsize = sampSizeVec,
#                                     replace = TRUE,
#                                     norm.votes = TRUE)
#                    }
# saveRDS(rf.fit2, paste0(outnam, ".rds"))

snow::stopCluster(cl)
##########

##### Boosted Regression Tree - No random effects #####
library(gbm)

outnam <- paste0("BRT_", rnam)
# Filter response table by variables chosen in HOTR_model_functions.r
bestvars_T <- fread("vars_to_keep_BRT_Oct28.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]
allcols <- c("response", varT)
indata <- train_dat[, ..allcols][, nlcd := as.factor(nlcd)][, nearType := as.factor(nearType)]

# Create model object
# the following settings are from multiple fitting steps in SAHM (which gets thru model fitting then dies)
brt.fit1 <- gbm(response~., data = indata, n.trees = 5000, interaction.depth = 20, shrinkage = 0.1049,
                bag.fraction = 0.75, cv.folds = 3, verbose = TRUE, class.stratify.cv = TRUE, n.cores = 20)

# Save model object
saveRDS(brt.fit1, paste0(outnam, ".rds"))
# Output some basic model evaluations while I have it loaded
# More can be found in Model_review_Nov.Rmd
summary(brt.fit1)

##########

##### Generalized Linear Mixed Model #####
library(lme4)

outnam <- paste0("GLMER_", rnam)
# Filter response table by variables chosen in HOTR_model_functions.r
bestvars_G <- fread("vars_to_keep_GLM_GAM.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_dat[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect

# Using package biomod2 just as a shortcut for initial formula
f_glm <- biomod2::makeFormula(respName = "response",
                     explVar = head(indata[, ..varG]),
                     type = "simple",
                     interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

# Create model object
glmm.fit1 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

# Save model object
saveRDS(glmm.fit1, paste0(outnam, ".rds"))
# Output some basic model evaluations while I have it loaded
# More can be found in Model_review_Nov.Rmd
summary(glmm.fit1)

##########
