library(stringr)
library(data.table)
library(dplyr)

# Note, if on Modeling Server,
# set the following from cmd.exe before starting:
# setx TMPDIR "E:\\MMF\\temp"

#ncores <- (parallel::detectCores()) - 1

pth <- "H:/HOTR_models/"
setwd(pth)
set.seed(9832)
rnam <- "Nov04_full"

train_dat <- fread("training_data_seed717.csv", header=TRUE, sep=",")
bestvars_T <- fread("vars_to_keep_RF_Oct22.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]

##### Random Forest ##### ##10/21/2019 This runs on Model Server##
library(randomForest)
library(doSNOW)

ntrees <- 2000
ncores <- 20 # needs to evenly divide the trees
outnam <- paste0("RF_", rnam)

### only need to run this section once ###
sub0 <- train_dat %>% dplyr::filter(response == 0) %>% sample_frac(0.1)
sub1 <- train_dat %>% dplyr::filter(response == 1) %>% sample_frac(0.1)
subtune <- as.data.table(bind_rows(sub1, sub0))

x <- tuneRF(subtune[, ..varT], y=subtune[, response],
            ntreeTry = 200,
            stepFactor = 1.5)
mt <- x[x[,2] == min(x[,2]),1]
rm(subtune, sub0, sub1)
### mt <- 4 ###

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

treeSubs <- ntrees/ncores
#indata <- train_dat[, response := as.factor(response)][, ..varT] #This will be a classification, not a regression
allcols <- c("response", varT)
indata <- train_dat[, ..allcols][, nlcd := as.factor(nlcd)][, nearType := as.factor(nearType)][, response := as.factor(response)]

rf.fit1 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata,
                                  y=indata[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  #strata = train_dat[, Grid_ID],
                                  #sampsize = sampSizeVec,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }

saveRDS(rf.fit1, paste0(outnam, ".rds"))
impvals <- importance(rf.fit1)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit1)

snow::stopCluster(cl)
##########

##### Boosted Regression Tree ##### ##10/28/2019 runs on Model Server##
library(gbm)

outnam <- paste0("BRT_", rnam)
bestvars_T <- fread("vars_to_keep_BRT_Oct28.csv", header=TRUE, sep=",")
varT <- bestvars_T$invar[bestvars_T$keep==1]
allcols <- c("response", varT)
indata <- train_dat[, ..allcols][, nlcd := as.factor(nlcd)][, nearType := as.factor(nearType)]

# the following settings are from multiple fitting steps in SAHM (which gets thru model fitting then dies)
brt.fit1 <- gbm(response~., data = indata, n.trees = 5000, interaction.depth = 20, shrinkage = 0.1049,
                bag.fraction = 0.75, cv.folds = 3, verbose = TRUE, class.stratify.cv = TRUE, n.cores = 20)

saveRDS(brt.fit1, paste0(outnam, ".rds"))
summary(brt.fit1)

##########

##### Generalized Linear Mixed Model ##### ##10/22/2019 This runs on Model Server##
library(lme4)

outnam <- paste0("GLMER_", rnam)
bestvars_G <- fread("vars_to_keep_GLM_GAM.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", "Grid_ID", varG)
indata <- train_dat[, Grid_ID := as.factor(Grid_ID)][, ..allcols] #Grid_ID = random effect
f_glm <- biomod2::makeFormula(respName = "response",
                     explVar = head(indata[, ..varG]),
                     type = "simple",
                     interaction.level = 0)
q <- str_c(c(as.character(f_glm[3]), "(1|Grid_ID)", "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm <- as.formula(paste("response ~", q, sep = " "))

glmm.fit1 <- glmer(f_glm, data = indata, family = "binomial", control = glmerControl(optimizer = "bobyqa"), nAGQ = 0)

saveRDS(glmm.fit1, paste0(outnam, ".rds"))
summary(glmm.fit1)

##########

##### Generalized Linear Model ##### ## Nov 4, 2019 runs on desktop in ~1 minute (!)
# no random effects
outnam <- paste0("GLM_", rnam)
bestvars_G <- fread("vars_to_keep_GLM_GAM.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", varG)
indata <- train_dat[, ..allcols]
f_glm2 <- biomod2::makeFormula(respName = "response",
                              explVar = head(indata[, ..varG]),
                              type = "simple",
                              interaction.level = 0)
q <- str_c(c(as.character(f_glm2[3]), "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
f_glm2 <- as.formula(paste("response ~", q, sep = " "))

glm.fit1 <- glm(f_glm2, data = indata, family = "binomial")

saveRDS(glm.fit1, paste0(outnam, ".rds"), version = 3)
summary(glm.fit1)

# ##### Generalized Additive Mixed Model ##### ##10/28/2019 - this runs for days and days and never ends##
# library(gamm4)
#
# outnam <- paste0("GAMM4_", rnam)
# # Don't try to smooth the binary variables
# varSm <- varG[str_which(varG, "nlcd\\.|nearType\\.|hab_non", negate = TRUE)]
# varOther <- varG[str_which(varG, "nlcd\\.|nearType\\.|hab_non")]
# f_gam <- biomod2::makeFormula(respName = 'response',
#                      explVar = head(indata[, ..varSm]),
#                      type = "s_smoother", interaction.level = 0)
# q <- str_c(c(as.character(f_gam[3]), varOther, "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
# q <- str_remove(q, "1 \\+ ")
# f_gam <- as.formula(paste("response ~", q, sep = " "))
#
# gamm.fit1 <- gamm4(f_gam, data = indata, random = ~(1|Grid_ID), family = "binomial")
#
# saveRDS(gamm.fit1, paste0(outnam, ".rds"))
# summary(gamm.fit1)
# ##########

##### Generalized Additive Mixed Model ##### ##10/29/2019 Runs on Model Server##
library(mgcv)
library(parallel)

outnam <- paste0("BAM_", rnam)
ncores <- 14 # do not exceed actual number of cores
cl <- parallel::makeCluster(ncores)

#No random effects version
f_gam <- as.formula("response ~ s(TRI, k=3) + s(DEM, k=3) + s(GDD5, k=3) +
                    s(ppt_sf, k=3) + s(ph, k=3) + s(TWI, k=3) + s(bd, k=3) +
                    s(ppt_ws, k=3) + s(clay, k=3) + s(om, k=3) + s(eastness, k=3) +
                    s(depth, k=3) + s(northness, k=3) + s(distToNon, k=3) +
                    s(nlcd_patch, k=3) + nlcd.Grassland + hab_non +
                    nlcd.Forest + nearType.DevelOpen + nearType.Wetland +
                    nearType.Forest + nearType.Cropland + nlcd.Developed + nearType.Other +
                    nlcd.Wetland + nlcd.Other + nlcd.Water + nlcd.DevelOpen +
                    nearType.Developed + nlcd.PastureHay + TWI:clay + GDD5:clay +
                    ppt_sf:clay + ppt_ws:clay")

bestvars_G <- fread("vars_to_keep_GLM_GAM.csv", header=TRUE, sep=",")
varG <- bestvars_G$invar[bestvars_G$keep==1]
allcols <- c("response", varG)
indata <- train_dat[, response := as.factor(response)][, ..allcols]
bip <- binomial(link = "logit")

gam_fit1 <- bam(f_gam, family = bip, data = indata, chunk.size = 10000,
                cluster = cl)

saveRDS(gam_fit1, paste0(outnam, ".rds"))
summary(gam_fit1)

### 11/04/2019 after running gam.check, decided to alter some k values in formula
f_gam2 <- as.formula("response ~ s(TRI, k = 3) + s(DEM, k = 5) + s(GDD5, k = 5) + s(ppt_sf,
    k = 5) + s(ph, k = 3) + s(TWI, k = 3) + s(bd, k = 5) + s(ppt_ws,
    k = 5) + s(clay, k = 10) + s(om, k = 3) + s(eastness, k = 3) +
    s(depth, k = 5) + s(northness, k = 3) + s(distToNon, k = 3) +
    s(nlcd_patch, k = 3) + nlcd.Grassland + hab_non + nlcd.Forest +
    nearType.DevelOpen + nearType.Wetland + nearType.Forest +
    nearType.Cropland + nlcd.Developed + nearType.Other + nlcd.Wetland +
    nlcd.Other + nlcd.Water + nlcd.DevelOpen + nearType.Developed +
    nlcd.PastureHay + TWI:clay + GDD5:clay + ppt_sf:clay + ppt_ws:clay")

gam_fit2 <- bam(f_gam2, family = bip, data = indata, chunk.size = 10000,
                cluster = cl)

saveRDS(gam_fit2, paste0(outnam, ".rds"))
summary(gam_fit2)

# #Random effects version ##runs out of memory##
# f_gam2 <- as.formula("response ~ s(TRI, k=3, bs='cr') + s(DEM, k=3, bs='cr') + s(GDD5, k=3, bs='cr') +
#                     s(ppt_sf, k=3, bs='cr') + s(ph, k=3, bs='cr') + s(TWI, k=3, bs='cr') + s(bd, k=3, bs='cr') +
#                     s(ppt_ws, k=3, bs='cr') + s(clay, k=3, bs='cr') + s(om, k=3, bs='cr') + s(eastness, k=3, bs='cr') +
#                     s(depth, k=3, bs='cr') + s(northness, k=3, bs='cr') + s(distToNon, k=3, bs='cr') +
#                     s(nlcd_patch, k=3, bs='cr') + nlcd.Grassland + hab_non +
#                     nlcd.Forest + nearType.DevelOpen + nearType.Wetland +
#                     nearType.Forest + nearType.Cropland + nlcd.Developed + nearType.Other +
#                     nlcd.Wetland + nlcd.Other + nlcd.Water + nlcd.DevelOpen +
#                     nearType.Developed + nlcd.PastureHay + TWI:clay + GDD5:clay +
#                     ppt_sf:clay + ppt_ws:clay + s(Grid_ID, bs = 're')")
#
# allcols <- c("response", "Grid_ID", varG)
# indata <- train_dat[, response := as.factor(response)][, Grid_ID := as.factor(Grid_ID)][, ..allcols]
# ncores <- 10 # ran out of memory at 14
# cl <- parallel::makeCluster(ncores)
#
# gam_fit2 <- bam(f_gam2, family = bip, data = indata, chunk.size = 45000,
#                 cluster = cl, gc.level = 2)
# outnam <- paste0("BAMre_", rnam)
# saveRDS(gam_fit2, paste0(outnam, ".rds"))
# summary(gam_fit2)

parallel::stopCluster(cl)
##########
