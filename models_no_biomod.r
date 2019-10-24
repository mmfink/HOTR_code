library(stringr)
library(data.table)
library(dplyr)

# Note, if on Modeling Server,
# set the following from cmd.exe before starting:
# setx TMPDIR "E:\\MMF\\temp"

#ncores <- (parallel::detectCores()) - 1

pth <- "H:/HOTR_models/"
setwd(pth)
set.seed(3782)
rnam <- "Oct22_full"

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
indata <- train_dat[, response := as.factor(response)][, ..varT] #This will be a classification, not a regression

rf.fit1 <- foreach(tree = rep(treeSubs,ncores),
                   .combine = randomForest::combine,
                   .packages = c("randomForest", "data.table"),
                   .multicombine = TRUE) %dopar% {
                     randomForest(indata,
                                  y=train_dat[, response],
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  #strata = factor(c(0,1)),
                                  #strata = train_dat[, Grid_ID],
                                  #sampsize = sampSizeVec,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                   }

saveRDS(rf.fit1, paste0(outnam, ".rds"), version = 3)
impvals <- importance(rf.fit1)
fwrite(impvals, file = paste0(outnam, "_varImp.txt"), row.names = TRUE)
varImpPlot(rf.fit1)

snow::stopCluster(cl)
##########

##### Boosted Regression Tree #####
library(gbm)


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

saveRDS(glmm.fit1, paste0(outnam, ".rds"), version = 3)
summary(glmm.fit1)

##########

##### Generalized Additive Mixed Model #####
library(gamm4)

outnam <- paste0("GAMM4_", rnam)
# Don't try to smooth the binary variables
varSm <- varG[str_which(varG, "nlcd\\.|nearType\\.|hab_non", negate = TRUE)]
varOther <- varG[str_which(varG, "nlcd\\.|nearType\\.|hab_non")]
f_gam <- biomod2::makeFormula(respName = 'response',
                     explVar = head(indata[, ..varSm]),
                     type = "s_smoother", interaction.level = 0)
q <- str_c(c(as.character(f_gam[3]), varOther, "TWI:clay", "GDD5:clay", "ppt_sf:clay", "ppt_ws:clay"), collapse = " + ")
q <- str_remove(q, "1 \\+ ")
f_gam <- as.formula(paste("response ~", q, sep = " "))

gamm.fit1 <- gamm4(f_gam, data = indata, random = ~(1|Grid_ID), family = "binomial")

saveRDS(gamm.fit1, paste0(outnam, ".rds"), version = 3)
summary(gamm.fit1)
##########

##### Evaluation Stats #####
library(ROCR)
library(PresenceAbsence) # for Kappa

rf.mod <- readRDS(paste0(pth, "RF_Oct21_full.rds"))
p_train <- as.vector(rf.mod$votes[,2])
pred_train <- prediction(p_train, rf.mod$y)
AUC <- performance(pred_train, measure = "auc")@y.values[[1]]
plot(performance(pred_train, "tpr", "fpr"), lwd= 3, col="red", main="RF_Oct21_full AUC")
sens <- mean(performance(pred_train, measure = "sens")@y.values[[1]])
spec <- mean(performance(pred_train, measure = "spec")@y.values[[1]])
TSS <- sens + spec - 1
##########

##### Raster Creation #####
library(raster)

tilebuilder <- function(raster, size = 1000, overlap = NULL, out = c('data.frame', 'list')){
  # From Rafael Wüest, 10/19/2019: sdm_parallel.html
  # 'size' is tile size in linear map units (e.g., size = 3000 is a tile 3,000m on a side)
  out <- match.arg(out)
  # get raster extents
  xmin <- xmin(raster)
  xmax <- xmax(raster)
  ymin <- ymin(raster)
  ymax <- ymax(raster)
  xmins <- c(seq(xmin ,xmax, by = size))
  ymins <- c(seq(ymin, ymax, by = size))
  exts <- expand.grid(xmin = xmins, ymin = ymins)
  exts$ymax <- exts$ymin + size
  exts$xmax <- exts$xmin + size
  # if overlapped tiles are requested, create new columns with buffered extents
  if (!is.null(overlap)){
    exts$yminb <- exts$ymin
    exts$xminb <- exts$xmin
    exts$ymaxb <- exts$ymax
    exts$xmaxb <- exts$xmax
    t1 <- (exts$ymin - overlap) >= ymin
    exts$yminb[t1] <- exts$ymin[t1] - overlap
    t2 <- (exts$xmin - overlap) >= xmin
    exts$xminb[t2] <- exts$xmin[t2]-overlap
    t3 <- (exts$ymax + overlap) <= ymax
    exts$ymaxb[t3] <- exts$ymax[t3] + overlap
    t4 <- (exts$xmax + overlap) <= xmax
    exts$xmaxb[t4] <- exts$xmax[t4] + overlap
  }
  exts$tile <- 1:nrow(exts)
  if (out == 'list') {
    lapply(exts$tile, function(x) extent(unlist(exts[x, c('xmin', 'xmax', 'ymin', 'ymax')])))
  } else {
    exts
  }
}

env_inputs <- fread("env_inputs_10022019.csv", header=TRUE, sep=",")
# also - "onehot_rasters.csv"
# then pair down list to only those inputs used in each model.
ncores <- (parallel::detectCores()) - 1
beginCluster(n=ncores)
layerStk <- raster::stack(env_inputs$raster, quick=TRUE)
names(layerStk) <- env_inputs$label
endCluster()
##########