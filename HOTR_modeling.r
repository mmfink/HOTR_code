library(stringr)
library(data.table)
library(dplyr)
library(raster)
library(doSNOW) # !!! Adjust to however Summit does parallelization !!!
library(biomod2)

#### Set some basic variables ####
set.seed(717)
cp <- 8 # number of cores/threads, adjust per platform

pth <- "H:/HOTR_models/"
setwd(pth)
source("HOTR_model_functions.r")

max_cor <- 0.7 # maximum allowable correlation between variables
test_pct <- 0.25 # proportion of response data to use for testing
file <- paste0(pth, "response.csv") # response table, from create_response_table.r
env_inputs <- paste0(pth, "env_inputs_08292019.csv") # List of raster environmental inputs
########

#### Read in and format the data ####
inputs <- read.table(env_inputs, header=TRUE, sep=",",
                     stringsAsFactors=FALSE)
respdata <- fread(file, header=TRUE, sep=",", nThread=cp)
respdata <- respdata %>% tidyr::drop_na() #Strip any records with NA

# randomly sample equal number of absence as there are presence
cntbl <- respdata %>% count(response)
numpts <- as.numeric(cntbl[2,2]) #number of presence pts
sub_abs <- respdata %>% dplyr::filter(response == 0) %>% sample_n(numpts)
train_data <- respdata %>% dplyr::filter(response == 1) %>% bind_rows(sub_abs)

# turn Grid_ID into a factor at this point
train_data <- as.data.table(train_data)
train_data[, Grid_ID := as.factor(Grid_ID)]

# Parse out Categorical Vars depending on Model Type
# Multi-categorical vars are 6 & 16, one-hot versions are 29 - 45
# X&Y are 46:47, moving those to the front, before 'response' [3] -> [5]
train_dat_G <- train_data %>% dplyr::select(c(1:2, 46:47, 3:5, 7:15, 17:45)) #GLM, GAM
train_dat_tree <- train_data %>% dplyr::select(c(1:2, 46:47, 3:28)) # RF, BRT
# Convert the categorical vars to factors
train_dat_tree[, c("nearType", "nlcd") := .(as.factor(nearType), as.factor(nlcd))]
########

#### Variable Selection ####
# Would prefer to do this interactively, but here's an automated way for now
# Start parallel cluster of workers
cl <- snow::makeCluster(cp, type = "SOCK")
registerDoSNOW(cl)

y <- train_dat_G[, response]
varcols_G <- names(train_dat_G[, 6:45])
vartbl <- devianceExp(train_dat_G, y, varcols_G)
bestvars_G <- numCorr(train_dat_G[, 6:45], max_cor, unlist(vartbl$dev_exp))

varcols_T <- names(train_dat_tree[, 6:30])
vartbl <- devianceExp(train_dat_tree, y, varcols_T)
bestvars_T <- numCorr(train_dat_tree[, c(6:7, 9:17, 19:30)], # can't use factors
                      max_cor, unlist(vartbl[c(-3, -13), 2]))
# Now have to put the factor vars back in as default keepers
bestvars_T <- bestvars_T %>%
  add_row(invar=unlist(vartbl[c(3,13), 1]),
          dev_exp=unlist(vartbl[c(3,13), 2]), keep=c(1,1)) %>%
  arrange(desc(dev_exp))

# save for posterity
write.csv(bestvars_G, paste0(pth, "vars_to_keep_GLM_GAM.csv"))
write.csv(bestvars_T, paste0(pth, "vars_to_keep_RF_BRT.csv"))

# # Final input var tables
# dat_G <- train_dat_G %>% dplyr::select(response, coords.x1, coords.x2, Grid_ID,
#                                        bestvars_G$invar[bestvars_G$keep==1])
# dat_T <- train_dat_tree %>% dplyr::select(response, coords.x1, coords.x2,
#                                           bestvars_T$invar[bestvars_T$keep==1])
# bmdat_G <- BIOMOD_FormatingData(resp.var = as.numeric(y),
#                                 expl.var = as.data.frame(dat_G[, 4:ncol(dat_G)]),
#                                 resp.xy = as.data.frame(dat_G[, 2:3]),
#                                 resp.name = "response")
# bmdat_T <- BIOMOD_FormatingData(resp.var = as.numeric(y),
#                                 expl.var = as.data.frame(dat_T[, 4:ncol(dat_T)]),
#                                 resp.xy = as.data.frame(dat_T[, 2:3]),
#                                 resp.name = "response")
########

#### GLM Responses ####
# Broken: biomod2 can't deal with glmm/glmer
# To be dealt with soonly. But not now.
########

#### Fitting the Models ####
# set model options
varG <- bestvars_G$invar[bestvars_G$keep==1]
varT <- bestvars_T$invar[bestvars_T$keep==1]
f_gam <- makeFormula(respName = 'response',
                     explVar = head(train_data[, ..varG]),
                     type = "s_smoother", interaction.level = 1)
q <- str_replace(as.character(f_gam[3]), "1", "s(Grid_ID,bs='re')")
f_gam <- as.formula(paste("response ~", q, sep = " "))
myBiomodOptions <- BIOMOD_ModelingOptions()
myBiomodOptions <- BIOMOD_ModelingOptions(
  GLM = list( type = 'quadratic',
              interaction.level = 1,
              myFormula = NULL,
              test = 'AIC',
              family = binomial(link = 'logit'),
              mustart = 0.5,
              control = glm.control(epsilon = 1e-07, maxit = 100, trace = FALSE) ),
  GAM = list( algo = 'BAM_mgcv',
              type = 's_smoother',
              k = 3,
              interaction.level = 1,
              myFormula = f_gam,
              family = binomial(link = 'logit'),
              method = 'REML',
              optimizer = c('outer','newton'),
              select = FALSE,
              knots = NULL,
              paraPen = NULL,
              control = list(nthreads = cp, irls.reg = 0, epsilon = 1e-07, maxit = 200, trace = FALSE
                             , mgcv.tol = 1e-07, mgcv.half = 15, rank.tol = 1.49011611938477e-08
                             , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                             , optim = list(factr=1e+07), newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
                             , outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, efs.lspmax = 15, efs.tol = 0.1
                             , keepData = FALSE, edge.correct = FALSE) ),
  RF = list( do.classif = TRUE,
             ntree = 2000,
             mtry = 'default',
             nodesize = 2,
             maxnodes = NULL),
  GBM = list( distribution = 'bernoulli',
              n.trees = 2000,
              interaction.depth = 3,
              n.minobsinnode = 5,
              shrinkage = 0.01,
              bag.fraction = 0.33,
              train.fraction = 1,
              cv.folds = 10,
              keep.data = FALSE,
              verbose = FALSE,
              perf.method = 'cv')
)

modopt <- list(3, (1 - test_pct)*100, c("ROC", "TSS", "POD", "KAPPA"))
models <- c("GAM", "GLM", "BRT", "RF")

whoknows <- foreach(m=models,
                    .packages = c("biomod2", "mgcv", "stats",
                                  "randomForest", "gbm", "dplyr, data.table"),
                    .export = c(varG, varT)) %dopar% {
  if("G" %in% m){invars <- varG
  } else {invars <- varT}
  modpar(m, dat=train_data, vars=invars, modopt_obj=myBiomodOptions, otheropts=modopt)
}

########
snow::stopCluster(cl)

#### Ensemble Model ####
# Don't want to run this until individual models are reviewed and tweaked as needed
# But it will look something like this

# myBiomodEM <- BIOMOD_EnsembleModeling(
#   modeling.output = myBiomodModelOut,
#   chosen.models = 'all',
#   em.by='all',
#   eval.metric = c('TSS'),
#   eval.metric.quality.threshold = c(0.7),
#   prob.mean = T,
#   prob.cv = T,
#   prob.ci = T,
#   prob.ci.alpha = 0.05,
#   prob.median = T,
#   committee.averaging = T,
#   prob.mean.weight = T,
#   prob.mean.weight.decay = 'proportional' )
########

#### Model Projections ####
# Projection over current conditions
# FIXME: will need to (somehow?) relate the one-hot vars back to the rasters. UGH.
# For now, just project BRT & RF
beginCluster(n=cp)
layerStk <- raster::stack(env_inputs$raster, quick=TRUE)
names(layerStk) <- env_inputs$label
endCluster()
# ***A foreach loop that uses get_built_models & BIOMOD_Projection??***
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = layerStk,
  proj.name = "current_test",
  selected.models = "all",
  compress = TRUE,
  build.clamping.mask = TRUE,
  do.stack = FALSE,
  keep.in.memory = FALSE,
  omit.na = FALSE,
  output.format = ".grd")
########