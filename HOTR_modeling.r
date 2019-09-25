library(biomod2)
library(stringr)
library(data.table)
library(dplyr)
library(raster)
library(doSNOW) # !!! Adjust to however Summit does parallelization !!!

#### Set some basic variables ####
set.seed(8301)
source("HOTR_model_functions.r")
pth <- "D:/GIS/Projects/HOTR/"

max_cor <- 0.7 # maximum allowable correlation between variables
test_pct <- 0.25 # proportion of response data to use for testing
file <- paste0(pth, "response_test.csv") # response table, from create_response_table.r
env_inputs <- paste0(pth, "env_inputs_08292019.csv") # List of raster environmental inputs
########

#### Read in and format the data ####
inputs <- read.table(env_inputs, header=TRUE, sep=",",
                     stringsAsFactors=FALSE)
respdata <- fread(file, header=TRUE, sep=",", nThread=6)
respdata <- respdata %>% tidyr::drop_na() #Strip any records with NA

# randomly sample equal number of absence as there are presence
cntbl <- respdata %>% count(response)
numpts <- as.numeric(cntbl[2,2]) #number of presence pts
sub_abs <- respdata %>% dplyr::filter(response == 0) %>% sample_n(numpts)
# get presence records in the same format
sub_pres <- respdata %>% dplyr::filter(response == 1)
train_data <- bind_rows(sub_pres, sub_abs)

# Parse out Categorical Vars depending on Model Type
# Multi-categorical vars are 6 & 16, one-hot versions are 29 - 45
# X&Y are 46:47, moving those to the front, before 'response' [3] -> [5]
train_dat_G <- train_data %>% dplyr::select(c(1:2, 46:47, 3:5, 7:15, 17:45)) #GLM, GAM
train_dat_tree <- train_data %>% dplyr::select(c(1:2, 46:47, 3:28)) # RF, BRT
########

#### Variable Selection ####
# Would prefer to do this interactively, but here's an automated way for now
# Start parallel cluster of workers
cl <- snow::makeCluster(8, type = "SOCK") # !!! Likewise !!!
registerDoSNOW(cl)

varcols_G <- names(train_dat_G[6:45])
y_G <- train_dat_G %>% pull(var = response)
vartbl <- devianceExp(train_dat_G, y_G, varcols_G)
bestvars_G <- numCorr(train_dat_G[, 6:45], max_cor, unlist(vartbl$dev_exp))

varcols_T <- names(train_dat_tree[6:30])
y_T <- train_dat_tree %>% pull(var = response)
vartbl <- devianceExp(train_dat_tree, y_T, varcols_T)
bestvars_T <- numCorr(train_dat_tree[, 6:30], max_cor, unlist(vartbl$dev_exp))

# save for posterity
write.csv(bestvars_G, paste0(pth, "vars_to_keep_GLM_GAM.csv"))
write.csv(bestvars_T, paste0(pth, "vars_to_keep_RF_BRT.csv"))

# Release the cluster so Biomod can do its thing
snow::stopCluster(cl)
########

#### GLM Responses ####
# To be dealt with soonly. But not now.
########

#### Fitting the Models ####
# set model options
myBiomodOptions <- BIOMOD_ModelingOptions()

########

