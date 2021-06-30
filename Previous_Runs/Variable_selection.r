##############################################################################
# Variable Selection
# Select a subset of all environmental variables that maximizes deviance
# explained and minimizes correlation with each other.
# This uses the response table already created in create_response_table.r
# Additionally, the functions contained in HOTR_model_functions.r are required
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
library(raster)
library(doSNOW)

#### Set some basic variables ####
set.seed(717)
cp <- 8 # number of cores/threads, adjust per platform
max_cor <- 0.7 # maximum allowable correlation between variables

pth <- "H:/HOTR_models/"
setwd(pth)
source("HOTR_model_functions.r")
file <- paste0(pth, "response.csv") # response table, from create_response_table.r
########

#### Read in and format the data ####
respdata <- fread(file, header=TRUE, sep=",")
respdata <- respdata %>% tidyr::drop_na() #Strip any records with NAs

# randomly sample equal number of absence as there are presence
cntbl <- respdata %>% count(response)
numpts <- as.numeric(cntbl[2,2]) #number of presence pts
sub_abs <- respdata %>% dplyr::filter(response == 0) %>% sample_n(numpts)
train_data <- respdata %>% dplyr::filter(response == 1) %>% bind_rows(sub_abs)
# Save this for future use
fwrite(train_data, file = "training_data_seed717.csv")
########

#### Variable Selection ####
# Would prefer to do this interactively, but here's an automated way for now
# Start parallel cluster of workers
cl <- snow::makeCluster(cp, type = "SOCK")
registerDoSNOW(cl)

# Parse out Categorical Vars depending on Model Type
# Multi-categorical vars are 6 & 16, one-hot versions are 29 - 45
# X&Y are 46:47, moving those to the front, before 'response' [3] -> [5]
train_dat_G <- train_data %>% dplyr::select(c(1:2, 46:47, 3:5, 7:15, 17:45)) #GLM, GAM
train_dat_tree <- train_data %>% dplyr::select(c(1:2, 46:47, 3:28)) # RF, BRT
# Convert the categorical vars to factors
train_dat_tree[, c("nearType", "nlcd") := .(as.factor(nearType), as.factor(nlcd))]

y <- train_dat_G[, response]
varcols_G <- names(train_dat_G[, 6:45])
vartbl <- devianceExp(train_dat_G, y, varcols_G)
bestvars_G <- numCorr(train_dat_G[, 6:45], max_cor, unlist(vartbl$dev_exp))

varcols_T <- names(train_dat_tree[, 6:30])
vartbl <- devianceExp(train_dat_tree, y, varcols_T)
bestvars_T <- numCorr(train_dat_tree[, c(6:7, 9:17, 19:30)], # can't use factors
                      max_cor, unlist(vartbl[c(-3, -13), 2]))

snow::stopCluster(cl)

# Now have to put the factor vars back in as default keepers
bestvars_T <- bestvars_T %>%
  add_row(invar=unlist(vartbl[c(3,13), 1]),
          dev_exp=unlist(vartbl[c(3,13), 2]), keep=c(1,1)) %>%
  arrange(desc(dev_exp))

# save for posterity
write.csv(bestvars_G, paste0(pth, "vars_to_keep_GLM_GAM.csv"))
write.csv(bestvars_T, paste0(pth, "vars_to_keep_RF_BRT.csv"))
########
