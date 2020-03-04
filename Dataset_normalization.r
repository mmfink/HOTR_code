##############################################################################
# Normalize all numeric environmental input values to 0-1
# Binary variables left as 0 or 1, northness/eastness left as -1 to +1
# Multi-categorical variables not changed
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 02/24/2020
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
library(rgdal)
library(doSNOW)

pth <- "H:/HOTR_models"
setwd(pth)
set.seed(655)

## First normalize the rasters ##
normem <- fread("env_inputs_to_normalize.csv", header=TRUE, sep=",")
beginCluster(16, type='SOCK')
for(r in 1:nrow(normem)){
  ras <- raster(normem$raster[r])
  outname <- str_c("H:/HOTR_models/PARC_norm/", normem$label[r], ".tif")
  mnval <- minValue(ras)
  mxval <- maxValue(ras)
  newras <- clusterR(ras, fun=function(x){(x-mnval)/(mxval-mnval)},
                     export=c('mnval', 'mxval'), filename=outname)
}
endCluster()

# BROKEN: didn't work for slope, it couldn't figure out the min/max vals
# so here it is done individually
mnval <- 0
mxval <- 8.7799577713013
oldslope <- raster("H:/HOTR_models/PARC/HOTR_Dinf_slope.tif")
outname <- "H:/HOTR_models/PARC_norm/slope.tif"
newras <- calc(oldslope, fun=function(x){(x-mnval)/(mxval-mnval)},
               filename=outname)
####

## Then build new response table with the normalized data ##
# Function to extract raster values
myextract <- function(inpts, inrow){
  rasname <- inrow$raster
  inlbl <- inrow$label
  inras <- raster(rasname)
  x <- extract(inras, inpts)
  names(x) <- inlbl
  return(x)
}

# Bring in the point data
prespts <- shapefile("D:/GIS/Projects/HOTR/BTPD_prespts.shp", integer64="allow.loss")
abspts <- shapefile("D:/GIS/Projects/HOTR/BTPD_abspts.shp", integer64="allow.loss")
prespts$response <- rep(1, times = length(prespts))
abspts$response <- rep(0, times = length(abspts))
allpts <- rbind(prespts, abspts)
allpts$ID <- row.names(allpts)

# read in Environmental Inputs
inputs <- fread("env_inputs_02202020.csv", header=TRUE, sep=",")

# Start parallel cluster of workers
cl <- snow::makeCluster(16, type = 'SOCK')
registerDoSNOW(cl)

# Run thru the inputs and extract values - this takes lots of time & memory
outvals <- foreach(r=1:nrow(inputs),
                   .packages = c("raster"),
                   .combine = "cbind") %dopar%
  myextract(allpts, inputs[r, 1:2])

outdt <- data.table(outvals)

# Name the outputs
for(n in names(outdt)){
  if(str_detect(n, "result")){
    z <- str_extract_all(n, "[:digit:]")
    if(nchar(z) > 1){
      z <- paste0(unlist(z), collapse = "")
    }
    idx <- as.numeric(z)
    nn <- inputs$label[idx]
    names(outdt)[idx] <- nn
  }
}

outdt <- tibble::rownames_to_column(outdt, var = "ID")
output <- sp::merge(allpts, outdt)
# Save response table
fwrite(as.data.frame(output), file = "response_norm_full.csv")

# Release the cluster
snow::stopCluster(cl)
####

## Finally, create new training dataset ##
# response table, created above
respdata <- fread(file.path(pth, "response_norm_full.csv"), header=TRUE, sep=",")
respdata <- respdata %>% tidyr::drop_na() #Strip any records with NAs

# randomly sample equal number of absence as there are presence
cntbl <- respdata %>% count(response)
numpts <- as.numeric(cntbl[2,2]) #number of presence pts
sub_abs <- respdata %>% dplyr::filter(response == 0) %>% sample_n(numpts)
train_data <- respdata %>% dplyr::filter(response == 1) %>% bind_rows(sub_abs)

# Save training dataset
fwrite(train_data, file = "training_data_seed655.csv")
########
