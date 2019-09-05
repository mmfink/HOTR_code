##############################################################################
# Extract raster values to pres/abs points
# Presence and Absence points already created in Generate_sample_points.r
# The output depends on the type of raster input ([N]umeric, [B]inary,
# [C]ategorical, [S]pecial). Type is specified in the inputs csv file.
#   Numeric and Binary (0, 1) are left as-is
#   Categorical inputs are transformed into one-hot encoding (dummy variables)
#   Special inputs use the value (integer) as a look-up to a Count attribute
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 09/03/2019
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

library(data.table)
library(dplyr)
library(raster)
library(rgdal)
library(foreign)
library(ade4)
library(stringr)
library(doSNOW)

# point data
prespts <- shapefile("D:/GIS/Projects/HOTR/BTPD_prespts.shp", integer64="allow.loss")
abspts <- shapefile("D:/GIS/Projects/HOTR/BTPD_abspts.shp", integer64="allow.loss")
allpts <- rbind(prespts, abspts)

# read in Environmental Inputs
inputs <- fread("env_inputs_08292019.csv", sep = ",", nThread = 2)

# Start parallel cluster of workers
cl <- snow::makeCluster(24, type = 'SOCK')
registerDoSNOW(cl)

# Function to extract correct values depending on type
myextract <- function(inpts, inrow){
  rasname <- inrow$raster
  inlbl <- inrow$label
  intype <- inrow$type
  if(intype == "S"){
    inras <- ratify(raster(rasname))
    dbfname <- str_c(rasname, ".vat.dbf")
    indbf <- read.dbf(dbfname)
    rat <- levels(inras)[[1]]
    rat <- left_join(rat, indbf, by=c("ID" = "Value"))
    levels(inras) <- rat
    a <- extract(inras, inpts)
    x <- factorValues(inras, a, att = "Count")
  } else {
    inras <- raster(rasname)
    x <- extract(inras, inpts)
  }
  names(x) <- inlbl
  return(x)
}

# Function to transform Categorical values to dummy variables
onehot <- function(lbl, df){
  key <- tribble(~IN, ~OUT,
                 "1", "Cropland",
                 "2", "Forest",
                 "3", "Water",
                 "4", "Developed",
                 "5", "Other",
                 "6", "Grassland",
                 "7", "Shrubland",
                 "8", "PastureHay",
                 "9", "DevelOpen",
                 "10", "Wetland")
  #there has *got* to be a cleaner way than this hideous pipe...
  onevec <- df %>%
    transmute_at(lbl, as.character) %>%
    rename_at(lbl, funs(str_replace(lbl, lbl, "IN"))) %>%
    left_join(key) %>%
    transmute(new = as.factor(OUT))
  dum <- acm.disjonctif(as.matrix(onevec))
  dum <- dum %>% setNames(paste0(lbl, names(.)))
  return(dum)
}

# Run thru the inputs and extract values - this takes lots of time & memory
outvals <- foreach(r=1:nrow(inputs),
                   .packages = c("raster", "stringr", "foreign", "dplyr"),
                   .combine = "cbind") %dopar%
  myextract(allpts, inputs[r, 1:3])

# Naming only works for the Special types for some reason,
# so now have to name the rest
for(n in names(outvals)){
  if(str_detect(n, "result")){
    z <- str_extract_all(n, "[:digit:]")
    if(nchar(z) > 1){
      z <- paste0(unlist(z), collapse = "")
    }
    idx <- as.numeric(z)
    nn <- inputs$label[idx]
    names(outvals)[idx] <- nn
  }
}

outvals <- tibble::rownames_to_column(outvals, var = "ID")
# Save in case of disaster with the next bit
fwrite(outvals, file = "D:/GIS/Projects/HOTR/outvals_raw.csv", nThread = 7)

# Make dummy variables and merge them back into the output
outv_onehot <- foreach(r=1:nrow(inputs),
                       .packages = c("ade4", "dplyr", "stringr"),
                       .combine = "cbind") %:%
  when(inputs$type[r] == "C") %dopar%
  onehot(inputs$label[r], outvals)

outvals <- bind_cols(outvals, outv_onehot)
output <- sp::merge(allpts, outvals, by.x="ptID", by.y="ID")

# Release the cluster
snow::stopCluster(cl)

fwrite(output, file = "D:/GIS/Projects/HOTR/response.csv", nThread = 16)
# Note that shapefile output truncates the fieldnames
shapefile(output, filename="D:/GIS/Projects/HOTR/allpts_response.shp", overwrite=TRUE)
