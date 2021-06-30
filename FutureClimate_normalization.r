##############################################################################
# Normalize future climate environmental input values to 0-1
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 02/04/2021
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
.libPaths("E:/MMF/R/win-library/3.6")
Sys.setenv(R_USER = "E:/MMF")
library(stringr)
library(data.table)
library(dplyr)
library(raster)
library(rgdal)
library(doSNOW)
library(ggplot2)

pth <- "H:/HOTR_models"
setwd(pth)
ncores <- 16 #modeling server

## Comparative Histograms ##
CurClim_sfppt <- raster("H:/HOTR_models/PARC_norm/ppt_sf.tif")
WWClim_sfppt <- raster("H:/HOTR_models/PARC_norm/WarmWet_ppt_sf.tif")
HDClim_sfppt <- raster("H:/HOTR_models/PARC_norm/HotDry_ppt_sf.tif")

cursamp <- sampleRandom(CurClim_sfppt, 25000)
wwsamp <- sampleRandom(WWClim_sfppt, 25000)
hdsamp <- sampleRandom(HDClim_sfppt, 25000)

dat <- data.table(cursamp, wwsamp, hdsamp)
names(dat) <- c("Current", "WarmWet", "HotDry")
longdat <- melt(dat, measure.vars = c(1:3), variable.name = "sfppt")

ggplot(data=longdat, aes(x=value)) +
  geom_density(aes(color=sfppt), lwd=1.5) +
  labs(x="normalized value") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  guides(color=guide_legend(title = "Summer+Fall Precip"))

# hist(CurClim_sfppt, maxpixels=10000, main="Current Summer-Fall Precip",
#      freq=FALSE)
# hist(WWClim_sfppt, maxpixels=10000, main="Warm-Wet Future Summer-Fall Precip",
#      freq=FALSE)
# hist(HDClim_sfppt, maxpixels=10000, main="Hot-Dry Future Summer-Fall Precip",
#      freq=FALSE)

CurClim_sfppt <- raster("H:/HOTR_models/PARC/gridMET_1994_2014_summer_fall_pr.tif")
WWClim_sfppt <- raster("H:/HOTR_models/PARC/IPSL_CM5A_LR_rcp45_MACA_2079_2099_summer_fall_pr.tif")
HDClim_sfppt <- raster("H:/HOTR_models/PARC/MIROC5_rcp85_MACA_2079_2099_summer_fall_pr.tif")

cursamp <- sampleRandom(CurClim_sfppt, 25000)
wwsamp <- sampleRandom(WWClim_sfppt, 25000)
hdsamp <- sampleRandom(HDClim_sfppt, 25000)

dat <- data.table(cursamp, wwsamp, hdsamp)
names(dat) <- c("Current", "WarmWet", "HotDry")
longdat <- melt(dat, measure.vars = c(1:3), variable.name = "sfppt")

ggplot(data=longdat, aes(x=value)) +
  geom_density(aes(color=sfppt), lwd=1.5) +
  labs(x="original value") +
  guides(color=guide_legend(title = "Summer+Fall Precip"))

## Normalize the rasters ##
CurClim_GDD <- raster("H:/HOTR_models/PARC/gridMET_1994_2014_annual_GDD5.tif")
gddmin <- minValue(CurClim_GDD)
gddmax <- maxValue(CurClim_GDD)
CurClim_sfppt <- raster("H:/HOTR_models/PARC/gridMET_1994_2014_summer_fall_pr.tif")
sfpptmin <- minValue(CurClim_sfppt)
sfpptmax <- maxValue(CurClim_sfppt)
CurClim_wsppt <- raster("H:/HOTR_models/PARC/gridMET_1994_2014_winter_spring_pr.tif")
wspptmin <- minValue(CurClim_wsppt)
wspptmax <- maxValue(CurClim_wsppt)

normem <- fread("clim_inputs_to_normalize2021.csv", header=TRUE, sep=",")
normem$min <- c(gddmin, gddmin, sfpptmin, sfpptmin, wspptmin, wspptmin)
normem$max <- c(gddmax, gddmax, sfpptmax, sfpptmax, wspptmax, wspptmax)

beginCluster(ncores, type='SOCK')
for(r in 1:nrow(normem)){
  ras <- raster(normem$raster[r])
  outname <- str_c("H:/HOTR_models/PARC_norm/", normem$label[r], ".tif")
  mnval <- normem$min[r]
  mxval <- normem$max[r]
  newras <- clusterR(ras, fun=function(x){(x-mnval)/(mxval-mnval)},
                     export=c('mnval', 'mxval'), filename=outname)
}
endCluster()






##################################
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
prespts <- shapefile("H:/HOTR_models/BTPD_prespts_Jan2021.shp", integer64="allow.loss")
abspts <- shapefile("H:/HOTR_models/BTPD_abspts_Jan2021.shp", integer64="allow.loss")
prespts$response <- rep(1, times = length(prespts))
abspts$response <- rep(0, times = length(abspts))
allpts <- rbind(prespts, abspts)
allpts$ID <- row.names(allpts)

# read in the normalized future climate inputs
inputs <- fread("futureclim_inputs_02042021.csv", header=TRUE, sep=",")

# Start parallel cluster of workers
cl <- snow::makeCluster(ncores, type = 'SOCK')
registerDoSNOW(cl)

# Run thru the inputs and extract values - this takes lots of time & memory
outvals <- foreach(r=1:nrow(inputs),
                   .packages = c("raster", "sp", "sf"),
                   .combine = "cbind") %dopar%
  myextract(allpts, inputs[r, 1:2])

outdt <- data.table(outvals)

# Release the cluster
snow::stopCluster(cl)

# Give columns sensible names
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
full_data <- sp::merge(allpts, outdt, by = "ID")

# Save response table
#full_data <- tidyr::drop_na(as.data.frame(full_data))
fwrite(as.data.frame(full_data), file = "response_futureClim_norm.csv")

####

## Replace current climate columns with each climate scenario
full_data <- fread("response_futureClim_norm.csv", header=TRUE, sep=",")
orig_resp <- fread("response_Jan2021.csv", header=TRUE, sep=",")
colorder <- names(orig_resp)
# HotDry
hotdry_resp <- orig_resp[, c(-9, -15, -16)][full_data[, c(1, 4:6)], on=.(ID=ID)]
setnames(hotdry_resp, c("hdGDD5", "hdppt_sf", "hdppt_ws"),
         c("GDD5", "ppt_sf", "ppt_ws"))
setcolorder(hotdry_resp, neworder = colorder)
fwrite(hotdry_resp, file = "HotDry_response.csv")

# WarmWet
warmwet_resp <- orig_resp[, c(-9, -15, -16)][full_data[, c(1, 7:9)], on=.(ID=ID)]
setnames(warmwet_resp, c("wwGDD5", "wwppt_sf", "wwppt_ws"),
         c("GDD5", "ppt_sf", "ppt_ws"))
setcolorder(warmwet_resp, neworder = colorder)
fwrite(warmwet_resp, file = "WarmWet_response.csv")
########
