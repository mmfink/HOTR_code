##############################################################################
# Depth-Weighted Aggregations of Polaris soils data
# (http://hydrology.cee.duke.edu/POLARIS)
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 06/26/2019
#
# Soils parameter units:
#    bd - bulk density, g/cm3
#    sand - sand percentage, %
#    clay - clay percentage, %
#    om - organic matter, log10(%)
#    ph - soil pH in H2O, N/A
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

library(raster)
library(rgdal)
library(doSNOW)
library(readr)

polaris_tiles <- read_csv("polaris_tiles.csv")
metrics <- c("clay", "sand", "bd", "om", "ph")
depths <- c("0_5", "5_15", "15_30", "30_60", "60_100")
wkdir <- "H:/Polaris/"

# Start parallel cluster of workers
cl <- snow::makeCluster(16, type = 'SOCK')
registerDoSNOW(cl)

# For each metric, make a raster stack of all the depths for each latlon tile
# Get the depth-weighted average for the stack, and
# Merge the average tiles into one raster and save as geotiff

dpave <- function(raslist, m){
  dpct <- c(0.05, 0.1, 0.15, 0.3, 0.4)
  trulist <- as.list(unlist(raslist))
  rastack <- stack(trulist)
  if(m == "om"){  #organic matter is in log10(%)
    tmpstack <- 10 ^ rastack
    tile <- sum(tmpstack * dpct)
  } else {
    tile <- sum(rastack * dpct)
  }
  return(tile)
}

for (m in metrics){
  prefix <- paste(wkdir, m, sep = "")

  x <- foreach(s=polaris_tiles$suffix) %:%
    foreach(d=depths, .combine = "c") %do%
      paste(prefix, d, s, sep = "_")

  y <- foreach(i=1:length(x), .packages = c("raster")) %dopar%
    dpave(x[i], m)

  y$filename <- paste(wkdir, m, "_dwave.tif", sep = "")
  y$format <- "GTiff"
  y$options <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=6", "BIGTIFF=YES")

  m <- do.call(merge, y)
}

# Release the cluster
snow::stopCluster(cl)
