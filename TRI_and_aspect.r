##############################################################################
# Terrain Ruggedness Index and calculation of aspect in terms of
# 'northness' and 'eastness'
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 07/10/2019
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

#### Set Variables ####
#Inputs
data.dir <- "H:/DEM_etc" #location of both inputs & outputs
dem <- "HOTR_DEM_filled.tif" #input elevation - use filled

#Outputs
out.north <- paste(data.dir, "HOTR_Northness.tif", sep = "/")
out.east <- paste(data.dir, "HOTR_Eastness.tif", sep = "/")
out.tri <- paste(data.dir, "HOTR_TRI.tif", sep = "/")
########

dem.ras <- raster(paste(data.dir, dem, sep = "/"))
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=6", "BIGTIFF=YES")

#### Aspect ####
aspect <- terrain(dem.ras, opt = "aspect", unit = "radians", neighbors = 8)
northness <- cos(aspect)
nf <- writeRaster(northness, filename = out.north, format = "GTiff", options = gopts)
eastness <- sin(aspect)
ef <- writeRaster(eastness, filename = out.east, format = "GTiff", options = gopts)

#### Terrain Ruggedness Index ####
tri <- terrain(dem.ras, opt = "TRI", filename = out.tri, format = "GTiff",
               options = gopts)