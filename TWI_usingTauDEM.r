##############################################################################
# Compute Topographic Wetness Index (TWI) using TauDEM
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 06/04/2019
#
# This file is a part of the Homes On The Range code repository.
#
# *** Requires installation of TauDEM v5 ***
# Which is Copyright (C) 2010-2015 David Tarboton, Utah State University
# http://hydrology.usu.edu/taudem/taudem5/
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

#### Set Variables ####
# Inputs
data.dir <- "Path/To/Folder" #location of both inputs & outputs
dem <- "HOTR_dem_clip.tif" #input elevation
pitmask <- "PitMask.tif" #use a mask to prevent pitremove from filling legit sinks
num.cores <- "20" #number of cores to use in processing

# Outputs
dem.filled <- paste(data.dir, "HOTR_DEM_filled.tif", sep = "/")
flowdir <- paste(data.dir, "HOTR_Dinf_flowdir.tif", sep = "/")
dem.slope <- paste(data.dir, "HOTR_Dinf_slope.tif", sep = "/")
area.upslope <- paste(data.dir, "HOTR_Dinf_accum.tif", sep = "/")
out.twi <- paste(data.dir, "HOTR_TWI.tif", sep = "/")
########

# Function to alter areas of zero slope, which cause NoData in TauDEM
fun <- function(x){
  x[x==0] <- 0.01 # using 0.01 on advice of David.Augustine@ARS.USDA.GOV
  return(x)
}

dem.name <- paste(data.dir, dem, sep = "/")

# Pitremove (fill the DEM) with a pitmask
pitm.name <- paste(data.dir, pitmask, sep = "/")
system(paste("mpiexec -n", num.cores, "pitremove -z", dem.name, "-depmask", pitm.name, "-fel", dem.filled))

# Flow direction and slope using D-infinity algorithm
system(paste("mpiexec -n", num.cores, "DinfFlowdir -ang", flowdir, "-slp", dem.slope, "-fel", dem.filled))

# Upslope accumulation using D-infinity algorithm
# Use the -nc flag to turn off checking for edge contamination.
# Prevents edge-caused NoData, however results for areas that would have been NoData
# will be *incorrect* to some degree because of incomplete upstream area.
system(paste("mpiexec -n", num.cores, "AreaDinf -ang", flowdir, "-sca", area.upslope, "-nc"))

# Adjust slope - prevents slope-caused NoData
slope.ras <- raster(dem.slope)
slope2 <- paste(data.dir, "HOTR_Dinf_slope2.tif", sep = "/") #won't let me overwrite file so have to make a second

beginCluster()
slope.mod <- clusterR(slope.ras, calc, args = list(fun=fun), filename = slope2)
endCluster()

#TWI
system(paste("mpiexec -n", num.cores, "TWI -slp", slope2, "-sca", area.upslope, "-twi", out.twi))
