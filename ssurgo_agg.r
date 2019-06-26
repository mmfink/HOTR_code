##############################################################################
# Depth and Area - Weighted Aggregations of SSURGO data
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 01/16/2019
#
# This file is a part of the Homes On The Range code repository.
# Adapted from
#https://casoilresource.lawr.ucdavis.edu/software/r-advanced-statistical-package/aggregating-ssurgo-data-r/
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

library(plyr)

# Function for calculating a weighted mean
# values passed in should be vectors of equal length
wt_mean <- function(property, weights){
  # compute thickness weighted mean, but only when we have enough data
  # otherwise return NA

  # Record which indices have data vs NA
  property.that.is.na <- which(is.na(property))
  property.that.is.not.na <- which(!is.na(property))

  if(length(property) - length(property.that.is.na) >= 1){
    prop.aggregated <-
      sum(weights[property.that.is.not.na] * property[property.that.is.not.na], na.rm=T) / sum(weights[property.that.is.not.na], na.rm=T)
  } else {
    prop.aggregated <- NA

  return(prop.aggregated)
  }
}

ph_mean <- function(property, weights){
  # weighted mean specific to pH
  property.that.is.na <- which(is.na(property))
  property.that.is.not.na <- which(!is.na(property))

  if(length(property) - length(property.that.is.na) >= 1){
    antilog.aggregated <-
      sum(weights[property.that.is.not.na] * (10^property[property.that.is.not.na]), na.rm=T) / sum(weights[property.that.is.not.na], na.rm=T)
    prop.aggregated <- log10(antilog.aggregated)
  } else {
    prop.aggregated <- NA
  }
}

profile_total <- function(property, thickness){
  # compute profile sum for properties such as AWS or return NA

  # Record which indices have data vs NA
  property.that.is.na <- which(is.na(property))
  property.that.is.not.na <- which(!is.na(property))

  if(length(property) - length(property.that.is.na) >= 1){
    prop.aggregated <-
      sum(thickness[property.that.is.not.na] * property[property.that.is.not.na], na.rm=T)
  } else {
    prop.aggregated <- NA
  }

  return(prop.aggregated)
}

# Function to perform horizon-thickness weighted aggregation
component_level_aggregation <- function(i, depth=100){

  # horizon thickness is our weighting vector
  if(depth > 0){
    hz.use <- which(i$hzdept_r < depth-2) #can't get exact depths, just in the neighborhood
  } else {
    hz.use <- which(i$hzdept_r)
  }

  # Figure out if we have any bedrock & compute min depth
  hz.br <- which((i$desgnmaster == 'R') & (!is.na(i$hzdept_r)))
  if(length(hz.br) == 0){
    depth_br <- NA_integer_
  } else {
    depth_br <- min(i$hzdept_r[hz.br], na.rm = T)
  }

  hz_thick <- i$hzdepb_r[hz.use] - i$hzdept_r[hz.use]

  # compute wt.mean aggregate values
  clay <- wt_mean(i$claytotal_r[hz.use], hz_thick)
  silt <- wt_mean(i$silttotal_r[hz.use], hz_thick)
  sand <- wt_mean(i$sandtotal_r[hz.use], hz_thick)
  orgm <- wt_mean(i$om_r[hz.use], hz_thick)
  ph <- ph_mean(i$ph1to1h2o_r[hz.use], hz_thick)

  # compute profile sum values (not used)
  #water_storage <- profile_total(i$awc_r[hz.use], hz_thick)

  # make a new dataframe out of the aggregate values
  d <- data.frame(cokey=unique(i$cokey), clay=clay, silt=silt, sand=sand, organic=orgm, pH=ph, depth_bedrock=depth_br)

  return(d)
}

mapunit_level_aggregation <- function(i){
  # component percentage is our weighting vector
  comppct <- i$comppct_r

  # wt. mean by component percent
  clay <- wt_mean(i$clay, comppct)
  silt <- wt_mean(i$silt, comppct)
  sand <- wt_mean(i$sand, comppct)
  orgm <- wt_mean(i$organic, comppct)
  ph <- ph_mean(i$pH, comppct)
  #water_storage <- wt_mean(i$water_storage, comppct)
  depth_br <- wt_mean(i$depth_bedrock, comppct)

  # make a new dataframe out of the aggregate values
  d <- data.frame(mukey=unique(i$mukey), clay=clay, silt=silt, sand=sand, organic=orgm, pH=ph, depth_bedrock=depth_br)

  return(d)
}
#------
# set a maximum depth to consider
maxdep <- 100 #use 0 for all possible depths

# load horizon and component data
chorizon <- read.csv('D:/GIS/Projects/HOTR/chorizon.csv')

# only keep some of the columns from the component & chtexturegrp tables
component <- read.csv('D:/GIS/Projects/HOTR/component.csv')[,c('mukey','cokey','comppct_r')]

# aggregate horizon data to the component level
chorizon.agg <- ddply(chorizon, .(cokey), .fun = component_level_aggregation, .progress='text')

# join up the aggregate chorizon data to the component table
comp.merged <- merge(component, chorizon.agg, by='cokey')

# aggregate component data to the map unit level
component.agg <- ddply(comp.merged, .(mukey), .fun=mapunit_level_aggregation, .progress='text')

# add the aggregated depth to bedrock b/c chorizon bedrock is mostly null
muagg <- read.csv('D:/GIS/Projects/HOTR/muaggatt.csv')[,c('mukey','brockdepmin')]
rock.data <- which(!is.na(muagg$brockdepmin))
final.df <- join(component.agg, muagg[rock.data,], by = 'mukey', match = 'first')

# save data back to CSV
write.csv(final.df, file='soil_data.csv', row.names=FALSE)

# save data back to CSV
write.csv(component.agg, file='soil_data.csv', row.names=FALSE)
