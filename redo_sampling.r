#############################################################################
# Re-do response table and the various datasets using the revised sampling
# cells and BTPD colony polygons, and make all data selections via sampling grids
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 11/19/2020
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
library(dplyr)
library(data.table)
library(raster)
library(dismo)
library(sf)
library(spatstat)
library(maptools)
library(doSNOW)
spatstat.options(checkpolygons=FALSE)
rasterOptions(tmpdir = "E:/MMF/temp")

pth <- "H:/HOTR_models"
setwd(pth)
set.seed(486)

st_erase <- function(x, y) {
  # Lower precision to allow this to complete sometime this century
  st_set_precision(x, units::set_units(1, cm))
  st_set_precision(y, units::set_units(1, cm))
  newst <- st_difference(x, st_union(st_combine(y)))
  return(newst)
}

mysample <- function(inpoly, samps){
  if(samps == 1){
    outpts <- as.data.frame(st_sample(inpoly, size = 1) %>%
                              st_sfc(), col.names = c("geometry"))
  } else {
    p_owin <- as(as(inpoly, "Spatial"), "owin")
    # Each point needs to be at least 60m from its neighbors
    outppp <- rSSI(60, n = samps, win = p_owin, giveup = 100)
    outtmp <- as(outppp, "SpatialPoints")
    outpts <- as.data.frame(st_as_sfc(outtmp), col.names = c("geometry"))
  }
  return(outpts)
}

# Function to extract raster values
myextract <- function(inpts, inrow){
  rasname <- inrow$raster
  inlbl <- inrow$label
  inras <- raster(rasname)
  x <- extract(inras, inpts)
  names(x) <- inlbl
  return(x)
}

### I) Re-do generation of sample points with the new polygon data ###

# shapefile of presence polygons. Already clipped to be within backpolys
pres_shp <- "pdog_polys_11172020.shp"
# shapefile of blocks where all points will be generated
back_shp <- "SamplingGrid2020_Nov16.shp"
# fieldname in backpolys that identifies each block
block_id <- "Grid_ID"

prespolys <- st_read(pres_shp, quiet = T)
backpolys <- st_read(back_shp, quiet = T)
# Only want to generate absence points in the censused blocks, not the
# grid cells containing ad hoc digitized colonies.
abspolys <- filter(backpolys, abs==1)

# Want every prespoly to get at least 1 point, polys >= 2 ha get 1 pt/ha
prespolys$Sampnum <- floor(prespolys$Shape_Area / 10000)
prespolys <- mutate(prespolys,
                    Sampnum = if_else(Sampnum == 0, 1, Sampnum))

# Want the number of absence points to be within the same order of magnitude
tot_pres <- sum(prespolys$Sampnum)
# 1.5x extra created to allow for removal of NAs and subsampling as necessary
tot_abs <- floor(tot_pres * 1.5)

# Do not want absence pts generated inside presence polys
abpoly_clip <- st_erase(abspolys, prespolys) #FIXME failed, used ArcGIS clip
# abpoly_clip <- st_read("abspoly_clip.shp", quiet = T)
abpoly_clip$Shape_Area <- st_area(abpoly_clip)
tot_area <- (sum(abpoly_clip$Shape_Area))
# Now that I've introduced units in the erase function,
# it wants to make Sampnum m^2, hence 'as.numeric'
abpoly_clip$Sampnum <- as.numeric(ceiling((abpoly_clip$Shape_Area / tot_area) * tot_abs))
abpoly_clip <- mutate(abpoly_clip,
                      Sampnum = if_else(Sampnum == 0, 1, Sampnum))

# *** Generate Presence Points ***
cl <- snow::makeCluster(16, type = "SOCK")
registerDoSNOW(cl)

prespts <- foreach(p=1:nrow(prespolys),
                   .combine=bind_rows,
                   .packages=c("sf", "spatstat", "maptools", "dplyr")) %dopar%
  mysample(prespolys$geometry[[p]], prespolys$Sampnum[[p]])

snow::stopCluster(cl)

prespts_out <- st_join(st_as_sf(prespts, crs = st_crs(backpolys)), st_as_sf(backpolys)) %>%
  dplyr::select(all_of(block_id))
st_write(prespts_out, "BTPD_prespts_Nov2020.shp")

# *** Generate Absence Points ***
cl <- snow::makeCluster(16, type = "SOCK")
registerDoSNOW(cl)

abspts <- foreach(p=1:nrow(abpoly_clip),
                  .combine=bind_rows,
                  .packages=c("sf", "spatstat", "maptools", "dplyr")) %dopar%
  mysample(abpoly_clip$geometry[[p]], abpoly_clip$Sampnum[[p]])

snow::stopCluster(cl)

abspts <- st_join(st_as_sf(abspts, crs = st_crs(backpolys)), st_as_sf(backpolys)) %>%
  dplyr::select(all_of(block_id))

# Remove any absence point that is within *500m* of a presence point.
##NOTE runs out of memory with raster::pointDistance
# thing <- pointDistance(abspts, prespts_out, lonlat = F)
# idx <- which(thing <= 500)
# ptsmod_idx <- unique(arrayInd(idx, dim(thing))[, 1])
# abspts_out <- abspts[-ptsmod_idx, ]
too_close <- st_is_within_distance(abspts, prespts_out, 500)
ptsmod_idx <- which(lengths(too_close) == 0)
abspts_out <- abspts[ptsmod_idx, ]
st_write(abspts_out, "BTPD_abspts_Nov2020.shp", driver = "ESRI Shapefile")

### II) Generate full response table ###
prespts_out$response <- rep(1, times = nrow(prespts_out))
abspts_out$response <- rep(0, times = nrow(abspts_out))
allpts <- as_Spatial(rbind(prespts_out, abspts_out))
allpts$ID <- seq_len(nrow(allpts))

# read in Environmental Inputs
inputs <- fread("env_inputs_04022020.csv", header=TRUE, sep=",")

cl <- snow::makeCluster(16, type = 'SOCK')
registerDoSNOW(cl)
# Run thru the inputs and extract values - this takes lots of time & memory
outvals <- foreach(r=1:nrow(inputs),
                   .packages = c("raster", "sp", "sf"),
                   .combine = "cbind") %dopar%
  myextract(as_Spatial(allpts), inputs[r, 1:2])

outdt <- data.table(outvals)

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
fwrite(as.data.frame(full_data), file = "response_Nov2020.csv")
full_data <- tidyr::drop_na(full_data)

### III) Split out into Training-Tweaking-Testing (70-15-15) datasets ###
test_pct <- 0.30 # proportion of response data to withhold from Training

# 1) Select all 'positive' Grid_IDs
# which means: From all pres pts, get unique Grid_IDs
pos_IDs <- as.character(unique(prespts_out$Grid_ID))

# 2) Remaining Grid_IDs are 'negative' cells
all_IDs <- (unique(full_data$Grid_ID))
posidx <- which(pos_IDs %in% all_IDs)
neg_IDs <- all_IDs[-posidx]

# 3) Select 70% of positive cells
train_pos <- sample(pos_IDs, length(pos_IDs)*(1-test_pct), replace = FALSE)
train_pospts <- dplyr::filter(full_data,
                              (Grid_ID %in% train_pos) &
                                (response==1))

# 4) Count how many absence points are in these cells,
#    then figure out how many more are needed to equal pres pts
absinpos <- nrow(dplyr::filter(full_data, (Grid_ID %in% train_pos) &
                                 (response==0)))
if(absinpos < nrow(train_pospts)){
  addabs <- nrow(train_pospts) - absinpos
  # 5) Randomly select 1 negative cell at a time, and add its
  #    absence points to the tally, until we have the right number
  newabs <- 0
  usecells <- vector("character")
  while(newabs < addabs){
    addcell <- sample(neg_IDs, 1)
    if(!addcell %in% usecells){
      usecells <- append(usecells, addcell)
      idx <- which((full_data$Grid_ID == addcell) &
                     (full_data$response==0))
      newabs <- newabs + length(idx)
    }
  }
}

train_neg <- append(train_pos, usecells)
train_negpts <- dplyr::filter(full_data,
                              (Grid_ID %in% train_neg) &
                                (response==0))

training_pts <- rbind(train_pospts, train_negpts)
fwrite(training_pts, file = "training_data11192020.csv")

# 6) Assign the 10-fold numbers and save as the Training data
folds <- 10
indata <- training_pts[, Grid_ID := as.factor(Grid_ID)]
kf <- kfold(levels(indata$Grid_ID), k=folds)
names(kf) <- levels(indata$Grid_ID)
kfdt <- data.table(Grid_ID=names(kf), KFsplit=kf, stringsAsFactors = TRUE)

# Add the Split column to the dataset
indata <- indata[kfdt, on=.(Grid_ID=Grid_ID)]
fwrite(indata, file = "training_kfold11192020.csv")

# Let's run some summary stats on the folds to check for weirdness.
library(ggplot2)
inroll <- rollup(indata, j = c(list(cnt=.N), lapply(.SD, mean)),
                 by = c("KFsplit", "response"), .SDcols = c(4:25))
subresp <- inroll[!is.na(response)]
setorder(subresp, KFsplit, response)
fwrite(subresp, file = "KFSplit_byResponse_Nov19.csv")

qplot(x=factor(KFsplit), data=indata, geom="bar", fill=factor(response),
      xlab="K-fold Split", ylab="Count",
      main="Presence/Absence per Split")

tograph <- melt(subresp, id.vars = c("KFsplit", "response"),
                measure.vars = c(4:25), variable.name = "metric")

qplot(factor(response), value, data=tograph,
      geom="boxplot", facets=.~factor(metric), color=factor(response),
      xlab="K-fold mean value spread for each metric by response")

# 7) Of the remaining 30% positive cells from 3), select half
#    and repeat steps 4 & 5.
remaining_pos <- pos_IDs[-which(train_pos %in% pos_IDs)]
remaining_neg <- neg_IDs[-which(train_neg %in% neg_IDs)]

twk_pos <- sample(remaining_pos, length(remaining_pos)/2, replace = FALSE)
twk_pospts <- dplyr::filter(full_data,
                              (Grid_ID %in% twk_pos) &
                                (response==1))

absinpos <- nrow(dplyr::filter(full_data, (Grid_ID %in% twk_pos) &
                                 (response==0)))
if(absinpos < nrow(twk_pospts)){
  addabs <- nrow(twk_pospts) - absinpos
  newabs <- 0
  usecells <- vector("character")
  while(newabs < addabs){
    addcell <- sample(remaining_neg, 1)
    if(!addcell %in% usecells){
      usecells <- append(usecells, addcell)
      idx <- which((full_data$Grid_ID == addcell) &
                     (full_data$response==0))
      newabs <- newabs + length(idx)
    }
  }
}

twk_neg <- append(twk_pos, usecells)
twk_negpts <- dplyr::filter(full_data,
                              (Grid_ID %in% twk_neg) &
                                (response==0))

# 8) Save this dataset as the Tweaking data
tweak_pts <- rbind(twk_pospts, twk_negpts)
fwrite(tweak_pts, file = "tweaking_data11192020.csv")

# 9) Repeat with the remaining 15% positive cells
test_pos <- remaining_pos[-which(twk_pos %in% remaining_pos)]
rest_neg <- remaining_neg[-which(twk_neg %in% remaining_neg)]

test_pospts <- dplyr::filter(full_data,
                              (Grid_ID %in% test_pos) &
                                (response==1))

absinpos <- nrow(dplyr::filter(full_data, (Grid_ID %in% test_pos) &
                                 (response==0)))

if(absinpos < nrow(test_pospts)){
  addabs <- nrow(test_pospts) - absinpos
  newabs <- 0
  usecells <- vector("character")
  while(newabs < addabs){
    addcell <- sample(rest_neg, 1)
    if(!addcell %in% usecells){
      usecells <- append(usecells, addcell)
      idx <- which((full_data$Grid_ID == addcell) &
                     (full_data$response==0))
      newabs <- newabs + length(idx)
    }
  }
}

test_neg <- append(test_pos, usecells)
test_negpts <- dplyr::filter(full_data,
                              (Grid_ID %in% test_neg) &
                                (response==0))

# 10) Save this dataset as the Testing data
test_pts <- rbind(test_pospts, test_negpts)
fwrite(test_pts, file = "testing_data11192020.csv")
