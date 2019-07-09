# Generate presence and absence points for HOTR
library(sf)
library(dplyr)
library(spatstat)
library(maptools)
library(doSNOW)
spatstat.options(checkpolygons=FALSE)

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

setwd('D:/Projects/HOTR')
# shapefile of presence polygons. Already clipped to be within backpolys
pres_shp <- "pdog_prespolys.shp"
# shapefile of blocks where all points will be generated
back_shp <- "GridCells_in_bnd.shp"
# fieldname in backpolys that identifies each block
block_id <- "Grid_ID"

prespolys <- st_read(pres_shp, quiet = T)
backpolys <- st_read(back_shp, quiet = T)

# We want every prespoly to get at least 1 point, polys >= 2 ha get 1 pt/ha
prespolys$Sampnum <- floor(prespolys$Shape_Area / 10000)
prespolys <- mutate(prespolys,
                    Sampnum = if_else(Sampnum == 0, 1, Sampnum))

# We want the number of absence points to be within the same order of magnitude
tot_pres <- sum(prespolys$Sampnum)
# 3x extra created to allow for subsampling as necessary
tot_abs <- (10 ^ (ceiling(log10(tot_pres))) * 3)

# Do not want absence pts generated in presence polys
bpoly_clip <- st_erase(backpolys, prespolys)
bpoly_clip$Shape_Area <- st_area(bpoly_clip)
tot_area <- (sum(bpoly_clip$Shape_Area))
# Now that I've introduced units in the erase function,
# it wants to make Sampnum m^2, hence 'as.numeric'
bpoly_clip$Sampnum <- as.numeric(ceiling((bpoly_clip$Shape_Area / tot_area) * tot_abs))
bpoly_clip <- mutate(bpoly_clip,
                     Sampnum = if_else(Sampnum == 0, 1, Sampnum))

# *** Generate Presence Points ***
# Start parallel cluster of workers
cl <- snow::makeCluster(16, type = "SOCK")
registerDoSNOW(cl)

prespts <- foreach(p=1:nrow(prespolys),
                   .combine=bind_rows,
                   .packages=c("sf", "spatstat", "maptools", "dplyr")) %dopar%
  mysample(prespolys$geometry[[p]], prespolys$Sampnum[[p]])

prespts_out <- st_join(st_as_sf(prespts, crs = st_crs(backpolys)), st_as_sf(backpolys)) %>%
  select(block_id)
st_write(prespts_out, "BTPD_prespts.shp")

# *** Generate Absence Points ***
abspts <- foreach(p=1:nrow(bpoly_clip),
                  .combine=bind_rows,
                  .packages=c("sf", "spatstat", "maptools", "dplyr")) %dopar%
  mysample(bpoly_clip$geometry[[p]], bpoly_clip$Sampnum[[p]])

abspts <- st_join(st_as_sf(abspts, crs = st_crs(backpolys)), st_as_sf(backpolys)) %>%
  select(block_id)
st_write(abspts, "BTPD_abspts.shp")

# Release the cluster
snow::stopCluster(cl)
