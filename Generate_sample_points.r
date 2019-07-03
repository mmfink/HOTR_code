# Generate presence and absence points for HOTR
library(sf)
library(dplyr)
library(spatstat)
library(maptools)
library(doSNOW)
spatstat.options(checkpolygons=FALSE)

setwd('D:/GIS/Projects/HOTR')
pres_shp <- "pdog_prespolys.shp" #shapefile of presence polygons. All are within (spatially) the backpolys
back_shp <- "GridCells_in_bnd.shp" #shapefile of blocks where all points will be generated
block_id <- "Grid_ID" #fieldname in backpolys that identifies each block

prespolys <- st_read(pres_shp, quiet = T)
backpolys <- st_read(back_shp, quiet = T)

# We want every prespoly to get at least 1 point, polys >= 2 ha get 1 pt/ha
prespolys$Sampnum <- floor(prespolys$Shape_Area / 10000)
prespolys <- mutate(prespolys,
                    Sampnum = if_else(Sampnum == 0, 1, Sampnum))

mysample <- function(inpoly, samps){
  if(samps == 1){
    outpts <- as.data.frame(st_sample(inpoly, size = 1) %>% st_sfc(), col.names = c("geometry"))
  } else {
    p_owin <- as(as(inpoly, "Spatial"), "owin")
    # Each point needs to be at least 60m from its neighbors
    outppp <- rSSI(60, n = samps, win = p_owin, giveup = 100)
    outtmp <- as(outppp, "SpatialPoints")
    outpts <- as(outtmp, "sf")
  }
  return(outpts)
}

# *** Generate Presence Points ***
# Start parallel cluster of workers
cl <- snow::makeCluster(5, type = "SOCK", outfile = "out.txt")
registerDoSNOW(cl)

prespts <- foreach(p=1:nrow(prespolys),
                   .combine=rbind,
                   .packages=c("sf", "spatstat", "maptools", "dplyr")) %dopar%
  mysample(prespolys$geometry[[p]], prespolys$Sampnum[[p]])

# for(p in 1:nrow(prespolys)){
#   if(p == 1){
#     # poly to point just to have a starting point file to append to
#     # will probably want to delete it from final output
#     prespts <- as.data.frame(prespolys$geometry[[p]] %>% st_centroid() %>% st_sfc(), col.names = c("geometry"))
#   }
#   newpts <- mysample(prespolys$geometry[[p]], prespolys$Sampnum[[p]])
#   prespts <- st_sf(bind_rows(prespts, newpts))
# }

prespts <- st_join(prespts, backpolys) %>% select(block_id)
st_write(prespts, "BTPD_prespts.shp")

# *** Generate Absence Points ***
# We want the number of absence points to be within the same order of magnitude
tot_pres <- sum(prespolys$Sampnum)
tot_abs <- 10 ^ (ceiling(log10(tot_pres)))

# Not sure how this will work as an owin, but don't want background pts generated in towns
bpoly_clip <- st_difference(backpolys, prespolys)
bpoly_clip$Shape_Area <- st_area(bpoly_clip)
tot_area <- sum(bpoly_clip$Shape_Area)
bpoly_clip$Sampnum <- ceiling((bpoly_clip$Shape_Area / tot_area) * tot_abs)
bpoly_clip <- mutate(bpoly_clip,
                     Sampnum = if_else(Sampnum == 0, 1, Sampnum))

abspts <- foreach(p=1:nrow(bpoly_clip),
                  .combine=st_sf(bind_rows()),
                  .packages=c("sf", "spatstat", "maptools", "dplyr")) %dopar%
  mysample(bpoly_clip$geometry[[p]], bpoly_clip$Sampnum[[p]])

abspts <- st_join(abspts, backpolys) %>% select(block_id)
st_write(abspts, "BTPD_abspts.shp")

# Release the cluster
snow::stopCluster(cl)
