##### Raster Creation #####
library(data.table)
library(raster)
library(rgdal)
library(doSNOW)

pth <- "H:/HOTR_models/"
setwd(pth)
rasterOptions(tmpdir = "E:/MMF/temp")
#Testing on a single model
BRT1 <- "BRT_Oct28_full.rds"
ncores <- (parallel::detectCores()) - 2
gopts <- c("COMPRESS=LZW", "TFW=YES", "NUM_THREADS=10", "BIGTIFF=YES")

some_inputs <- fread("env_inputs_10022019.csv", header = TRUE, sep = ",")
add_inputs <- fread("onehot_rasters.csv", header = TRUE, sep = ",")
env_inputs <- rbind(some_inputs[, .(raster, label)], add_inputs[, .(raster, label)])

# pair down list to only those inputs used in each model.
mod <- readRDS(paste0(pth,BRT1))
subdf <- env_inputs[label %in% mod$var.names]

# make raster stack
layerStk <- raster::stack(subdf$raster, quick=TRUE)
names(layerStk) <- subdf$label

bz <- blockSize(layerStk) #5,075 blocks of 15 rows
slab <- 203 #25 hunks of 203 blocks

out <- raster(env_inputs$raster[1]) #whyyyy do I have to load a raster to write another raster?
out <- writeStart(out, filename = "BRT_Oct28_full_prob.tif", format="GTiff", overwrite=TRUE, options = gopts)

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)

for(j in seq.int(1, 5076, by=203)){
  oy <- foreach(i=j:(j+(slab-1)), .packages = c("raster", "gbm"),
                .combine = "rbind") %dopar% {
                  chunk <- getValues(layerStk, bz$row[i], bz$nrows[i])
                  pred <- gbm::predict.gbm(mod, as.data.frame(chunk),
                                           n.trees=4998, type="response", single.tree = TRUE)
                }
  for(r in j:(j+(slab-1))){
    out <- writeValues(out, oy, bz$row[r])
  }
}

snow::stopCluster(cl)

out <- writeStop(out)
