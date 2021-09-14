setwd("E:/Dropbox/TRIMpaper_2021/shapefiles")
library(sf)
library(dplyr)
species.1 <- read.table("specieslist.tsv", header = F, sep = "\t")[,1]



freshwater1 <- read_sf(dsn = ".", layer = "FW_FISH_PART1")
freshwater1 <- freshwater1 %>% filter(binomial %in% species.1)

freshwater2 <- read_sf(dsn = ".", layer = "FW_FISH_PART2")
freshwater2 <- freshwater2 %>% filter(binomial %in% species.1)

marine1 <- read_sf(dsn = ".", layer = "MARINEFISH_PART1")
marine1 <- marine1 %>% filter(binomial %in% species.1)

marine2 <- read_sf(dsn = ".", layer = "MARINEFISH_PART2")
marine2 <- marine2 %>% filter(binomial %in% species.1)

marine3 <- read_sf(dsn = ".", layer = "MARINEFISH_PART3")
marine3 <- marine3 %>% filter(binomial %in% species.1)

sharkrays <- read_sf(dsn = ".", layer = "SHARKS_RAYS_CHIMAERAS")
sharkrays <- sharkrays %>% filter(binomial %in% species.1)

sharkrays <- read_sf(dsn = ".", layer = "SHARKS_RAYS_CHIMAERAS")
sharkrays <- sharkrays %>% filter(binomial %in% species.1)

hagfish <- read_sf(dsn = ".", layer = "HAGFISH")
hagfish <- hagfish %>% filter(binomial %in% species.1)

marine <- rbind(marine1, marine2, marine3)
freshwater <- rbind(freshwater1, freshwater2)

rm(freshwater1, freshwater2, marine1, marine2, marine3)
cat(unique(marine$binomial), unique(freshwater$binomial), unique(sharkrays$binomial)
    
spat.dat.all <- unique(rbind(freshwater, marine, hagfish, sharkrays))

## lamprey data is present in the database, but not available for download apparently.
# 80 species present, of 140 possible

specieslist.iucn <- sort(unique(spat.dat.all$binomial))

minmax <- c()
for (i in specieslist.iucn){
  tmpdat <- spat.dat.all %>% filter(binomial == i)
  minmax <- rbind(minmax,st_bbox(tmpdat))
}

minmax <- as.data.frame(unlist(minmax))
minmax$species <- as.character(specieslist.iucn)

minmax <- minmax[,c(5,2,4,1,3)]



tmpa <- spat.dat.all$binomial
tmpb <- spat.dat.all$marine
tmpc <- spat.dat.all$freshwater

tmpd <- as.data.frame(t(unlist(rbind(tmpa, tmpb, tmpc))), stringsAsFactors = F, row.names = F)
colnames(tmpd) <- c("species","mfish", "ffish")
tmpd <- unique(tmpd)

minmax <- merge(minmax, tmpd)
write.table(minmax, file = "iucn.data.tsv")