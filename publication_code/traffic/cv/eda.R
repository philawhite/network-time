library(ggmap)
library(ggplot2)
# library(maps)
# library(mapdata)
# library(sf)
# library(lubridate)
# library(scales)
library(rgdal)

dat = read.csv("grid_data_215.csv")[,-1]
dat$loc_ind = dat$grid_loc

dist_mat = as.matrix(rdist(dat$milepoint))
time_mat = as.matrix(rdist(dat$time_numeric))/7

cord.dec = SpatialPoints(cbind(dat$long, dat$lat), proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:2850"))
cord.UTM@coords[,1] = (cord.UTM@coords[,1] - min(cord.UTM@coords[,1]))/1000
cord.UTM@coords[,2] = (cord.UTM@coords[,2] - min(cord.UTM@coords[,2]))/1000

euc_dist = as.matrix(rdist(cbind(cord.UTM@coords[,1],cord.UTM@coords[,2])))*0.62137119224


scord.dec = SpatialPoints(unique(cbind(dat$long, dat$lat)), proj4string = CRS("+proj=longlat"))
scord.UTM <- spTransform(scord.dec, CRS("+init=epsg:2850"))
scord.UTM@coords[,1] = (scord.UTM@coords[,1] - min(scord.UTM@coords[,1]))/1000
scord.UTM@coords[,2] = (scord.UTM@coords[,2] - min(scord.UTM@coords[,2]))/1000

small_euc =  as.matrix(rdist(cbind(scord.UTM@coords[,1],scord.UTM@coords[,2])))*0.62137119224

small_gc = geosphere::distm(as.matrix(unique(dat[,c("long","lat")])),fun = distCosine)/1000*0.62137119224
small_net = as.matrix(rdist(unique(dat$milepoint)))

plot(small_net,small_euc)

dat$x = cord.UTM@coords[,1]
dat$y = cord.UTM@coords[,2]


ggplot(dat) +
  geom_point(aes(long,lat,col = as.factor(milepoint)))


