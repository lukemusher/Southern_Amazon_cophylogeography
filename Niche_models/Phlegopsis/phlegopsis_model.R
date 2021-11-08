## Author: Peter J. Galante
## This loads occurrence records, creates a background region, sets up, and performs model tuning and selection.
## Then models are projected to time periods and thresholded using the ETSS. Finally they are projected in raw output for use in "timeSlicePlot.R"
setwd("~/AMNH/Roosevelt_Project/Niche_models/Phlegopsis/")
library(spThin);library(ENMeval);library(rgdal);library(rgeos)
bio<-stack(list.files(path = '~/AMNH/GIS_layers/Worldclim_bioclimatic/Current_wc2.1_30s_bio/', pattern = '\\.tif$', full.names = T))
##optimize function
optimize <- function(res) {
  ###Remove any candidate model which has an AUC less than 0.51= models with no discrimination
  opt.auc <- res[res$train.AUC >= 0.5,]
  ###Remove any candidates which have no parameters
  no.param <- opt.auc[opt.auc$parameters > 1,]
  ###Remove any candidates where the AIC score was NA (too many parameters)
  noAICNA<- no.param[!is.na(no.param$parameters),]
  #noAICNA<- no.param[which(!is.na(no.param$AICc)),]
  ###Remove any models which have an OR of zero
  noOR0 <- noAICNA[noAICNA$avg.test.or10pct != 0,]
  ###Order the remaining list by lowest OR then highest AUC, sequentially
  ordered<-noOR0[with(noOR0, order(avg.test.or10pct, -avg.test.AUC)), ]
  ###Grab the settings of that first model (the optimal model)
  ordered[1,]
}
######################################################################
#####################  phlegopsis  ##################################
######################################################################
## phlegopsis - thin points and select dataset
s1<-read.csv('./phleg_thinned_thin1.csv')[,2:3]
s2<-read.csv('./phleg_thinned_thin2.csv')[,2:3]
s3<-read.csv('./phleg_thinned_thin3.csv')[,2:3]
s4<-read.csv('./phleg_thinned_thin4.csv')[,2:3]
s5<-read.csv('./phleg_thinned_thin5.csv')[,2:3]

length(is.na(extract(bio[[1]], s1))[is.na(extract(bio[[1]], s1))==T]) # 0
length(is.na(extract(bio[[1]], s2))[is.na(extract(bio[[1]], s2))==T]) # 0
length(is.na(extract(bio[[1]], s3))[is.na(extract(bio[[1]], s3))==T]) # 0
length(is.na(extract(bio[[1]], s4))[is.na(extract(bio[[1]], s4))==T]) # 0
length(is.na(extract(bio[[1]], s5))[is.na(extract(bio[[1]], s5))==T]) # 0

# load points
randomThin <- sample(1:5,1) # 4
phlegopsis<-read.csv('./phleg_thinned_thin1.csv')[,2:3]
# Generate background data
mcp <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}
##  Set a 10km buffer distance
buff.dist <- 5000
## Creating A MCP background for thinned phlegopsis data
ca.MCP<-mcp(phlegopsis)
crs(ca.MCP)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ca.MCP<-spTransform(ca.MCP, CRSobj = crs('+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'))
ca.shp <- gBuffer(ca.MCP, width = buff.dist)
ca.shp<-spTransform(ca.shp, CRSobj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ca.env.crp<-crop(bio, ca.shp)
ca.env<-mask(ca.env.crp, ca.shp)
plot(ca.env[[1]])
ca.backg<-randomPoints(ca.env, 10000)
write.csv(ca.backg, 'phlegopsisBackg.csv', row.names = F)
caback<-read.csv('phlegopsisBackg.csv')


## Tune in ENMeval

#Sys.setenv((NOAWT=TRUE)) #only for certain macs before loading rJava
library(rJava)
options(java.parameters = "-Xmx12g")
library(dismo)
library(raster)
library(maptools)
library(GISTools)

cares<-ENMevaluate(occ = phlegopsis, env = ca.env, 
                   bg.coords = caback, RMvalues = seq(1,6,0.5), 
                   fc = c("L", "LQ", "H", "LQH", "LQHP"),
                   method = 'block', rasterPreds = T, parallel = F, 
                   algorithm = "maxent.jar", numCores = 1)

write.csv(cares@results, "phlegopsisResults.csv", row.names=F)
Phleg<-read.csv('phlegopsisResults.csv')
PhlegOpt<-optimize(Phleg) # H_1
#dir.create("phlegopsisModel")
PhlegMod<-maxent(x = ca.env,
               p = phlegopsis,
               a = caback,
               path = "phlegopsisModel",
               args=c(
                 'betamultiplier=1',
                 'linear=false',
                 'quadratic=false',
                 'product=false',
                 'threshold=false',
                 'hinge=true',
                 'threads=2',
                 'responsecurves=true',
                 'jackknife=true',
                 'askoverwrite=false'
               )
)
saveRDS(PhlegMod, "phlegopsisModel/phlegopsisMaxent.RDS")
PhlegMod<-readRDS("phlegopsisModel/phlegopsisMaxent.RDS")

## Predict current to full extent
bio<-stack(list.files(path = '~/AMNH/GIS_layers/Worldclim_bioclimatic/Current_wc2.1_30s_bio/', pattern = '\\.tif$', full.names = T))
colnames(phlegopsis) <- c("x","y")
e1 <- c(extent(phlegopsis)[1]-2, extent(phlegopsis)[2]+2, extent(phlegopsis)[3]-2, extent(phlegopsis)[4]+2)
bio <- crop(bio, e1)

PhlegCurFull<-predict(
  object = PhlegMod,
  x = bio,
  filename = "./phlegopsisModel/CurrentPhlegFull.tif",
  na.rm=T,
  format = "GTiff",
  overwrite=T,
  args = 'cloglog'
)

LGM <- stack(list.files(path = "~/AMNH/GIS_layers/LGM_ccsm4/", pattern = "\\.tif$", full.names = T))
LGM <- crop(LGM, e1)

PhlegCurFull<-predict(
  object = PhlegMod,
  x = LGM,
  filename = "./phlegopsisModel/LGMPhlegFull.tif",
  na.rm=T,
  format = "GTiff",
  overwrite=T,
  args = 'cloglog'
)


