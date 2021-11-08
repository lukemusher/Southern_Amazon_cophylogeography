#=============================================================================================================#
# Script created by Glenn Seeholzer, contact at seeholzer.glenn@gmail.com
# Script created in version R 3.3.2 
# This script: runs univariate and multivariate matrix regressions testing for isolation-by-adaptation, 
#			isolation-by-environment, and isolation-by-distance in Cranioleuca antisiensis. Where 
#			geographic distance is used in place of dispersal distance.  
#			Runs a commonality analysis assesing the relative important of these processes in shaping 
#			genetic structure of C. antisiensis
#			Produces Figure S8 & Table S2
# Usage notes: run line by line
#=============================================================================================================#
rm(list=ls())
#set working directory
library(plyr)
library(StAMPP)
library(gdistance)
library(ecodist)
library(maptools)
library(hierfstat)
library(pcadapt)
library(tidyverse)
source('~/AMNH/Roosevelt_Project/Data/CAs/supporting.functions/FUN.add.alpha.R', chdir = TRUE)
source('~/AMNH/Roosevelt_Project/Data/CAs/supporting.functions/MMRR.R', chdir = F)
setwd('~/AMNH/Roosevelt_Project/Data/CAs/Hypocnemus/')

# ####RAW MAXENT output CURRENT
r = raster('~/AMNH/Roosevelt_Project/Niche_models/Hypocnemis/HypocModel/CurrentHypocFull.tif')
#rnew = r
#tr = transition(rnew,mean,directions=16)
#trC = geoCorrection(tr,'c',scl=T)
#save(trC,file='cost.trC.raw.Rdata')
get(load('./cost.trC.raw.Rdata'))

#load coordinates and other files
data = read.csv('Indivs.Table.Hypocnemus.csv',stringsAsFactors=F)
#loc <- as.matrix(data[,c('long','lat')])
#dcost = costDistance(trC,loc)
#cost.raw = as.matrix(dcost)
#write.table(cost.raw,file='cost.raw.txt',row.names=T,col.names=T,quote=F,sep='\t')
cost.raw<-read.table('cost.raw.txt')

#average TSNE distances
setwd('~/AMNH/Roosevelt_Project/Data/CAs/Hypocnemus/H_st_TSNE/')
files <- list.files(path="./", pattern="*.csv", full.names=F, recursive=FALSE)

tbl <- sapply(files, read_csv, simplify=FALSE)

t<-data.frame(matrix(0, nrow = 60, ncol = 60))

for (i in 1:length(files)){
  d<-data.frame(as.matrix(dist(as.data.frame(tbl[i])[,2:3])))
  t<-d+t
}

TSNE.raw<-t/length(files)
write.csv(TSNE.raw,file="./TSNE.raw.csv")

#load coordinates and other files
setwd('~/AMNH/Roosevelt_Project/Data/CAs/Hypocnemus/')

df <- data.frame(pops = data$Code, ID = data$Code, lat = data$lat, lon = data$long)
df <- df[order(df$ID),]
head(df$ID)


bar = read.delim('mbar.hypoc.txt')

detach(package:tidyr)
detach("package:tidyverse", unload=TRUE)


####	ENV
####	import bioclim and enviro rasters
files <- list.files(path="../BioClim/", pattern='tif', full.names=TRUE)
BIOCLIM = stack(files) #create raster stack of BIOCLIM rasters
envdata = data.frame(pop=data[,'Code'],as.data.frame(extract(BIOCLIM, data[,c('long','lat')]))) #extract BIOCLIM data for each individual
tmp = aggregate(. ~ pop, data= envdata,mean,na.rm=T) #take mean for each locality
my.prc = prcomp(tmp[,-1], center=T, scale = T, na.omit=T) #do PCA
pc.values = predict(my.prc)#[,1:ncol(tmp[,-1])]
pop.env.raw = data.frame(pc.values)
rownames(pop.env.raw) = tmp[,1]

####	DISTANCE MATRIX
#tmp = aggregate(. ~ Code, data= data[,c('Code','long','lat')],mean,na.rm=T)
#pop.loc = as.matrix(tmp[,c('long','lat')])
#rownames(pop.loc) = tmp[,'Code']
#DIST.raw = spDists(pop.loc,longlat=TRUE) ; colnames(DIST.raw) = rownames(DIST.raw) = rownames(pop.loc)

#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------ Create distance matrices of population data for predictors ------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#

# Libraries
library(yhat) # Commonality analysis

#-transformation: functions for formatting predictor matrices and standardizing the values so coefficients comparable as beta-weights
standardize = function(x){tmp=(x-mean(x,na.rm=T))/sd(x,na.rm=T);diag(tmp)=0;return(tmp)}

prep.standardize <- function(x){x=data.matrix(x)
x=x[lower.tri(x, diag = FALSE)]
x=(x-mean(x))/sqrt(var(x)) 
x}

prep <- function(x){x=data.matrix(x)
x=x[lower.tri(x, diag = FALSE)]
x}

################################################################################################
####	Raw Distance Matrices
################################################################################################
mTSNE = as.matrix(TSNE.raw)
mTSNE[upper.tri(mTSNE)] = mTSNE[lower.tri(mTSNE)]
mDIST = as.matrix(cost.raw)
mbar = as.matrix(bar)
mEnv = as.matrix(dist(pop.env.raw))

################################################################################################
####	Lists of distance matrices formatted for Multiple Matrix Regression on Distance Matrices (MMRR, Wang 2013)
################################################################################################
#MMRR to get regression coefficients and significance values for multivariate model
#data for MMRR

#Do Commonality Analysis to get Unique, Common, and Total contribution of each predictor variable
#for univariate and regressions
predictors = list(DIST=mDIST,ENV=mEnv,RIV=mbar)
predictors = lapply(predictors,standardize)
x = predictors
resBootstrap.list = list()
fit = MMRR(Y = mTSNE, X = x, nperm=1000)

################################################################################################
####	Distance matrices formatted for Commonality Analysis
################################################################################################
####	for commonality analysis
TSNE = as.numeric(prep.standardize(mTSNE))
DIST = as.numeric(prep.standardize(mDIST))
ENV = as.numeric(prep.standardize(mEnv))
RIV = as.numeric(prep(mbar))

data=data.frame(TSNE = TSNE, DIST = DIST, ENV = ENV, RIV = RIV)
ca = regr(lm(TSNE ~ DIST + ENV + RIV,data=data))

############################
############################
#Commonality analysis For Body Mass with confidence intervals
############################
############################
#bootstrap procedure modified from code in supplementary of Prunier et al. (2015), which was based on methods in Peterman et al. (2014)

#	Prunier JG, Colyn M, Legendre X, Nimon KF, Flamand MC (2015) Multicollinearity in spatial genetics: Separating the wheat from the chaff using commonality analyses. Molecular Ecology, 24, 263–283.
#	Peterman WE, Connette GM, Semlitsch RD, Eggert LS (2014) Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology, 23, 2402–2413.

nperm=1000
n.predictors = 3
ncombos = ((2^n.predictors)-1)
boot=matrix(data = 0, nrow = nperm, ncol = ncombos)
resBootstrap=data.frame(n=rep(0,ncombos),o=rep(0,ncombos),l=rep(0,ncombos),u=rep(0,ncombos),p=rep(0,ncombos))
n=ncol(mTSNE)
sn=0.9*n

for (i in 1:nperm){
  rarray = sort(sample(n,sn,replace=F))
  mmTSNE = mTSNE[rarray,rarray][lower.tri(mTSNE[rarray,rarray],diag=F)]
  mmDIST = prep(mDIST[rarray,rarray])
  mmENV = prep(mEnv[rarray,rarray])
  mmRIV = prep(mbar[rarray,rarray])
  
  comm=regr(lm(mmTSNE ~ mmDIST + mmENV + mmRIV))
  boot[i,]=comm$Commonality_Data$CC[c(1:ncombos),1]
  print(i)
}

for (i in 1:ncombos){
  q=quantile(boot[,i], c(.025,.975))
  resBootstrap[i,1]=i
  resBootstrap[i,3]=q[1]
  resBootstrap[i,4]=q[2]
}

resBootstrap[,2]=ca$Commonality_Data$CC[c(1:ncombos),1]
resBootstrap[,5]=ca$Commonality_Data$CC[c(1:ncombos),2]
rownames(resBootstrap) = rownames(ca$Commonality_Data$CC)[1:ncombos]
resBootstrap[,c('o','l','u')] = resBootstrap[,c('o','l','u')]

new.rownames = gsub('\\s|[?!Unique$|Common$|to$|and$]','',rownames(resBootstrap))
rownames(resBootstrap) = new.rownames

write.csv(resBootstrap,file = "resBootstrap_TSNE.csv")

pdf('Hypocnemis_CA_TSNE.pdf',width=8,height=5, bg='transparent')
tmp = resBootstrap
tmp = round(tmp,3)
tmp$n = rev(tmp$n)

par(mar=c(4,2,4,20))
xlim = c(min(tmp[,3]), max(tmp[,4]))
plot(tmp[,2],tmp[,1],xlim=xlim,font=5,lab=c(ncombos, 7, 1),xaxt="n",yaxt="n",cex.lab=1,xlab='',ylab='',cex=1.5,col="steelblue3", pch = 16) 
arrows(tmp[,3],tmp[,1],tmp[,4],tmp[,1],code=3,length=0.05,angle=90,lwd=1)
abline(v=0,lty=2)
axis(1,cex.axis=1)
mtext('Correlation Coefficient',1,line=2.5,cex=1.2)
#yaxis
x=cbind(1:15,rev(rownames(tmp)))
axis(4,at=x[,1],labels=x[,2],cex.axis=1,las=2)
mtext('Predictor Sets',2,cex=1.2,line=.5)
#upper xaxis
seq(0,round(max(tmp[,4]),2),by=round(max(tmp[,4]),2)/5)
at = seq(0,round(max(tmp[,4]),2),by=round(max(tmp[,4]),2)/5)
labels= round(at/sum(tmp[,'o']),2)
axis(3,at=at,labels=labels,cex.axis=1)
mtext('% Total',3,line=2.5,cex=1.2)
#coefficients
mtext('Coefficient',4,at = 7.6,las=2,line=12,adj=.5)
axis(4,line=12,at=tmp[,1],labels=tmp[,'o'],tick=F,cex.axis=1,las=2,hadj=1)
mtext('% Total',4,at = 7.6,las=2,line=17,adj=.5)
axis(4,line=17,at=tmp[,1],labels=tmp[,'p'],tick=F,cex.axis=1,las=2,hadj=1)

mtext('Total',4,at = 0.2,las=2,line=8,adj=.5)
mtext(sum(tmp[,'o']),4,at = 0.2,las=2,line=12,adj=.5)
mtext(100,4,at = 0.2,las=2,line=17,adj=.5)

dev.off()
