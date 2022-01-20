#=============================================================================================================#
# Script created by Luke Musher, contact at ljm357@drexel.edu
# Script created in version R 3.3.2 
# Adapted from scripts from Glenn Seeholzer (Seeholzer and Brumfield 2017) as well as
#     from Prunier et al. 2014 
# This script: runs multivariate logistic regressions testing for isolation-by-distance, 
#			isolation-by-environment, and isolation-by-barrier in Amazonian birds Where 
#			dispersal distance (least cost-path) is used in place of geographic (euclidiean) distance.  
#			Runs a binomial GLM and commonality analysis that subsamples and removes duplicate localities
#     assesing the relative important of these processes in shaping 
#			genetic structure of Willisornis
#
# Bootstrap procedure modified from code in supplementary of Prunier et al. (2015), which was based on methods in Peterman et al. (2014)
#
#	Prunier JG, Colyn M, Legendre X, Nimon KF, Flamand MC (2015) Multicollinearity in spatial genetics: Separating the wheat from the chaff using commonality analyses. Molecular Ecology, 24, 263–283.
#	Peterman WE, Connette GM, Semlitsch RD, Eggert LS (2014) Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology, 23, 2402–2413.
#
# Usage notes: run line by line
#=============================================================================================================#
rm(list=ls())
#set working directory
library(fmsb)
library(descr) # addon for Commonality analysis
library(plyr)
library(StAMPP)
library(gdistance)
library(ecodist)
library(maptools)
library(hierfstat)
library(pcadapt)
source("~/AMNH/Roosevelt_Project/Data/CAs/supporting.functions/CA_logistic.R")
source("~/AMNH/Roosevelt_Project/Data/CAs/supporting.functions/CA_semiStdCoef.R")

setwd('~/AMNH/Roosevelt_Project/Data/CAs/Willisornis/')

DXY.raw<-read.delim("W_po_distances_12Jan2022.csv")

# ####RAW MAXENT output CURRENT, generate dispersal distance matrix
# r = raster('~/AMNH/Roosevelt_Project/Niche_models/Malacoptila/MalacoptilaModel/CurrentMalaFull.tif')
# rnew = r
#tr = transition(rnew,mean,directions=16)
#trC = geoCorrection(tr,'c',scl=T)
#save(trC,file='cost.trC.raw.Rdata')
#get(load('./cost.trC.raw.Rdata'))

#load coordinates and other files
data = read.csv('Indivs.Table.Willisornis.csv',stringsAsFactors=F)
#loc <- as.matrix(data[,c('long','lat')])
#dcost = costDistance(trC,loc)
#cost.raw = as.matrix(dcost)
#write.table(cost.raw,file='cost.raw.txt',row.names=T,col.names=T,quote=F,sep='\t')
cost.raw<-read.table('cost.raw.txt')

#load barrier data
bar = read.delim('mbar.willisornis.txt')

####	ENV
####	import bioclim and enviro rasters
files <- list.files(path="../BioClim/", pattern='tif', full.names=TRUE)
BIOCLIM = stack(files) #create raster stack of BIOCLIM rasters
envdata = data.frame(pop=data[,'Code'],as.data.frame(extract(BIOCLIM, data[,c('long','lat')]))) #extract BIOCLIM data for each individual
tmp = aggregate(. ~ pop, data= envdata,mean,na.rm=T) #take mean for each locality
my.prc = prcomp(tmp[,-1], center=T, scale = T, na.omit=T) #do DXY
pc.values = predict(my.prc)#[,1:ncol(tmp[,-1])]
pop.env.raw = data.frame(pc.values)
rownames(pop.env.raw) = tmp[,1]

# ####	Generate geographic DISTANCE MATRIX
# tmp = aggregate(. ~ Code, data= data[,c('Code','long','lat')],mean,na.rm=T)
# pop.loc = as.matrix(tmp[,c('long','lat')])
# rownames(pop.loc) = tmp[,'Code']
# DIST.raw = spDists(pop.loc,longlat=TRUE) ; colnames(DIST.raw) = rownames(DIST.raw) = rownames(pop.loc)

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
prep <- function(x){ x=data.matrix(x)
x=x[upper.tri(x, diag = FALSE)]
x=(x-mean(x))/sqrt(var(x)) }

lg <- function(x){ exp(x)/(1 + exp(x)) }

################################################################################################
####	Raw Distance Matrices
################################################################################################
mDXY = as.matrix(dist(DXY.raw))
mDIST = as.matrix(cost.raw)
mbar = as.matrix(bar)
mEnv = as.matrix(dist(pop.env.raw))

################################################################################################
####	Distance matrices formatted for Commonality Analysis
################################################################################################
####	z-transformation for commonality analysis
DXY = as.numeric(prep(mDXY))

# identify end of first "valley" in the kernal density
# this is the threshold above which we define successful divergence
# below this threshold we consider values equivalent to failure to diverge

hist(DXY, prob=T, col="skyblue1",breaks=75); rug(DXY)
lines(density(DXY, adj=.25),lwd=3)
abline(v=-0.65,col="red",lwd=2,lty=2)
DXY[DXY> -0.65]<-1
DXY[DXY< -0.65]<-0

### z-transform
DIST = as.numeric(prep(mDIST))
ENV = as.numeric(prep(mEnv))
RIV = as.numeric(prep(mbar))

# LRDM and CA
temp2=data.frame(DXY=DXY,DIST=DIST,ENV=ENV,RIV=RIV)
ca=cc4log(temp2,"DXY",list("DIST","ENV","RIV"), "N") 
ca

glmt=glm(DXY ~ DIST+ENV+RIV, family=binomial("logit"),data=temp2)
summary(glmt)
pseudoR=NagelkerkeR2(glmt)
pseudoR$R2

# Plot
minTOT=min(min(ENV),min(RIV),min(DIST))
maxTOT=max(max(ENV),max(RIV),max(DIST))
Zscore=seq(minTOT,maxTOT,0.1)
p1=rep(mean(ENV),length(Zscore));
p2=rep(mean(RIV),length(Zscore));
p0=rep(mean(DIST),length(Zscore));
for (i in c(1,2,3)){
  if(i==1){Newdata = data.frame(DIST=Zscore,ENV=p1,RIV=p2)
  prGLM <- predict(glmt,newdata=Newdata,type="link",se=TRUE)
  plot(Zscore,lg(prGLM$fit),type="l",ylim=c(0,1),col="black",lty=1,lwd=3)
  par(new=T)}
  else if(i==2){Newdata = data.frame(DIST=p0,ENV=Zscore,RIV=p2)
  prGLM <- predict(glmt,newdata=Newdata,type="link",se=TRUE)
  plot(Zscore,lg(prGLM$fit),type="l",ylim=c(0,1),col="darkred",lty=1,lwd=3)
  par(new=T)}
  else if(i==3){Newdata = data.frame(DIST=p0,ENV=p1,RIV=Zscore)
  prGLM <- predict(glmt,newdata=Newdata,type="link",se=TRUE)
  plot(Zscore,lg(prGLM$fit),type="l",ylim=c(0,1),col="steelblue3",lty=1,lwd=3)
  par(new=T)}
}	

stdBeta=rbind(0,0,0)
stdBeta2=rbind(0,0,0)
for (i in c(2:4)){
  Newdata = data.frame(DIST=DIST,ENV=ENV,RIV=RIV)
  prGLM <- predict(glmt,newdata=Newdata,type="response",se=TRUE)
  mprob=mean(prGLM$fit)
  mprob2=0.5
  stdBeta[i-1]=(1/(1+exp(-(log(mprob/(1-mprob))+0.5*glmt$coefficients[i]))))-(1/(1+exp(-(log(mprob/(1-mprob))-0.5*glmt$coefficients[i]))))
  stdBeta2[i-1]=(1/(1+exp(-(log(mprob2/(1-mprob2))+0.5*glmt$coefficients[i]))))-(1/(1+exp(-(log(mprob2/(1-mprob2))-0.5*glmt$coefficients[i]))))
}
stdBeta
exp(stdBeta)

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
n=ncol(mDXY)
sn=0.9*n

for (i in 1:nperm){
  rarray = sort(sample(n,sn,replace=F))
  mmDXY = prep(mDXY[rarray,rarray])
  mmDIST = prep(mDIST[rarray,rarray])
  mmENV = prep(mEnv[rarray,rarray])
  mmRIV = prep(mbar[rarray,rarray])
  mmDXY[mmDXY> -0.65]<-1
  mmDXY[mmDXY< -0.65]<-0
  prout=data.frame(DXY=mmDXY,DIST=mmDIST,ENV=mmENV,RIV=mmRIV)
  comm=cc4log(prout,"DXY",list("DIST","ENV","RIV"), "N") 
  boot[i,]=comm$CC[c(1:ncombos),1]
  print(i)
}

for (i in 1:ncombos){
  q=quantile(boot[,i], c(.025,.975))
  resBootstrap[i,1]=i
  resBootstrap[i,3]=q[1]
  resBootstrap[i,4]=q[2]
}

resBootstrap[,2]=ca$CC[c(1:ncombos),1]
resBootstrap[,5]=ca$CC[c(1:ncombos),2]
rownames(resBootstrap) = rownames(ca$CC)[1:ncombos]
resBootstrap[,c('o','l','u')] = resBootstrap[,c('o','l','u')]

new.rownames = gsub('\\s|[?!Unique$|Common$|to$|and$]','',rownames(resBootstrap))
rownames(resBootstrap) = new.rownames

write.csv(resBootstrap,file = "resBootstrap_DXY.csv")

pdf('Willis_CA_logistic_DXY.pdf',width=8,height=5, bg='transparent')
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
