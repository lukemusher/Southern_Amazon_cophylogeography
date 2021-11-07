require("ape")
library(phytools)
library(ape)
library(phangorn)
library(ips)
library(parallel)

#must have exactly same genes and trees--if there are 20 genes there must be 20 genes 
#and named appropriately so you know you have the right tree with each gene
#ideally you will have uces and gts named exactly the same e.g. uce-10.nexus uce-10.tre
#the easiest way to do this is make a list of the genes, and use a for loop to change
# the names of the RAxML outputs to just the locus name
#e.g. in linux/unix: for i in `cat names.genes`; do mv RAxML*$i* $i.tre; done
#runs slowly on datasets with lots of tips
#Start fresh R to maximize efficiency

setwd("~/AMNH/Roosevelt_Project/Data/genetrees/willisornis/PreFiltered/")

files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees <- list.files(path="./",pattern="*.tre",full.names=F,recursive=FALSE)

outgroup<-"reference"

pis<-function(x){
  x<-read.nexus.data(x)
  len<-length(x[[1]])
  site<-c()
  counts=0
  for(j in 1:length(x[[1]])){
    for(i in 1:length(x)){
      site[i]<-(x[[i]][j])
    }
    if(!("n" %in% unique(site))&!("-" %in% unique(site))){
      if (length(unique(site))>=2){
        counts=counts+1
      }
    }
    else if(!("n" %in% unique(site))){
      if (length(unique(site))>=3){
        counts=counts+1
      }
    }
    else if(!("-" %in% unique(site))){
      if (length(unique(site))>=3){
        counts=counts+1
      }
    }
    else if (length(unique(site))>4){
      counts=counts+1
    }
    else counts=counts
  }
  return(counts)
}

site_count<-function(x){
  x<-read.nexus.data(x)
  len<-length(x[[1]])
  return(len)
}

is.clocklike<-function(dna,tree){
  print(dna)
  x<-read.nexus.data(dna)
  x<-as.DNAbin(x)
  X<-as.phyDat(x)
  t<-c()
  counts=0
  tree<-read.tree(tree)
  for(i in 1:length(tree$tip.label)){
    if(!(tree$tip.label[i] %in% labels(X))){
      counts=counts+1
      t[counts]<-tree$tip.label[i]
    }
  }
  t<-t[!is.na(t)]
  tree1<-drop.tip(tree,tree$tip.label[na.omit(match(t,tree$tip.label))])
  tree1<-force.ultrametric(tree1, method="nnls")
  fit<-pml(tree1, X, model="GTRG")   #do likelihood estimation of the distance tree 
  chr<-tryCatch(chronopl(phy=root(tree1,outgroup), lambda=1.0, age.min=1, CV=T),error=function(e) NULL)
  fit2<-tryCatch(pml(chr, X, model="GTRG"), error=function(e) "NA")   #do likelihood estimation of the distance tree 
  if(fit2!="NA"){
    LR=(2*(fit$logLik-fit2$logLik))
  }
  if(fit2=="NA"){
    LR=NULL
  }
  return(LR)
}

get.likelihood<-function(dna,tree){
  x<-read.nexus.data(dna)
  x<-as.DNAbin(x)
  X<-as.phyDat(x)
  t<-c()
  counts=0
  tree<-read.tree(tree)
  for(i in 1:length(tree$tip.label)){
    if(!(tree$tip.label[i] %in% labels(X))){
      counts=counts+1
      t[counts]<-tree$tip.label[i]
    }
  }
  t<-t[!is.na(t)]
  tree1<-drop.tip(tree,tree$tip.label[na.omit(match(t,tree$tip.label))])
  fit<-pml(tree1, X, model="GTRG")   #do likelihood estimation of the distance tree 
  return(fit$logLik)
}

no_cores <- 3
cl <- makeCluster(no_cores)
clusterExport(cl, list("read.nexus.data","pis","site_count","is.clocklike"))
var.sites<-parLapply(cl,files, pis)
save(var.sites,file = "../var.sites.Rdata")
seq.length<-parLapply(cl,files,site_count)
save(seq.length,file = "../seq.length.Rdata")
stopCluster(cl)

LR<-mcmapply(is.clocklike,files,trees,mc.cores = 3)
save(LR,file = "../LR.Rdata")

lh<-mcmapply(get.likelihood,files,trees,mc.cores=3)
save(lh,file = "../likelihoods.Rdata")

LR<-get(load("../LR.Rdata"))
var.sites<-get(load("../var.sites.Rdata"))
seq.length<-get(load("../seq.length.Rdata"))
lh<-get(load("../likelihoods.Rdata"))

list2vector<-function(x){
  y<-c()
  for (i in 1:length(x)){
    if (x[i]!="NULL"){
      y[i]<-x[[i]]
    }
  }
  return(y)
}

var.sites.v<-list2vector(var.sites)
seq.length.v<-list2vector(seq.length)
var.frac.v<-var.sites.v/seq.length.v
LR.v<-list2vector(LR)
lh.v<-list2vector(lh)

df<-as.data.frame(cbind(files,trees,lh.v,seq.length.v,var.frac.v,var.sites.v,LR.v))
colnames(df)<-c("alignment","tree","logLikelihood","seq.length","var.frac","var.sites","LR")
write.csv(df,file="../Phleg_gt_metrics.csv")
