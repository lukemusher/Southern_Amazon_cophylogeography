library(bayestestR)
library(dplyr)
library(ggplot2)

setwd("~/AMNH/Roosevelt_Project/Data/GPhocs/")

minVal <- 5E-10
maxVal <- 7E-09

m<-rnorm(n = 1000000, mean = 4.6E-09,sd = 1E-9)
m <- pmax(minVal, m)
m <- pmin(maxVal, m)
#m<-runif(n=1000000, min = 1E-9, max = 7E-9)
g<-rnorm(n=1000000, mean = 2,sd = 0.25)

par(mfrow=c(1,1))
hist(m,col='black', main = "sampled mutation rates")
hist(g,col='black', main = "sampled generation times")

tau_to_time<-function(tau){
  t<-tau*2/(4.6E-09)
  return(t)
}

theta_to_Ne<-function(theta){
  Ne<-sample(theta, replace = T,size = 1)/(4*(sample(m,replace = T,size=1)))
  return(Ne)
}

tau_to_time_rand<-function(tau){
  t<-(sample(tau, replace = T,size = 1)*sample(g,replace = T, size=1))/(sample(m,replace = T,size=1))
  return(t)
}

tau_to_time_rand_mut<-function(tau){
  t<-sample(tau, replace = T,size = 1)*2/(sample(m,replace = T,size=1))
  return(t)
}

#######read data

gcy<-read.csv("G_cy_var_mcmc.csv",header=T)
#l<-length(gcy_var$tau_root)
#plot(1:l,gcy$Full.ld.ln)
#plot(1:l,tau_to_time(gcy$tau_root))

hst<-read.csv("H_st_var_mcmc.csv",header=T)
#l<-length(hst$tau_root)
#plot(1:l,hst$Full.ld.ln)
#plot(1:l,tau_to_time(hst$tau_root))

mru<-read.csv("M_ru_var_mcmc.csv",header=T)
#l<-length(mru$tau_root)
#plot(1:l,mru$Full.ld.ln)
#plot(1:l,tau_to_time(mru$tau_root))

tae<-read.csv("T_ae_var_mcmc.csv",header=T)
#l<-length(tae$tau_root)
#plot(1:l,tae$Full.ld.ln)
#plot(1:l,tau_to_time(tae$tau_root))

pni<-read.csv("P_ni_var_mcmc.csv",header=T)
#l<-length(pni$tau_root)
#plot(1:l,pni$Full.ld.ln)
#plot(1:l,tau_to_time(pni$tau_root))

wpo<-read.csv("W_po_var_mcmc.csv",header=T)
#l<-length(wpo$tau_root)
#plot(1:l,wpo$Full.ld.ln)
#plot(1:l,tau_to_time(wpo$tau_root))

#estimate uncertainty from mut + gen time
gcy_root<-replicate(n = 100000,expr = tau_to_time_rand(gcy$tau_root))
gcy_tap<-replicate(n = 100000,expr = tau_to_time_rand(gcy$tau_Para_MaRo))
gcy_pur<-replicate(n = 100000,expr = tau_to_time_rand(gcy$tau_Inam_Puru))
gcy_arro<-replicate(n = 100000,expr = tau_to_time_rand(gcy$tau_Para_MaRo_ArSuTa))
#plot(density(gcy_root))

mru_root<-replicate(n = 100000,expr = tau_to_time_rand(mru$tau_root))
mru_tap<-replicate(n = 100000,expr = tau_to_time_rand(mru$tau_Para_MaArTa))
mru_mad<-replicate(n = 100000,expr = tau_to_time_rand(mru$tau_Puru_JiMa))
mru_pur<-replicate(n = 100000,expr = tau_to_time_rand(mru$tau_Inam_Puru_JiMa))
#plot(density(mru_root))

hst_root<-replicate(n = 100000,expr = tau_to_time_rand(hst$tau_root))
hst_tap<-replicate(n = 100000,expr = tau_to_time_rand(hst$tau_Para_Mach))
hst_arro<-replicate(n = 100000,expr = tau_to_time_rand(hst$tau_Para_Mach_RoArSuTa))
#plot(density(hst_root))

tae_root<-replicate(n = 100000,expr = tau_to_time_rand(tae$tau_root))
tae_mad<-replicate(n = 100000,expr = tau_to_time_rand(tae$tau_InPu_JiMaArSuTa))
tae_arro<-replicate(n = 100000,expr = tau_to_time_rand(tae$tau_JiMaRoArSuTa))
#plot(density(tae_root))

pni_root<-replicate(n = 100000,expr = tau_to_time_rand(pni$tau_root))
pni_tap<-replicate(n = 100000,expr = tau_to_time_rand(pni$tau_Para_Rondonia))
pni_pur<-replicate(n = 100000,expr = tau_to_time_rand(pni$tau_InPu))
pni_arro<-replicate(n = 100000,expr = tau_to_time_rand(pni$tau_Rondonia))
#plot(density(pni_root))

wpo_root<-replicate(n = 100000,expr = tau_to_time_rand(wpo$tau_root))
wpo_tap<-replicate(n = 100000,expr = tau_to_time_rand(wpo$tau_Para_SuTa))
wpo_mad<-replicate(n = 100000,expr = tau_to_time_rand(wpo$tau_InPu_Rondonia))
wpo_sut<-replicate(n = 100000,expr = tau_to_time_rand(wpo$tau_Para_SuTa))
#plot(density(wpo_root))


#estimate uncertainty from mut only
gcy_root_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(gcy$tau_root))
gcy_tap_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(gcy$tau_Para_MaRo))
gcy_pur_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(gcy$tau_Inam_Puru))
gcy_arro_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(gcy$tau_Para_MaRo_ArSuTa))
#plot(density(gcy_root))

mru_root_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(mru$tau_root))
mru_tap_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(mru$tau_Para_MaArTa))
mru_mad_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(mru$tau_Puru_JiMa))
mru_pur_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(mru$tau_Inam_Puru_JiMa))
#plot(density(mru_root))

hst_root_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(hst$tau_root))
hst_tap_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(hst$tau_Para_Mach))
hst_arro_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(hst$tau_Para_Mach_RoArSuTa))
#plot(density(hst_root))

tae_root_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(tae$tau_root))
tae_mad_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(tae$tau_InPu_JiMaArSuTa))
tae_arro_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(tae$tau_JiMaRoArSuTa))
#plot(density(tae_root))

pni_root_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(pni$tau_root))
pni_tap_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(pni$tau_Para_Rondonia))
pni_pur_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(pni$tau_InPu))
pni_arro_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(pni$tau_Rondonia))
#plot(density(pni_root))

wpo_root_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(wpo$tau_root))
wpo_tap_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(wpo$tau_Para_SuTa))
wpo_mad_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(wpo$tau_InPu_Rondonia))
wpo_sut_mut<-replicate(n = 100000,expr = tau_to_time_rand_mut(wpo$tau_Para_SuTa))
#plot(density(wpo_root))

############################################
#####################Ne#####################
############################################

par(mfrow=c(1,1))

#plot Ne assuming 4.6E-9
plot(c(1,13),c(0,0),ylim = c(0,1.25E6), type = 'n', main = "θ by population", xlab = "species", ylab = "τ", xaxt='n')
axis(1, at=c(2,4,6,8,10,12), labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1,10,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1.4,mean(theta_to_Ne(gcy$theta_Inam)),pch = 16)
c<-ci(theta_to_Ne(gcy$theta_Inam),method = "HDI",ci = 0.95)
arrows(1.4, c[[2]],1.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(1.7,mean(theta_to_Ne(gcy$theta_Puru)),pch = 16)
c<-ci(theta_to_Ne(gcy$theta_Puru),method = "HDI",ci = 0.95)
arrows(1.7, c[[2]],1.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(2,mean(theta_to_Ne(gcy$theta_MaRo)),pch = 16)
c<-ci(theta_to_Ne(gcy$theta_MaRo),method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(2.3,mean(theta_to_Ne(gcy$theta_ArSuTa)),pch = 16)
c<-ci(theta_to_Ne(gcy$theta_ArSuTa),method = "HDI",ci = 0.95)
arrows(2.3, c[[2]],2.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(2.7,mean(theta_to_Ne(gcy$theta_Para)),pch = 16)
c<-ci(theta_to_Ne(gcy$theta_Para),method = "HDI",ci = 0.95)
arrows(2.7, c[[2]],2.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(3.4,mean(theta_to_Ne(mru$theta_Inam)),pch = 16)
c<-ci(theta_to_Ne(mru$theta_Inam),method = "HDI",ci = 0.95)
arrows(3.4, c[[2]],3.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(3.7,mean(theta_to_Ne(mru$theta_Puru)),pch = 16)
c<-ci(theta_to_Ne(mru$theta_Puru),method = "HDI",ci = 0.95)
arrows(3.7, c[[2]],3.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(4,mean(theta_to_Ne(mru$theta_JiMa)),pch = 16)
c<-ci(theta_to_Ne(mru$theta_JiMa),method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(4.3,mean(theta_to_Ne(mru$theta_MaArTa)),pch = 16)
c<-ci(theta_to_Ne(mru$theta_MaArTa),method = "HDI",ci = 0.95)
arrows(4.3, c[[2]],4.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(4.6,mean(theta_to_Ne(mru$theta_Para)),pch = 16)
c<-ci(theta_to_Ne(mru$theta_Para),method = "HDI",ci = 0.95)
arrows(4.6, c[[2]],4.6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(5.4,mean(theta_to_Ne(hst$theta_Puru)),pch = 16)
c<-ci(theta_to_Ne(hst$theta_Puru),method = "HDI",ci = 0.95)
arrows(5.4, c[[2]],5.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(5.7,mean(theta_to_Ne(hst$theta_JiGu)),pch = 16)
c<-ci(theta_to_Ne(hst$theta_JiGu),method = "HDI",ci = 0.95)
arrows(5.7, c[[2]],5.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(6,mean(theta_to_Ne(hst$theta_Mach)),pch = 16)
c<-ci(theta_to_Ne(hst$theta_Mach),method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(6.3,mean(theta_to_Ne(hst$theta_RoArSuTa)),pch = 16)
c<-ci(theta_to_Ne(hst$theta_RoArSuTa),method = "HDI",ci = 0.95)
arrows(6.3, c[[2]],6.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(6.6,mean(theta_to_Ne(hst$theta_Para)),pch = 16)
c<-ci(theta_to_Ne(hst$theta_Para),method = "HDI",ci = 0.95)
arrows(6.6, c[[2]],6.6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(7.4,mean(theta_to_Ne(tae$theta_InPu)),pch = 16)
c<-ci(theta_to_Ne(tae$theta_InPu),method = "HDI",ci = 0.95)
arrows(7.4, c[[2]],7.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(7.7,mean(theta_to_Ne(tae$theta_JiMaRo)),pch = 16)
c<-ci(theta_to_Ne(tae$theta_JiMaRo),method = "HDI",ci = 0.95)
arrows(7.7, c[[2]],7.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(8,mean(theta_to_Ne(tae$theta_ArSuTa)),pch = 16)
c<-ci(theta_to_Ne(tae$theta_ArSuTa),method = "HDI",ci = 0.95)
arrows(8, c[[2]],8, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(8.3,mean(theta_to_Ne(tae$theta_Para)),pch = 16)
c<-ci(theta_to_Ne(tae$theta_Para),method = "HDI",ci = 0.95)
arrows(8.3, c[[2]],8.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(9.4,mean(theta_to_Ne(pni$theta_Inam)),pch = 16)
c<-ci(theta_to_Ne(pni$theta_Inam),method = "HDI",ci = 0.95)
arrows(9.4, c[[2]],9.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(9.7,mean(theta_to_Ne(pni$theta_Puru)),pch = 16)
c<-ci(theta_to_Ne(pni$theta_Puru),method = "HDI",ci = 0.95)
arrows(9.7, c[[2]],9.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(10,mean(theta_to_Ne(pni$theta_JiMaRo)),pch = 16)
c<-ci(theta_to_Ne(pni$theta_JiMaRo),method = "HDI",ci = 0.95)
arrows(10, c[[2]],10, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(10.3,mean(theta_to_Ne(pni$theta_ArSuTa)),pch = 16)
c<-ci(theta_to_Ne(pni$theta_ArSuTa),method = "HDI",ci = 0.95)
arrows(10.3, c[[2]],10.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(10.6,mean(theta_to_Ne(pni$theta_Para)),pch = 16)
c<-ci(theta_to_Ne(pni$theta_Para),method = "HDI",ci = 0.95)
arrows(10.6, c[[2]],10.6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(11.4,mean(theta_to_Ne(wpo$theta_InPu)),pch = 16)
c<-ci(theta_to_Ne(wpo$theta_InPu),method = "HDI",ci = 0.95)
arrows(11.4, c[[2]],11.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(11.7,mean(theta_to_Ne(wpo$theta_Rondonia)),pch = 16)
c<-ci(theta_to_Ne(wpo$theta_Rondonia),method = "HDI",ci = 0.95)
arrows(11.7, c[[2]],11.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(12,mean(theta_to_Ne(wpo$theta_SuTa)),pch = 16)
c<-ci(theta_to_Ne(wpo$theta_SuTa),method = "HDI",ci = 0.95)
arrows(12, c[[2]],12, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

points(12.3,mean(theta_to_Ne(wpo$theta_Para)),pch = 16)
c<-ci(theta_to_Ne(wpo$theta_Para),method = "HDI",ci = 0.95)
arrows(12.3, c[[2]],12.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#######################
par(mfrow=c(2,2))
#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,0.009), type = 'n', main = "τ(Root)", xlab = "species", ylab = "τ", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1,10,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(gcy$tau_root),pch = 16)
c<-ci(gcy$tau_root,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru$tau_root),pch = 16)
c<-ci(mru$tau_root,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst$tau_root),pch = 16)
c<-ci(hst$tau_root,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae$tau_root),pch = 16)
c<-ci(tae$tau_root,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni$tau_root),pch = 16)
c<-ci(pni$tau_root,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo$tau_root),pch = 16)
c<-ci(wpo$tau_root,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,4E6), type = 'n', main = "µ = 4.6E-09, G = 2", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,500000,1000000,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(tau_to_time(gcy$tau_root)),pch = 16)
c<-ci(tau_to_time(gcy$tau_root),method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(tau_to_time(mru$tau_root)),pch = 16)
c<-ci(tau_to_time(mru$tau_root),method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(tau_to_time(hst$tau_root)),pch = 16)
c<-ci(tau_to_time(hst$tau_root),method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tau_to_time(tae$tau_root)),pch = 16)
c<-ci(tau_to_time(tae$tau_root),method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(tau_to_time(pni$tau_root)),pch = 16)
c<-ci(tau_to_time(pni$tau_root),method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(tau_to_time(wpo$tau_root)),pch = 16)
c<-ci(tau_to_time(wpo$tau_root),method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)


#################################
#################################
#plot root times for each species with uncertainty for mut only
plot(c(1,6),c(0,0),ylim = c(0,5E6), type = 'n', main = "Uncertainty in µ", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,750000,1500000,2250000), labels=c(0,0.75,1.5,2.25))

#Galbula
points(1,mean(gcy_root_mut),pch = 16)
c<-ci(gcy_root_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_root_mut),pch = 16)
c<-ci(mru_root_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_root_mut),pch = 16)
c<-ci(hst_root_mut,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_root_mut),pch = 16)
c<-ci(tae_root_mut,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_root_mut),pch = 16)
c<-ci(pni_root_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo_root_mut),pch = 16)
c<-ci(wpo_root_mut,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species with uncertainty
plot(c(1,6),c(0,0),ylim = c(0,4E6), type = 'n', main = "Uncertainty in µ and G", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1500000,3000000,4500000), labels=c(0,1.5,3,4.5))

#Galbula
points(1,mean(gcy_root),pch = 16)
c<-ci(gcy_root,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_root),pch = 16)
c<-ci(mru_root,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_root),pch = 16)
c<-ci(hst_root,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_root),pch = 16)
c<-ci(tae_root,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_root),pch = 16)
c<-ci(pni_root,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo_root),pch = 16)
c<-ci(wpo_root,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

par(mfrow=c(4,4))

#################################
#################################
#################################
########TAPAJOS BARRIER##########
#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,0.0045), type = 'n', main = "τ(Tapajós)", xlab = "species", ylab = "τ", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1,10,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(gcy$tau_Para_MaRo),pch = 16)
c<-ci(gcy$tau_Para_MaRo,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru$tau_Para_MaArTa),pch = 16)
c<-ci(mru$tau_Para_MaArTa,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst$tau_Para_Mach),pch = 16)
c<-ci(hst$tau_Para_Mach,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae$tau_root),pch = 16)
c<-ci(tae$tau_root,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni$tau_Para_Rondonia),pch = 16)
c<-ci(pni$tau_Para_Rondonia,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo$tau_Para_SuTa),pch = 16)
c<-ci(wpo$tau_Para_SuTa,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "µ = 4.6E-9, G = 2", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,500000,1000000,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(tau_to_time(gcy$tau_Para_MaRo)),pch = 16)
c<-ci(tau_to_time(gcy$tau_Para_MaRo),method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(tau_to_time(mru$tau_Para_MaArTa)),pch = 16)
c<-ci(tau_to_time(mru$tau_Para_MaArTa),method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(tau_to_time(hst$tau_Para_Mach)),pch = 16)
c<-ci(tau_to_time(hst$tau_Para_Mach),method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tau_to_time(tae$tau_root)),pch = 16)
c<-ci(tau_to_time(tae$tau_root),method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(tau_to_time(pni$tau_Para_Rondonia)),pch = 16)
c<-ci(tau_to_time(pni$tau_Para_Rondonia),method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(tau_to_time(wpo$tau_Para_SuTa)),pch = 16)
c<-ci(tau_to_time(wpo$tau_Para_SuTa),method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)


#################################
#################################
#plot root times for each species with uncertainty for mut only
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "Uncertainty in µ", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,750000,1500000,2250000), labels=c(0,0.75,1.5,2.25))

#Galbula
points(1,mean(gcy_tap_mut),pch = 16)
c<-ci(gcy_tap_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_tap_mut),pch = 16)
c<-ci(mru_tap_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_tap_mut),pch = 16)
c<-ci(hst_tap_mut,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_root_mut),pch = 16)
c<-ci(tae_root_mut,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_tap_mut),pch = 16)
c<-ci(pni_tap_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo_tap_mut),pch = 16)
c<-ci(wpo_tap_mut,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species with uncertainty
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "Uncertainty in µ and G", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1500000,3000000,4500000), labels=c(0,1.5,3,4.5))

#Galbula
points(1,mean(gcy_tap),pch = 16)
c<-ci(gcy_tap,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_tap),pch = 16)
c<-ci(mru_tap,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_tap),pch = 16)
c<-ci(hst_tap,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_root),pch = 16)
c<-ci(tae_root,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_tap),pch = 16)
c<-ci(pni_tap,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo_tap),pch = 16)
c<-ci(wpo_tap,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#############################
#############################
#############################
#####MADEIRA BARRIER#########
#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,0.0045), type = 'n', main = "τ(Madeira)", xlab = "species", ylab = "τ", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1,10,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(gcy$tau_root),pch = 16)
c<-ci(gcy$tau_root,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru$tau_Puru_JiMa),pch = 16)
c<-ci(mru$tau_Puru_JiMa,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst$tau_root),pch = 16)
c<-ci(hst$tau_root,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae$tau_InPu_JiMaArSuTa),pch = 16)
c<-ci(tae$tau_InPu_JiMaArSuTa,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni$tau_root),pch = 16)
c<-ci(pni$tau_root,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo$tau_InPu_Rondonia),pch = 16)
c<-ci(wpo$tau_InPu_Rondonia,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "µ = 4.6E-9, G = 2", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,500000,1000000,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(tau_to_time(gcy$tau_root)),pch = 16)
c<-ci(tau_to_time(gcy$tau_root),method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(tau_to_time(mru$tau_Puru_JiMa)),pch = 16)
c<-ci(tau_to_time(mru$tau_Puru_JiMa),method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(tau_to_time(hst$tau_root)),pch = 16)
c<-ci(tau_to_time(hst$tau_root),method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tau_to_time(tae$tau_InPu_JiMaArSuTa)),pch = 16)
c<-ci(tau_to_time(tae$tau_InPu_JiMaArSuTa),method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(tau_to_time(pni$tau_root)),pch = 16)
c<-ci(tau_to_time(pni$tau_root),method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(tau_to_time(wpo$tau_InPu_Rondonia)),pch = 16)
c<-ci(tau_to_time(wpo$tau_InPu_Rondonia),method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)


#################################
#################################
#plot root times for each species with uncertainty for mut only
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "Uncertainty in µ", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,750000,1500000,2250000), labels=c(0,0.75,1.5,2.25))

#Galbula
points(1,mean(gcy_root_mut),pch = 16)
c<-ci(gcy_root_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_mad_mut),pch = 16)
c<-ci(mru_mad_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_root_mut),pch = 16)
c<-ci(hst_root_mut,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_mad_mut),pch = 16)
c<-ci(tae_mad_mut,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_root_mut),pch = 16)
c<-ci(pni_root_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo_mad_mut),pch = 16)
c<-ci(wpo_mad_mut,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species with uncertainty
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "Uncertainty in µ and G", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1500000,3000000,4500000), labels=c(0,1.5,3,4.5))

#Galbula
points(1,mean(gcy_root),pch = 16)
c<-ci(gcy_root,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_mad),pch = 16)
c<-ci(mru_mad,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_root),pch = 16)
c<-ci(hst_root,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_mad),pch = 16)
c<-ci(tae_mad,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_root),pch = 16)
c<-ci(pni_root,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo_mad),pch = 16)
c<-ci(wpo_mad,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#############################
#############################
#############################
#####PURUS BARRIER#########
#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,0.0045), type = 'n', main = "τ(Purus)", xlab = "species", ylab = "τ", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1,10,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(gcy$tau_Inam_Puru),pch = 16)
c<-ci(gcy$tau_Inam_Puru,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru$tau_Inam_Puru_JiMa),pch = 16)
c<-ci(mru$tau_Inam_Puru_JiMa,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
#points(3,mean(hst$tau_root),pch = 16)
#c<-ci(hst$tau_root,method = "HDI",ci = 0.95)
#arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
#points(4,mean(tae$tau_InPu_JiMaArSuTa),pch = 16)
#c<-ci(tae$tau_InPu_JiMaArSuTa,method = "HDI",ci = 0.95)
#arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni$tau_InPu),pch = 16)
c<-ci(pni$tau_InPu,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(wpo$tau_InPu_Rondonia),pch = 16)
#c<-ci(wpo$tau_InPu_Rondonia,method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "µ = 4.6E-9, G = 2", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,500000,1000000,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(tau_to_time(gcy$tau_Inam_Puru)),pch = 16)
c<-ci(tau_to_time(gcy$tau_Inam_Puru),method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(tau_to_time(mru$tau_Inam_Puru_JiMa)),pch = 16)
c<-ci(tau_to_time(mru$tau_Inam_Puru_JiMa),method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
#points(3,mean(tau_to_time(hst$tau_root)),pch = 16)
#c<-ci(tau_to_time(hst$tau_root),method = "HDI",ci = 0.95)
#arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
#points(4,mean(tau_to_time(tae$tau_InPu_JiMaArSuTa)),pch = 16)
#c<-ci(tau_to_time(tae$tau_InPu_JiMaArSuTa),method = "HDI",ci = 0.95)
#arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(tau_to_time(pni$tau_InPu)),pch = 16)
c<-ci(tau_to_time(pni$tau_InPu),method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(tau_to_time(wpo$tau_InPu_Rondonia)),pch = 16)
#c<-ci(tau_to_time(wpo$tau_InPu_Rondonia),method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)


#################################
#################################
#plot root times for each species with uncertainty for mut only
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "Uncertainty in µ", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,750000,1500000,2250000), labels=c(0,0.75,1.5,2.25))

#Galbula
points(1,mean(gcy_pur_mut),pch = 16)
c<-ci(gcy_pur_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_pur_mut),pch = 16)
c<-ci(mru_pur_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
#points(3,mean(hst_root_mut),pch = 16)
#c<-ci(hst_root_mut,method = "HDI",ci = 0.95)
#arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
#points(4,mean(tae_mad_mut),pch = 16)
#c<-ci(tae_mad_mut,method = "HDI",ci = 0.95)
#arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_pur_mut),pch = 16)
c<-ci(pni_pur_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(wpo_mad_mut),pch = 16)
#c<-ci(wpo_mad_mut,method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species with uncertainty
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "Uncertainty in µ and G", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1500000,3000000,4500000), labels=c(0,1.5,3,4.5))

#Galbula
points(1,mean(gcy_pur),pch = 16)
c<-ci(gcy_pur,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_pur),pch = 16)
c<-ci(mru_pur,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
#points(3,mean(hst_root),pch = 16)
#c<-ci(hst_root,method = "HDI",ci = 0.95)
#arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
#points(4,mean(tae_mad),pch = 16)
#c<-ci(tae_mad,method = "HDI",ci = 0.95)
#arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_pur),pch = 16)
c<-ci(pni_pur,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(wpo_mad),pch = 16)
#c<-ci(wpo_mad,method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)


#############################
#############################
#############################
#####ARIPUANA-ROOSEVELT BARRIER#########
#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,0.0045), type = 'n', main = "τ(Aripuanã)", xlab = "species", ylab = "τ", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1,10,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(gcy$tau_Para_MaRo_ArSuTa),pch = 16)
c<-ci(gcy$tau_Para_MaRo_ArSuTa,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru$tau_root),pch = 16)
c<-ci(mru$tau_root,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst$tau_Para_Mach_RoArSuTa),pch = 16)
c<-ci(hst$tau_Para_Mach_RoArSuTa,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae$tau_JiMaRoArSuTa),pch = 16)
c<-ci(tae$tau_JiMaRoArSuTa,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni$tau_Rondonia),pch = 16)
c<-ci(pni$tau_Rondonia,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(wpo$tau_InPu_Rondonia),pch = 16)
#c<-ci(wpo$tau_InPu_Rondonia,method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species without uncertainty
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "µ = 4.6E-9, G = 2", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,500000,1000000,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(tau_to_time(gcy$tau_Para_MaRo_ArSuTa)),pch = 16)
c<-ci(tau_to_time(gcy$tau_Para_MaRo_ArSuTa),method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(tau_to_time(mru$tau_root)),pch = 16)
c<-ci(tau_to_time(mru$tau_root),method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(tau_to_time(hst$tau_Para_Mach_RoArSuTa)),pch = 16)
c<-ci(tau_to_time(hst$tau_Para_Mach_RoArSuTa),method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tau_to_time(tae$tau_JiMaRoArSuTa)),pch = 16)
c<-ci(tau_to_time(tae$tau_JiMaRoArSuTa),method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(tau_to_time(pni$tau_Rondonia)),pch = 16)
c<-ci(tau_to_time(pni$tau_Rondonia),method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(tau_to_time(wpo$tau_InPu_Rondonia)),pch = 16)
#c<-ci(tau_to_time(wpo$tau_InPu_Rondonia),method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)


#################################
#################################
#plot root times for each species with uncertainty for mut only
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "Uncertainty in µ", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,750000,1500000,2250000), labels=c(0,0.75,1.5,2.25))

#Galbula
points(1,mean(gcy_arro_mut),pch = 16)
c<-ci(gcy_arro_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_root_mut),pch = 16)
c<-ci(mru_root_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_arro_mut),pch = 16)
c<-ci(hst_arro_mut,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_arro_mut),pch = 16)
c<-ci(tae_arro_mut,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_arro_mut),pch = 16)
c<-ci(pni_arro_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(wpo_mad_mut),pch = 16)
#c<-ci(wpo_mad_mut,method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#################################
#################################
#plot root times for each species with uncertainty
plot(c(1,6),c(0,0),ylim = c(0,3E6), type = 'n', main = "Uncertainty in µ and G", xlab = "species", ylab = "age(Ma)", xaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1500000,3000000,4500000), labels=c(0,1.5,3,4.5))

#Galbula
points(1,mean(gcy_arro),pch = 16)
c<-ci(gcy_arro,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_root),pch = 16)
c<-ci(mru_root,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_arro),pch = 16)
c<-ci(hst_arro,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_arro),pch = 16)
c<-ci(tae_arro,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_arro),pch = 16)
c<-ci(pni_arro,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(wpo_mad),pch = 16)
#c<-ci(wpo_mad,method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#####################
######6 plots########
#####################

#################################
#################################
par(mfrow=c(2,3))
#hist(m,col='black', main = "sampled mutation rates", xlab='µ')
#abline(v = mean(m), col='red',lwd=2, lty=2)
#arrows(mean(m)-1E-9,1E5,mean(m)+1E-9,1E5,code = 3, lty = 1, col = 'red', angle=90, length=0.05, lwd=2)
###################################################
########### Crown Div Times With mtDNA ############
#### mtDNA times from Silva et al 2019 SciAdv #####
###################################################

#par(mfrow=c(1,1))
#plot root times for each species with uncertainty for mut only
plot(c(0.5,6.5),c(0,0),ylim = c(0,3.5E6), type = 'n', main = "Crown Divergence Times", xlab = "species", ylab = "age(Ma)", xaxt='n',yaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
axis(2, at=c(0,1000000,2000000,3000000), labels=c(0,1.0,2.0,3.0))

#Galbula
points(0.90,mean(gcy_root_mut),pch = 16)
c<-ci(gcy_root_mut,method = "HDI",ci = 0.95)
arrows(0.90, c[[2]],0.90, c[[3]],code = 3, angle=90, length=0.025, lwd=1)
#mtdna
points(1.10,690000,pch = 16,col='red')
arrows(1.10, 360000,1.10, 970000,code = 3, angle=90, length=0.025, lwd=1,col='red')

#Malacoptila
points(1.90,mean(mru_root_mut),pch = 16)
c<-ci(mru_root_mut,method = "HDI",ci = 0.95)
arrows(1.90, c[[2]],1.90, c[[3]],code = 3, angle=90, length=0.025, lwd=1)
#mtdna
points(2.10,2480000,pch = 16, col="red")
arrows(2.10, 1720000,2.10, 3160000,code = 3, angle=90, length=0.025, lwd=1,col='red')

#Hypocnemis
points(2.90,mean(hst_root_mut),pch = 16)
c<-ci(hst_root_mut,method = "HDI",ci = 0.95)
arrows(2.90, c[[2]],2.90, c[[3]],code = 3, angle=90, length=0.025, lwd=1)
#mtdna
points(3.10,1580000,pch = 16,col='red')
arrows(3.10, 780000,3.10, 2670000,code = 3, angle=90, length=0.025, lwd=1,col='red')

#Thamnophilus
points(3.90,mean(tae_root_mut),pch = 16)
c<-ci(tae_root_mut,method = "HDI",ci = 0.95)
arrows(3.90, c[[2]],3.90, c[[3]],code = 3, angle=90, length=0.025, lwd=1)
#mtdna
points(4.10,1530000,pch = 16, col='red')
arrows(4.10, 940000,4.10, 2050000,code = 3, angle=90, length=0.025, lwd=1,col='red')

#Phlegopsis
points(4.90,mean(pni_root_mut),pch = 16)
c<-ci(pni_root_mut,method = "HDI",ci = 0.95)
arrows(4.90, c[[2]],4.90, c[[3]],code = 3, angle=90, length=0.025, lwd=1)
#mtdna
points(5.10,760000,pch = 16, col='red')
arrows(5.10, 540000, 5.10, 1000000,code = 3, angle=90, length=0.025, lwd=1,col='red')

#Willisornis
points(5.90,mean(wpo_root_mut),pch = 16)
c<-ci(wpo_root_mut,method = "HDI",ci = 0.95)
arrows(5.90, c[[2]],5.90, c[[3]],code = 3, angle=90, length=0.025, lwd=1)
#mtdna
points(6.10,1300000,pch = 16, col='red')
arrows(6.10, 550000,6.10, 2310000,code = 3, angle=90, length=0.025, lwd=1,col='red')

#plot root times for each species with uncertainty for mut only
plot(c(1,6),c(0,0),ylim = c(0,1.5E6), type = 'n', main = "Tapajos", xlab = "species", ylab = "age(Ma)", xaxt='n',yaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
axis(2, at=c(0,500000,1000000,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
points(1,mean(gcy_tap_mut),pch = 16)
c<-ci(gcy_tap_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_tap_mut),pch = 16)
c<-ci(mru_tap_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_tap_mut),pch = 16)
c<-ci(hst_tap_mut,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_root_mut),pch = 16)
c<-ci(tae_root_mut,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_tap_mut),pch = 16)
c<-ci(pni_tap_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo_tap_mut),pch = 16)
c<-ci(wpo_tap_mut,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#plot root times for each species with uncertainty for mut only
plot(c(1,6),c(0,0),ylim = c(0,2E6), type = 'n', main = "Madeira", xlab = "species", ylab = "age(Ma)", xaxt='n', yaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
axis(2, at=c(0,500000,1000000,1500000,2000000), labels=c(0,0.5,1.0,1.5,2.0))

#Galbula
points(1,mean(gcy_root_mut),pch = 16)
c<-ci(gcy_root_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_mad_mut),pch = 16)
c<-ci(mru_mad_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_root_mut),pch = 16)
c<-ci(hst_root_mut,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_mad_mut),pch = 16)
c<-ci(tae_mad_mut,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_root_mut),pch = 16)
c<-ci(pni_root_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
points(6,mean(wpo_mad_mut),pch = 16)
c<-ci(wpo_mad_mut,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#plot root times for each species with uncertainty for mut only
plot(c(1,6),c(0,0),ylim = c(0,1.1E6), type = 'n', main = "Aripuanã-Roosevelt", xlab = "species", ylab = "age(Ma)", xaxt='n',yaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
axis(2, at=c(0,250000,500000,750000,1000000), labels=c(0,0.25,0.5,0.75,1.0))

#Galbula
points(1,mean(gcy_arro_mut),pch = 16)
c<-ci(gcy_arro_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_root_mut),pch = 16)
c<-ci(mru_root_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
points(3,mean(hst_arro_mut),pch = 16)
c<-ci(hst_arro_mut,method = "HDI",ci = 0.95)
arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
points(4,mean(tae_arro_mut),pch = 16)
c<-ci(tae_arro_mut,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_arro_mut),pch = 16)
c<-ci(pni_arro_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(wpo_sut_mut),pch = 16)
#c<-ci(wpo_sut_mut,method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

plot(c(1,6),c(0,0),ylim = c(0,1E6), type = 'n', main = "Purus", xlab = "species", ylab = "age(Ma)", xaxt='n', yaxt='n')
axis(1, at=1:6, labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
axis(2, at=c(0,250000,500000,750000,1000000), labels=c(0,0.25,0.5,0.75,1.0))

#Galbula
points(1,mean(gcy_pur_mut),pch = 16)
c<-ci(gcy_pur_mut,method = "HDI",ci = 0.95)
arrows(1, c[[2]],1, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
points(2,mean(mru_pur_mut),pch = 16)
c<-ci(mru_pur_mut,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
#points(3,mean(hst_root_mut),pch = 16)
#c<-ci(hst_root_mut,method = "HDI",ci = 0.95)
#arrows(3, c[[2]],3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
#points(4,mean(tae_mad_mut),pch = 16)
#c<-ci(tae_mad_mut,method = "HDI",ci = 0.95)
#arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
points(5,mean(pni_pur_mut),pch = 16)
c<-ci(pni_pur_mut,method = "HDI",ci = 0.95)
arrows(5, c[[2]],5, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
#points(6,mean(wpo_mad_mut),pch = 16)
#c<-ci(wpo_mad_mut,method = "HDI",ci = 0.95)
#arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#plot Ne 
plot(c(1,13),c(0,0),ylim = c(0,1.6E6), type = 'n', main = "θ by population", xlab = "species", ylab = "τ", xaxt='n')
axis(1, at=c(2,4,6,8,10,12), labels=c("Gal.","Mal.","Hyp.","Tha.","Phl.","Wil."))
#axis(2, at=c(0,1,10,1500000), labels=c(0,0.5,1.0,1.5))

#Galbula
gcyne1<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_Inam))
points(1.4,mean(gcyne1),pch = 16)
c<-ci(gcyne1,method = "HDI",ci = 0.95)
arrows(1.4, c[[2]],1.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

gcyne2<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_Puru))
points(1.7,mean(gcyne2),pch = 16)
c<-ci(gcyne2,method = "HDI",ci = 0.95)
arrows(1.7, c[[2]],1.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

gcyne3<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_MaRo))
points(2,mean(gcyne3),pch = 16)
c<-ci(gcyne3,method = "HDI",ci = 0.95)
arrows(2, c[[2]],2, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

gcyne4<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_ArSuTa))
points(2.3,mean(gcyne4),pch = 16)
c<-ci(gcyne4,method = "HDI",ci = 0.95)
arrows(2.3, c[[2]],2.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

gcyne5<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_Para))
points(2.7,mean(gcyne5),pch = 16)
c<-ci(gcyne5,method = "HDI",ci = 0.95)
arrows(2.7, c[[2]],2.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Malacoptila
mrune1<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_Inam))
points(3.4,mean(mrune1),pch = 16)
c<-ci(mrune1,method = "HDI",ci = 0.95)
arrows(3.4, c[[2]],3.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

mrune2<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_Puru))
points(3.7,mean(mrune2),pch = 16)
c<-ci(mrune2,method = "HDI",ci = 0.95)
arrows(3.7, c[[2]],3.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

mrune3<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_JiMa))
points(4,mean(mrune3),pch = 16)
c<-ci(mrune3,method = "HDI",ci = 0.95)
arrows(4, c[[2]],4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

mrune4<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_MaArTa))
points(4.3,mean(mrune4),pch = 16)
c<-ci(mrune4,method = "HDI",ci = 0.95)
arrows(4.3, c[[2]],4.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

mrune5<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_Para))
points(4.6,mean(mrune5),pch = 16)
c<-ci(mrune5,method = "HDI",ci = 0.95)
arrows(4.6, c[[2]],4.6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Hypocnemis
hstne1<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_Puru))
points(5.4,mean(hstne1),pch = 16)
c<-ci(hstne1,method = "HDI",ci = 0.95)
arrows(5.4, c[[2]],5.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

hstne2<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_JiGu))
points(5.7,mean(hstne2),pch = 16)
c<-ci(hstne2,method = "HDI",ci = 0.95)
arrows(5.7, c[[2]],5.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

hstne3<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_Mach))
points(6,mean(hstne3),pch = 16)
c<-ci(hstne3,method = "HDI",ci = 0.95)
arrows(6, c[[2]],6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

hstne4<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_RoArSuTa))
points(6.3,mean(hstne4),pch = 16)
c<-ci(hstne4,method = "HDI",ci = 0.95)
arrows(6.3, c[[2]],6.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

hstne5<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_Para))
points(6.6,mean(hstne5),pch = 16)
c<-ci(hstne5,method = "HDI",ci = 0.95)
arrows(6.6, c[[2]],6.6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Thamnophilus
taene1<-replicate(n = 100000,expr = theta_to_Ne(tae$theta_InPu))
points(7.4,mean(taene1),pch = 16)
c<-ci(taene1,method = "HDI",ci = 0.95)
arrows(7.4, c[[2]],7.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

taene2<-replicate(n = 100000,expr = theta_to_Ne(tae$theta_JiMaRo))
points(7.7,mean(taene2),pch = 16)
c<-ci(taene2,method = "HDI",ci = 0.95)
arrows(7.7, c[[2]],7.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

taene3<-replicate(n = 100000,expr = theta_to_Ne(tae$theta_ArSuTa))
points(8,mean(taene3),pch = 16)
c<-ci(taene3,method = "HDI",ci = 0.95)
arrows(8, c[[2]],8, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

taene4<-replicate(n = 100000,expr = theta_to_Ne(tae$theta_Para))
points(8.3,mean(taene4),pch = 16)
c<-ci(taene4,method = "HDI",ci = 0.95)
arrows(8.3, c[[2]],8.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Phlegopsis
pnine1<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_Inam))
points(9.4,mean(pnine1),pch = 16)
c<-ci(pnine1,method = "HDI",ci = 0.95)
arrows(9.4, c[[2]],9.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

pnine2<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_Puru))
points(9.7,mean(pnine2),pch = 16)
c<-ci(pnine2,method = "HDI",ci = 0.95)
arrows(9.7, c[[2]],9.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

pnine3<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_JiMaRo))
points(10,mean(pnine3),pch = 16)
c<-ci(pnine3,method = "HDI",ci = 0.95)
arrows(10, c[[2]],10, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

pnine4<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_ArSuTa))
points(10.3,mean(pnine4),pch = 16)
c<-ci(pnine4,method = "HDI",ci = 0.95)
arrows(10.3, c[[2]],10.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

pnine5<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_Para))
points(10.6,mean(pnine5),pch = 16)
c<-ci(pnine5,method = "HDI",ci = 0.95)
arrows(10.6, c[[2]],10.6, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

#Willisornis
wpone1<-replicate(n = 100000,expr = theta_to_Ne(wpo$theta_InPu))
points(11.4,mean(wpone1),pch = 16)
c<-ci(wpone1,method = "HDI",ci = 0.95)
arrows(11.4, c[[2]],11.4, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

wpone2<-replicate(n = 100000,expr = theta_to_Ne(wpo$theta_Rondonia))
points(11.7,mean(wpone2),pch = 16)
c<-ci(wpone2,method = "HDI",ci = 0.95)
arrows(11.7, c[[2]],11.7, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

wpone3<-replicate(n = 100000,expr = theta_to_Ne(wpo$theta_SuTa))
points(12,mean(wpone3),pch = 16)
c<-ci(wpone3,method = "HDI",ci = 0.95)
arrows(12, c[[2]],12, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

wpone4<-replicate(n = 100000,expr = theta_to_Ne(wpo$theta_Para))
points(12.3,mean(wpone4),pch = 16)
c<-ci(wpone4,method = "HDI",ci = 0.95)
arrows(12.3, c[[2]],12.3, c[[3]],code = 3, angle=90, length=0.025, lwd=1)

###################################################
##############       MIGRATION      ###############
###################################################

#Galbula:
mean(gcy$m_Inam..Puru) #[1] 27.64018
ci(gcy$m_Inam..Puru,ci=0.95)#[2.75, 42.72]
mean(gcy$m_Puru..Inam)#[1] 30.65931
ci(gcy$m_Puru..Inam,ci=0.95)#[18.89, 46.63]

mean(gcy$m_Puru..MaRo) #[1] 5.621368
ci(gcy$m_Puru..MaRo,ci=0.95) #[2.98, 8.83]
mean(gcy$m_MaRo..Puru) #[1] 28.78329
ci(gcy$m_MaRo..Puru,ci=0.95)#[20.93, 38.20]

#Malacoptila
mean(mru$m_Puru_JiMa..MaArTa)#[1] 60.98444
ci(mru$m_Puru_JiMa..MaArTa,ci=0.95)#[1.13, 350.88]
mean(mru$m_MaArTa..Puru_JiMa)#[1] 57.35782
ci(mru$m_MaArTa..Puru_JiMa,ci=0.95)#[11.18, 194.78]

mean(mru$m_JiMa..MaArTa)#[1] 1.476423
ci(mru$m_JiMa..MaArTa,ci=0.95)#[0.00, 5.88]
mean(mru$m_MaArTa..JiMa)#[1] 8.70897
ci(mru$m_MaArTa..JiMa,ci=0.95)#[5.10, 12.61]

#Hypocnemis
mean(hst$m_Puru..JiGu)#[1] 0.4422647
ci(hst$m_Puru..JiGu,ci=0.95)#[0.00, 3.32]
mean(hst$m_JiGu..Puru)#[1] 13.88097
ci(hst$m_JiGu..Puru,ci=0.95)#[5.90, 22.21]

mean(hst$m_Puru..Mach)#[1] 0.7751072
ci(hst$m_Puru..Mach,ci=0.95)#[0.01, 2.55]
mean(hst$m_JiGu..Mach)#[1] 5.595002
ci(hst$m_JiGu..Mach,ci=0.95)#[2.95, 8.70]

mean(hst$m_JiGu..RoArSuTa) #[1] 3.066259
ci(hst$m_JiGu..RoArSuTa,ci=0.95) #[1.23, 5.37]
mean(hst$m_RoArSuTa..JiGu) #[1] 0.8458656
ci(hst$m_RoArSuTa..JiGu,ci=0.95) #[0.00, 6.35]

mean(hst$m_JiGu..Para) #[1] 1.376127
ci(hst$m_JiGu..Para,ci=0.95) #[0.04, 3.76]
mean(hst$m_Para..JiGu) #[1] 0.7703442
ci(hst$m_Para..JiGu,ci=0.95) #[0.00, 3.76]

#Thamnophilus
mean(tae$m_InPu..JiMaRo)#[1] 2.814413
ci(tae$m_InPu..JiMaRo,ci=0.95)#[0.00, 17.35]
mean(tae$m_JiMaRo..InPu)#[1] 18.0037
ci(tae$m_JiMaRo..InPu,ci=0.95)#[0.40, 32.75]

mean(tae$m_InPu..ArSuTa)#[1] 0.8233518
ci(tae$m_InPu..ArSuTa,ci=0.95)#[0.00, 6.94]
mean(tae$m_ArSuTa..InPu)#[1] 4.826842
ci(tae$m_ArSuTa..InPu,ci=0.95)#[0.00, 19.49]

mean(tae$m_JiMaRoArSuTa..Para)#[1] 42.91648
ci(tae$m_JiMaRoArSuTa..Para,ci=0.95)#[23.18, 73.85]
mean(tae$m_Para..JiMaRoArSuTa)#[1] 24.63901
ci(tae$m_Para..JiMaRoArSuTa,ci=0.95)#[13.72, 43.73]

#Phlegopsis
mean(pni$m_Puru..Rondonia)#[1] 133.1478
ci(pni$m_Puru..Rondonia,ci=0.95)#[52.36, 313.58]
mean(pni$m_Rondonia..Puru)#[1] 771.321
ci(pni$m_Rondonia..Puru,ci=0.95)#[332.53, 1903.02]

mean(pni$m_Puru..Para) #[1] 13.69924
ci(pni$m_Puru..Para,ci=0.95)#[7.11, 21.81]
mean(pni$m_Para..Puru)#[1] 12.49582
ci(pni$m_Para..Puru,ci=0.95)#[4.03, 21.00]

mean(pni$m_JiMaRo..ArSuTa)#[1] 0.4816192
ci(pni$m_JiMaRo..ArSuTa,ci=0.95)#[0.00, 5.16]
mean(pni$m_ArSuTa..JiMaRo)#[1] 2.198343
ci(pni$m_ArSuTa..JiMaRo,ci=0.95)#[0.00, 22.09]

mean(pni$m_ArSuTa..Para)#[1] 26.51823
ci(pni$m_ArSuTa..Para,ci=0.95)#[10.04, 47.62]
mean(pni$m_Para..ArSuTa)#[1] 31.83947
ci(pni$m_Para..ArSuTa,ci=0.95)#[18.04, 47.36]

#Willisornis
mean(wpo$m_InPu..Para)#[1] 1.615682
ci(wpo$m_InPu..Para,ci=0.95)#[0.91, 2.42]
mean(wpo$m_Para..InPu)#[1] 1.951523
ci(wpo$m_Para..InPu,ci=0.95)#[1.28, 2.70]

mean(wpo$m_Rondonia..SuTa)#[1] 4.061343
ci(wpo$m_Rondonia..SuTa,ci=0.95)#[3.03, 5.20]
mean(wpo$m_SuTa..Rondonia)#[1] 1.739683
ci(wpo$m_SuTa..Rondonia,ci=0.95)#[1.09, 2.51]

###########ancestral Ne#########
#Galbula
gcyne1<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_Inam_Puru))
mean(gcyne1) #350176.1
ci(gcyne1,method = "HDI",ci = 0.95)#[176071.43, 566763.95]

gcyne2<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_Para_MaRo_ArSuTa))
mean(gcyne2) #233044
ci(gcyne2,method = "HDI",ci = 0.95)#[135357-358724]

gcyne3<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_root))
mean(gcyne3) #170759.5
ci(gcyne3,method = "HDI",ci = 0.95)#[102857.14, 254217.73]

gcyne4<-replicate(n = 100000,expr = theta_to_Ne(gcy$theta_Para_MaRo))
mean(gcyne4) #18415.25
ci(gcyne4,method = "HDI",ci = 0.95)#[4999.70, 36548.03]


mrune1<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_root))
mean(mrune1) #105169.5
ci(mrune1,method = "HDI",ci = 0.95)#[62852.11, 159211.74]

mrune2<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_Inam_Puru_JiMa))
mean(mrune2) #77029.58
ci(mrune2,method = "HDI",ci = 0.95)#[40964.13, 123017.90]

mrune2<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_Inam_Puru_JiMa))
mean(mrune2) #77029.58
ci(mrune2,method = "HDI",ci = 0.95)#[40964.13, 123017.90]

mrune3<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_Puru_JiMa))
mean(mrune3) #95059.58
ci(mrune3,method = "HDI",ci = 0.95)#[46403.35, 155727.57]

mrune4<-replicate(n = 100000,expr = theta_to_Ne(mru$theta_Para_MaArTa))
mean(mrune4) #112902.7
ci(mrune4,method = "HDI",ci = 0.95)#[65357.14, 173367.35]

hstne1<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_root))
mean(hstne1) #150497.6
ci(hstne1,method = "HDI",ci = 0.95)#[88571.43, 229219.14]

hstne2<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_Para_Rondonia))
mean(hstne2) #352442.1
ci(hstne2,method = "HDI",ci = 0.95)#[199583.75, 546154.71]

hstne3<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_Para_Mach_RoArSuTa))
mean(hstne3) #472954.8
ci(hstne3,method = "HDI",ci = 0.95)#[273489.54, 727502.51]

hstne4<-replicate(n = 100000,expr = theta_to_Ne(hst$theta_Para_Mach))
mean(hstne4) #197998.6
ci(hstne4,method = "HDI",ci = 0.95)#[19707.17, 430387.62]

taene1<-replicate(n = 100000,expr = theta_to_Ne(tae$theta_root))
mean(taene1) #144746.9
ci(taene1,method = "HDI",ci = 0.95)#[85349.88, 220559.61]

taene2<-replicate(n = 100000,expr = theta_to_Ne(tae$theta_InPu_JiMaArSuTa))
mean(taene2) #501731
ci(taene2,method = "HDI",ci = 0.95)#[262073.16, 794785.29]

taene3<-replicate(n = 100000,expr = theta_to_Ne(tae$theta_JiMaRoArSuTa))
mean(taene3) #712391.2
ci(taene3,method = "HDI",ci = 0.95)#[352572.49, 1151244.08]

pnine1<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_root))
mean(pnine1) #213951.8
ci(pnine1,method = "HDI",ci = 0.95)#[129273.60, 319280.31]

pnine2<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_Para_Rondonia))
mean(pnine2) #180850.6
ci(pnine2,method = "HDI",ci = 0.95)#[106040.29, 275344.16]

pnine3<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_Rondonia))
mean(pnine3) #137006.7
ci(pnine3,method = "HDI",ci = 0.95)#[79539.65, 208237.87]

pnine4<-replicate(n = 100000,expr = theta_to_Ne(pni$theta_InPu))
mean(pnine4) #421874.5
ci(pnine4,method = "HDI",ci = 0.95)#[249987.66, 635332.07]

wpone1<-replicate(n = 100000,expr = theta_to_Ne(wpo$theta_root))
mean(wpone1) #152158.3
ci(wpone1,method = "HDI",ci = 0.95)#[90357.14, 228755.94]

wpone2<-replicate(n = 100000,expr = theta_to_Ne(wpo$theta_InPu_Rondonia))
mean(wpone2) #278629.2
ci(wpone2,method = "HDI",ci = 0.95)#[165714.29, 419840.39]

wpone3<-replicate(n = 100000,expr = theta_to_Ne(wpo$theta_Para_SuTa))
mean(wpone3) #200760.1
ci(wpone3,method = "HDI",ci = 0.95)#[118927.54, 302428.84]


