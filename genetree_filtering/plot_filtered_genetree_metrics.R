rm(list=ls())
library(gmodels)
setwd("~/AMNH/Roosevelt_Project/Data/genetrees/willisornis/")

tab<-read.csv("Willis_gt_metrics.csv",header = T)
attach(tab)
names(tab)

var.sites<-as.numeric(var.sites)
LR<-as.vector(LR)

plot(var.sites,LR)
abline(lm(LR~var.sites),col="red")
summary(lm(LR~var.sites))

plot(var.frac,LR)
abline(lm(LR~var.frac),col="red")
summary(lm(LR~var.frac))

plot(var.sites,logLikelihood)
abline(lm(logLikelihood~var.sites),col="red")
summary(lm(logLikelihood~var.sites))

plot(var.frac,logLikelihood)
abline(lm(logLikelihood~var.frac),col="red")
summary(lm(logLikelihood~var.frac))

plot(LR,logLikelihood)

mod<-glm(LR~var.sites*logLikelihood*var.frac)
summary(mod)
cooksd<-cooks.distance(mod)
#plot(cooksd, pch="+", cex=1, main="Influential Obs by Cooks distance")
#abline(h=4*mean(cooksd, na.rm=T), col="red")
#points(x=ifelse(cooksd>4*mean(cooksd, na.rm=T), names(cooksd),""), col="red")
#text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>3*mean(cooksd, na.rm=T), names(cooksd),""),col="red")

influential <- as.numeric(names(cooksd)[cooksd > 4*mean(cooksd, na.rm=T)])
tab[influential, ]
influential
x<-c(1:length(LR))
filtered_tab<-tab[! x %in% influential,]
write.csv(filtered_tab,file="filtered_tab.csv")

filtered_tab<-read.csv("filtered_tab.csv")

plot(filtered_tab$seq.length,filtered_tab$var.sites)
plot(filtered_tab$seq.length,filtered_tab$var.frac)
plot(filtered_tab$var.frac,filtered_tab$var.sites)

plot(filtered_tab$var.sites,filtered_tab$LR)
plot(filtered_tab$var.sites,filtered_tab$logLikelihood)
plot(filtered_tab$LR,filtered_tab$logLikelihood)

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

conf<-confidence_interval(filtered_tab$var.sites,0.975)
conf

influential_tab<-tab[influential,]
write.csv(influential_tab,file="influential_tab.csv")

plot(influential_tab$var.sites,influential_tab$LR)
plot(influential_tab$var.sites,influential_tab$logLikelihood)
plot(influential_tab$LR,influential_tab$logLikelihood)


setwd('PreFiltered/')

