rm(list=ls())

library(PopGenome)

setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Thamnophilus/')
# Load the data
snp <- readData("./vcf/", format="VCF")

# This is complex object, with several slots
get.sum.data(snp)
show.slots(snp)

# You can access the different "slots" by using the "@" sign:
snp@n.sites
get.individuals(snp)

# Set populations
Inambari <- c("T_ae_T14205_In","T_ae_A7968_In","T_ae_A8087_In","T_ae_28095_In","T_ae_75520_In")
Purus<-c("T_ae_82491_pu","T_ae_82606_pu","T_ae_82608_pu","T_ae_A21562_pu","T_ae_A24377_pu","T_ae_A2720_pu","T_ae_A2833_pu","T_ae_A3559_pu","T_ae_A439_pu","T_ae_T12313_pu","T_ae_T13219_pu")
Jiparana<-c("T_ae_T3237_jigu","T_ae_T3260_jigu","T_ae_T3376_jigu","T_ae_A317_jigu","T_ae_A324_jigu")
Machado<-c("T_ae_A458_ma","T_ae_A509_ma","T_ae_A522_ma","T_ae_J249_ma","T_ae_J252_ma","T_ae_J261_ma","T_ae_J298_ma","T_ae_T13278_ma","T_ae_T2166_ma","T_ae_T2207_ma","T_ae_T4355_ma")
Roosevelt<-c("T_ae_J319_roar","T_ae_J419_roar","T_ae_J678_roar","T_ae_J683_roar")
Aripuana<-c("T_ae_J598_arsu","T_ae_J616_arsu","T_ae_86147_arsu","T_ae_86229_arsu","T_ae_86569_arsu","T_ae_80508_arsu","T_ae_80716_arsu","T_ae_81278_arsu","T_ae_81290_arsu","T_ae_81338_arsu")
Sucunduri<-c("T_ae_T10229_suta","T_ae_T14601_suta","T_ae_T16625_suta","T_ae_T16704_suta","T_ae_T19436_suta","T_ae_T24591_suta","T_ae_T742_suta","T_ae_T7941_suta","T_ae_85274_suta","T_ae_A15069_suta","T_ae_A16041_suta")
Para<-c("T_ae_A15279_pa","T_ae_T10679_pa","T_ae_T13575_pa","T_ae_T17307_pa","T_ae_T19420_pa","T_ae_T8268_pa")

snp  <- set.populations(snp, list(Inambari,Purus,Jiparana,Machado,Roosevelt,Aripuana,Sucunduri,Para))

snp@populations # check if it worked

# Diversities and FST (by scaffold)
snp <- F_ST.stats(snp) # this does the calculations and 
snp@Pi

snp <- linkage.stats(snp)
snp@Wall.B

snp <- neutrality.stats(snp)
snp@Tajima.D

snp@theta_Watterson

pops<-c("Inambari","Purus","Jiparana","Machado","Roosevelt","Aripuana","Sucunduri","Para")
pi<-snp@Pi[1,]
theta<-snp@theta_Watterson[1,]
B<-snp@Wall.B[1,]
Taj.D<-snp@Tajima.D[1,]

gen_stats<-as.data.frame(cbind(pops,pi,theta,B,Taj.D))
gen_stats

write.table(gen_stats, file="stats.txt", sep="\t")

snp@nuc.F_ST.pairwise

####sliding window stats###
snp.slide <- sliding.window.transform(snp,width=10, jump=5,type=1,whole.data=TRUE)
snp.slide <- neutrality.stats(snp.slide)
snp.slide <- F_ST.stats(snp.slide)
snp.slide <- linkage.stats(snp.slide)

snp.slide@nuc.F_ST.pairwise
write.csv(snp.slide@nuc.F_ST.pairwise,file="pairwise.fst.csv")

Purus<-as.numeric(snp.slide@nuc.F_ST.pairwise[1,])
Madeira.Upper1<-as.numeric(snp.slide@nuc.F_ST.pairwise[2,])
Madeira.Upper2<-as.numeric(snp.slide@nuc.F_ST.pairwise[8,])
Jiparana<-as.numeric(snp.slide@nuc.F_ST.pairwise[14,])
Roosevelt<-as.numeric(snp.slide@nuc.F_ST.pairwise[19,])
Aripuana.Lower<-as.numeric(snp.slide@nuc.F_ST.pairwise[20,])
Tapajos.Upper<-as.numeric(snp.slide@nuc.F_ST.pairwise[22,])
Aripuana.Upper<-as.numeric(snp.slide@nuc.F_ST.pairwise[23,])
Sucunduri<-as.numeric(snp.slide@nuc.F_ST.pairwise[26,])
Tapajos.Lower<-as.numeric(snp.slide@nuc.F_ST.pairwise[28,])

par(mfrow=c(1,1))

pw.fst.rivs<-as.data.frame(cbind(Purus,Madeira.Upper1,Madeira.Upper2,Jiparana,Roosevelt,Aripuana.Upper,Aripuana.Lower,Sucunduri,Tapajos.Upper,Tapajos.Lower))

riv.fst.tab<-data.frame()

counts=0

for (i in 1:length(names(pw.fst.rivs))){
  for (j in 1:length(pw.fst.rivs[,1])){
    counts=counts+1
    riv.fst.tab[counts,1]<-names(pw.fst.rivs)[i]
    riv.fst.tab[counts,2]<-pw.fst.rivs[j,i]
  }
}

riv.fst.tab

library(ggplot2)
p <- ggplot(riv.fst.tab, aes(x=V1, y=V2),reorder(names(pw.fst.rivs))) + 
  geom_boxplot() +
  xlab("") +
  ylab("FST") + 
  geom_jitter(shape=1, position=position_jitter(0.05))
p
