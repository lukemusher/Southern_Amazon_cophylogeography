rm(list=ls())

library(PopGenome)

setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Phlegopsis/')
# Load the data
snp <- readData("./vcf/", format="VCF")

# This is complex object, with several slots
get.sum.data(snp)
show.slots(snp)

# You can access the different "slots" by using the "@" sign:
snp@n.sites
get.individuals(snp)

# Set populations
Inambari <- c("P_ni_A7862_In","P_ni_A7911_In","P_ni_A7928_In","P_ni_T6243_In")
Purus<-c("P_ni_T4404_pu","P_ni_T3609_pu","P_ni_T3611_pu","P_ni_T3817_pu","P_ni_T4043_pu","P_ni_T4051_pu","P_ni_T4313_pu","P_ni_T5850_pu","P_ni_T15938_pu","P_ni_T5940_pu","P_ni_T5974_pu","P_ni_80034_pu")
Jiparana<-c("P_ni_T15863_jigu","P_ni_T15868_jigu","P_ni_T15871_jigu","P_ni_T22153_jigu","P_ni_T3261_jigu","P_ni_A3255_jigu")
Machado<-c("P_ni_J434_ma","P_ni_J461_ma","P_ni_J462_ma","P_ni_J485_ma","P_ni_T369_ma","P_ni_T443_ma","P_ni_T467_ma","P_ni_A2418_ma","P_ni_A542_ma","P_ni_J210_ma","P_ni_J227_ma","P_ni_J260_ma")
Roosevelt<-c("P_ni_J363_roar","P_ni_J371_roar","P_ni_J373_roar","P_ni_J381_roar","P_ni_J385_roar","P_ni_J389_roar","P_ni_J417_roar","P_ni_J684_roar","P_ni_J724_roar","P_ni_J361_roar")
Aripuana<-c("P_ni_J551_arsu","P_ni_J602_arsu","P_ni_J603_arsu","P_ni_J614_arsu","P_ni_J617_arsu","P_ni_80555_arsu","P_ni_80684_arsu","P_ni_80802_arsu","P_ni_80874_arsu","P_ni_85430_arsu","P_ni_86072_arsu")
Sucunduri<-c("P_ni_T11888_suta","P_ni_T14543_suta","P_ni_T10204_suta","P_ni_T10967_suta","P_ni_T16698_suta","P_ni_T9076_suta","P_ni_A15120_suta","P_ni_85721_suta","P_ni_77876_suta","P_ni_78155_suta")
Para<-c("P_ni_T10673_pa","P_ni_T10940_pa","P_ni_T11193_pa","P_ni_T11222_pa","P_ni_T12345_pa","P_ni_T12854_pa","P_ni_T1642_pa","P_ni_T18703_pa","P_ni_A14342_pa","P_ni_A15277_pa","P_ni_A7066_pa")

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


PopGenome####sliding window stats###
snp.slide <- sliding.window.transform(snp,width=50, jump=50,type=1,whole.data=TRUE)
snp.slide <- neutrality.stats(snp.slide)
snp.slide <- F_ST.stats(snp.slide)
snp.slide <- linkage.stats(snp.slide)

snp.slide@Tajima.D

boxplot(snp.slide@Pi)

############THREE POPS##############
# Load the data
snp <- readData("./vcf/", format="VCF")

# This is complex object, with several slots
get.sum.data(snp)
show.slots(snp)

# You can access the different "slots" by using the "@" sign:
snp@n.sites

# Set populations
Inambari <- c("G_cy_T12385_In",	"G_cy_T23196_In", "G_cy_T310_In", "G_cy_T7636_In","G_cy_82508_pu", "G_cy_T12275_pu", "G_cy_T12279_pu", "G_cy_T12392_pu", "G_cy_T13184_pu", "G_cy_T26228_pu", "G_cy_T26229_pu", "G_cy_T26252_pu")
Rondonia_S<-c("G_cy_T3343_jigu", "G_cy_T3384_jigu", "G_cy_T3385_jigu","G_cy_J296_ma", "G_cy_J477_ma", "G_cy_J773_ma", "G_cy_T13251_ma", "G_cy_T363_ma", "G_cy_T364_ma","G_cy_J691_roar", "G_cy_J694_roar")
Rondonia_N<-c("G_cy_80582_arsu", "G_cy_80701_arsu", "G_cy_80801_arsu", "G_cy_80826_arsu", "G_cy_81108_arsu", "G_cy_81118_arsu", "G_cy_85499_arsu", "G_cy_85678_arsu","G_cy_85356_suta", "G_cy_86297_suta", "G_cy_86321_suta", "G_cy_86458_suta", "G_cy_86478_suta", "G_cy_T14558_suta", "G_cy_T16693_suta", "G_cy_T18563_suta", "G_cy_T18620_suta")
Para<-c("G_cy_T10897_pa", "G_cy_T11062_pa", "G_cy_T16771_pa", "G_cy_T1705_pa", "G_cy_T18744_pa", "G_cy_T19429_pa", "G_cy_T19520_pa", "G_cy_T19765_pa", "G_cy_T2497_pa", "G_cy_T6579_pa", "G_cy_T9133_pa")

snp  <- set.populations(snp, list(Inambari,Rondonia_S,Rondonia_N,Para))

snp@populations # check if it worked

# Diversities and FST (by scaffold)
snp <- F_ST.stats(snp) # this does the calculations and 
snp@Pi

snp <- linkage.stats(snp)
snp@Wall.B

snp <- neutrality.stats(snp)
snp@Tajima.D
snp@
  
  ####sliding window stats###
  snp.slide <- sliding.window.transform(snp,width=50, jump=50,type=1,whole.data=TRUE)
snp.slide <- neutrality.stats(snp.slide)
snp.slide <- F_ST.stats(snp.slide)
snp.slide <- F_ST.stats.2(snp.slide)
snp.slide <- linkage.stats(snp.slide)
snp.slide <- detail.stats(snp.slide)

snp.slide@Tajima.D

boxplot(snp.slide@Pi)
boxplot(snp.slide@theta_Watterson)

