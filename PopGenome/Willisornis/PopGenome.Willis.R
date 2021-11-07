rm(list=ls())

library(PopGenome)

setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Willisornis/')
# Load the data
snp <- readData("./vcf/", format="VCF")

# This is complex object, with several slots
get.sum.data(snp)
show.slots(snp)

# You can access the different "slots" by using the "@" sign:
snp@n.sites
get.individuals(snp)

# Set populations
Inam<- c("W_po_A815_In", "W_po_T23485_In","W_po_T3645_In","W_po_T4097_In")
Puru<- c("W_po_A424_pu","W_po_A7736_pu","W_po_A101_pu","W_po_T12400_pu",  "W_po_T13148_pu", "W_po_A9697_pu","W_po_A10701_pu","W_po_A1381_pu","W_po_A265_pu", "W_po_A2725_pu", "W_po_A2749_pu", "W_po_A2759_pu")
JiGu<-c("W_po_A308_jigu", "W_po_T7625_jigu", "W_po_T7627_jigu","W_po_T15888_jigu", "W_po_T15889_jigu","W_po_A3264_jigu", "W_po_A334_jigu", "W_po_A368_jigu")
Mach<-c("W_po_A465_ma","W_po_T2190_ma",  "W_po_T3157_ma",  "W_po_T388_ma",  "W_po_A882_ma", "W_po_J439_ma",  "W_po_A894_ma", "W_po_A896_ma", "W_po_A898_ma", "W_po_A902_ma","W_po_A472_ma","W_po_A24413_ma","W_po_A4235_ma","W_po_J163_ma", "W_po_J192_ma", "W_po_J202_ma", "W_po_J209_ma", "W_po_J222_ma", "W_po_J226_ma", "W_po_J302_ma")
Roar<-c("W_po_J335_roar", "W_po_J340_roar", "W_po_J391_roar", "W_po_J421_roar","W_po_J644_roar", "W_po_J766_roar", "W_po_J768_roar")
ArSu<-c("W_po_A891_arsu","W_po_80798_arsu", "W_po_81036_arsu", "W_po_85474_arsu", "W_po_86085_arsu","W_po_J561_arsu", "W_po_J563_arsu", "W_po_J610_arsu")
SuTa<-c("W_po_A16298_suta","W_po_T11887_suta","W_po_81265_suta", "W_po_T10153_suta", "W_po_T10234_suta")
Para<-c( "W_po_T12793_pa", "W_po_T7132_pa",  "W_po_T9346_pa","W_po_T10694_pa","W_po_25484_pa","W_po_A7346_pa","W_po_A12284_pa",  "W_po_A14569_pa", "W_po_A16227_pa")

snp  <- set.populations(snp, list(Inam,Puru,JiGu,Mach,Roar,ArSu,SuTa,Para))

snp@populations # check if it worked

# Diversities and FST (by scaffold)
snp <- F_ST.stats(snp) # this does the calculations and 
snp@Pi

snp <- linkage.stats(snp)
snp@Wall.B

snp <- neutrality.stats(snp)
snp@Tajima.D
snp@pi
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

