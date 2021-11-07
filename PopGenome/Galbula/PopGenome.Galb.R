rm(list=ls())

library(PopGenome)

setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Galbula/')
# Load the data
snp <- readData("./vcf/", format="VCF")

# This is complex object, with several slots
get.sum.data(snp)
show.slots(snp)

# You can access the different "slots" by using the "@" sign:
snp@n.sites
get.individuals(snp)

# Set populations
Inambari <- c("G_cy_T12385_In",	"G_cy_T23196_In", "G_cy_T310_In", "G_cy_T7636_In")
Purus<-c("G_cy_82508_pu", "G_cy_T12275_pu", "G_cy_T12279_pu", "G_cy_T12392_pu", "G_cy_T13184_pu", "G_cy_T26228_pu", "G_cy_T26229_pu", "G_cy_T26252_pu")
Jiparana<-c("G_cy_T3343_jigu", "G_cy_T3384_jigu", "G_cy_T3385_jigu")
Machado<-c("G_cy_J296_ma", "G_cy_J477_ma", "G_cy_J773_ma", "G_cy_T13251_ma", "G_cy_T363_ma", "G_cy_T364_ma")
Roosevelt<-c("G_cy_J691_roar", "G_cy_J694_roar")
Aripuana<-c("G_cy_80582_arsu", "G_cy_80701_arsu", "G_cy_80801_arsu", "G_cy_80826_arsu", "G_cy_81108_arsu", "G_cy_81118_arsu", "G_cy_85499_arsu", "G_cy_85678_arsu")
Sucunduri<-c("G_cy_85356_suta", "G_cy_86297_suta", "G_cy_86321_suta", "G_cy_86458_suta", "G_cy_86478_suta", "G_cy_T14558_suta", "G_cy_T16693_suta", "G_cy_T18563_suta", "G_cy_T18620_suta")
Para<-c("G_cy_T10897_pa", "G_cy_T11062_pa", "G_cy_T16771_pa", "G_cy_T1705_pa", "G_cy_T18744_pa", "G_cy_T19429_pa", "G_cy_T19520_pa", "G_cy_T19765_pa", "G_cy_T2497_pa", "G_cy_T6579_pa", "G_cy_T9133_pa")

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
####sliding window stats###
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

