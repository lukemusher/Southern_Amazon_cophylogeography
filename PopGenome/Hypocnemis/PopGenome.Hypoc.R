rm(list=ls())

library(PopGenome)

setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Hypocnemis/')
# Load the data
snp <- readData("./vcf/", format="VCF")

# This is complex object, with several slots
get.sum.data(snp)
show.slots(snp)

# You can access the different "slots" by using the "@" sign:
snp@n.sites
get.individuals(snp)

# Set populations
Inambari <- c()
Purus<-c("H_pe_A24346_pu")
Jiparana<-c("H_oc_A311_jigu","H_oc_T15847_jigu","H_oc_T355_jigu")
Machado<-c("H_ro_A547_ma","H_ro_J213_ma","H_ro_J305_ma","H_ro_J508_ma","H_ro_J774_ma","H_ro_J775_ma","H_ro_J796_ma","H_ro_T1842_ma","H_ro_T366_ma","H_ro_T385_ma","H_ro_T471_ma")
Roosevelt<-c("H_st_J364_roar","H_st_J368_roar","H_st_J370_roar","H_st_J374_roar","H_st_J408_roar","H_st_J621_roar","H_st_J664_roar","H_st_J665_roar","H_st_J711_roar","H_st_J762_roar","H_st_J765_roar")
Aripuana<-c("H_ro_A551_ma","H_ro_A409_ma","H_ro_A410_ma","H_ro_A521_ma","H_st_A272_arsu","H_st_A273_arsu","H_st_J525_arsu","H_st_J530_arsu","H_st_J536_arsu","H_st_J572_arsu","H_st_77860_arsu","H_st_78249_arsu","H_st_80727_arsu","H_st_80800_arsu","H_st_81143_arsu","H_st_85680_arsu","H_st_85970_arsu")
Sucunduri<-c("H_st_A9955_pa","H_st_86405_suta","H_st_A14546_suta","H_st_A4899_suta","H_st_T10207_suta","H_st_T11900_suta","H_st_T12194_suta","H_st_T24564_suta","H_st_T7114_suta","H_st_81279_suta")
Para<-c("H_st_A11597_pa","H_st_A15208_pa","H_st_A15210_pa","H_st_A16571_pa","H_st_A7597_pa","H_st_T16744_pa","H_st_T17858_pa")

snp  <- set.populations(snp, list(Purus,Jiparana,Machado,Roosevelt,Aripuana,Sucunduri,Para))

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
Inambari <- c("G_sp_T23299_In","G_sp_T16241_In","G_sp_T4085_In","G_sp_A059_In","G_sp_A7903_In","G_sp_A870_In","G_sp_A875_In","G_sp_T12393_pu","G_sp_T13157_pu","G_sp_T13704_pu","G_sp_A9702_pu","G_sp_82591_pu","G_sp_82593_pu","G_sp_82595_pu","G_sp_A10327_pu","G_sp_A10333_pu","G_sp_A1495_pu","G_sp_A16885_pu","G_sp_A191_pu","G_sp_A2721_pu","G_sp_A2755_pu","G_sp_A3189_pu","G_sp_A879_pu","G_sp_A3020_jigu","G_sp_A3022_jigu","G_sp_A3027_jigu","G_sp_A329_jigu","G_sp_T3205_jigu","G_sp_T3219_jigu","G_sp_T15825_jigu","G_sp_T22166_jigu")
Ron_par<-c("G_sp_A2378_ma","G_sp_T2107_ma","G_sp_A906_ma","G_sp_A2401_ma","G_sp_A24267_ma","G_sp_A4227_ma","G_sp_A461_ma","G_sp_A553_ma","G_sp_J454_ma","G_sp_J458_ma","G_sp_J482_ma","G_sp_J484_ma","G_sp_J486_ma","G_sp_T420_ma","G_sp_T460_ma","G_sp_T799_ma","G_sp_J337_roar","G_sp_J343_roar","G_sp_A895_arsu","G_sp_J355_roar","G_sp_J388_roar","G_sp_J645_roar","G_sp_J675_roar","G_sp_J720_roar")
Aripuana<-c("G_sp_T16289_arsu","G_sp_A845_arsu","G_sp_J532_arsu","G_sp_J538_arsu","G_sp_J564_arsu","G_sp_80451_arsu","G_sp_80639_arsu","G_sp_81132_arsu","G_sp_81852_arsu","G_sp_81861_arsu","G_sp_86271_arsu")
Sucunduri<-c("G_sp_T14563_suta","G_sp_A8463_suta","G_sp_A890_suta","G_sp_T1532_suta","G_sp_T10159_suta","G_sp_T16644_suta","G_sp_T19297_suta","G_sp_T10238_suta","G_sp_T11811_suta","G_sp_A14324_suta","G_sp_A15178_suta","G_sp_A5172_suta")
Para<-c("G_sp_A14339_pa","G_sp_A14414_pa","G_sp_A7053_pa","G_sp_T19382_pa","G_sp_A9937_pa","G_sp_T11190_pa","G_sp_T11232_pa","G_sp_T24586_pa")

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

