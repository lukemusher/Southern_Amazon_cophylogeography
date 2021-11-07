rm(list=ls())

library(PopGenome)

setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Thamnophilus/')
# Load the data

snp <- readData(path = './vcf/',include.unknown = T,format = 'VCF',big.data = T)

# You can access the different "slots" by using the "@" sign:
snp@n.sites

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

####sliding window stats###
snp.slide <- sliding.window.transform(snp,width=5000, jump=2500,type=1,whole.data=TRUE)
snp.slide <- F_ST.stats(snp.slide, mode = "haplotype")
write.csv(snp.slide@Pi,file="pi.csv")

Inambari<-as.numeric(snp.slide@Pi[,1])
Purus<-as.numeric(snp.slide@Pi[,2])
Jiparana<-as.numeric(snp.slide@Pi[,3])
Machado<-as.numeric(snp.slide@Pi[,4])
Roosevelt<-as.numeric(snp.slide@Pi[,5])
Aripuana<-as.numeric(snp.slide@Pi[,6])
Sucunduri<-as.numeric(snp.slide@Pi[,7])
Para<-as.numeric(snp.slide@Pi[,8])

pi.rivs<-as.data.frame(cbind(Inambari,Purus,Jiparana,Machado,Aripuana,Roosevelt,Sucunduri,Para))

riv.pi.tab<-data.frame()

counts=0

for (i in 1:length(names(pi.rivs))){
  for (j in 1:length(pi.rivs[,1])){
    counts=counts+1
    riv.pi.tab[counts,1]<-names(pi.rivs)[i]
    riv.pi.tab[counts,2]<-pi.rivs[j,i]
    riv.pi.tab[counts,3]<-as.character("Thamnophilus")
  }
}

riv.pi.tab

setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Galbula/')
# Load the data
snp <- readData(path = './vcf/',include.unknown = T,format = 'VCF',big.data = T)

# This is complex object, with several slots
get.sum.data(snp)

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

####sliding window stats###
snp.slide <- sliding.window.transform(snp,width=5000, jump=2500,type=1,whole.data=TRUE)
snp.slide <- F_ST.stats(snp.slide, mode = "haplotype")

write.csv(snp.slide@Pi,file="pi.csv")

Inambari<-as.numeric(snp.slide@Pi[,1])
Purus<-as.numeric(snp.slide@Pi[,2])
Jiparana<-as.numeric(snp.slide@Pi[,3])
Machado<-as.numeric(snp.slide@Pi[,4])
Roosevelt<-as.numeric(snp.slide@Pi[,5])
Aripuana<-as.numeric(snp.slide@Pi[,6])
Sucunduri<-as.numeric(snp.slide@Pi[,7])
Para<-as.numeric(snp.slide@Pi[,8])

pi.rivs<-as.data.frame(cbind(Inambari,Purus,Jiparana,Machado,Aripuana,Roosevelt,Sucunduri,Para))

for (i in 1:length(names(pi.rivs))){
  for (j in 1:length(pi.rivs[,1])){
    counts=counts+1
    riv.pi.tab[counts,1]<-names(pi.rivs)[i]
    riv.pi.tab[counts,2]<-pi.rivs[j,i]
    riv.pi.tab[counts,3]<-as.character("Galbula")
  }
}


setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Phlegopsis/')
# Load the data
snp <- readData(path = './vcf/',include.unknown = T,format = 'VCF',big.data = T)

# This is complex object, with several slots
get.sum.data(snp)

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

####sliding window stats###
snp.slide <- sliding.window.transform(snp,width=5000, jump=2500,type=1,whole.data=TRUE)
snp.slide <- F_ST.stats(snp.slide, mode = "haplotype")

write.csv(snp.slide@Pi,file="pi.csv")

Inambari<-as.numeric(snp.slide@Pi[,1])
Purus<-as.numeric(snp.slide@Pi[,2])
Jiparana<-as.numeric(snp.slide@Pi[,3])
Machado<-as.numeric(snp.slide@Pi[,4])
Roosevelt<-as.numeric(snp.slide@Pi[,5])
Aripuana<-as.numeric(snp.slide@Pi[,6])
Sucunduri<-as.numeric(snp.slide@Pi[,7])
Para<-as.numeric(snp.slide@Pi[,8])

pi.rivs<-as.data.frame(cbind(Inambari,Purus,Jiparana,Machado,Aripuana,Roosevelt,Sucunduri,Para))

for (i in 1:length(names(pi.rivs))){
  for (j in 1:length(pi.rivs[,1])){
    counts=counts+1
    riv.pi.tab[counts,1]<-names(pi.rivs)[i]
    riv.pi.tab[counts,2]<-pi.rivs[j,i]
    riv.pi.tab[counts,3]<-as.character("Phlegopsis")
  }
}


setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Willisornis/')
# Load the data
snp <- readData(path = './vcf/',include.unknown = T,format = 'VCF',big.data = T)

# Set populations
Inambari <- c("W_po_A815_In","W_po_A815_In.2","W_po_T12362_In","W_po_T12362_In.2","W_po_T23485_In","W_po_T23485_In.2","W_po_T3645_In","W_po_T3645_In.2","W_po_A16828_In","W_po_A16828_In.2","W_po_T4097_In","W_po_T4097_In.2")
Purus<-c("W_po_A7736_pu","W_po_A7736_pu.2","W_po_A9697_pu","W_po_A9697_pu.2","W_po_A424_pu","W_po_A424_pu.2","W_po_A101_pu","W_po_A101_pu.2","W_po_A10701_pu","W_po_A10701_pu.2","W_po_A1381_pu","W_po_A1381_pu.2","W_po_A2725_pu","W_po_A2725_pu.2","W_po_A2749_pu","W_po_A2749_pu.2","W_po_A2759_pu","W_po_A2759_pu.2","W_po_T12400_pu","W_po_T12400_pu.2","W_po_T13148_pu","W_po_T13148_pu.2")
Jiparana<-c("W_po_T7625_jigu","W_po_T7625_jigu.2","W_po_T7627_jigu","W_po_T7627_jigu.2","W_po_A308_jigu","W_po_A308_jigu.2","W_po_A3264_jigu","W_po_A3264_jigu.2","W_po_A334_jigu","W_po_A334_jigu.2","W_po_A368_jigu","W_po_A368_jigu.2","W_po_T15888_jigu","W_po_T15888_jigu.2","W_po_T15889_jigu","W_po_T15889_jigu.2","W_po_A265_pu","W_po_A265_pu.2")
Machado<-c("W_po_T2190_ma","W_po_T2190_ma.2","W_po_T3157_ma","W_po_T3157_ma.2","W_po_T388_ma","W_po_T388_ma.2","W_po_T6424_ma","W_po_T6424_ma.2","W_po_A24413_ma","W_po_A24413_ma.2","W_po_A4235_ma","W_po_A4235_ma.2","W_po_A465_ma","W_po_A465_ma.2","W_po_A472_ma","W_po_A472_ma.2","W_po_A882_ma","W_po_A882_ma.2","W_po_A894_ma","W_po_A894_ma.2","W_po_A896_ma","W_po_A896_ma.2","W_po_A898_ma","W_po_A898_ma.2","W_po_A902_ma","W_po_A902_ma.2","W_po_J163_ma","W_po_J163_ma.2","W_po_J192_ma","W_po_J192_ma.2","W_po_J202_ma","W_po_J202_ma.2","W_po_J209_ma","W_po_J209_ma.2","W_po_J222_ma","W_po_J222_ma.2","W_po_J226_ma","W_po_J226_ma.2","W_po_J302_ma","W_po_J302_ma.2","W_po_J439_ma","W_po_J439_ma.2","W_po_J448_ma","W_po_J448_ma.2","W_po_J453_ma","W_po_J453_ma.2")
Roosevelt<-c("W_po_J335_roar","W_po_J335_roar.2","W_po_J340_roar","W_po_J340_roar.2","W_po_J391_roar","W_po_J391_roar.2","W_po_J421_roar","W_po_J421_roar.2","W_po_J644_roar","W_po_J644_roar.2","W_po_J766_roar","W_po_J766_roar.2","W_po_J768_roar","W_po_J768_roar.2")
Aripuana<-c("W_po_J561_arsu","W_po_J561_arsu.2","W_po_J563_arsu","W_po_J563_arsu.2","W_po_J610_arsu","W_po_J610_arsu.2","W_po_80798_arsu","W_po_80798_arsu.2","W_po_81036_arsu","W_po_81036_arsu.2","W_po_85474_arsu","W_po_85474_arsu.2","W_po_86085_arsu","W_po_86085_arsu.2","W_po_A891_arsu","W_po_A891_arsu.2")
Sucunduri<-c("W_po_T10153_suta","W_po_T10153_suta.2","W_po_T10234_suta","W_po_T10234_suta.2","W_po_T11887_suta","W_po_T11887_suta.2","W_po_T14615_suta","W_po_T14615_suta.2","W_po_85827_suta","W_po_85827_suta.2","W_po_81265_suta","W_po_81265_suta.2","W_po_A12265_suta","W_po_A12265_suta.2","W_po_A16298_suta","W_po_A16298_suta.2")
Para<-c("W_po_T10694_pa","W_po_T10694_pa.2","W_po_A9868_pa","W_po_A9868_pa.2","W_po_T12793_pa","W_po_T12793_pa.2","W_po_T7132_pa","W_po_T7132_pa.2","W_po_T9346_pa","W_po_T9346_pa.2","W_po_A7346_pa","W_po_A7346_pa.2","W_po_A14569_pa","W_po_A14569_pa.2","W_po_A16227_pa","W_po_A16227_pa.2","W_po_A12284_pa","W_po_A12284_pa.2","W_po_25484_pa","W_po_25484_pa.2")

snp  <- set.populations(snp, list(Inambari,Purus,Jiparana,Machado,Roosevelt,Aripuana,Sucunduri,Para))

####sliding window stats###
snp.slide <- sliding.window.transform(snp,width=5000, jump=2500,type=1,whole.data=TRUE)
snp.slide <- F_ST.stats(snp.slide, mode = "haplotype")

write.csv(snp.slide@Pi,file="pi.csv")

Inambari<-as.numeric(snp.slide@Pi[,1])
Purus<-as.numeric(snp.slide@Pi[,2])
Jiparana<-as.numeric(snp.slide@Pi[,3])
Machado<-as.numeric(snp.slide@Pi[,4])
Roosevelt<-as.numeric(snp.slide@Pi[,5])
Aripuana<-as.numeric(snp.slide@Pi[,6])
Sucunduri<-as.numeric(snp.slide@Pi[,7])
Para<-as.numeric(snp.slide@Pi[,8])

pi.rivs<-as.data.frame(cbind(Inambari,Purus,Jiparana,Machado,Aripuana,Roosevelt,Sucunduri,Para))

for (i in 1:length(names(pi.rivs))){
  for (j in 1:length(pi.rivs[,1])){
    counts=counts+1
    riv.pi.tab[counts,1]<-names(pi.rivs)[i]
    riv.pi.tab[counts,2]<-pi.rivs[j,i]
    riv.pi.tab[counts,3]<-as.character("Willisornis")
  }
}

setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Malacoptila/')
# Load the data
snp <- readData("./vcf/", format="VCF",include.unknown = T,big.data = T)

# Set populations
Inambari <- c("M_ru_A7875_In","M_ru_T14456_In","M_ru_T23245_In","M_ru_T23246_In","M_ru_T23416_In","M_ru_T23478_In")
Purus<-c("M_ru_A10311_pu","M_ru_A10329_pu","M_ru_A1380_pu","M_ru_A2741_pu","M_ru_A436_pu","M_ru_A440_pu")
Jiparana<-c("M_ru_T3228_jigu","M_ru_T7634_jigu","M_ru_A207_jigu")
Machado2<-c("M_ru_T368_ma","M_ru_T381_ma","M_ru_T476_ma","M_ru_T494_ma")
MachadoS<-c("M_ru_A474_ma","M_ru_J265_ma","M_ru_T3164_ma","M_ru_T13253_ma")
Roosevelt<-c("M_ru_J640_roar","M_ru_J676_roar")
Aripuana<-c("M_ru_80819_arsu","M_ru_81347_arsu","M_ru_85426_arsu")
Sucunduri<-c("M_ru_A9903_pa","M_ru_T14532_suta","M_ru_T14622_suta","M_ru_T753_suta","M_ru_T11780_suta","M_ru_T10184_suta","M_ru_A5487_suta","M_ru_A11834_suta","M_ru_A15176_suta","M_ru_77750_suta","M_ru_78182_suta","M_ru_85919_suta")
Para<-c("M_ru_T1649_pa","M_ru_T16553_pa","M_ru_T19782_pa","M_ru_T6500_pa","M_ru_T6577_pa","M_ru_A9235_pa","M_ru_A16195_pa","M_ru_T11079_pa","M_ru_T11238_pa","M_ru_T12541_pa")

snp  <- set.populations(snp, list(Inambari,Purus,Jiparana,Machado2,MachadoS,Roosevelt,Aripuana,Sucunduri,Para))

####sliding window stats###
snp.slide <- sliding.window.transform(snp,width=5000, jump=2500,type=1,whole.data=TRUE)
snp.slide <- F_ST.stats(snp.slide, mode = "haplotype")

write.csv(snp.slide@Pi,file="pi.csv")

Inambari<-as.numeric(snp.slide@Pi[,1])
Purus<-as.numeric(snp.slide@Pi[,2])
Jiparana<-as.numeric(snp.slide@Pi[,3])
Machado2<-as.numeric(snp.slide@Pi[,4])
Machado<-as.numeric(snp.slide@Pi[,5])
Roosevelt<-as.numeric(snp.slide@Pi[,6])
Aripuana<-as.numeric(snp.slide@Pi[,7])
Sucunduri<-as.numeric(snp.slide@Pi[,8])
Para<-as.numeric(snp.slide@Pi[,9])

pi.rivs<-as.data.frame(cbind(Inambari,Purus,Jiparana,Machado2,Machado,Roosevelt,Aripuana,Sucunduri,Para))

for (i in 1:length(names(pi.rivs))){
  for (j in 1:length(pi.rivs[,1])){
    counts=counts+1
    riv.pi.tab[counts,1]<-names(pi.rivs)[i]
    riv.pi.tab[counts,2]<-pi.rivs[j,i]
    riv.pi.tab[counts,3]<-as.character("Malacoptila")
  }
}


setwd('~/AMNH/Roosevelt_Project/Data/PopGenome/Hypocnemis/')

snp <- readData(path = './vcf/',include.unknown = T,format = 'VCF',big.data = T)

Purus<-c("H_pe_A24346_pu")
Jiparana<-c("H_oc_A311_jigu","H_oc_T15847_jigu","H_oc_T355_jigu")
Machado<-c("H_ro_A547_ma","H_ro_J213_ma","H_ro_J305_ma","H_ro_J508_ma","H_ro_J774_ma","H_ro_J775_ma","H_ro_J796_ma","H_ro_T1842_ma","H_ro_T366_ma","H_ro_T385_ma","H_ro_T471_ma")
Roosevelt<-c("H_st_J364_roar","H_st_J368_roar","H_st_J370_roar","H_st_J374_roar","H_st_J408_roar","H_st_J621_roar","H_st_J664_roar","H_st_J665_roar","H_st_J711_roar","H_st_J762_roar","H_st_J765_roar")
Aripuana<-c("H_ro_A551_ma","H_ro_A409_ma","H_ro_A410_ma","H_ro_A521_ma","H_st_A272_arsu","H_st_A273_arsu","H_st_J525_arsu","H_st_J530_arsu","H_st_J536_arsu","H_st_J572_arsu","H_st_77860_arsu","H_st_78249_arsu","H_st_80727_arsu","H_st_80800_arsu","H_st_81143_arsu","H_st_85680_arsu","H_st_85970_arsu")
Sucunduri<-c("H_st_A9955_pa","H_st_86405_suta","H_st_A14546_suta","H_st_A4899_suta","H_st_T10207_suta","H_st_T11900_suta","H_st_T12194_suta","H_st_T24564_suta","H_st_T7114_suta","H_st_81279_suta")
Para<-c("H_st_A11597_pa","H_st_A15208_pa","H_st_A15210_pa","H_st_A16571_pa","H_st_A7597_pa","H_st_T16744_pa","H_st_T17858_pa")

snp  <- set.populations(snp, list(Purus,Jiparana,Machado,Roosevelt,Aripuana,Sucunduri,Para))

snp@populations # check if it worked


####sliding window stats###
snp.slide <- sliding.window.transform(snp,width=5000, jump=2500,type=1,whole.data=TRUE)
snp.slide <- F_ST.stats(snp.slide, mode = "haplotype")

write.csv(snp.slide@Pi,file="pi.csv")

Purus<-as.numeric(snp.slide@Pi[,1])
Jiparana<-as.numeric(snp.slide@Pi[,2])
Machado<-as.numeric(snp.slide@Pi[,3])
Roosevelt<-as.numeric(snp.slide@Pi[,4])
Aripuana<-as.numeric(snp.slide@Pi[,5])
Sucunduri<-as.numeric(snp.slide@Pi[,6])
Para<-as.numeric(snp.slide@Pi[,7])

pi.rivs<-as.data.frame(cbind(Purus,Jiparana,Machado,Roosevelt,Aripuana,Sucunduri,Para))

for (i in 1:length(names(pi.rivs))){
  for (j in 1:length(pi.rivs[,1])){
    counts=counts+1
    riv.pi.tab[counts,1]<-names(pi.rivs)[i]
    riv.pi.tab[counts,2]<-pi.rivs[j,i]
    riv.pi.tab[counts,3]<-as.character("Hypocnemis")
  }
}

###Finalize and save matrix and get fst/1-fst

colnames(riv.pi.tab)<-c("River","Pi","Species")
write.csv(riv.pi.tab,file="~/AMNH/Roosevelt_Project/Data/PopGenome/Pairwise.FSTs.all.csv")

tab<-read.csv("~/AMNH/Roosevelt_Project/Data/PopGenome/Pairwise.FSTs.all.csv",header=T)
attach(tab)

library(ggplot2)
library(wesanderson)

cols<-c(wes_palette(name="Darjeeling1"),wes_palette(name="Darjeeling2"))

p <- ggplot(riv.pi.tab, aes(x=River, y=Pi,fill=Species),reorder(names(pi.rivs))) + 
  geom_boxplot(outlier.shape=NA) +
  xlab("") +
  ylab("Pi") +
#  ylim(0,2) +
  scale_fill_manual(values = cols) +
  ggtitle("Haplotype Diversity by Population")
p

p2 <- ggplot(riv.pi.tab, aes(x=Species, y=Pi,fill=River),reorder(names(pi.rivs))) + 
  geom_boxplot(outlier.shape=NA) +
  xlab("") +
  ylab("Pi") +
  #  ylim(0,2) +
  scale_fill_manual(values = cols) +
  ggtitle("Haplotype Diversity by Population")
p2

pdf('../pi.boxplot.pdf',width=13,height=4.5, bg='transparent')
p
dev.off()

pdf('../pi.boxplot2.pdf',width=13,height=4.5, bg='transparent')
p2
dev.off()

####Thamnophilus inset
tham<-tab[Species=="Thamnophilus",]
names(tham)
tham[River=="Inambari",]

p3 <- ggplot(tham, aes(x=Species, y=Pi,fill=River)) + 
  geom_boxplot(outlier.shape=NA) +
  xlab("") +
  ylab("Pi") +
  #  ylim(0,2) +
  scale_fill_manual(values = cols) +
  ggtitle("Thamnophilus")
p3

pdf('../pi.boxplot.tham.pdf',width=4.5,height=3, bg='transparent')
p3
dev.off()
