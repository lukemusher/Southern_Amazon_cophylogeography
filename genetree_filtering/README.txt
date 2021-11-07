###Loci are min 75percent individuals represented
###Loci are min 400bp long
###Loci trees are rooted using the reference genome sequence
###Loci are then filtered using a modified version of the protocol from Musher et al 2018 
###For the filtering protocol
###We estimate 
###(1) the number of variablse sites, 
###(2) the likelihood ratio of clocklike to non-clocklike models (phytools)
###(3) the likelihood of the gene tree given the sequence alignment (ape)
###(4) the proportion of variable sites (variable sites/total sites)

##These data are used to generate a glm where genetree loglikelihood is a function of the other three variables
###Outliers and extreme values in the model are discarded by using cooks distance (anything over 2*mean cooks distance)
###We then show that #var sites is the best predictor of better clocklikeness and likelihood
