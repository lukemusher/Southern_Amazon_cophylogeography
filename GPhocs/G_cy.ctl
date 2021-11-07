GENERAL-INFO-START

	seq-file            G_cy_2000_rand.gphocs
	trace-file          G_cy_var_mcmc.log				
#	num-loci            5000
	locus-mut-rate      VAR 1.0
	
	burn-in	          50000
	mcmc-iterations	  500000
	iterations-per-log  10
	logs-per-line       10
#	start-mig           10000

	find-finetunes		FALSE
	finetune-coal-time	0.01		
	finetune-mig-time	0.3		
	finetune-theta		0.04
	finetune-mig-rate	0.02
	finetune-tau		0.0000008
	finetune-mixing		0.003
	finetune-locus-rate 0.3
	
	tau-theta-print		1.0
	tau-theta-alpha		3.0			# for STD/mean ratio of 100%
	tau-theta-beta		1000.0		# for mean of 1e-4

	mig-rate-print		1.0
	mig-rate-alpha		0.002
	mig-rate-beta		0.00001

GENERAL-INFO-END

CURRENT-POPS-START	

	POP-START
		name		Inam
		samples		G_cy_T12385_In d	G_cy_T23196_In d	G_cy_T7636_In d	G_cy_T3384_jigu d	G_cy_T3385_jigu d
	POP-END

	POP-START
		name		Puru
		samples		G_cy_82508_pu d	G_cy_T12275_pu d	G_cy_T12279_pu d	G_cy_T12392_pu d	G_cy_T26252_pu d
			POP-END

	POP-START
		name		MaRo
		samples		G_cy_J296_ma d	G_cy_J691_roar d	G_cy_T13251_ma d	G_cy_T363_ma d	G_cy_T364_ma d
	POP-END

	POP-START
		name		ArSuTa
		samples		G_cy_80801_arsu d	G_cy_80826_arsu d	G_cy_81118_arsu d	G_cy_85499_arsu d	G_cy_T18563_suta d
	POP-END

	POP-START
		name		Para
		samples		G_cy_T10897_pa d	G_cy_T11062_pa d	G_cy_T1705_pa d	G_cy_T19429_pa d	G_cy_T19765_pa d
	POP-END
	

CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			Para_MaRo
		children		Para		MaRo
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			Para_MaRo_ArSuTa
		children		Para_MaRo		ArSuTa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			Inam_Puru
		children		Inam		Puru
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			root
		children		Inam_Puru	Para_MaRo_ArSuTa
		tau-initial	0.003
#		tau-alpha	3.0
#		tau-beta		1000.0	
#		theta-alpha	5.0
#		theta-beta		1000.0	
		finetune-tau			0.00000286
	POP-END

ANCESTRAL-POPS-END

MIG-BANDS-START	
	BAND-START		
       source  MaRo
       target  Puru
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Puru
       target  Inam
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Puru
       target  MaRo
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Inam
       target  Puru
       mig-rate-print 0.1
	BAND-END

MIG-BANDS-END
