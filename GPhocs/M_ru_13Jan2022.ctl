GENERAL-INFO-START

	seq-file            M_ru_2000_rand.gphocs
	trace-file          M_ru_var_mcmc_13Jan2022.log				
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
	mig-rate-alpha		1.2
	mig-rate-beta		0.01

GENERAL-INFO-END

CURRENT-POPS-START	

	POP-START
		name		Inam
		samples		M_ru_A7875_In d	M_ru_T14456_In d	M_ru_T23246_In d	M_ru_T23416_In d	M_ru_T23478_In d
	POP-END

	POP-START
		name		Puru
		samples		M_ru_A10311_pu d	M_ru_A1380_pu d	M_ru_A2741_pu d	M_ru_A436_pu d	M_ru_A440_pu d
			POP-END

	POP-START
		name		JiMa
		samples		M_ru_A207_jigu d	M_ru_T7634_jigu d	M_ru_T13253_ma d	M_ru_A474_ma d	M_ru_J265_ma d
	POP-END

	POP-START
		name		MaArTa
		samples		M_ru_77750_suta d	M_ru_81347_arsu d	M_ru_85426_arsu d	M_ru_J640_roar d	M_ru_T14622_suta d	
	POP-END

	POP-START
		name		Para
		samples		M_ru_A9235_pa d	M_ru_A9903_pa d	M_ru_T11238_pa d	M_ru_T1649_pa d	M_ru_T16553_pa d	
	POP-END
	

CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			Para_MaArTa
		children		Para		MaArTa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			Puru_JiMa
		children		Puru		JiMa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			Inam_Puru_JiMa
		children		Inam		Puru_JiMa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			root
		children		Inam_Puru_JiMa	Para_MaArTa
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
       source  MaArTa
       target  JiMa
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  JiMa
       target  MaArTa
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  MaArTa
       target  Puru_JiMa
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Puru_JiMa
       target  MaArTa
       mig-rate-print 0.1
	BAND-END

MIG-BANDS-END

