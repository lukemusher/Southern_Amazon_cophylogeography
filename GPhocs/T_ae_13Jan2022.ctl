GENERAL-INFO-START

	seq-file            T_ae_2000_rand.gphocs
	trace-file          T_ae_var_mcmc_13Jan2022.log
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
		name		InPu
		samples		T_ae_A2833_pu d	T_ae_A439_pu d	T_ae_A8087_In d	T_ae_T12313_pu d	T_ae_T13219_pu d
	POP-END

	POP-START
		name		JiMaRo
		samples		T_ae_A317_jigu d	T_ae_A509_ma d	T_ae_J319_roar d	T_ae_J419_roar d	T_ae_J678_roar d
			POP-END

	POP-START
		name		ArSuTa
		samples		T_ae_80508_arsu d	T_ae_80716_arsu d	T_ae_81278_arsu d	T_ae_85274_suta d	T_ae_86147_arsu d
	POP-END

	POP-START
		name		Para
		samples		T_ae_A15279_pa d	T_ae_T10679_pa d	T_ae_T13575_pa d	T_ae_T24591_suta d	T_ae_T8268_pa d
	POP-END
	

CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			JiMaRoArSuTa
		children		JiMaRo		ArSuTa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			InPu_JiMaArSuTa
		children		InPu		JiMaRoArSuTa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			root
		children		InPu_JiMaArSuTa		Para
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
       source  Para
       target  JiMaRoArSuTa
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  JiMaRoArSuTa
       target  Para
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  InPu
       target  JiMaRo
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  JiMaRo
       target  InPu
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  InPu
       target  ArSuTa
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  ArSuTa
       target  InPu
       mig-rate-print 0.1
	BAND-END

MIG-BANDS-END
