GENERAL-INFO-START

	seq-file            H_st_2000_rand.gphocs
	trace-file          H_st_var_mcmc.log
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
		name		Puru
		samples		H_pe_A24346_pu d
	POP-END

	POP-START
		name		JiGu
		samples		H_oc_T15847_jigu d	H_oc_A311_jigu d
			POP-END

	POP-START
		name		Mach
		samples		H_ro_A551_ma d	H_ro_J508_ma d	H_ro_J774_ma d	H_ro_J775_ma d	H_ro_T366_ma d
	POP-END

	POP-START
		name		RoArSuTa
		samples		H_st_81143_arsu d	H_st_A272_arsu d	H_st_J664_roar d	H_st_J665_roar d	H_st_T7114_suta d
	POP-END

	POP-START
		name		Para
		samples		H_st_A11597_pa d	H_st_A15210_pa d	H_st_A16571_pa d	H_st_T16744_pa d
	POP-END
	

CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			Para_Mach
		children		Para		Mach
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			Para_Mach_RoArSuTa
		children		Para_Mach		RoArSuTa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			Para_Rondonia
		children		Para_Mach_RoArSuTa		JiGu
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			root
		children		Para_Rondonia	Puru
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
       source  JiGu
       target  RoArSuTa
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  RoArSuTa
       target  JiGu
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  JiGu
       target  Para
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Para
       target  JiGu
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  JiGu
       target  Puru
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Puru
       target  JiGu
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  JiGu
       target  Mach
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Mach
       target  JiGu
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Puru
       target  Mach
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Mach
       target  Puru
       mig-rate-print 0.1
	BAND-END

MIG-BANDS-END
