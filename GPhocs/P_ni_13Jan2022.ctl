GENERAL-INFO-START

	seq-file            P_ni_2000_rand.gphocs
	trace-file          P_ni_var_mcmc_13Jan2022.log
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
		samples		P_ni_A7862_In d	P_ni_A7911_In d	P_ni_A7928_In d
	POP-END

	POP-START
		name		Puru
		samples		P_ni_80034_pu d	P_ni_T22153_jigu d	P_ni_T3817_pu d	P_ni_T4404_pu d	P_ni_T5974_pu d
			POP-END

	POP-START
		name		Para
		samples		P_ni_A14342_pa d	P_ni_T10673_pa d	P_ni_T11222_pa d	P_ni_T12854_pa d	P_ni_T1642_pa d
	POP-END

	POP-START
		name		JiMaRo
		samples		P_ni_J210_ma d	P_ni_J361_roar d	P_ni_J363_roar d	P_ni_J371_roar d	P_ni_T3261_jigu d
	POP-END

	POP-START
		name		ArSuTa
		samples		P_ni_80684_arsu d	P_ni_86072_arsu d	P_ni_J602_arsu d	P_ni_T14543_suta d	P_ni_T9076_suta d
	POP-END

CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			InPu
		children		Inam		Puru
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			Rondonia
		children		JiMaRo		ArSuTa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			Para_Rondonia
		children		Rondonia		Para
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			root
		children		Para_Rondonia		InPu
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
       target  Puru
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Puru
       target  Para
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Rondonia
       target  Puru
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Puru
       target  Rondonia
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  ArSuTa
       target  Para
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Para
       target  ArSuTa
       mig-rate-print 0.1
	BAND-END


	BAND-START		
       source  JiMaRo
       target  ArSuTa
       mig-rate-print 0.1
	BAND-END


	BAND-START		
       source  ArSuTa
       target  JiMaRo
       mig-rate-print 0.1
	BAND-END

MIG-BANDS-END
