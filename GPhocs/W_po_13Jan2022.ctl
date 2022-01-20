GENERAL-INFO-START

	seq-file            W_po_2000_rand.gphocs
	trace-file          W_po_var_mcmc_13Jan2022.log				
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
		samples		W_po_A2725_pu d	W_po_A424_pu d	W_po_A815_In d	W_po_T12362_In d	W_po_T13148_pu d
	POP-END

	POP-START
		name		Rondonia
		samples		W_po_A3264_jigu d	W_po_A472_ma d	W_po_A882_ma d	W_po_T15889_jigu d	W_po_T2190_ma d
			POP-END

	POP-START
		name		SuTa
		samples		W_po_81265_suta d	W_po_A16298_suta d	W_po_T10153_suta d	W_po_T10234_suta d	W_po_T11887_suta d
	POP-END

	POP-START
		name		Para
		samples		W_po_A12284_pa d	W_po_A14569_pa d	W_po_T12793_pa d	W_po_T7132_pa d	W_po_T9346_pa d
	POP-END

CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			Para_SuTa
		children		Para		SuTa
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			InPu_Rondonia
		children		InPu		Rondonia
		tau-initial	0.001
#		tau-alpha		1.0	
#		tau-beta		50000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END

	POP-START
		name			root
		children		InPu_Rondonia	Para_SuTa
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
       source  InPu
       target  Para
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Para
       target  InPu
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  SuTa
       target  Rondonia
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Rondonia
       target  SuTa
       mig-rate-print 0.1
	BAND-END

MIG-BANDS-END
