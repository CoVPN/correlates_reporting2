all: cop_mediation report


cop_mediation:
ifeq ($(TRIAL),$(filter $(TRIAL), vat08_combined))
	ifeq ($(stage),1) 
		echo "stage 1"
	else
		echo "stage 2"
		Rscript code/cop_mediation_sanofi_stage2.R FR
		Rscript code/cop_mediation_sanofi_stage2.R Day43
	endif

else ifeq ($(TRIAL),$(filter $(TRIAL), nextgen_mock))
	Rscript code/cop_mediation_nextgen_mock.R D31toM12_nextgen_mock_sera
	#Rscript code/cop_mediation_nextgen_mock.R D31toM12_nextgen_mock_tcell

endif


report:
	Rscript code/make_tables.R
