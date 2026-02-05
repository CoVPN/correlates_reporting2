all: cop_mediation

cop_mediation:
ifeq ($(TRIAL),$(filter $(TRIAL), vat08_combined))
    ifeq ($(stage),1) 
	echo "stage 1"
    else
	echo "stage 2"
	Rscript code/mediation_report_sanofi_stage2.R FR
	Rscript code/mediation_report_sanofi_stage2.R Day43
    endif
endif

