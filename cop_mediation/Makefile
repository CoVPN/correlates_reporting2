## all              : run mediation analysis
all:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/mediation.R D29
	Rscript code/mediation.R D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_na_mock janssen_la_mock janssen_sa_mock))
	Rscript code/mediation.R D29
	Rscript code/mediation.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_real janssen_pooled_realADCP janssen_na_real janssen_na_realADCP janssen_la_real janssen_la_realADCP janssen_sa_real janssen_sa_realADCP))
	Rscript code/mediation.R D29IncludeNotMolecConfirmed
	Rscript code/mediation.R D29IncludeNotMolecConfirmedstart1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_partA  janssen_la_partA ))
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmed
	Rscript code/cor_coxph.R D29SevereIncludeNotMolecConfirmed
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_partAsenior janssen_pooled_partAnonsenior janssen_na_partA janssen_na_partAsenior janssen_na_partAnonsenior janssen_la_partAsenior janssen_la_partAnonsenior janssen_sa_partA  janssen_sa_partAnonsenior)) # janssen_sa_partAsenior has too few cases
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmed
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705 hvtn705V1V2 hvtn705second hvtn705secondprimary))
	Rscript code/mediation.R D210
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19))
	Rscript code/mediation.R D35
else ifeq ($(TRIAL),$(filter $(TRIAL), azd1222))
	Rscript code/mediation.R D29
	Rscript code/mediation.R D57
	Rscript code/mediation.R D29start28
else ifeq ($(TRIAL),$(filter $(TRIAL), azd1222_bAb))
	Rscript code/cor_coxph.R D29
	Rscript code/cor_coxph.R D57
	Rscript code/cor_coxph.R D29start28
endif

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
# .PHONY: help params data report results
