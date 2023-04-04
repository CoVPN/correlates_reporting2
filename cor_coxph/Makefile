## all           : clean saved results and get results to populate report
all: clean cor_coxph 

## clean         : delete saved bootstrap results
clean:
	rm -rf output/$(TRIAL)/*

## cor_coxph     : create tables and figures
cor_coxph:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/cor_coxph.R D29
	Rscript code/cor_coxph.R D57
	Rscript code/cor_coxph.R D57over29
	Rscript code/cor_coxph.R D29start28
	Rscript code/cor_graphical_time_to_event.R D29D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_na_mock janssen_la_mock janssen_sa_mock))
	Rscript code/cor_coxph.R D29
	Rscript code/cor_coxph.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_EUA janssen_na_EUA janssen_la_EUA janssen_sa_EUA))
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmed
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmedstart1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_partA  janssen_la_partA ))
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmed
	Rscript code/cor_coxph.R D29SevereIncludeNotMolecConfirmed
	Rscript code/cor_coxph.R D29ModerateIncludeNotMolecConfirmed
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_na_partA janssen_sa_partA))
	# not enough cases for severe endpoints
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmed
	Rscript code/cor_coxph.R D29ModerateIncludeNotMolecConfirmed
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_partAsenior janssen_pooled_partAnonsenior janssen_na_partAsenior janssen_na_partAnonsenior janssen_la_partAsenior janssen_la_partAnonsenior janssen_sa_partAnonsenior)) 
	# janssen_sa_partAsenior not included for having too few cases
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmed
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705 hvtn705V1V2 hvtn705second hvtn705secondprimary))
	Rscript code/cor_coxph.R D210
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19))
	Rscript code/cor_coxph.R D35
else ifeq ($(TRIAL),$(filter $(TRIAL), azd1222))
	Rscript code/cor_coxph.R D29
	Rscript code/cor_coxph.R D57
	Rscript code/cor_coxph.R D57over29
	Rscript code/cor_coxph.R D29start28
	Rscript code/cor_graphical_time_to_event.R D29D57
else ifeq ($(TRIAL),$(filter $(TRIAL), azd1222_bAb)) ## for some reason if this is merged with azd1222, azd1222 won't run
	Rscript code/cor_coxph.R D29
	Rscript code/cor_coxph.R D57
	Rscript code/cor_coxph.R D57over29
	Rscript code/cor_coxph.R D29start28
	Rscript code/cor_graphical_time_to_event.R D29D57
else ifeq ($(TRIAL),$(filter $(TRIAL), vat08m_naive vat08m_nonnaive))
	Rscript code/cor_coxph.R D22
	Rscript code/cor_coxph.R D43
	Rscript code/cor_coxph.R D22omi
	Rscript code/cor_coxph.R D43omi
	Rscript code/cor_graphical_time_to_event.R D22D43
else ifeq ($(TRIAL),$(filter $(TRIAL), profiscov))
	Rscript code/cor_coxph.R D43
	Rscript code/cor_coxph.R D43start48
	Rscript code/cor_coxph.R D91
else ifeq ($(TRIAL),$(filter $(TRIAL), profiscov_lvmn))
	Rscript code/cor_coxph.R D43
endif

collate:
	Rscript code/cor_coxph_collate_meta_VE.R
	Rscript -e "rmarkdown::render('collate_ensemble.Rmd', output_file='covpn_correlates_cor_coxph_janssen_collated.pdf')"

meta:
	Rscript code/cor_coxph_meta_VE_curves.R
	Rscript -e "rmarkdown::render('meta_cve.Rmd', output_file='covpn_correlates_cor_coxph_1meta.pdf')"


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph
