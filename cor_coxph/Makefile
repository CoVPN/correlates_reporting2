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
	Rscript code/cor_graphical_time_to_event.R D29D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_na_mock janssen_la_mock janssen_sa_mock))
	Rscript code/cor_coxph.R D29
	Rscript code/cor_coxph.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_realbAb janssen_pooled_realPsV janssen_pooled_realADCP janssen_na_realbAb janssen_na_realPsV janssen_na_realADCP janssen_la_realbAb janssen_la_realPsV janssen_la_realADCP janssen_sa_realbAb janssen_sa_realPsV janssen_sa_realADCP))
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmed
	Rscript code/cor_coxph.R D29IncludeNotMolecConfirmedstart1
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705 hvtn705V1V2))
	Rscript code/cor_coxph.R D210
endif


## report        : create a standalone report
meta_report:
	Rscript code/cor_coxph_collate.R
	Rscript code/cor_coxph_meta_VE_curves.R
	Rscript -e "rmarkdown::render('meta_cve.Rmd', output_file='covpn_correlates_cor_coxph_meta_cve.pdf')"
	Rscript -e "rmarkdown::render('collate_ensemble.Rmd', output_file='covpn_correlates_cor_coxph_ensemble_collated.pdf')"


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph
