## all           : clean saved results and get results to populate report
all: clean cor_poisson report

report:
	bash _build_chapter.sh

	Rscript -e "rmarkdown::render('code/posthoc_analyses.Rmd',           output_file='../posthoc_analyses_$$(date +%Y%m%d).pdf')"

	#Rscript -e "rmarkdown::render('code/comparative_immunogenicity.Rmd', output_file='../comparative_immunogenicity_$$(date +%Y%m%d).pdf')"
	#Rscript -e "rmarkdown::render('code/diproperm_run.Rmd',              output_file='../diproperm_run_$$(date +%Y%m%d).pdf')"
	#Rscript code/Plotting_Main_YutongJin.R


## clean         : delete saved bootstrap results
clean:
	rm -rf output/$(TRIAL)/*

## cor_poisson     : create tables and figures
cor_poisson:
ifeq ($(TRIAL),$(filter $(TRIAL), iliad_ib202p))
	Rscript code/cor_poisson_iliad_ib202p.R C0iliad_ib202p_PPAI
	Rscript code/cor_poisson_iliad_ib202p.R C0iliad_ib202p_PPAI_C14only
	#Rscript code/cor_poisson_iliad_ib202p.R C0iliad_ib202p
	#Rscript code/cor_poisson_iliad_ib202p.R C0iliad_ib202p_C14only
endif


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_poisson
