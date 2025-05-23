## all              : twophase_plots + demog_plots
all: clean immuno_graphical

## twophase_plots   : rcdfs, boxplots, scatter plots of assays saved in output/
#twophase_plots: code/descriptive_graphics_two_phase_plots.R \
#  params plot_data make_assay_lim
#	Rscript code/descriptive_graphics_two_phase_plots.R

## demog_plots      : rcdfs, boxplots by demographics saved in output/demographics, 
## violin_pair_plots      : violin, pair plots saved in output, code/descriptive_graphics_violin_pair.R
## code/descriptive_graphics_demographic_strata.R, code/descriptive_graphics_violin_pair.R
immuno_graphical: 
ifeq ($(TRIAL),$(filter $(TRIAL), nextgen_mock vat08_combined janssen_partA_VL))
	Rscript code/descriptive_graphics_data_preprocess.R
	Rscript code/make_assay_lim.R
	Rscript code/descriptive_graphics_two_phase_plots.R
	Rscript code/descriptive_graphics_violin_pair.R
## Rscript code/descriptive_graphics_violin_pair_extra.R
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_partA))
	Rscript code/descriptive_graphics_data_preprocess.R
	Rscript code/make_assay_lim.R
	Rscript code/descriptive_graphics_two_phase_plots.R
	#Rscript code/descriptive_graphics_demographic_strata.R
	Rscript code/descriptive_graphics_violin_pair.R
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19_stage2))
	Rscript code/descriptive_graphics_data_preprocess.R
	Rscript code/make_assay_lim.R
	Rscript code/descriptive_graphics_two_phase_plots.R
else
	Rscript code/descriptive_graphics_data_preprocess.R
	Rscript code/make_assay_lim.R
	Rscript code/descriptive_graphics_two_phase_plots.R
	Rscript code/descriptive_graphics_demographic_strata.R
endif	

## plot_data        : relabel and clean data sets for plotting
# creates clean_data/long_twophase_data.rds and clean_data/twophase_data.rds
#plot_data: code/descriptive_graphics_data_preprocess.R params
#	Rscript code/descriptive_graphics_data_preprocess.R

# parameters that may vary across runs
#params: code/params.R

# define the axis limits for the assays
#make_assay_lim: params plot_data code/make_assay_lim.R 
#	Rscript code/make_assay_lim.R

## clean            : delete contents of data_clean/ and figs/
clean:
	rm -f data_clean/*.rds
	rm -f output/$(TRIAL)/*.pdf
	rm -f output/$(TRIAL)/demographics/*.pdf
	rm -f Rplots.pdf

## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean immuno_graphical
