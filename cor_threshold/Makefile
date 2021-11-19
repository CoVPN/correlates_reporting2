## all              : run threshold analysis and get plots
all: clean params data threshold plot

allbutclean: params data threshold plot

## threshold        : run threshold analysis 
threshold:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/Run_Threshold_analysis.R D29
	Rscript code/Run_Threshold_analysis.R D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_pooled_real janssen_na_mock janssen_na_real janssen_la_mock janssen_la_real janssen_sa_mock janssen_sa_real))
	Rscript code/Run_Threshold_analysis.R D29
	Rscript code/Run_Threshold_analysis.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705))
	Rscript code/Run_Threshold_analysis.R D210
endif

## data             : prepare data for analysis
data: 
	Rscript code/clean_data.R
	
## plot             : plot results
plot:  
	Rscript code/plotting.R

# parameters that may vary across runs
params:  
	Rscript code/params.R

## clean            : delete contents of data_clean/ and figs/
clean:
	rm -f data_clean/*.csv
	rm -f data_clean/Thresholds_by_marker/*.csv
	rm -f figs/simultaneous_CI/*.pdf
	rm -f figs/pointwise_CI/*.pdf
	rm -f figs/pointwise_CI/*.rds
	rm -f figs/simultaneous_CI/*.rds

## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"

## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('report.Rmd', output_file='cor_threshold_report_$(TRIAL).pdf')"


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean plot params threshold data
