## all              : run threshold analysis and get plots
all: clean threshold plot

allbutclean: threshold plot

## threshold        : run threshold analysis 
threshold:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
    Rscript code/params.R D29
  Rscript code/clean_data.R D29
    Rscript code/Run_Threshold_analysis.R D29
    
    Rscript code/params.R D57
    Rscript code/clean_data.R D57
    Rscript code/Run_Threshold_analysis.R D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_trial_EUA janssen_pooled_mock janssen_pooled_EUA janssen_na_mock janssen_na_EUA janssen_la_mock janssen_la_EUA janssen_sa_mock janssen_sa_EUA))
    Rscript code/params.R D29
    Rscript code/clean_data.R D29
    Rscript code/Run_Threshold_analysis.R D29
    
    Rscript code/params.R D29start1
    Rscript code/clean_data.R D29start1
    Rscript code/Run_Threshold_analysis.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705))
    Rscript code/params.R D210
    Rscript code/clean_data.R D210
    Rscript code/Run_Threshold_analysis.R D210
endif

## data             : prepare data for analysis
data: 
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
    Rscript code/clean_data.R D29
    Rscript code/clean_data.R D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_trial_EUA janssen_pooled_mock janssen_pooled_EUA janssen_na_mock janssen_na_EUA janssen_la_mock janssen_la_EUA janssen_sa_mock janssen_sa_EUA))
    Rscript code/clean_data.R D29
    Rscript code/clean_data.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705))
    Rscript code/clean_data.R D210
endif
     
    
## plot             : plot results
plot:  
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
    Rscript code/params.R D29
    Rscript code/plotting.R D29
    Rscript code/params.R D57
    Rscript code/plotting.R D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_trial_EUA janssen_pooled_mock janssen_pooled_EUA janssen_na_mock janssen_na_EUA janssen_la_mock janssen_la_EUA janssen_sa_mock janssen_sa_EUA))
    Rscript code/params.R D29
    Rscript code/plotting.R D29
    Rscript code/params.R D29start1
    Rscript code/plotting.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705))
    Rscript code/params.R D210
    Rscript code/plotting.R D210
endif

# parameters that may vary across runs
params:  
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
    Rscript code/params.R D29
    Rscript code/params.R D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_trial_EUA janssen_pooled_mock janssen_pooled_EUA janssen_na_mock janssen_na_EUA janssen_la_mock janssen_la_EUA janssen_sa_mock janssen_sa_EUA))
    Rscript code/params.R D29
    Rscript code/params.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705))
    Rscript code/params.R D210
endif


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
