## all              : run threshold analysis and get plots
all: clean threshold report

allbutclean: threshold report

## clean            : delete contents of data_clean/ and figs/
clean:
	rm -f data_clean/*.csv
	rm -f data_clean/Thresholds_by_marker/*.csv
	rm -f figs/simultaneous_CI/*.pdf
	rm -f figs/pointwise_CI/*.pdf
	rm -f figs/pointwise_CI/*.rds
	rm -f figs/simultaneous_CI/*.rds

threshold:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/clean_data.R 						D29
	Rscript code/Run_Threshold_analysis.R D29
	Rscript code/plotting.R 							D29

	Rscript code/clean_data.R 						D57
	Rscript code/Run_Threshold_analysis.R D57
	Rscript code/plotting.R 							D57
	
else ifeq ($(TRIAL),$(filter $(TRIAL), moderna_boost))
	Rscript code/clean_data.R 						BD29
	Rscript code/Run_Threshold_analysis.R BD29
	Rscript code/plotting.R 							BD29
	
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_trial_EUA janssen_pooled_mock janssen_pooled_EUA janssen_na_mock janssen_na_EUA janssen_la_mock janssen_la_EUA janssen_sa_mock janssen_sa_EUA))
	Rscript code/clean_data.R 						D29
	Rscript code/Run_Threshold_analysis.R D29
	Rscript code/plotting.R 							D29
	
	Rscript code/clean_data.R 						D29start1
	Rscript code/Run_Threshold_analysis.R D29start1
	Rscript code/plotting.R 							D29start1
	
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_partA janssen_la_partA))
	Rscript code/clean_data.R 						D29IncludeNotMolecConfirmed
	Rscript code/Run_Threshold_analysis.R D29IncludeNotMolecConfirmed
	Rscript code/plotting.R 							D29IncludeNotMolecConfirmed

	Rscript code/clean_data.R 						D29SevereIncludeNotMolecConfirmed
	Rscript code/Run_Threshold_analysis.R D29SevereIncludeNotMolecConfirmed
	Rscript code/plotting.R 							D29SevereIncludeNotMolecConfirmed

	Rscript code/clean_data.R 						D29ModerateIncludeNotMolecConfirmed
	Rscript code/Run_Threshold_analysis.R D29ModerateIncludeNotMolecConfirmed
	Rscript code/plotting.R 							D29ModerateIncludeNotMolecConfirmed
	
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_na_partA janssen_sa_partA))
	Rscript code/clean_data.R 						D29IncludeNotMolecConfirmed
	Rscript code/Run_Threshold_analysis.R D29IncludeNotMolecConfirmed
	Rscript code/plotting.R 							D29IncludeNotMolecConfirmed
	
	Rscript code/clean_data.R 						D29ModerateIncludeNotMolecConfirmed
	Rscript code/Run_Threshold_analysis.R D29ModerateIncludeNotMolecConfirmed
	Rscript code/plotting.R 							D29ModerateIncludeNotMolecConfirmed

else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705second))
	Rscript code/clean_data.R 						D210
	Rscript code/Run_Threshold_analysis.R D210
	Rscript code/plotting.R 							D210

endif


report:
	bash _build_chapter.sh


## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean plot params threshold data
