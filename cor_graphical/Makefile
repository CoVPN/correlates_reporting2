## all            : clear current results and make figures
all: clean figs

## figs           : clean data and make figures
figs: 
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/cor_data_preprocess.R D29D57 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D29D57
	
	Rscript code/cor_data_preprocess_2.R D29 # Brian's
	Rscript code/CoR_assay_graphics.R D29
	Rscript code/CoR_wrcdf_with_VE_lines.R D29
	
	Rscript code/cor_data_preprocess_2.R D57 # Brian's
	Rscript code/CoR_assay_graphics.R D57
	Rscript code/CoR_wrcdf_with_VE_lines.R D57

else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_na_mock janssen_la_mock janssen_sa_mock))
	Rscript code/cor_data_preprocess.R D29 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D29
	
	Rscript code/cor_data_preprocess.R D29start1 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D29start1
	
	Rscript code/cor_data_preprocess_2.R D29 # Brian's
	Rscript code/CoR_assay_graphics.R D29
	Rscript code/CoR_wrcdf_with_VE_lines.R D29
	
	Rscript code/cor_data_preprocess_2.R D29start1 # Brian's
	Rscript code/CoR_assay_graphics.R D29start1
	Rscript code/CoR_wrcdf_with_VE_lines.R D29start1
	
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_realbAb janssen_pooled_realADCP janssen_pooled_realPsV janssen_na_real janssen_la_real janssen_sa_real))
	Rscript code/cor_data_preprocess.R D29IncludeNotMolecConfirmed # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D29IncludeNotMolecConfirmed
	
	Rscript code/cor_data_preprocess.R D29IncludeNotMolecConfirmedstart1 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D29IncludeNotMolecConfirmedstart1
	
	Rscript code/cor_data_preprocess_2.R D29IncludeNotMolecConfirmed # Brian's
	Rscript code/CoR_assay_graphics.R D29IncludeNotMolecConfirmed
	Rscript code/CoR_wrcdf_with_VE_lines.R D29IncludeNotMolecConfirmed
	
	Rscript code/cor_data_preprocess_2.R D29IncludeNotMolecConfirmedstart1 # Brian's
	Rscript code/CoR_assay_graphics.R D29IncludeNotMolecConfirmedstart1
	Rscript code/CoR_wrcdf_with_VE_lines.R D29IncludeNotMolecConfirmedstart1
	
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19))
	Rscript code/cor_data_preprocess.R D35 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D35
	
	Rscript code/cor_data_preprocess_2.R D35 # Brian's
	Rscript code/CoR_assay_graphics.R D35
	Rscript code/CoR_wrcdf_with_VE_lines.R D35
	
endif

## clean          : clean up data_clean and output folder
clean:
	rm -f data_clean/* && rm -rf output/$(TRIAL)/*

## report        : create a standalone report
report:

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: figs clean all
