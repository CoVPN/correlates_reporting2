## all            : clear current results and make figures
all: clean figs

## figs           : clean data and make figures
figs: 
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real azd1222 azd1222_bAb))
	Rscript code/cor_data_preprocess.R D29D57 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D29D57

	Rscript code/cor_data_preprocess_2.R D29 # Brian's
	Rscript code/CoR_assay_graphics.R D29
	Rscript code/CoR_wrcdf_with_VE_lines.R D29

	Rscript code/cor_data_preprocess_2.R D57 # Brian's
	Rscript code/CoR_assay_graphics.R D57
	Rscript code/CoR_wrcdf_with_VE_lines.R D57

else ifeq ($(TRIAL),$(filter $(TRIAL), vat08m_naive vat08m_nonnaive))
	Rscript code/cor_data_preprocess.R D22D43 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D22D43
	Rscript code/cor_data_preprocess.R D22D43omi # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D22D43omi

	Rscript code/cor_data_preprocess_2.R D22 # Brian's
	Rscript code/CoR_assay_graphics.R D22
	Rscript code/CoR_wrcdf_with_VE_lines.R D22
	Rscript code/cor_data_preprocess_2.R D43 # Brian's
	Rscript code/CoR_assay_graphics.R D43
	Rscript code/CoR_wrcdf_with_VE_lines.R D43

	Rscript code/cor_data_preprocess_2.R D22omi # Brian's
	Rscript code/CoR_assay_graphics.R D22omi
	Rscript code/CoR_wrcdf_with_VE_lines.R D22omi
	Rscript code/cor_data_preprocess_2.R D43omi # Brian's
	Rscript code/CoR_assay_graphics.R D43omi
	Rscript code/CoR_wrcdf_with_VE_lines.R D43omi

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

else ifeq ($(TRIAL), $(filter $(TRIAL), janssen_pooled_EUA janssen_na_EUA janssen_la_EUA janssen_sa_EUA))
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
	
else ifeq ($(TRIAL), $(filter $(TRIAL), janssen_pooled_partA janssen_na_partA janssen_la_partA janssen_sa_partA janssen_pooled_partAsenior janssen_na_partAsenior janssen_la_partAsenior  janssen_pooled_partAnonsenior janssen_na_partAnonsenior janssen_la_partAnonsenior janssen_sa_partAnonsenior))
	Rscript code/cor_data_preprocess.R D29IncludeNotMolecConfirmed # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D29IncludeNotMolecConfirmed

	Rscript code/cor_data_preprocess_2.R D29IncludeNotMolecConfirmed # Brian's
	Rscript code/CoR_assay_graphics.R D29IncludeNotMolecConfirmed
	Rscript code/CoR_wrcdf_with_VE_lines.R D29IncludeNotMolecConfirmed
	
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19))
	Rscript code/cor_data_preprocess.R D35 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D35

	Rscript code/cor_data_preprocess_2.R D35 # Brian's
	Rscript code/CoR_assay_graphics.R D35
	Rscript code/CoR_wrcdf_with_VE_lines.R D35

else ifeq ($(TRIAL),$(filter $(TRIAL), profiscov))
	Rscript code/cor_data_preprocess.R D43D91 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D43D91

	Rscript code/cor_data_preprocess_2.R D43 # Brian's
	Rscript code/CoR_assay_graphics.R D43
	Rscript code/CoR_wrcdf_with_VE_lines.R D43

	Rscript code/cor_data_preprocess_2.R D91 # Brian's
	Rscript code/CoR_assay_graphics.R D91
	Rscript code/CoR_wrcdf_with_VE_lines.R D91

else ifeq ($(TRIAL),$(filter $(TRIAL), profiscov_lvmn))
	Rscript code/cor_data_preprocess.R D43 # Yiwen's
	Rscript code/cor_graphics_violin_scatter.R D43

	Rscript code/cor_data_preprocess_2.R D43 # Brian's
	Rscript code/CoR_assay_graphics.R D43
	Rscript code/CoR_wrcdf_with_VE_lines.R D43
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
