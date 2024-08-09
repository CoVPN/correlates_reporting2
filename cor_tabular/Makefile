## all              : tables_output
all: tables_output

## tables_output    : make summary tables and save data to output/
tables_output: 
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/make_table_all.R D29
	Rscript code/make_table_all.R D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_na_mock janssen_la_mock janssen_sa_mock))
	Rscript code/make_table_all.R D29
	Rscript code/make_table_all.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_EUA janssen_na_EUA janssen_la_EUA janssen_sa_EUA))
	Rscript code/make_table_all.R D29IncludeNotMolecConfirmed
	Rscript code/make_table_all.R D29IncludeNotMolecConfirmedstart1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_partA  janssen_la_partA))
	Rscript code/make_table_all.R D29IncludeNotMolecConfirmed
	Rscript code/make_table_all.R D29ModerateIncludeNotMolecConfirmed
	Rscript code/make_table_all.R D29SevereIncludeNotMolecConfirmed
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_partAsenior janssen_pooled_partAnonsenior janssen_na_partA janssen_na_partAsenior janssen_na_partAnonsenior janssen_la_partAsenior janssen_la_partAnonsenior janssen_sa_partA  janssen_sa_partAnonsenior)) # janssen_sa_partAsenior has too few cases
	Rscript code/make_table_all.R D29IncludeNotMolecConfirmed
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705 hvtn705V1V2 hvtn705second hvtn705secondprimary ))
	Rscript code/make_table_all.R D210
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19))
	Rscript code/make_table_all.R D35
else ifeq ($(TRIAL),$(filter $(TRIAL), azd1222 azd1222_bAb))
	Rscript code/make_table_all.R D29
	Rscript code/make_table_all.R D57
else ifeq ($(TRIAL),$(filter $(TRIAL), vat08m_naive vat08m_nonnaive))
	Rscript code/make_table_all.R D22
	Rscript code/make_table_all.R D43
	Rscript code/make_table_all.R D22omi
	Rscript code/make_table_all.R D43omi
else ifeq ($(TRIAL),$(filter $(TRIAL), profiscov))
	Rscript code/make_table_all.R D43
	Rscript code/make_table_all.R D91
else ifeq ($(TRIAL),$(filter $(TRIAL), profiscov_lvmn))
	Rscript code/make_table_all.R D43	
else ifeq ($(TRIAL),$(filter $(TRIAL), id27hpv))
	Rscript code/make_table_all_hpv.R M18	
	Rscript code/make_table_all_hpv.R M18sus
else ifeq ($(TRIAL),$(filter $(TRIAL), id27hpvnAb))
	Rscript code/make_table_all_hpv.R M18nAb	
	Rscript code/make_table_all_hpv.R M18susnAb	
else ifeq ($(TRIAL),$(filter $(TRIAL), covail))
	Rscript code/make_table_all_covail.R D15to91
	Rscript code/make_table_all_covail.R D92to181
	Rscript code/make_table_all_covail.R D15to181
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_partA_VL))
	Rscript code/make_table_all_VL.R D29VLvariant
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19_stage2))
	Rscript code/make_table_all_stg2.R D35prevent19_stage2_delta 
	Rscript code/make_table_all_stg2.R D35prevent19_stage2_severe
else ifeq ($(TRIAL),$(filter $(TRIAL), azd1222_stage2))
	Rscript code/make_table_all_stg2_azd.R D57azd1222_stage2_delta_nAb 
	Rscript code/make_table_all_stg2_azd.R D57azd1222_stage2_delta_bAb 
	Rscript code/make_table_all_stg2_azd.R D57azd1222_stage2_severe_nAb
	Rscript code/make_table_all_stg2_azd.R D57azd1222_stage2_severe_bAb
else ifeq ($(TRIAL),$(filter $(TRIAL), nvx_uk302))
	Rscript code/make_table_all.R D35nvx_uk302
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19nvx))
	Rscript code/make_table_all.R D35prevent19nvx
else ifeq ($(TRIAL),$(filter $(TRIAL), vat08_combined))
	Rscript code/make_table_all_sanofi.R D22vat08_combined_M6_bAb
	Rscript code/make_table_all_sanofi.R D43vat08_combined_M6_bAb
	Rscript code/make_table_all_sanofi.R D22vat08_combined_M6_nAb
	Rscript code/make_table_all_sanofi.R D43vat08_combined_M6_nAb
endif

## clean            : delete contents of data_clean/ and output/ directories
clean:
	rm -f data_clean/*.Rdata
	rm -f output/*.Rdata
	rm -f output/*.pdf

## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"

## report        : create a standalone report
report:

## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean style tables_output
