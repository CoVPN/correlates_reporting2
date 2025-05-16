## all              : tables_output
all: tables_output

## tables_output    : make summary tables and save data to output/
tables_output: 

ifeq ($(TRIAL),$(filter $(TRIAL), prevent19_stage2))
	Rscript code/make_table_all_stg2.R D35prevent19_stage2_delta 
	Rscript code/make_table_all_stg2.R D35prevent19_stage2_severe 
else ifeq ($(TRIAL),$(filter $(TRIAL), nextgen_mock))
	Rscript code/make_table_all_stg2.R D31toM12_nextgen_mock 
	Rscript code/make_table_all_stg2.R D31toM12_nextgen_mock_tcell
else code/make_table_all.R clean_data
	Rscript code/make_table_all.R
endif

## clean_data       : clean up and reformat input data, and save to data_clean/
clean_data: code/make_clean_data.R set_params
	Rscript code/make_clean_data.R

## set_params       : set up parameters (e.g., lloq) and modify included tables
set_params: code/make_parameters.R
	Rscript code/make_parameters.R

## clean            : delete contents of data_clean/ and output/ directories
clean:
	rm -f data_clean/*.Rdata
	rm -f output/*.Rdata
	rm -f output/*.pdf

## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"

## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean style tables_output
