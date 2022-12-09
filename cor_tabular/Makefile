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
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_EUA janssen_na_EUA janssen_la_EUA janssen_sa_EUA ))
    Rscript code/make_table_all.R D29IncludeNotMolecConfirmed
    Rscript code/make_table_all.R D29IncludeNotMolecConfirmedstart1
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
