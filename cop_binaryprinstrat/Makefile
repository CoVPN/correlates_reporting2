## all           : clean saved results and get results to populate report
all: clean cop_binaryprinstrat 

## clean         : delete saved results
clean:
    rm -rf output/$(TRIAL)/*

## cop_binaryprinstrat     : create tables and figures
cop_binaryprinstrat:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
    Rscript code/binary_correlate_VE_analysis.R D29 1
    Rscript code/binary_correlate_VE_analysis.R D57 1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_na_mock janssen_la_mock janssen_sa_mock))
    Rscript code/binary_correlate_VE_analysis.R D29 2
    Rscript code/binary_correlate_VE_analysis.R D29start1 2
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_EUA janssen_na_EUA janssen_la_EUA janssen_sa_EUA))
    Rscript code/binary_correlate_VE_analysis.R D29IncludeNotMolecConfirmed 2
    Rscript code/binary_correlate_VE_analysis.R D29IncludeNotMolecConfirmedstart1 2
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705 hvtn705V1V2))
    Rscript code/binary_correlate_VE_analysis.R D210 1
endif


## report        : create a standalone report
report:
    Rscript -e "rmarkdown::render('cop_binaryprinstrat_report.Rmd', output_file='cop_binaryprinstrat_report_$(TRIAL).pdf')"


# type 'make help' to show all make commands
help: Makefile
    @sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cop_binaryprinstrat
