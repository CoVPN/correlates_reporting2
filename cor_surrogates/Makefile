## all           : clean saved results and get results to populate report
all: run_cvsl_varsets get_vimp createRDAfiles_fromSLobjects tables_figures

## clean         : delete all SL results
clean:
	rm -f output/*
	rm -f figs/*

run_cvsl_varsets:
	Rscript code/run_cvsl_varsets.R D210 36
	Rscript code/run_cvsl_varsets.R D210 37
	Rscript code/run_cvsl_varsets.R D210 38
	Rscript code/run_cvsl_varsets.R D210 39
	Rscript code/run_cvsl_varsets.R D210 40
	Rscript code/run_cvsl_varsets.R D210 41
	Rscript code/run_cvsl_varsets.R D210 42
	Rscript code/run_cvsl_varsets.R D210 43
	Rscript code/run_cvsl_varsets.R D210 44
	Rscript code/run_cvsl_varsets.R D210 45
	Rscript code/run_cvsl_varsets.R D210 46
	Rscript code/run_cvsl_varsets.R D210 47

get_vimp:
	Rscript code/get_vimp.R

## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects: run_cvsl_varsets
	Rscript code/createRDAfiles_fromSLobjects.R

## tables_figures   : create tables and figures
tables_figures: createRDAfiles_fromSLobjects
	Rscript code/tables_figures.R

## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<
