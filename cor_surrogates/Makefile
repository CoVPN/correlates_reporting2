## all           : clean saved results and get results to populate report
all: clean run_cvsl_varsets createRDAfiles_fromSLobjects tables_figures 

## clean         : delete all SL results
clean:
	rm -f output/*
	rm -f figs/*

run_cvsl_varsets: 
	Rscript code/run_cvsl_varsets.R D210 1
	Rscript code/run_cvsl_varsets.R D210 2
	Rscript code/run_cvsl_varsets.R D210 3
	Rscript code/run_cvsl_varsets.R D210 4
	Rscript code/run_cvsl_varsets.R D210 5
	Rscript code/run_cvsl_varsets.R D210 6
	Rscript code/run_cvsl_varsets.R D210 7
	Rscript code/run_cvsl_varsets.R D210 8
	Rscript code/run_cvsl_varsets.R D210 9
	Rscript code/run_cvsl_varsets.R D210 10
	Rscript code/run_cvsl_varsets.R D210 11
	Rscript code/run_cvsl_varsets.R D210 12
	Rscript code/run_cvsl_varsets.R D210 13
	Rscript code/run_cvsl_varsets.R D210 14
	Rscript code/run_cvsl_varsets.R D210 15


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
