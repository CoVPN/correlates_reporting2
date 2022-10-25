#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "********* Please provide a module name, e.g. cor_coxph, as an argument."
    exit
fi


echo $TRIAL
#Rscript -e "bookdown::clean_book(TRUE)"
if [ $1 = "immuno_graphical" ] || [ $1 = "immuno_tabular" ]; then
	Rscript -e "bookdown::render_book(input = 'index_immuno.Rmd', output_file = 'covpn_correlates_$1_$TRIAL.pdf', config_file = '_bookdown_$1.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"
else 
	Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_$1_$TRIAL.pdf', config_file = '_bookdown_$1.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"	
fi

rm -f *.log _report_cor/*.aux _report_cor/*.log _report_cor/*.lof _report_cor/*.tex _report_cor/*.aux _report_cor/*.toc _report_cor/*.lot _report_cor/*.out
