export module=cor_logistic
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_${module}_$TRIAL.pdf', config_file = '_bookdown_$module.yml', output_format = bookdown::pdf_document2(toc_depth=4), quiet=TRUE)"	
cd $module
