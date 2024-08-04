export module=cor_coxph
cd ..
if [$TRIAL = "vat08_combined"]; then

  if [ $stage = 1 ]; then
    echo "stage 1"
  	Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_${module}_vat08_stage1.pdf', config_file = '_bookdown_${module}.yml', output_format = bookdown::pdf_document2(toc_depth=5), quiet=TRUE)"	
  else
    echo "stage 2"
  	Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_${module}_vat08_stage2.pdf', config_file = '_bookdown_${module}.yml', output_format = bookdown::pdf_document2(toc_depth=5), quiet=TRUE)"	
  fi
  cd $module

else 

	Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_${module}_${TRIAL}.pdf',   config_file = '_bookdown_${module}.yml', output_format = bookdown::pdf_document2(toc_depth=5), quiet=TRUE)"	

fi