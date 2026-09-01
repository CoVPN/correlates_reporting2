export module=cor_coxph

Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = '${module}_${TRIAL}_$(date +%Y%m%d).pdf',   config_file = '_bookdown_${module}.yml', output_format = bookdown::pdf_document2(toc_depth=5), quiet=TRUE)"	

