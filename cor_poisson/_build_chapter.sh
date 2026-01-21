export module=cor_poisson
Rscript -e "setwd('..'); bookdown::render_book(input = 'index_cor.Rmd', output_file = 'correlates_${module}_${TRIAL}_$(date +%Y%m%d).pdf',   config_file = '_bookdown_${module}.yml', output_format = bookdown::pdf_document2(toc_depth=5), quiet=TRUE)"


