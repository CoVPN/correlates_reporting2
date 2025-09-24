# Sanofi Stage 2 CoP manuscript NP threshold code 

Set up note: 

- Assume that we work with R 4.4.2.
- Assume that we have followed the steps in using_renv_for_reproducibility.Rmd at the root level. 
- To install package dependencies, run the following commands once:
  ```{r}
  renv::install("isotone")
  renv::install("tlverse/tmle3")
  renv::install("jpspeng/npthreshold")
  renv::snapshot()
  ```
- Later, when you start R in this project folder, renv will activate and provide the right package versions for use.



To reproduce the report, follow these steps:

- Change the directory of the data file in this line in sanofi_stage2_npthreshold_report_2025June25.Rmd: 
  ```{r}
  df <- read.csv("vat08_combined_data_processed_20250321.csv")
  ```
  It is better to replace this line with the following if working on the SCHARP file system:
  ```{r}
  config.reporting <- config::get(config = "vat08_combined", file="../../config.yml") 
  dat<-read.csv(config.reporting$data_cleaned)
  ```
- To render the R markdown files, run the following commands in a bash shell:
  ```{bash}
  Rscript -e "rmarkdown::render('sanofi_stage2_npthreshold_report_2025June25.Rmd')"
  ```


For any questions, please email jpspeng@uw.edu. 