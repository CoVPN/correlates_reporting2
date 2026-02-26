# Coding for Reproducibility

1.	**Portability**. 
    -   Avoid using absolute paths because they will break if someone else downloads the code and tries to run it. 
    -	If the code sources a utility functions file, make sure that file is part of the code base, either at the project level or at the module level.
  	-	If the code writes results to subdirectories, make sure the subdirectories exist through code. For example, the following command creates a folder named 'output' in the current directory (nothing happens if the folder already exists). The recursive option is needed if the path is more than one level deep.
    ```
    dir.create("output", showWarnings = FALSE, recursive=TRUE)
    ```    
2.	Use the package **renv** to manage R system and package versions. See the section below for details.
3.	Every project- or module-level README.md should have a **Reproducibility section**, which ideally should include 4-part bash commands such as:
    ```{bash}
    # a) obtaining the code
    wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/sanofi_stage2_R4.0.zip
    unzip sanofi_stage2_R4.0.zip
    cd correlates_reporting2-sanofi_stage2_R4.0
    
    # b) restore R package dependencies
    R
        Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
        renv::restore()
    
    # c) edit config.yml so that the data_cleaned field uder vat08_combined points to a local copy of vat08_combined_data_processed_20250417.csv
    
    # d) generate report pdf
    export TRIAL=vat08_combined
    export stage=2    
    cd cor_coxph
    make 
    ```
4.  Including an **appendix** at the end of the report Rmd to show the code commit and file name (ANALYSIS_READY_DATA_FILE_NAME should be replaced with your file name):
    ````
    # Appendix
    ```{r, echo = FALSE, message = FALSE, warning = FALSE, results='asis'}
    commit_hash <- system("git rev-parse HEAD", intern = TRUE)
    git_url <- sub("\\.git$", paste0("/commits/", commit_hash), system("git remote get-url origin", intern = TRUE))
    cat("This report was built with code from [", sprintf("**[this commit](%s)**", git_url), "] and data from [", sprintf("**[this file](%s)**", ANALYSIS_READY_DATA_FILE_NAME), "].", sep="")
    ```    
    ````
5.  Include a **date string** in the report file name to show the date on which report was produced (e.g., filename ending ‘20251023’ for October 23, 2025).  This can be done via, e.g.
    ```
    Rscript -e "rmarkdown::render('cor_threshold_barda_mock.Rmd', output_file='cor_threshold_barda_mock_$(date +%Y%m%d).pdf')"
    ```
6.	**Running** reports. There are several options:
    -	If there are only a few Rmd files to be rendered, an Rscript call is sufficient. This often happens at the project level, e.g., https://github.com/CoVPN/correlates_reporting2/blob/master/cor_threshold/sanofi_stage2/README.md. 
    -	Use a Makefile or a bash script to run analyses and generate reports. This often happens at the module level, e.g., https://github.com/CoVPN/correlates_reporting2/blob/master/cor_coxph/Makefile
    -	If high performance cluster/slurm is used and there are dependencies between steps, multiple scripts may be needed, e.g., https://github.com/CoVPN/correlates_reporting2/blob/master/cop_exposureproximal/ensemble_severe_correlates/README.md


## Using renv

renv can be used at one of three levels: repo-level, module-level, and project level. We recommend the module-level.

See https://hvtn-sdmc.github.io/reproducibility-core/using_renv.html#/title-slide for more details on using renv.






