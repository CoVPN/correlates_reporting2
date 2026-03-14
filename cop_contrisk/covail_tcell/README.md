# COVAIL T cell Correlates Controlled Risk Curves and Other Figures 

## Reproducibility

This project uses a project-level renv.lock. 


```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/covail_tcell.zip
unzip covail_tcell.zip
cd correlates_reporting2-covail_tcell


# b) restore R package dependencies
# Assume that we have R 4.2.2 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # install.packages("renv") # if needed
    renv::restore()


# c) edit all R scripts so that references to CD4+T_dt-summary_COMPASS.csv and covail_data_processed_20250612.csv point to a local copy.


# d) generate report pdf
Rscript Figure_3_Panel_A.R
Rscript Figure_3_Panel_B.R
Rscript Figure_4_FS_markers.R
Rscript Figure_4_primary_T_cell_markers.R
Rscript Figure_4_secondary_markers.R
Rscript Figure_5.R
Rscript Figure_6.R
```

