# Sanofi Stage 2 CoP manuscript NP threshold code 

This folder contains code for the non-parametric threshold analysis in the Sanofi Stage 2 correlates of protection manuscript. 

The report is the HTML file. 

The code to create the report is in the Rmd file. 

To run the code, be sure that you have all the package dependencies including the working version of the `npthreshold` packages: 

```{r}
devtools::install_github("jpspeng/npthreshold")
```
Also, be sure to change the directory of the data file in this line: 

```{r}
df <- read.csv("vat08_combined_data_processed_20250321.csv")
```

For any questions, please email jpspeng@uw.edu. 