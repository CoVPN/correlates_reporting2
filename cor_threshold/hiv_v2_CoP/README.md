# HIV V2 Antibody CoP manuscript NP threshold code 

This folder contains code for the non-parametric threshold analysis in the HIV V2 Antibody correlates of protection manuscript. 

The report is the HTML file. 

The code to create the report is in the Rmd file. 

To run the code, be sure that you have all the package dependencies including the working version of the `npthreshold` packages: 

```{r}
devtools::install_github("jpspeng/npthreshold")
```
Also, be sure to change the directory of the data file in this line: 

```{r}
df_hvtn702 <- read.csv("data/controlledVEdata702.csv")
```

For any questions, please email jpspeng@uw.edu. 