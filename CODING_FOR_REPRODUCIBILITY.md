# Coding for Reproducibility

1.	**Portability**. 
    -   Avoid using absolute paths because they will break if someone else downloads the code and tries to run it. 
    -	If the code sources a utility functions file, make sure that file is part of the code base, either at the project level or at the module level.
  	-	If the code writes results to subdirectories, make sure the subdirectories exist through code. For example, the following command creates a folder named 'output' in the current directory (nothing happens if the folder already exists). The recursive option is needed if the path is more than one level deep.
    ```
    dir.create("output", showWarnings = FALSE, recursive=TRUE)
    ```    
2.	Use the package **renv** to manage R system and package versions. See the section below for details.
3.	**Running** reports. There are several options:
    -	If there are only a few Rmd files to be rendered, an Rscript call is sufficient. This often happens at the project level, e.g., https://github.com/CoVPN/correlates_reporting2/blob/master/cor_threshold/sanofi_stage2/README.md. 
    -	Use a Makefile or a bash script to run analyses and generate reports. This often happens at the module level, e.g., https://github.com/CoVPN/correlates_reporting2/blob/master/cor_coxph/Makefile
    -	If high performance cluster/slurm is used and there are dependencies between steps, multiple scripts may be needed, e.g., https://github.com/CoVPN/correlates_reporting2/blob/master/cop_exposureproximal/ensemble_severe_correlates/README.md
4.	Every project- or module-level README.md should have a **Reproducibility section**. 
    -	Project-level example: https://github.com/CoVPN/correlates_reporting2/blob/master/cor_threshold/sanofi_stage2/README.md
    -	Module-level example: https://github.com/CoVPN/correlates_reporting2/blob/master/cor_coxph/README.md
5.  Including an **appendix** at the end of the report Rmd to show the code commit and file name (ANALYSIS_READY_DATA_FILE_NAME should be replaced with your file name):
    ````
    # Appendix
    ```{r, echo = FALSE, message = FALSE, warning = FALSE, results='asis'}
    commit_hash <- system("git rev-parse HEAD", intern = TRUE)
    git_url <- sub("\\.git$", paste0("/commits/", commit_hash), system("git remote get-url origin", intern = TRUE))
    cat("This report was built with ", sprintf("**[code](%s)**", git_url), " and ", sprintf("**[data](%s)**", ANALYSIS_READY_DATA_FILE_NAME), ".", sep="")
    ```    
    ````
6.  Include a **date string** in the report file name to show the date on which report was produced (e.g., filename ending ‘20251023’ for October 23, 2025).  This can be done via, e.g.
    ```
    Rscript -e "rmarkdown::render('cor_threshold_barda_mock.Rmd', output_file='cor_threshold_barda_mock_$(date +%Y%m%d).pdf')"
    ```

## Using renv for Reproducibility

renv can be used at one of three levels: repo-level, module-level, and project level. The following instructions show how to use renv when starting a new project. Don't follow these instructions if you are trying to reproduce results for an existing project. Instead, follow the instructions in the project (e.g., https://github.com/CoVPN/correlates_reporting2/blob/master/cor_threshold/sanofi_stage2/README.md) or the module that contains the project (e.g, https://github.com/CoVPN/correlates_reporting2/blob/master/cor_coxph/README.md). 


### Setting up renv at the project-level for a new project

This is the most reproducible way of using renv because each project/manuscript has its own renv.lock.

To start, make sure there are two files in the project folder:
- .gitignore  This file tells git which files/folders to ignore. It should contain the following lines:
```
.html
.pdf
.Rhistory
.RData
.Rproj.user
renv/
```
- README.md  This file should have a section titled Reproducibility, which details how to reproduce the reports.

Open a new R console in the project folder in a terminal and run the following commands. Note that we use renv 1.1.5 (or 0.13.2 in earlier code)
```{r}
# if the following does not work, it is because 1.1.5 is the current release. Then try install.packages("renv")
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz",
  repos = NULL,
  type = "source"
)

packageVersion("renv")
```

Run the following R command at the project level to initialize:

```{r}
renv::init()
```

The call to init creates three files/folders:

- renv.lock  This is a manifest file that contains the R packages and their versions used in the project. If the project is not empty, init() looks for package dependencies in the files, installs them, and writes the dependencies to renv.lock. If there are many package dependencies, this step can take a while.
- renv/  This folder contains the local library of packages used in the project. In fact, to reduce redundancy, it contains links to the packages installed in your global R library.
- .Rprofile  This file is sourced when the project is opened and activates the renv environment It is helpful to edit the file and add a line to show the current directory so that its content is as follows:
  ```{r}
  source("renv/activate.R")
  print(getwd())
  ```

At this point, you can close the R console.
```{r}
q()
```

Now if you start a new R console in the project folder, you should see a message from renv, indicating that the renv environment is activated.

For an example of a project using project-level renv, check out cor_threshold/sanofi_stage2 (https://github.com/CoVPN/correlates_reporting2/blob/master/cor_threshold/sanofi_stage2/README.md).

### Setting up renv at the module-level for a new module

This is not currently done. There seems to be little advantages of using renv at the module level.

### Setting up renv at the repo-level

This has already been done for the repo. For an example of using the repo-level renv, check out cor_coxph (https://github.com/CoVPN/correlates_reporting2/blob/master/cor_coxph/README.md).




### General tips for using renv

Use renv::install and not remotes or DevTools to install or update pacakges from CRAN, github, or a url, e.g.,
- renv::install("github_id/package_name")
- renv::install("survival")
- renv::install("https://cran.r-project.org/src/contrib/Archive/kyotil/kyotil_2024.7-31.tar.gz")

renv::snapshot() updates renv.lock after new or updated packages are installed.

If you encounter warnings like: curl: (22) The requested URL returned error: 403 rate limit exceeded. The solution is to first authenticate with your personal github access token by running the following in R:
```{r}
Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxx")
```

If something inexplicable goes wrong, check the $HOME directory to see if there are .Rprofile and renv/. If yes, delete them and try again. 

