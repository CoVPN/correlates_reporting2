# Using renv for Reproducibility

renv can be used at one of three levels: repo-level, module-level, and project level. For example, cor_threshold is a module-level directory and under that, sanofi_stage2 is a project-level directory. 

## Setting up renv at the project-level for a new project

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

Open a new R console in the project folder in a terminal and run the following commands. Note that we use renv 0.13.2, which uses renv/activate.R, instead of newer versions because of some errors with the newer versions. (If in a slurm env, load an appropriate R module and a CMmake module. The latter is needed to install some packages, e.g., nloptr, lme4.
)
```{r}
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/renv/renv_0.13.2.tar.gz",
  repos = NULL,
  type = "source"
)

packageVersion("renv")  # should show ‘0.13.2’

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

## Setting up renv at the module-level for a new module

This is not currently done. There seems to be little advantages of using renv at the module level.

## Setting up renv at the repo-level

This has already been done for the repo. For an example of using the repo-level renv, check out cor_coxph (https://github.com/CoVPN/correlates_reporting2/blob/master/cor_coxph/README.md).




## General tips for using renv

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

