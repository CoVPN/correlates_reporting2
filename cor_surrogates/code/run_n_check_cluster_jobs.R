# Sys.setenv(TRIAL = "hvtn705second")
# Sys.setenv(TRIAL = "moderna_real")
# Sys.setenv(TRIAL = "janssen_pooled_partA")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("summarise", "dplyr")

# Get varset_names.csv 
varsets <- read.csv(paste0("output/", Sys.getenv("TRIAL"), "/varset_names.csv"), stringsAsFactors = FALSE) 

# Run all sbatch jobs through R!
system(paste0("export TRIAL=", Sys.getenv("TRIAL")))
for(varset_number in 1:nrow(varsets)){
  system(paste("sbatch code/submit_cluster_job.sh", varset_number))
}

# Run all sbatch jobs through R!
system(paste0("export TRIAL=", Sys.getenv("TRIAL")))
for(varset_number in 1:nrow(varsets)){
  system(paste("sbatch code/submit_VIMPcluster_job.sh", varset_number))
}

################################################################################
# Make checks on which jobs failed due to random renv or python issues, and resubmit them! 
# Keep making this check, until all jobs are run successfully!
cron_add(command = cron_rscript("code/cron_job.R"), frequency = '/3 * * * *', at = Sys.time() , id = 'test_linux_run', description = "testing linux scheduler")

################################################################################
# If re-submitting job, then check output files that have been newly created!
# If all jobs have finished running (i.e. there are no jobs running as listed by "squeue -al -u UID")
# then run the varset_no that have OutputFileStatus as "old"!
outputFiles <- file.info(paste0("output/", list.files(path = "output", pattern = ".rds"))) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "filename") %>%
  filter(str_detect(filename, "CVSLfits_vacc_EventIndPrimaryD57_")) %>%
  separate(col = filename, into = c("left", "right"), sep = "_EventIndPrimaryD57_") %>%
  rename(filename = right) %>%
  mutate(filename = gsub(".rds", "\\1", filename)) %>%
  select(filename, mtime, size, uname) 

approx.DateTime.of.job.start <- as.POSIXct(readline(prompt="Enter approximate time of job start: \n FORMAT should be YYYY-MM-DD HH:MM:SS"))
#Enter: 2022-02-23 21:00:00

varsets %>% left_join(outputFiles, by = c("varset_names" = "filename")) %>%
  arrange(mtime) %>%
  mutate(OutputFileStatus = ifelse(mtime > approx.DateTime.of.job.start, "new", "old")) %>%
  arrange(varset_no)

################################################################################
# Get list of jobs that are running on cluster nodes in dataframe
jobsRunning <- system("squeue -al -u bborate", intern = TRUE)[-1] 
jobsRunning_df <- data.frame(JOBID = character(),
                             PARTITION = character(),
                             NAME = character(),
                             USER = character(),
                             STATE = character(),
                             TIME = character(),
                             TIME_LIMI = character(),
                             NODES = character(),
                             NODELIST_REASON = character())

parse_string_to_columns <- function(dat, vec){
  vec <- strsplit(vec, "\\s+")[[1]]
  vec <- vec[vec != ""]
  if(grepl("_", vec[1], fixed = TRUE)){
    vecdat <- data.frame(t(vec)) 
    names(vecdat) <- names(dat)
    dat <- bind_rows(dat, vecdat)
  }
  return(dat)
} 

for(i in 2:length(jobsRunning)){
  jobsRunning_df <- parse_string_to_columns(jobsRunning_df, jobsRunning[i])
}

jobsRunning_df <- jobsRunning_df %>% 
  separate(col = JOBID, into = c("JOBID", "JOBID_varset"), sep = "_") %>%
  mutate(JOBID_varset = as.numeric(as.character(JOBID_varset))) %>%
  arrange(JOBID_varset)

jobsRunning_df
