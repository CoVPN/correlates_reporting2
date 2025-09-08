# Sys.setenv(TRIAL = "hvtn705second")
# Sys.setenv(TRIAL = "moderna_real")
# Sys.setenv(TRIAL = "janssen_pooled_partA")
# Sys.setenv(TRIAL = "covail_tcell")
# COR = "D15to91covail_tcell"
# Dataset-level metadata (export TRIAL=covail_tcell) 
# Objective-level metadata (export COR=D15to91covail_tcell)
# Sys.setenv(TRIAL = "covail_tcell")
# COR = "D15to91covail_tcell"
# COR = "D15to181covail_tcell"

Sys.setenv(TRIAL = "covail_xassays")
# COR = "D15to91covail_xassays"
COR = "D15to181covail_xassays"
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("summarise", "dplyr")

# sbatch_input_vars <- data.frame(
#   COR = c("D15to91covail_xassays", "D15to91covail_xassays", "D15to91covail_xassays", "D15to91covail_xassays", "D15to181covail_xassays", "D15to181covail_xassays", "D15to181covail_xassays", "D15to181covail_xassays", "D15to181covail_xassays", "D15to181covail_xassays"),
#   non_naive = c("FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE"),
#   trt_arms = c("1dosemRNA", "1dosemRNA", "1dose", "1dose", "1dosemRNA", "1dosemRNA", "1dose", "1dose", "1dosemRNA", "1dose"),
#   bRiskFactors_includes_insert.stage.info = c("FALSE", "TRUE", "FALSE", "TRUE", "FALSE", "TRUE", "FALSE", "TRUE", "NA", "NA"),
#   stringsAsFactors = FALSE
# )
# 
# # Example parameter grid
# params_df <- expand.grid(
#   COR = c("D15to91covail_xassays", "D15to181covail_xassays"),
#   non_naive = c("TRUE", "FALSE"),
#   trt_arms = c("1dosemRNA", "1dose"),
#   bRiskFactors_includes_insert_stage_info = c("TRUE", "FALSE"),
#   varset_number = 1:length(varset_names)
# )
# 
# # Add job_id column starting from 1
# params_df$job_id <- 1:nrow(params_df)
# 
# # Save to CSV
# write.csv(params_df, "job_params.csv", row.names = FALSE)

##################################################

# job_table <- data.frame()
# 
# for (i in 1:nrow(sbatch_input_vars)) {
#   # You define this per combination (or load from some setup script)
#   #source("code/cor_surrogates_setup.R", local = TRUE)  # defines varset_names based on sbatch_input_vars[i, ]
#   
#   n_varsets <- length(varset_names)
#   combo <- sbatch_input_vars[i, ]
#   
#   combo_expanded <- combo[rep(1, n_varsets), ]
#   combo_expanded$varset_number <- 1:n_varsets
#   job_table <- rbind(job_table, combo_expanded)
# }
# 
# # Add job_id column
# job_table$job_id <- 1:nrow(job_table)
# 
# # Save for cluster jobs
# write.csv(job_table, "job_table.csv", row.names = FALSE)

# 
# for (i in 1:nrow(sbatch_input_vars)) {
#   # Extract parameters for current row
#   COR <- sbatch_input_vars$COR[i]
#   non_naive <- sbatch_input_vars$non_naive[i]
#   trt_arms <- sbatch_input_vars$trt_arms[i]
#   bRiskFactors_includes_insert.stage.info <- sbatch_input_vars$bRiskFactors_includes_insert.stage.info[i]
#   
#   job_id = 1
#   # common setup for CV super learners and variable importance
#   source(here::here("code", "cor_surrogates_setup.R"))
#   
#   for (varset_number in 1) {
#     # env_vars <- paste0(
#     #   "COR=", COR, ",",
#     #   "non_naive=", non_naive, ",",
#     #   "trt_arms=", trt_arms, ",",
#     #   "bRiskFactors_includes_insert.stage.info=", bRiskFactors_includes_insert.stage.info
#     # )
#     
#     cmd <- paste("sbatch code/submit_cluster_job.sh ", 
#                  varset_number,
#                  shQuote(COR), 
#                  shQuote(non_naive), 
#                  shQuote(trt_arms), 
#                  shQuote(bRiskFactors_includes_insert.stage.info))
#     system(cmd)
#   }
# }
# 

job_id = 1
# common setup for CV super learners and variable importance
source(here::here("code", "cor_surrogates_setup.R"))

# Get varset_names.csv 
varsets <- read.csv(here(file_path, "varset_names.csv"), stringsAsFactors = FALSE) 

# Run all sbatch jobs through R!
system(paste0("export TRIAL=", Sys.getenv("TRIAL")))
for(varset_number in 1:nrow(varsets)){
  system(paste("sbatch code/submit_cluster_job.sh", varset_number))
}

# 
# for(varset_number in c(27, 35, 38, 104)){
#   system(paste("sbatch code/submit_cluster_job.sh", varset_number))
# }

##########################################################
# CHECK IF ALL JOBS ARE COMPLETED SUCCESSFULLY based off expected job ids every 2 mins!
# Set your working directory to where the output files are stored
# file_path = "/home/bborate/forked_correlates_directory/correlates_reporting2/cor_surrogates/output/covail_xassays/nonnaive_1dosemRNA_allcases_briskscore"

expected_ids <- 1:nrow(varsets)
check_interval <- 60  # seconds (i.e., 1 minute)

repeat {
  # List all .rds files
  rds_files <- list.files(file_path, pattern = "\\.rds$", full.names = TRUE)
  
  # Extract job IDs from filenames using regex
  actual_ids <- as.numeric(gsub(".*?_([0-9]+)_.*\\.rds", "\\1", rds_files))
  actual_ids <- actual_ids[!is.na(actual_ids)]  # remove unmatched entries
  
  # Find which jobs are still missing
  missing_ids <- setdiff(expected_ids, actual_ids)
  
  if (length(missing_ids) == 0) {
    message("ALL jobs finished successfully.")
    
    # Source your next script
    source(here::here("code", "createRDAfiles_fromSLobjects.R"))
    source(here::here("code", "tables_figures.R"))
    
    break  # Exit the loop
  } else {
    message(paste("Waiting... MISSING job IDs:", paste(missing_ids, collapse = ", ")))
    
    # Sleep for 3 minutes before re-checking
    Sys.sleep(check_interval)
  }
}
#########################################################
# RUN ONLY THE JOBS THAT RAN UNSUCCESSFULLY !

# for(varset_number in job_ids <- c(19, 28, 70, 78, 91)){
#   system(paste("sbatch code/submit_cluster_job.sh", varset_number))
# }
# For non-naive, 1dosemRNA, missing jobs are:  1   3   4  13  19  20  21  26  28  31  33  35  37  39  40  41  42  43  45  51  60  61  65  70  73  78  79  81  83  85  88  90  91  95  96  98 100 102 103 105 106 107 109
# For non-naive, 1dose, missing jobs are:  1  27  61  

# source(here::here("code", "createRDAfiles_fromSLobjects.R"))
# source(here::here("code", "tables_figures.R"))
#########################################################
# KNIT THE PDF REPORT !
# Run the following in Rstudio
# Make sure the path is correct
file_path = "/home/bborate/forked_correlates_directory/correlates_reporting2/cor_surrogates/output/covail_xassays/naive_1dosemRNA_allcases_bRiskFactors.insert.stage"
Sys.setenv(file_path = file_path)
Sys.setenv(TRIAL = "covail_xassays")

rmarkdown::render("report.Rmd", 
                  output_format = "pdf_document",
                  paste0(basename(file_path), ".pdf"))
# Move the file
file.rename(from = paste0(basename(file_path), ".pdf"), 
            to = paste0(here("cor_surrogates", "covail_xassays_multivariable_SuperLearner_reports"), "/", basename(file_path), ".pdf")) 


# # Run all sbatch jobs through R!
# system(paste0("export TRIAL=", Sys.getenv("TRIAL")))
# for(varset_number in 1:nrow(varsets)){
#   system(paste("sbatch code/submit_VIMPcluster_job.sh", varset_number))
# }

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
