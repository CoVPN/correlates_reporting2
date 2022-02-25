# Make checks on which jobs failed due to random renv or python issues, and resubmit them! 
# Keep making this check, until all jobs are run successfully!
slurm.out.files <- list.files(path = ".", pattern = "(slurm).*\\.out$")
for(i in length(slurm.out.files)){
  if(system(paste("gawk 'END {print}'", slurm.out.files[i]), intern = TRUE) == "Execution halted"){
    file.remove(slurm.out.files[i])
    varset_number = gsub(".out", "", strsplit(slurm.out.files[i], "_")[[1]][2])
    system(paste("sbatch code/submit_cluster_job.sh", varset_number))
  }
}
