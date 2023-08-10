#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(stringr)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
if (!dir.exists(paste0(save.results.to, "/demographics")))  dir.create(paste0(save.results.to, "/demographics"))

# Define age cutoff based on trial
age.cutoff <- ifelse(study_name %in% c("ENSEMBLE", "MockENSEMBLE"), 60, 65)

trt.labels <- c("Placebo", "Vaccine")
bstatus.labels <- c("Baseline Neg", "Baseline Pos")
bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")


all_assays <- c("bindSpike", "bindSpike_B.1.1.7", "bindSpike_B.1.351", "bindSpike_P.1", 
                "bindRBD", "bindRBD_B.1.1.7", "bindRBD_B.1.351", "bindRBD_P.1", 
                "bindN",
                "pseudoneutid50", "pseudoneutid80", "pseudoneutid50sa", "pseudoneutid50la",
                "liveneutmn50")
bAb_assays <- c("bindSpike", "bindSpike_B.1.1.7", "bindSpike_B.1.351", "bindSpike_P.1", 
                "bindRBD", "bindRBD_B.1.1.7", "bindRBD_B.1.351", "bindRBD_P.1", 
                "bindN")
nAb_assays <- c("pseudoneutid50", "pseudoneutid80", "pseudoneutid50sa", "pseudoneutid50la")
live_assays <- c("liveneutmn50")
times <- c("B", paste0("Day", config$timepoints), paste0("Delta", config$timepoints, "overB"))

# Depends on the Incoming data
if(include_bindN && grepl("bind", assays) && !grepl("bindN", assays) && !grepl("janssen_.+partA.*", attr(config,"config"))){
  assay_immuno <- all_assays[all_assays %in% c(assays, "bindN")]
  labels.assays.all <- c("Binding Antibody to N", labels.assays)
  names(labels.assays.all)[1] <- "bindN"
  labels.assays <- labels.assays.all[assay_immuno]
} else {
  assay_immuno <- assays
}

labels.assays.short <- c("Anti Spike IgG (BAU/ml)", 
                         "Anti Spike B.1.1.7 IgG (BAU/ml)", 
                         "Anti Spike B.1.351 IgG (BAU/ml)", 
                         "Anti Spike P.1 IgG (BAU/ml)", 
                         "Anti RBD IgG (BAU/ml)", 
                         "Anti RBD B.1.1.7 IgG (BAU/ml)", 
                         "Anti RBD B.1.351 IgG (BAU/ml)", 
                         "Anti RBD P.1 IgG (BAU/ml)", 
                         "Anti N IgG (BAU/ml)",
                         "Pseudovirus-nAb ID50 (IU50/ml)", 
                         "Pseudovirus-nAb ID80 (IU50/ml)", 
                         "Pseudovirus-nAb ID50 (SA) (IU50/ml)", 
                         "Pseudovirus-nAb ID50 (LA) (IU50/ml)", 
                         "Live Virus-mnAb ID50 (IU50/ml)",
                         "Phagocytic Score")
names(labels.assays.short) <- c("bindSpike", "bindSpike_B.1.1.7", "bindSpike_B.1.351", "bindSpike_P.1", 
                                "bindRBD", "bindRBD_B.1.1.7", "bindRBD_B.1.351", "bindRBD_P.1", 
                                "bindN",
                                "pseudoneutid50", "pseudoneutid80", 
                                "pseudoneutid50sa", 
                                "pseudoneutid50la", 
                                "liveneutmn50",
                                "ADCP")




# axis labeling
labels.axis <- outer(
  rep("", length(times)),
  labels.assays.short[assay_immuno],
  "%.%"
)
labels.axis <- as.data.frame(labels.axis)
rownames(labels.axis) <- times


## redefine the labels used in the immuno_graphical report

labels.title <- outer(
  labels.assays[assay_immuno],
  ": " %.%
    labels.time[names(labels.time) %in% times],
  paste0
) 


labels.title <- as.data.frame(labels.title)
colnames(labels.title) <- times
# NOTE: hacky solution to deal with changes in the number of markers
rownames(labels.title) <- assay_immuno
labels.title <- as.data.frame(t(labels.title))


labels.title2 <- apply(labels.title, c(1, 2), function(st) {
  str_replace(st, ":", "\n")
})

# creating short and long labels
labels.assays.short <- labels.axis[1, ]
labels.assays.long <- labels.title
