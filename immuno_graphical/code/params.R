#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(stringr)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
if (!dir.exists(paste0(save.results.to, "/demographics")))  dir.create(paste0(save.results.to, "/demographics"))
config <- config::get(config = Sys.getenv("TRIAL"))

# Define age cutoff based on trial
age.cutoff <- ifelse(study_name %in% c("ENSEMBLE", "MockENSEMBLE", "VAT08"), 60, 65)

trt.labels <- c("Placebo", "Vaccine")
if (study_name !="VAT08"){
  bstatus.labels <- c("Baseline Neg", "Baseline Pos")
  bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")
  bstatus.labels.3 <- c("baseline negative", "baseline positive")
} else {
  bstatus.labels <-  c("Naive", "Non-naive")
  bstatus.labels.2 <- bstatus.labels.3 <- c("naive", "non-naive")
}

bAb_assays <- assays[grepl("bind", assays)]
nAb_assays <- assays[grepl("pseudoneut", assays)]
live_assays <- assays[grepl("liveneut", assays)]
names(labels.assays.short) <- assays

times <- c("B", paste0("Day", config$timepoints), paste0("Delta", config$timepoints, "overB"))

# Depends on the incoming data
if(include_bindN && grepl("bind", assays) && !grepl("bindN", assays) && !grepl("janssen_.+partA.*", attr(config,"config"))){
  assay_immuno <- assays[assays %in% c(assays, "bindN")]
  labels.assays.all <- c("Binding Antibody to N", labels.assays)
  names(labels.assays.all)[1] <- "bindN"
  labels.assays <- labels.assays.all[assay_immuno]
} else {
  assay_immuno <- assays
}

# rerun the definitions in _common.R below based on assay_immuno
# axis labeling
labels.axis <- outer(rep("", length(times)), labels.assays.short[assays], "%.%")
labels.axis <- as.data.frame(labels.axis)
rownames(labels.axis) <- times

# short labels
labels.assays.short <- labels.axis[1, ]

# title labeling
labels.title <- outer(labels.assays[assays], ": " %.% labels.time, paste0)
labels.title <- as.data.frame(labels.title)
colnames(labels.title) <- times
# NOTE: hacky solution to deal with changes in the number of markers
rownames(labels.title)[seq_along(assays)] <- assays
labels.title <- as.data.frame(t(labels.title))




# new set of titles in two lines
labels.title2 <- apply(labels.title, c(1, 2), function(st) {
  str_replace(st, ":", "\n")
})


if (is.null(uloqs)) {
  llox_labels=assay_metadata$llox_label; names(llox_labels)=assays
  lloqs=assay_metadata$lloq; names(lloqs)=assays
  uloqs=assay_metadata$uloq; names(uloqs)=assays
  lods=assay_metadata$lod; names(lods)=assays
  lloxs=ifelse(llox_labels=="lloq", lloqs, lods)
  lloxs=ifelse(llox_labels=="pos", assay_metadata$pos.cutoff, lloxs)
}
