#Sys.setenv(TRIAL = "covail")
Sys.setenv(DESCRIPTIVE= 1)

#-----------------------------------------------
#setwd("~/correlates_reporting2_covail/cor_graphical")
#here::i_am("./code/cor_data_preprocess_covail.R")

# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

source(here::here("..", "_common.R"))

#colnames(assay_metadata) = gsub("X[.]+","", colnames(assay_metadata))
source(here::here("code/cor_graphics_functions_covail.R"))
#-----------------------------------------------

times=c("B","Day15","Day29","Day91","Delta15overB")
#times=c("BD1","BD29","DD1","DeltaBD29overBD1","DeltaDD1overBD1")


#remove pseudoneutid50Duke_BA.2.12.1
assays <- subset(assay_metadata,assay_metadata$panel=="id50")$assay

uloqs=assay_metadata$uloq; names(uloqs)=assays
pos.cutoffs=assay_metadata$pos.cutoff; names(pos.cutoffs)=assays
loqs=assay_metadata$lloq; names(loqs)=assays

#dat$EventIndPrimary=dat$EventIndOmicronBD29
#dat$EventTimePrimary=dat$EventTimeOmicronBD29
#dat.cp <- dat

#dat.mock <- read.csv("/trials/covpn/COVAILcorrelates/analysis/correlates/adata/covail_data_processed_20240105.csv")

#remove pseudoneutid50Duke_BA.2.12.1
dat.cp <- dat.mock%>%
  select(-contains(c('BA.2.12.1')))%>%
  mutate(TrtB_moderna=ifelse(!is.na(TrtB) & arm %in% c(1,2,4:6), TrtB, NA),
         TrtB_pfizer=ifelse(!is.na(TrtB) & arm %in% c(7:9,12,16,17), TrtB, NA))

# assign values above the uloq to the uloq
# for (a in assays){
#   for (t in c("B","Day15","Day29","Day91","Delta15overB")){
#     if (paste0(t, a) %in% colnames(dat.cp)){
#       dat.cp[, paste0(t, a)] = ifelse(dat.cp[, paste0(t, a)] > log10(uloqs[a]), log10(uloqs[a]), dat.cp[, paste0(t, a)])
#     } else {print(paste0(t, a, " doesn't exist"))}
#   }
# }

# recalculate Delta15overB after the censoring above in order to calculate the response rate for cases 
# for (a in assays){
#   if (paste0("Day15", a) %in% colnames(dat.cp) & paste0("B", a) %in% colnames(dat.cp)){
#     dat.cp[, paste0("Delta15overB", a)] = dat.cp[, paste0("Day15", a)] - dat.cp[, paste0("B", a)]
#   } else {print(paste0("B", a, "or Day15", a, " doesn't exist"))}
# }



immuno_vars <- colnames(dat.cp)[grepl("^Bpseudoneutid50|^Day15|^Day29|^Day91|^Delta", colnames(dat.cp))]

# variable selection
dat <- dat.cp %>% select(Ptid, 
                         treatment_actual,TrtmRNA,TrtonedosemRNA,arm,TrtA,TrtB,TrtC,TrtB_moderna,TrtB_pfizer,
                         stage,
                         naive,
                         Perprotocol,Immunemarkerset,ImmunemarkersetD92toD181,
                         COVIDtimeD22toD91,COVIDIndD22toD91,COVIDtimeD92toD181,COVIDIndD92toD181,COVIDtimeD22toD181,COVIDIndD22toD181,COVIDtimeD22toend,COVIDIndD22toend,
                         AsympInfectIndD15to29, AsympInfectIndD30to91,AsympInfectIndD15to91,AsympInfectIndD92to181,AsympInfectIndD182to271,AsympInfectIndD15to181,
                         ph1.D15, ph1.D92, wt.D15,
                         all_of(immuno_vars))

# create case vs non-case indicators
# Case = COVID-19 endpoint in the interval [≥ 7 days post BD29 and ≥ Dec 1 2021, May 2022 data base lock date]. As described in the appendix the COVID-19 endpoint is documented to be Omicron BA.1 if possible whereas for some non-naive COVID-19 endpoints there was not lineage data available to document the case to be Omicron BA.1.
# Non-case = Did not acquire COVID-19 (of any strain) in the interval [BD1, data base lock date].


dat <- dat %>%
  filter(ph1.D15==1)

dat.2 <- dat

write.csv(dat.2, file = here::here("data_clean", "cor_data_plot3.csv"), row.names=F)
saveRDS(dat.2, file = here::here("data_clean", "cor_data_plot3.rds"))

## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat %>%
  replicate(length(assays),., simplify = FALSE) %>%
  bind_rows()

dat.long.assay_value.names <- times
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assays),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (tt in seq_along(times)) {
  dat_mock_col_names <- paste(times[tt], assays, sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[tt]] <- unlist(lapply(
    # B, Day15, Day29, Day91, Delta15overB
    dat_mock_col_names,

    function(nn) {
      if (nn %in% colnames(dat)) {
        dat[, nn]
      } else {
        rep(NA, nrow(dat))
      }
    }
  ))
}

dat.long.assay_value$assay <- rep(assays, each = nrow(dat))
dat.long <- cbind(dat.long.subject_level, dat.long.assay_value)

## change the labels of the factors for plot labels
dat.long$TrtA <- factor(dat.long$TrtA, levels = c(0, 1), labels = c("Pfizer vaccines", "Moderna vaccines"))
dat.long$TrtB <- factor(dat.long$TrtB, levels = c(0, 1), labels = c("mRNA Omicron-Containing", "mRNA Prototype"))
dat.long$TrtB_moderna <- factor(dat.long$TrtB_moderna, levels = c(0, 1), labels = c("mRNA Omicron-Containing", "mRNA Prototype"))
dat.long$TrtB_pfizer <- factor(dat.long$TrtB_pfizer, levels = c(0, 1), labels = c("mRNA Omicron-Containing", "mRNA Prototype"))
dat.long$TrtC <- factor(dat.long$TrtC, levels = c(0, 1), labels = c("mRNA Monovalent", "mRNA Bivalent"))
dat.long$assay <- factor(dat.long$assay, levels = assays, labels = assays)
dat.long$naive <- factor(dat.long$naive, levels = c(0, 1), labels = c("non-Naive", "Naive"))

# add label = LLoQ, uloq values to show in the plot
dat.long$LLoD = with(dat.long, log10(lods[as.character(assay)]))
dat.long$pos.cutoffs = with(dat.long, log10(pos.cutoffs[as.character(assay)]))
#dat.long$LLoQ = with(dat.long, log10(loqs[as.character(assay)]))
dat.long$lb = with(dat.long,"LoD")
dat.long$lbval = with(dat.long,LLoD)

dat.long$ULoQ = with(dat.long, log10(uloqs[as.character(assay)]))
dat.long$lb2 = "ULoQ" #with(dat.long, ifelse(grepl("bind", assay), "ULoQ", ""))
dat.long$lbval2 =  dat.long$ULoQ #with(dat.long, ifelse(grepl("bind", assay), ULoQ, -99))

# Here, only filter based on ph2.BD29==1. Filtering by ph2.DD1 will occur downstream,
# since it should only happen for DD1-related figures.
# two timepoints study: ph2.tinterm
dat.long.cor.subset <- dat.long 

# write.csv(dat.long.cor.subset, file = here::here("data_clean", "scatter_rug_data_plot5.csv"), row.names=F)
# saveRDS(dat.long.cor.subset, file = here::here("data_clean", "scatter_rug_data_plot5.rds"))

# long to longer format by time
dat.longer.cor.subset <- dat.long.cor.subset %>%
  tidyr::pivot_longer(cols = all_of(times), names_to = "time", values_to = "value")

# phase 2 filters: 
#    include +++, ++- at BD1  for Omicron Cases
#    include +++, ++- at BD29 for Omicron Cases
#    include +++      at DD1  for Omicron Cases
#    include ++ at both BD1 and BD29 for Non-Cases
# if(length(times[!grepl("Delta", times)]) > 2) {
#     
#   dat.longer.cor.subset <- dat.longer.cor.subset %>% 
#     filter(!(time == "DD1" & cohort_event!="Omicron Cases"))  # only keep those +++ at DD1
#   # table(dat.longer.cor.subset$cohort_event, dat.longer.cor.subset$time, dat.longer.cor.subset$assay)
# }

# define response rates
resp <- getResponder(dat.cp, post_times = times[times %in% c("B","Day15")], 
               assays=assays, pos.cutoffs = pos.cutoffs)

resp_by_time_assay <- resp[, c("Ptid", colnames(resp)[grepl("Resp", colnames(resp))])] %>%
  tidyr::pivot_longer(!Ptid, names_to = "category", values_to = "response")

# create Non-cases, D22-D91 Cases, D92-D181 Cases, D22-D181 Cases indicators
dat.longer.cor.subset <- dat.longer.cor.subset %>%
  mutate(cohort_event = case_when(COVIDIndD22toD91==1 ~ "D22-D91 Cases",
                                  COVIDIndD92toD181==1 ~ "D92-D181 Cases",
                                  COVIDIndD22toD181==0 ~ "Non-Cases")
         ) %>%
    rbind(dat.longer.cor.subset %>%
            mutate(cohort_event = case_when(COVIDIndD22toD181==1 ~ "D22-D181 Cases"))%>%
            filter(cohort_event=="D22-D181 Cases"))%>%
  mutate(cohort_event=factor(cohort_event,levels = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases")))


# 1127 unique ids, 3 timepoints, 6 assays
dat.longer.cor.subset <- dat.longer.cor.subset %>%
  mutate(category=paste0(time, assay, "Resp")) %>%
  select(Ptid, time, assay, category, arm, TrtmRNA,TrtonedosemRNA,TrtA,TrtB,TrtC,TrtB_moderna,TrtB_pfizer,
         stage,
         naive,
         value, wt.D15,
         lbval,
         lbval2,
         lb,
         lb2,
         cohort_event) %>%
  left_join(resp_by_time_assay, by=c("Ptid", "category"))%>%
  filter(time %in% c("B","Day15","Delta15overB"))

# subsets for violin/line plots
#### figure specific data prep
# 1. define response rate:
# 2. make subsample datasets such that the violin plot only shows <= 100 non-case data points


#### for figures 1.1, non-case vs case, 13 one-dose mRNA booster arms pooled, by B, Day 15, Delta15overB
groupby_vars1.1=c("cohort_event", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.1 <- get_desc_by_group(subset(dat.longer.cor.subset,dat.longer.cor.subset$TrtonedosemRNA==1), groupby_vars1.1)
write.csv(dat.longer.cor.subset.plot1.1, file = here::here("data_clean", "longer_cor_data_plot1.1.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.1, file = here::here("data_clean", "longer_cor_data_plot1.1.rds"))

#### for figures 1.1.2, non-case vs case, 13 one-dose mRNA booster arms pooled (naive vs. non-naive), by B, Day 15, Delta15overB
groupby_vars1.1.2=c("cohort_event", "naive", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.1.2 <- get_desc_by_group(subset(dat.longer.cor.subset,dat.longer.cor.subset$TrtonedosemRNA==1), groupby_vars1.1.2)
write.csv(dat.longer.cor.subset.plot1.1.2, file = here::here("data_clean", "longer_cor_data_plot1.1.2.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.1.2, file = here::here("data_clean", "longer_cor_data_plot1.1.2.rds"))

#adhoc Peter: I wonder if you could look up a few numbers for me.
#-	Number of the 13 one-dose mRNA ppts that are naïve at baseline
#-	Number of the 13 one-dose mRNA ppts that are non-naïve at baseline.

length(unique(subset(dat.longer.cor.subset.plot1.1.2,dat.longer.cor.subset.plot1.1.2$naive=="Naive"& time=="B")$Ptid))
#629

length(unique(subset(dat.longer.cor.subset.plot1.1.2,dat.longer.cor.subset.plot1.1.2$naive=="non-Naive"& time=="B")$Ptid))
#356

#### for figures 1.2, non-case vs case, 3 Sanofi booster arms pooled, by B, Day 15, Delta15overB
groupby_vars1.2=c("cohort_event", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.2 <- get_desc_by_group(subset(dat.longer.cor.subset,dat.longer.cor.subset$stage==3), groupby_vars1.2)
write.csv(dat.longer.cor.subset.plot1.2, file = here::here("data_clean", "longer_cor_data_plot1.2.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.2, file = here::here("data_clean", "longer_cor_data_plot1.2.rds"))

#adhoc:	Number of the 3 Sanofi recomb at baseline
length(unique(subset(dat.longer.cor.subset.plot1.2,dat.longer.cor.subset.plot1.2$time=="B")$Ptid))
#142


#### for figures 1.3, non-case vs case, mRNA Moderna vs. mRNA Pfizer, by B, Day 15, Delta15overB
groupby_vars1.3=c("cohort_event", "TrtA", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.3 <- get_desc_by_group(subset(dat.longer.cor.subset,!is.na(dat.longer.cor.subset$TrtA)), groupby_vars1.3)
write.csv(dat.longer.cor.subset.plot1.3, file = here::here("data_clean", "longer_cor_data_plot1.3.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.3, file = here::here("data_clean", "longer_cor_data_plot1.3.rds"))


#### for figures 1.4, non-case vs case, mRNA Prototype vs. mRNA Omicron-Containing, by B, Day 15, Delta15overB
groupby_vars1.4=c("cohort_event", "TrtB", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.4 <- get_desc_by_group(subset(dat.longer.cor.subset,!is.na(dat.longer.cor.subset$TrtB)), groupby_vars1.4)
write.csv(dat.longer.cor.subset.plot1.4, file = here::here("data_clean", "longer_cor_data_plot1.4.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.4, file = here::here("data_clean", "longer_cor_data_plot1.4.rds"))


#### for figures 1.4.1, non-case vs case, mRNA Prototype vs. mRNA Omicron-Containing (Moderna), by B, Day 15, Delta15overB
groupby_vars1.4.1=c("cohort_event", "TrtB_moderna", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.4.1 <- get_desc_by_group(subset(dat.longer.cor.subset,!is.na(dat.longer.cor.subset$TrtB_moderna)), groupby_vars1.4.1)
write.csv(dat.longer.cor.subset.plot1.4.1, file = here::here("data_clean", "longer_cor_data_plot1.4.1.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.4.1, file = here::here("data_clean", "longer_cor_data_plot1.4.1.rds"))

#### for figures 1.4.2, non-case vs case, mRNA Prototype vs. mRNA Omicron-Containing (Pfzier), by B, Day 15, Delta15overB
groupby_vars1.4.2=c("cohort_event", "TrtB_pfizer", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.4.2 <- get_desc_by_group(subset(dat.longer.cor.subset,!is.na(dat.longer.cor.subset$TrtB_pfizer)), groupby_vars1.4.2)
write.csv(dat.longer.cor.subset.plot1.4.2, file = here::here("data_clean", "longer_cor_data_plot1.4.2.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.4.2, file = here::here("data_clean", "longer_cor_data_plot1.4.2.rds"))


#adhoc:	
#Number of the Stage 2 Pfizer Omicron at baseline
length(unique(subset(dat.longer.cor.subset.plot1.4.2,dat.longer.cor.subset.plot1.4.2$stage==2 & TrtB_pfizer=="mRNA Omicron-Containing"& time=="B")$Ptid))
#154
#Number of the Stage 2 Pfizer Prototype at baseline
length(unique(subset(dat.longer.cor.subset.plot1.4.2,dat.longer.cor.subset.plot1.4.2$stage==2 & TrtB_pfizer=="mRNA Prototype"& time=="B")$Ptid))
#47


#### for figures 1.5, non-case vs case, mRNA Bivalent vs. mRNA Monovalent, by B, Day 15, Delta15overB
groupby_vars1.5=c("cohort_event", "TrtC", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.5 <- get_desc_by_group(subset(dat.longer.cor.subset,!is.na(dat.longer.cor.subset$TrtC)), groupby_vars1.5)
write.csv(dat.longer.cor.subset.plot1.5, file = here::here("data_clean", "longer_cor_data_plot1.5.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.5, file = here::here("data_clean", "longer_cor_data_plot1.5.rds"))



