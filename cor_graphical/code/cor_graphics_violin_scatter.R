#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
Args <- commandArgs(trailingOnly=TRUE)
COR=Args[1]
if (grepl("IncludeNotMolecConfirmed", COR)) {incNotMol <- "IncludeNotMolecConfirmed"
} else {incNotMol <- ""}
#-----------------------------------------------

source(here::here("code", "cor_process_function.R"))
source(here::here("code", "cor_violin_scatter_function.R"))
source(here::here("..", "_common.R"))
# load parameters
source(here::here("code", "params.R"))

library(scales)
library(tidyverse)
library(here)
library(cowplot) # ggsave2()
library(gridExtra)
library(grid)

## load data 
longer_cor_data <- readRDS(here("data_clean", "longer_cor_data.rds"))
longer_cor_data_plot1 <- readRDS(here("data_clean", "longer_cor_data_plot1.rds")) # at level of trt and assay
plot.25sample1 <- readRDS(here("data_clean", "plot.25sample1.rds"))
if (study_name=="IARCHPV") {
  longer_cor_data_plot1.2 <- readRDS(here("data_clean", "longer_cor_data_plot1.2.rds")) # at level of trt, assay and 
  plot.25sample1.2 <- readRDS(here("data_clean", "plot.25sample1.2.rds"))
}
if (study_name!="IARCHPV") { # IARCHPV doesn't have high risk variable
  longer_cor_data_plot3 <- readRDS(here("data_clean", "longer_cor_data_plot3.rds"))
  plot.25sample3 <- readRDS(here("data_clean", "plot.25sample3.rds"))
}
if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & COR=="D29variant") {
  longer_cor_data_plot_variant <- readRDS(here("data_clean", "longer_cor_data_plot_variant.rds"))
}

### variables for looping
plots <- assays
bstatus <- as.character(unique(longer_cor_data$Bserostatus))

if (study_name=="IARCHPV") { trt = c("pooled", trt.labels)
} else { trt = trt.labels } # add pooled arm as the first arm for IARCHPV
if (study_name=="IARCHPV") {
  labels.assays.short.clean_ = unlist(str_extract_all(gsub("Score", "(Score)", labels.assays.short), "\\([^()]+\\)")) # Get the parenthesis and what is inside
  labels.assays.short.clean = substring(labels.assays.short.clean_, 2, nchar(labels.assays.short.clean_)-1) # Remove parenthesis
  labels.assays.names = names(labels.assays)
  labels.assays = paste(gsub("HPV", "HPV ", gsub("Binding Antibody |to L1, L2 " , "", labels.assays)), "titers")
  names(labels.assays) = labels.assays.names
} else {labels.assays.short.clean = labels.assays.short}
plots_ytitles <- labels.assays.short.clean
plots_titles <- labels.assays[names(labels.assays) %in% names(labels.assays.short)]
timesls <- list(labels.time[(names(labels.time) %in% times) & !grepl("fold-rise", labels.time)][-1], 
                labels.time[(names(labels.time) %in% times) & !grepl("fold-rise", labels.time)],
                labels.time[(names(labels.time) %in% times) & grepl("fold-rise over D1", labels.time)])
if (do.fold.change.overB==0) {timesls[[3]]<-NULL}
if (study_name=="IARCHPV") {timesls[[1]]<-NULL}# single timepoint without baseline

# x-axis need a wrapped version of label
case_grp1 = ifelse(study_name=="PROFISCOV", "Early Post-Peak Cases", "Intercurrent Cases")
case_grp1_wrap = ifelse(study_name=="PROFISCOV", "Early\nPost-Peak\nCases", "Intercurrent\nCases")
case_grp2 = ifelse(study_name=="PROFISCOV", "Late Post-Peak Cases", 
                   ifelse(study_name=="IARCHPV", "Any HPV Cases", "Post-Peak Cases"))
case_grp2_wrap = ifelse(study_name=="PROFISCOV", "Late\nPost-Peak\nCases", 
                        ifelse(study_name=="IARCHPV", "Any HPV Cases", "Post-Peak\nCases"))
ctrl <- ifelse(study_name=="IARCHPV", "Controls", "Non-Cases")

cohort_event_lb = c(case_grp1, case_grp2, ctrl)
names(cohort_event_lb) = c(case_grp1, case_grp2, ctrl)

## common variables with in loop
min_max_plot <- longer_cor_data %>% group_by(assay) %>% summarise(min=min(value, na.rm=T), max=max(value, na.rm=T))
mins <- min_max_plot$min
names(mins) <- min_max_plot$assay
maxs <- min_max_plot$max
names(maxs) <- min_max_plot$assay


# labels
if (length(timepoints)==1){

  x_lb <- c(if (study_name!="IARCHPV") "Day 1", if (study_name!="IARCHPV") paste0("Day ", tpeak), #"Day 2-14\nCases", paste0("Day 15-", 28+tpeaklag, "\nCases"), 
            if (study_name=="IARCHPV") times, 
            case_grp2_wrap,
            ctrl)
  names(x_lb) <- c(if (study_name!="IARCHPV") "Day 1", if (study_name!="IARCHPV") paste0("Day ", tpeak), #"Day 2-14 Cases", paste0("Day 15-", 28+tpeaklag, " Cases"), 
                   if (study_name=="IARCHPV") times, 
                   case_grp2,
                   ctrl)
  
  col_val <- c(#"#FF5EBF", "#0AB7C9", 
    "#FF6F1B", "#810094")
  names(col_val) <- c(#"Day 2-14 Cases", paste0("Day 15-", 28+tpeaklag, " Cases"), 
    case_grp2,
    ctrl)
  
  shp_val <- c(#18, 16, 
    17, 15)
  names(shp_val) <- c(#"Day 2-14 Cases", paste0("Day 15-", 28+tpeaklag, " Cases"), 
    case_grp2,
    ctrl)
  
} else if (COR != "D29variant") {
  x_lb <- c("Day 1", paste0("Day ", tinterm), paste0("Day ", tpeak), 
            paste0("D",tinterm, "\nfold-rise\nover D1"), paste0("D",tpeak, "\nfold-rise\nover D1"), 
            case_grp1_wrap, 
            case_grp2_wrap,
            ctrl)
  
  names(x_lb) <- c("Day 1", paste0("Day ", tinterm), paste0("Day ", tpeak), 
                   paste0("D",tinterm, " fold-rise over D1"), paste0("D",tpeak, " fold-rise over D1"), 
                   case_grp1, 
                   case_grp2,
                   ctrl)
  
  col_val <- c("#0AB7C9","#FF6F1B","#810094")
  names(col_val) <- c(case_grp1, 
                      case_grp2,
                      ctrl)
  
  shp_val <- c(16, 17, 15)
  names(shp_val) <- c(case_grp1, 
                      case_grp2,
                      ctrl)
}

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))


if (COR != "D29variant") {
  #### Figure 1. violin+box plot, case vs non-case, (Day 1), Day 29, and Day 57 if exists
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      if (nrow(subset(longer_cor_data_plot1, Bserostatus == bstatus[j])) == 0 ) next
      for (k in 1:length(trt)) {
        for (t in 1:length(timesls)) { # v1 and v2
          for (case_set in c(#if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") "severe", 
            # comment out on 4/7/2023 because only ENSEMBLE partA primary manuscript needs to be looped through "sev" 
            "Perprotocol")){
          
            y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
            y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]) + 0.25 + ifelse(log10(uloqs[plots[i]])==maxs[plots[i]], 0.1, 0))
            rate.y.pos <- max(y.lim)
            
            font_index <- ifelse(study_name=="IARCHPV", 1.2, 1)
            ll.cex <- 8.16 * font_index
            prop.cex <- 7 * font_index
            
            if (study_name!="IARCHPV"){
              p <- violin_box_plot(dat=subset(longer_cor_data_plot1, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                    dat.sample=subset(plot.25sample1, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                    ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                    facetby=vars(cohort_event),
                                    ylim=y.lim,
                                    ybreaks=y.breaks,
                                    prop.cex=prop.cex,
                                    ll.cex=ll.cex,
                                    group.num=length(levels(longer_cor_data_plot1$cohort_event)),
                                    rate.y.pos=rate.y.pos,
                                    n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                    xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                    )
              g <- grid.arrange(p, bottom = textGrob("All data points for cases are shown. Non-Case data points are shown for all eligible participants or for a random sample of 100 eligible participants, whichever is larger", x = 1, hjust = 1, gp = gpar(fontsize = 15)))
              file_name <- paste0("Linebox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", gsub("-", "_", trt[k]), "_", gsub(" ","",bstatus[j]), "_", if(case_set=="severe") "severe_", "v",t,"_", study_name, ".pdf")
              suppressWarnings(ggsave2(plot = g, filename = paste0(save.results.to, file_name), width = 16, height = 11))
            }
            
            for (casetype in c("Any HPV", unique(longer_cor_data_plot1$persistentindicator[!is.na(longer_cor_data_plot1$persistentindicator)]))) {# loop through any cases and specific breakthrough cases
              
              # skip if specific case type doesn't match with plot type
              if (casetype!="Any HPV" & !grepl(gsub(" ","", casetype),  plots[i])) next
              
              # change label for case type
              if (casetype=="Any HPV") {longer_cor_data_plot1_ = longer_cor_data_plot1
              } else {
                longer_cor_data_plot1_ = longer_cor_data_plot1.2 %>% 
                  filter(as.character(cohort_event) %in% c("Controls", casetype)) %>%
                  mutate(cohort_event = factor(cohort_event,
                                               levels = c(casetype, "Controls")))
                }
              x_lb_ = gsub("Any HPV", casetype, x_lb)
              names(x_lb_) = gsub("Any HPV", casetype, names(x_lb_))
              col_val_ = col_val
              names(col_val_) = NULL #gsub("Any HPV", casetype, names(col_val_))
              shp_val_ = shp_val
              names(shp_val_) = NULL #gsub("Any HPV", casetype, names(shp_val_))
              cohort_event_lb_ = gsub("Any HPV", casetype, cohort_event_lb)
              names(cohort_event_lb_) = gsub("Any HPV", casetype, cohort_event_lb_)
              
              casetype = gsub(" ", "", casetype)
              
              # when k=1 (trt=="pooled") so all trt arms are selected
              p <- violin_box_plot(dat=       subset(longer_cor_data_plot1_, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                    dat.sample=subset(longer_cor_data_plot1_, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                    ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                    x="cohort_event",
                                    xtitle="Cohort Event",
                                    facetby=vars(time),
                                    ylim=y.lim,
                                    type="noline",
                                    ybreaks=y.breaks,
                                    prop.cex=prop.cex,
                                    ll.cex=ll.cex,
                                    pt.size=1.5,
                                    group.num=length(timesls[[t]]),
                                    rate.y.pos=rate.y.pos,
                                    axis.text.x.cex=25 * font_index * 1.6,
                                    axis.text.y.cex=25 * font_index,
                                    n_rate="N_RespRate",
                                    xlabel=x_lb_,
                                    col=col_val_,
                                    shape=shp_val_,
                                    col_lb=cohort_event_lb_,
                                    shp_lb=cohort_event_lb_,
                                    global.size=25 * 1.6
                                    )
              
              file_name <- paste0("Violinbox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", gsub("-", "_", trt[k]), if(bstatus[j]!="") "_", gsub(" ","",bstatus[j]), 
                                  if(casetype!="") "_", casetype, "_", if(case_set=="severe") "severe_", "v",t, "_", study_name, ".pdf")
              suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 11))
            }
          }
        }
      }
    }
  }
  
  if (study_name!="IARCHPV"){
  #### Figure 2. violin + box plot, case vs non-case, (Day 1), Day 29, and Day 57 if exists, by Age, HighRisk, Sex, Race and Ethnic group
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(timesls)) {
          for (case_set in c(#if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") "severe", 
            "Perprotocol")){
            # comment out on 4/7/2023 because only ENSEMBLE partA primary manuscript needs to be looped through "sev" 
            for (s in c("age_geq_65_label","highrisk_label","sex_label","minority_label","Dich_RaceEthnic")) {
              
              groupby_vars2 <- c("Trt", "Bserostatus", "cohort_event", "time", "assay", s)
              
              # define response rate
              longer_cor_data_plot2 <- get_resp_by_group(longer_cor_data, groupby_vars2)
                
              if(s=="Dich_RaceEthnic"){
                longer_cor_data_plot2 <- subset(longer_cor_data_plot2, Dich_RaceEthnic %in% c("Hispanic or Latino","Not Hispanic or Latino"))
              }
              
              longer_cor_data_plot2 <- longer_cor_data_plot2 %>%
                mutate(N_RespRate = ifelse(grepl("Day", time), N_RespRate, ""),
                       lb = ifelse(grepl("Day", time), lb, ""),
                       lbval = ifelse(grepl("Day", time), lbval, NA),
                       lb2 = ifelse(grepl("Day", time), lb2, ""),
                       lbval2 = ifelse(grepl("Day", time), lbval2, NA)) # set fold-rise resp to ""
    
              # make subsample
              plot.25sample2 <- get_sample_by_group(longer_cor_data_plot2, groupby_vars2)
    
              y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
              y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]) + 1)
              rate.y.pos <- max(y.lim)
              
              ll.cex <- 7.5
              prop.cex <- 6.6
              
              p <- violin_box_plot(dat=subset(longer_cor_data_plot2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                    dat.sample=subset(plot.25sample2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1),  
                                    ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                    facetby=as.formula(paste("~",s,"+cohort_event")),
                                    ylim=y.lim,
                                    ybreaks=y.breaks,
                                    prop.cex=prop.cex,
                                    ll.cex=ll.cex,
                                    rate.y.pos=rate.y.pos,
                                    group.num=length(levels(longer_cor_data_plot2$cohort_event)),
                                    n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                    xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                    )
              g <- grid.arrange(p, bottom = textGrob("All data points for cases are shown. Non-Case data points are shown for all eligible participants or for a random sample of 100 eligible participants, whichever is larger", x = 1, hjust = 1, gp = gpar(fontsize = 15)))
              s1 <- ifelse(s=="age_geq_65_label", "Age", ifelse(s=="highrisk_label", "Risk", ifelse(s=="sex_label","Sex", ifelse(s=="minority_label","RaceEthnic", ifelse(s=="Dich_RaceEthnic","Dich_RaceEthnic",NA)))))
              file_name <- paste0("Linebox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", s1, "_", if(case_set=="severe") "severe_", "v", t,"_", study_name, ".pdf")
              suppressWarnings(ggsave2(plot = g, filename = paste0(save.results.to, file_name), width = 16, height = 11))
              
              p <- violin_box_plot(dat=       subset(longer_cor_data_plot2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                    dat.sample=subset(longer_cor_data_plot2, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                    ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                    x="cohort_event",
                                    xtitle="Cohort Event",
                                    facetby=as.formula(paste("~",s,"+time")),
                                    ylim=y.lim,
                                    type="noline",
                                    ybreaks=y.breaks,
                                    prop.cex=prop.cex,
                                    ll.cex=ll.cex,
                                    pt.size=1.5,
                                    group.num=length(timesls[[t]]),
                                    rate.y.pos=rate.y.pos,
                                    axis.text.x.cex=20,
                                    n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                    xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                    )
              file_name <- paste0("Violinbox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_", s1, "_", if(case_set=="severe") "severe_", "v", t,"_", study_name, ".pdf")
              suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 11))
              
            }
          }
        }
      }
    }
  }
  
  
  #### Figure 3. violin + box plot, case vs non-case, (Day 1), Day 29, and Day 57 if exists, by if Age >=65 and if at risk
  for (i in 1:length(plots)) {
    for (j in 1:length(bstatus)) {
      for (k in 1:length(trt)) {
        for (t in 1:length(timesls)) {
          for (case_set in c(#if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") "severe", 
            "Perprotocol")){
            # comment out on 4/7/2023 because only ENSEMBLE partA primary manuscript needs to be looped through "sev" 
    
            y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
            y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]) + 1.5)
            rate.y.pos <- max(y.lim)
            
            prop.cex <- 6.9
            ll.cex <- 7.5
            
            p <- violin_box_plot(dat=subset(longer_cor_data_plot3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                  dat.sample=subset(plot.25sample3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1),
                                  ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                  facetby=as.formula("age_risk_label~cohort_event"),
                                  ylim=y.lim,
                                  ybreaks=y.breaks,
                                  facetopt = "grid",
                                  prop.cex=prop.cex,
                                  ll.cex=ll.cex,
                                  rate.y.pos=rate.y.pos,
                                  group.num=length(levels(longer_cor_data_plot3$cohort_event)),
                                  n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                  xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                  )
            g <- grid.arrange(p, bottom = textGrob("All data points for cases are shown. Non-Case data points are shown for all eligible participants or for a random sample of 100 eligible participants, whichever is larger", x = 1, hjust = 1, gp = gpar(fontsize = 15)))
            file_name <- paste0("Linebox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_Age_Risk_", if(case_set=="severe") "severe_", "v", t,"_", study_name, ".pdf")
            suppressWarnings(ggsave2(plot = g, filename = paste0(save.results.to, file_name), width = 16, height = 13.5))
            
            p <- violin_box_plot(dat=       subset(longer_cor_data_plot3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                  dat.sample=subset(longer_cor_data_plot3, assay==plots[i] & Bserostatus==bstatus[j] & Trt==trt[k] & !is.na(value) & time %in% unlist(timesls[t]) & eval(as.name(case_set))==1), 
                                  ytitle=plots_ytitles[i],toptitle=plots_titles[i],
                                  x="cohort_event",
                                  xtitle="Cohort Event",
                                  facetby=as.formula("age_risk_label~time"),
                                  ylim=y.lim,
                                  type="noline",
                                  ybreaks=y.breaks,
                                  facetopt = "grid",
                                  prop.cex=prop.cex,
                                  ll.cex=ll.cex,
                                  pt.size=1.5,
                                  rate.y.pos=rate.y.pos,
                                  group.num=length(timesls[[t]]),
                                  axis.text.x.cex=20,
                                  n_rate=paste0("N_RespRate", if(case_set=="severe") "_severe"),
                                  xlabel=gsub("\nCases", ifelse(case_set=="severe", "\nSevere\nCases", "\nCases"), x_lb)
                                  )
            file_name <- paste0("Violinbox_", gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])), "_", trt[k], "_", gsub(" ","",bstatus[j]), "_Age_Risk_", if(case_set=="severe") "severe_", "v", t,"_", study_name, ".pdf")
            suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 13.5))
          }
        }
      }
    }
  }
    
    
  #### Figure 4. Scatter plot, assay vs. age in years, case vs non-case, (Day 1), Day 29, and Day 57 if exists
  for (i in 1:length(plots)) {
    for (d in 1:length(timesls[[2]])) { # Day 1, Day 29, Day 57
      for (c in c("Vaccine","all")) {
        for (case_set in c(#if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") "severe", 
          # comment out on 4/7/2023 because only ENSEMBLE partA primary manuscript needs to be looped through "sev" 
          "Perprotocol")){
          
          ds.tmp <- subset(longer_cor_data, assay==plots[i] & time==timesls[[2]][d] & eval(as.name(case_set))==1)
          ds.tmp$size <- with(ds.tmp, ifelse(cohort_event == "Non-Cases", 2.5, 4))
          
          if (timesls[[2]][d]==tail(timesls[[2]], n=1)) {ds.tmp <- ds.tmp %>% 
            filter(!(time==tail(timesls[[2]], n=1) & cohort_event %in% c(case_grp1#,"Day 2-14 Cases", "Day 15-29 Cases", "Day 15-35 Cases"
                                                                         ))) %>%
            mutate(cohort_event = factor(cohort_event, levels = tail(levels(cohort_event), 2)))
          }
          
          y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
          y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
          
          # subset for vaccine arm
          if (c=="Vaccine"){ds.tmp <- subset(ds.tmp, Trt=="Vaccine")}
          
          p <- scatter_plot(
              dat = ds.tmp,
              x="Age", 
              y="value", 
              xtitle="Age (years)",
              ytitle=plots_ytitles[i],
              title=paste0(plots_titles[i],": ",timesls[[2]][d]),
              colby="cohort_event", 
              shaby="cohort_event",
              col=col_val,
              shape=shp_val,
              col_lb=gsub(" Cases", ifelse(case_set=="severe", " Severe Cases", " Cases"), cohort_event_lb),
              shp_lb=gsub(" Cases", ifelse(case_set=="severe", " Severe Cases", " Cases"), cohort_event_lb),
              facetby=as.formula("~Bserostatus+Trt"),
              y.lim=c(floor(mins[plots[i]]), ceiling(maxs[plots[i]])), 
              y.breaks=seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]])),
              x.breaks=seq(from=18, to=86, by=17)
          )

          file_name <- paste0("scatter_",gsub("bind","",gsub("pseudoneut","pnAb_",plots[i])),"_",c,"_",gsub(" ","",timesls[[2]][d]),"_", if(case_set=="severe") "severe_", study_name, ".pdf")
          suppressMessages(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 12.5, height = 11))
        }
      }
    }
  }
  
  #### Figure 5. Scatter plot, assay vs. days since Day 29/Day 1, cases only, 1 panel per assay
  for (i in 1:length(plots)) {
    for (c in c("Vaccine","all")) {
      for (case_set in c(#if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") "severe", 
        "Perprotocol")){
        # comment out on 4/7/2023 because only ENSEMBLE partA primary manuscript needs to be looped through "sev" 
        if(length(timepoints)==1) {
          timesince <- labels.time[(names(labels.time) %in% times) & !grepl("fold-rise", labels.time)] 
        } else {timesince <- labels.time[(names(labels.time) %in% times) & !grepl("fold-rise", labels.time)][-1]}
        
        ds.tmp <- longer_cor_data %>%
          filter(assay==plots[i]) %>%
          filter(eval(as.name(case_set))==1) %>%
          filter(!(time==timesince[2] & cohort_event %in% c(case_grp1#,"Day 2-14 Cases", "Day 15-29 Cases", "Day 15-35 Cases"
                                                            ))) %>% 
          # case only and remove "intercurrent" cases from the last timepoint
          filter(time %in% timesince) %>%
          filter(!cohort_event == "Non-Cases") %>% 
          mutate(cohort_event = factor(cohort_event, levels = head(levels(cohort_event), -1)))
        
        xvar <- ifelse(length(timepoints)>1, paste0("EventTimePrimaryD", tinterm),
                       ifelse(incNotMol=="IncludeNotMolecConfirmed", gsub(tpeak, "1", config.cor$EventTimePrimary), "EventTimePrimaryD1"))
        xlb <- ifelse(length(timepoints)>1, paste0("Days Since the Day ", tinterm," Visit"), "Days Since the Day 1 Visit")
        y.breaks <- seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
        y.lim <- c(floor(mins[plots[i]]), ceiling(maxs[plots[i]]))
        x.breaks <- seq(from=0, to=max(ds.tmp[, xvar], na.rm=T), by=floor(max(ds.tmp[, xvar], na.rm=T)/5))
        x.lim <- c(min(ds.tmp[, xvar], na.rm=T), max(ds.tmp[, xvar], na.rm=T))
        
        # subset for vaccine baseline neg arm
        if (c=="Vaccine"){ds.tmp <- subset(ds.tmp, Trt=="Vaccine")}
      
        p <- ggplot(ds.tmp, aes(x = !!as.name(xvar), y = value, group = time)) + 
          facet_wrap(~Bserostatus+Trt, nrow = 1) + 
          geom_point(alpha = 1, aes(color = cohort_event, shape = cohort_event), size = 4) + 
          geom_line(aes(group = Ptid)) + 
          scale_y_continuous(limits=y.lim, breaks=y.breaks, labels=math_format(10^.x)) +
          scale_x_continuous(limits= x.lim, breaks = x.breaks) +
          scale_color_manual(values = col_val, labels = gsub(" Cases", ifelse(case_set=="severe", " Severe Cases", " Cases"), cohort_event_lb), drop=F) +
          scale_shape_manual(values = shp_val, labels = gsub(" Cases", ifelse(case_set=="severe", " Severe Cases", " Cases"), cohort_event_lb), drop=F) +
          #guides(color = guide_legend(nrow=1)) +
          labs(title = paste0(plots_titles[i], ": ", paste(timesince, collapse=" and ")), x = xlb, y = plots_ytitles[i],
               color="Category", shape="Category") +
          theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
                panel.grid = element_blank(),
                legend.title = element_text(size=22),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(size=ifelse(c=="Vaccine", 27, 19)))
        
        file_name <- paste0("scatter_daysince_",gsub("bind", "", gsub("pseudoneut","pnAb_",plots[i])), "_", c, "_", if(case_set=="severe") "severe_", study_name, ".pdf")
        suppressMessages(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 12.5, height = 11))
      }
    }
  }
  }
}

if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & COR=="D29variant") {
 # Latin America, 5 PsV markers, baseline negative, vaccine
  assay_la <- c("pseudoneutid50","pseudoneutid50_Zeta","pseudoneutid50_Mu","pseudoneutid50_Gamma","pseudoneutid50_Lambda")
  longer_cor_data_plot_variant_la <- longer_cor_data_plot_variant %>%
    filter(assay %in% assay_la & Region == 1 & Bserostatus=="Baseline Neg" & Trt=="Vaccine" & !is.na(value) & time == "Day 29") %>%
    mutate(assay = factor(assay, levels = assay_la, labels = subset(assay_metadata, assay %in% assay_la)$assay_label))
  
  p <- violin_box_plot(dat=longer_cor_data_plot_variant_la, 
                       dat.sample=longer_cor_data_plot_variant_la,
                       ytitle="Pseudovirus-nAb Levels (AU50/ml)",toptitle="Pseudovirus-nAb Levels at Day 29, in Latin America",
                       x="cohort_event",
                       xtitle="Cohort",
                       facetby=as.formula(paste("~","assay")),
                       ylim=c(0, 3.1),
                       type="noline",
                       ybreaks=c(0, 1, 2, 3),
                       prop.cex=6.5,
                       ll.cex=6.5,
                       pt.size=1.5,
                       group.num=2,
                       rate.y.pos=3,
                       axis.text.x.cex=20,
                       col=c("#FF6F1B","#0AB7C9"),
                       shape=c(17, 17),
                       colby="cohort_event",
                       shaby="cohort_event",
                       col_lb=c("Post-Peak Cases","Non-Cases"),
                       shp_lb=c("Post-Peak Cases","Non-Cases"),
                       n_rate="N_RespRate",
                       xlabel=c("Post-Peak Cases","Non-Cases")
  )
  file_name <- "Violinbox_Day29_vaccine_bseroneg_NAb_LA.pdf"
  suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 11))
  
  # Southern America, 3 PsV markers, baseline negative, vaccine
  assay_sa <- c("pseudoneutid50","pseudoneutid50_Delta","pseudoneutid50_Beta")
  longer_cor_data_plot_variant_sa <- longer_cor_data_plot_variant %>%
    filter(assay %in% assay_sa & Region == 2 & Bserostatus=="Baseline Neg" & Trt=="Vaccine" & !is.na(value) & time == "Day 29") %>%
    mutate(assay = factor(assay, levels = assay_sa, subset(assay_metadata, assay %in% assay_sa)$assay_label))
  
  p <- violin_box_plot(dat=longer_cor_data_plot_variant_sa, 
                       dat.sample=longer_cor_data_plot_variant_sa,
                       ytitle="Pseudovirus-nAb Levels (AU50/ml)",toptitle="Pseudovirus-nAb Levels at Day 29, in South Africa",
                       x="cohort_event",
                       xtitle="Cohort",
                       facetby=as.formula(paste("~","assay")),
                       ylim=c(0, 3.1),
                       type="noline",
                       ybreaks=c(0, 1, 2, 3),
                       prop.cex=6.5,
                       ll.cex=6.5,
                       pt.size=1.5,
                       group.num=2,
                       rate.y.pos=3,
                       axis.text.x.cex=20,
                       col=c("#FF6F1B","#0AB7C9"),
                       shape=c(17, 17),
                       colby="cohort_event", 
                       shaby="cohort_event",
                       col_lb=c("Post-Peak Cases","Non-Cases"),
                       shp_lb=c("Post-Peak Cases","Non-Cases"),
                       n_rate="N_RespRate",
                       xlabel=c("Post-Peak Cases","Non-Cases")
  )
  file_name <- "Violinbox_Day29_vaccine_bseroneg_NAb_SA.pdf"
  suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 11))
  
  
  # Northern America, 1 PsV marker, baseline negative, vaccine
  assay_na <- c("pseudoneutid50")
  longer_cor_data_plot_variant_na <- longer_cor_data_plot_variant %>%
    filter(assay %in% assay_na & Region == 0 & Bserostatus=="Baseline Neg" & Trt=="Vaccine" & !is.na(value) & time == "Day 29") %>%
    mutate(assay = factor(assay, levels = assay_na, subset(assay_metadata, assay %in% assay_na)$assay_label))
  
  p <- violin_box_plot(dat=longer_cor_data_plot_variant_na, 
                       dat.sample=longer_cor_data_plot_variant_na,
                       ytitle="Pseudovirus-nAb Levels (AU50/ml)",toptitle="Pseudovirus-nAb Levels at Day 29, in United States",
                       x="cohort_event",
                       xtitle="Cohort",
                       facetby=as.formula(paste("~","assay")),
                       ylim=c(0, 3.1),
                       type="noline",
                       ybreaks=c(0, 1, 2, 3),
                       prop.cex=9,
                       ll.cex=9,
                       pt.size=3,
                       group.num=2,
                       rate.y.pos=3,
                       axis.text.x.cex=25,
                       col=c("#FF6F1B","#0AB7C9"),
                       shape=c(17, 17),
                       colby="cohort_event", 
                       shaby="cohort_event",
                       col_lb=c("Post-Peak Cases","Non-Cases"),
                       shp_lb=c("Post-Peak Cases","Non-Cases"),
                       n_rate="N_RespRate",
                       xlabel=c("Post-Peak Cases","Non-Cases")
  )
  file_name <- "Violinbox_Day29_vaccine_bseroneg_NAb_US.pdf"
  suppressWarnings(ggsave2(plot = p, filename = paste0(save.results.to, file_name), width = 16, height = 11))
  
}
