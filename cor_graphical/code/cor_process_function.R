#-----------------------------------------------
# obligatory to append to the top of each script
#-----------------------------------------------

library(tidyverse)

#' Generate response calls and fold-rise indicator for endpoints:
#' Responders at each pre-defined timepoint are defined as participants who
#' had baseline values below the LLOD with detectable ID80 neutralization titer
#' above the assay LLOD, or as participants with baseline values above
#' the LLOD with a 4-fold increase in neutralizing antibody titer.
#'
#' @param data The dataframe with the endpoint of interest at baseline and
#'  post-baseline.
#' @param bl The variable of endpoint value at baseline
#' @param post The variable of endpoint value at post-baseline
#' @param folds Folds to generate fold-rise indicator
#' @param cutoff LLOD or LLOQ
#' @param respoderFR Fold-rise used for response call positivity criteria
#'
#' @return Response calls and fold-rise indicator for \code{bl}: \emph{bl}Resp,
#'  \emph{bl}FR\emph{folds}
getResponder <- function(data,
                         cutoff.name, 
                         times=times, 
                         assays=assays, 
                         #folds=c(2, 4),
                         grtns=c(2, 4),
                         responderFR = 4,
                         pos.cutoffs = pos.cutoffs) {
  
  #cutoff <- get(paste0(cutoff.name, "s"))
  for (i in times){
    for (j in assays){
      post <- paste0(i, j)
      bl <- paste0("B", j)
      delta <- paste0("Delta", gsub("Day", "", i), "overB", j) 
      
      # add these two rows for assays only done at certain timepoints for stage 1 studies, but not baseline
      if (!bl %in% colnames(data) & !grepl("stage2", COR)) {data[, bl] <- 0}
      if (!delta %in% colnames(data) & !grepl("stage2", COR)) {data[, delta] <- data[, post] - data[, bl]}
      
      if (bl %in% colnames(data)) {data[, bl] <- pmin(data[, bl], log10(uloqs[j]))}
      data[, post] <- pmin(data[, post], log10(uloqs[j]))
      #data[, delta] <- ifelse(10^data[, post] < cutoff[j], log10(cutoff[j]/2), data[, post])-ifelse(10^data[, bl] < cutoff[j], log10(cutoff[j]/2), data[, bl])
      
      #for (k in folds){
      #  data[, paste0(post, k, cutoff.name)] <- as.numeric(10^data[, post] >= k*cutoff[j])
      #}
      
      if (delta %in% colnames(data)) {# this delta only exists for stage 1 studies
        for (k in grtns){
          data[, paste0(post, "FR", k)] <- as.numeric(10^data[, delta] >= k)
        }
      }
      
      # if (!is.na(pos.cutoffs[j])
      if ((grepl("bind|ACE2", j) & !study_name %in% c("VAT08", "NextGen_Mock")) | COR == "D29VLvariant" | grepl("stage2", COR)) {
        data[, paste0(post, "Resp")] <- as.numeric(data[, post] > log10(pos.cutoffs[j]))
        if (bl %in% colnames(data)) {data[, paste0(bl, "Resp")] <- as.numeric(data[, bl] > log10(pos.cutoffs[j]))}
      } else if (study_name == "NextGen_Mock") {
        
        # 1 if baseline < lloq, post-baseline >= 4 * lloq
        # if lloq <= Baseline < uloq, post-baseline >= 4 * baseline
        data[, paste0(bl, "Resp")] <- as.numeric(data[, bl] > log10(pos.cutoffs[j]))
        data[, paste0(post, "Resp")] <- as.numeric(
          (data[, bl] < log10(pos.cutoffs[j]) & data[, post] > 4 * log10(pos.cutoffs[j])) |
            (data[, bl] >= log10(pos.cutoffs[j]) & data[, bl] < log10(uloqs[j]) & as.numeric(10^data[, post]/10^data[, bl] >= responderFR)) |
            data[, bl] < log10(uloqs[j]) & data[, post] >= log10(uloqs[j]) & as.numeric(uloqs[j]/10^data[, bl] < responderFR) ) 
        
        over_uloq <- which(!is.na(data[, bl]) & data[, bl] >= log10(uloqs[j]))
        data[over_uloq, paste0(bl, "Resp")] <- NA
        data[over_uloq, paste0(post, "Resp")] <- NA
        
      } else {
        data[, paste0(post, "Resp")] <- as.numeric(
          (data[, bl] < log10(pos.cutoffs[j]) & data[, post] > log10(pos.cutoffs[j])) |
            (data[, bl] >= log10(pos.cutoffs[j]) & data[, paste0(post, "FR", responderFR)] == 1))
        data[, paste0(bl, "Resp")] <- as.numeric(data[, bl] > log10(pos.cutoffs[j]))
      }
    }
  }
  return(data)
}

# a function to define response rate by group
get_resp_by_group <- function(dat=dat, group=group){
  
  if(length(timepoints)==1) {
    dat[which(dat$time=="Day 1"), "wt"] = 1
    dat[which(dat$time!="Day 1"), "wt"] = dat[which(dat$time!="Day 1"), config.cor$wt] # wt.D29 or wt.D29start1
    # special case for Janssen_partA_VL, ancestral assays use the wt.D29 instead of wt.D29variant
    if ("wt.D29" %in% colnames(dat)) {dat[which(dat$time!="Day 1" & dat$assay %in% c("bindSpike","pseudoneutid50")), "wt"]=dat[which(dat$time!="Day 1" & dat$assay %in% c("bindSpike","pseudoneutid50")),"wt.D29"]}
  } else if (study_name=="VAT08"){
    dat[which(dat$time=="Day 1"), "wt"]=1 # for intercurrent cases, we don't need to adjust for the weight because all of them are from the same stratum
   
    dat[which(dat$time=="Day 22" & grepl("bind", dat$assay)), "wt"]=dat[which(dat$time=="Day 22" & grepl("bind", dat$assay)), "wt.D22.bAb"]
    dat[which(dat$time=="Day 22" & grepl("pseudoneutid", dat$assay)), "wt"]=dat[which(dat$time=="Day 22" & grepl("pseudoneutid", dat$assay)), "wt.D22.nAb"]
    
    dat[which(dat$time=="Day 43" & grepl("bind", dat$assay)), "wt"]=dat[which(dat$time=="Day 43" & grepl("bind", dat$assay)), "wt.D43.bAb"]
    dat[which(dat$time=="Day 43" & grepl("pseudoneutid", dat$assay)), "wt"]=dat[which(dat$time=="Day 43" & grepl("pseudoneutid", dat$assay)), "wt.D43.nAb"]
  } else if (study_name == "NextGen_Mock") {
    dat$wt = dat[, "ph2.AB.trackA"]  # ccIAS, # initial, 5 timepoints, on Track A ccIAS-PBMC
    dat$wt2 = dat[, config.cor$wt] # final, 3 timepoints, on whole ccIAS-PBMC
    
  } else {
    dat[which(dat$time=="Day 1"), "wt"] = 1 # for intercurrent cases, we don't need to adjust for the weight because all of them are from the same stratum
    dat[which(dat$time==paste0("Day ",timepoints[1])), "wt"] = dat[which(dat$time==paste0("Day ",timepoints[1])), paste0("wt.D",timepoints[1])]
    dat[which(dat$time==paste0("Day ",timepoints[length(timepoints)])), "wt"] = dat[which(dat$time==paste0("Day ",timepoints[length(timepoints)])), paste0("wt.D",timepoints[length(timepoints)])]
  }
  
  complete <- complete.cases(dat[, group])
  
  dat_resp_by_group <-
    dat %>% filter(complete==1) %>%
    group_by_at(group) %>%
    mutate(counts = ifelse(study_name == "NextGen_Mock", sum(!is.na(response) & as.numeric(Track == "A")), n()),
           #counts_severe = sum(severe, na.rm=T),
           # comment out on 4/7/2023 because only ENSEMBLE partA primary manuscript needs to be looped through "sev" 
           num = ifelse(study_name == "NextGen_Mock", sum(response * wt * as.numeric(Track == "A"), na.rm=T), sum(response * wt, na.rm=T)), 
           #num_severe = sum(response * wt & severe==1, na.rm=T),
           denom = ifelse(study_name == "NextGen_Mock", sum(wt * as.numeric(Track == "A"), na.rm=T), sum(wt, na.rm=T)),
           #denom_severe = sum(wt & severe==1, na.rm=T),
           # test N_RespRate = paste0(counts, "\n", sum(response, na.rm=T), ",", round(sum(response, na.rm=T)/counts*100, 1), ",", LLoD),
           N_RespRate = ifelse(!grepl("Delta|fold", time) && !is.na(pos.cutoffs), paste0(counts, "\n",round(num/denom*100, 1),"%"), ""),
           #N_RespRate_severe = paste0(counts_severe, "\n",round(num_severe/denom_severe*100, 1),"%"),
           min = min(value),
           q1 = quantile(value, 0.25, na.rm=T),
           median = median(value, na.rm=T),
           q3 = quantile(value, 0.75, na.rm=T),
           max= max(value)) %>%
    mutate(N_RespRate = ifelse(grepl("mdw", assay), "", N_RespRate))
  
  if (study_name == "NextGen_Mock") {
    dat$wt2 <- dat[[config.cor$wt]]
  }
  
  if (study_name == "NextGen_Mock") {
    dat_resp_by_group2 <-
      dat %>% filter(complete == 1 & grepl("bind|pseudo", assay)) %>%
      filter(.data[[config.cor$ph2]] == 1) %>% # condition for the whole RIS for bAb/nAb and RIS-PBMC for ICS
      group_by_at(group) %>%
      mutate(counts = sum(!is.na(response)),
             num = sum(response * wt2, na.rm=T),
             denom = sum(wt2, na.rm=T),
             #N_RespRate = paste0(counts, "\n",round(num/denom*100, 1),"%"),
             N_RespRate = ifelse(!grepl("Delta|fold", time) && !is.na(pos.cutoffs), paste0(counts, "\n", round(num/denom*100, 1),"%"), ""), # RespRate at Delta timepoints will be ""
             min = min(value, na.rm=T),
             q1 = quantile(value, 0.25, na.rm=T),
             median = median(value, na.rm=T),
             q3 = quantile(value, 0.75, na.rm=T),
             max= max(value, na.rm=T)) %>%
      #bind_rows(dat %>% filter(complete == 1 & grepl("T4|T8", assay)) %>%
      #            filter(.data[[config.cor$ph2]] == 1) %>% # condition for the whole RIS for bAb/nAb and RIS-PBMC for ICS
      #            group_by_at(group) %>%
      #            mutate(counts = sum(!is.na(response)),
      #                   num = sum(response * wt, na.rm=T),
      #                   denom = sum(wt, na.rm=T),
      #                   N_RespRate = ifelse(!grepl("Delta", time) && !is.na(pos.cutoffs), paste0(counts, "\n", round(num/denom*100, 1),"%"), ""), # RespRate at Delta timepoints will be ""
      #                   min = min(value, na.rm=T),
      #                   q1 = quantile(value, 0.25, na.rm=T),
      #                   median = median(value, na.rm=T),
      #                   q3 = quantile(value, 0.75, na.rm=T),
      #                   max= max(value, na.rm=T))) %>%
      mutate(N_RespRate = ifelse(grepl("mdw", assay), "", N_RespRate))
  }
  
  if (exists("dat_resp_by_group2")) {
    return(list(dat_resp_by_group = dat_resp_by_group, dat_resp_by_group2 = dat_resp_by_group2))
  } else {
    return(dat_resp_by_group)
  }
}

# a function to get 100 non-case sample by group for points & lines in line plot 
get_sample_by_group <- function(dat.plot, group){
  
  dat.sample <- dat.plot %>%
    group_by_at(group) %>%
    sample_n((ifelse(n()>=100 & cohort_event=="Non-Cases", 100, n())), replace=F) %>% filter(time==paste0(substr(times[2],1,3)," ", substr(times[2],4,5))) %>% # eg time=="Day 29"
    ungroup() %>%
    select(c("Ptid", group[!group %in% "time"])) %>%
    inner_join(dat.plot, by=c("Ptid", group[!group %in% "time"]))
  return(dat.sample)
} 
