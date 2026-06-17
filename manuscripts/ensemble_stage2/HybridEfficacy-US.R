set.seed(206)
if(!exists("spath")) {
  spath <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2"
}
source(paste0(spath, "/code/Methods-CausalCumulativeIncidence.R"))
if(!exists("data.path")) {
  data.path <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"
}
ensemble.dat <- read.csv(data.path) # read.csv(paste0(spath, "/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"))
ensemble.dat <- ensemble.dat[ensemble.dat$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat$History %in% c("History a", "History b") &
                             ensemble.dat$Region == 0,]
ensemble.dat <- ensemble.dat[ensemble.dat$BOOSTDT < as.Date("2022-05-22"),]
ensemble.dat$White <- with(ensemble.dat, ifelse(Black == 0 & Asian == 0 & NatAmer == 0 & PacIsl == 0 &
                                                  Multiracial == 0 & Notreported == 0 & Unknown == 0, 1, 0))
ensemble.dat$OtherRace <- with(ensemble.dat, ifelse(Notreported == 1 | Unknown == 1, 1, 0))

# cohort x history
ensemble.dat$DelayedVaccine <- ifelse(ensemble.dat$Cohort == "Cohort 1" & ensemble.dat$History == "History a", 1, 0)
ensemble.dat$ImmediateVaccine <- ifelse(ensemble.dat$Cohort == "Cohort 1" & ensemble.dat$History == "History b", 1, 0)
ensemble.dat$DelayedHybrid <- ifelse(ensemble.dat$Cohort == "Cohort 2" & ensemble.dat$History == "History a", 1, 0)
ensemble.dat$ImmediateHybrid <- ifelse(ensemble.dat$Cohort == "Cohort 2" & ensemble.dat$History == "History b", 1, 0)

ensemble.dat$group <- character(nrow(ensemble.dat))
ensemble.dat$group <- ifelse(ensemble.dat$ImmediateVaccine == 0,
                             ensemble.dat$group, rep("Immediate Vaccine", nrow(ensemble.dat)))
ensemble.dat$group <- ifelse(ensemble.dat$DelayedVaccine == 0,
                             ensemble.dat$group, rep("Delayed Vaccine", nrow(ensemble.dat)))
ensemble.dat$group <- ifelse(ensemble.dat$ImmediateHybrid == 0,
                             ensemble.dat$group, rep("Immediate Hybrid", nrow(ensemble.dat)))
ensemble.dat$group <- ifelse(ensemble.dat$DelayedHybrid == 0,
                             ensemble.dat$group, rep("Delayed Hybrid", nrow(ensemble.dat)))

ensemble.dat$BMIDich <- ifelse(ensemble.dat$BMI %in% c(1, 2), 0, 1)
ensemble.dat$BMIDich[is.na(ensemble.dat$BMI)] <- NA

ensemble.dat$RegionUS <- ifelse(ensemble.dat$Region == 0, 1, 0)
ensemble.dat$RegionLA <- ifelse(ensemble.dat$Region == 1, 1, 0)
ensemble.dat$RegionRSA <- ifelse(ensemble.dat$Region == 2, 1, 0)

ensemble.dat$HistoryA <- ifelse(ensemble.dat$History == "History a", 1, 0)

# Event variables for efficacy analyses
ensemble.dat$AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         0, ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - 13)
ensemble.dat$AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc)

ensemble.dat$AdjModEventIndPrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         0, ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         0, ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - 13)
ensemble.dat$AdjModEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD14_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModEventTimePrimaryBD14_CensorAtOutsideVacc)

ensemble.dat$AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         0, ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         1, ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - 13)
ensemble.dat$AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc)

ensemble.dat$AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         0, ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         1, ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - 13)
ensemble.dat$AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc)

ensemble.dat$AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         0, ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc < 14,
         1, ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - 13)
ensemble.dat$AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc)

# calendar time
ensemble.dat$CalendarDateDose2 <- Date(nrow(ensemble.dat)) 
for(i in 1:nrow(ensemble.dat)) {
  if(ensemble.dat$group[i] %in% c("Immediate Hybrid", "Immediate Vaccine")) {
    ensemble.dat$CalendarDateDose2[i] <- as.Date(ensemble.dat$TRTSDT[i]) + ensemble.dat$NumberdaysD1toYear1booster[i]
  } else {
    ensemble.dat$CalendarDateDose2[i] <- as.Date(ensemble.dat$CROSSDT[i]) + ensemble.dat$NumberdaysD1toYear1booster[i]
  }
}

ensemble.dat$CalendarDateDose2 <- as.Date(ensemble.dat$BOOSTDT)
ensemble.dat$CalendarDateDose2Early <- ifelse(ensemble.dat$CalendarDateDose2 < as.Date("2021-12-01"), 1, 0)



#######
A <- with(ensemble.dat, ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
W <- with(ensemble.dat, cbind(Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex), HistoryA,
                              CalendarDateDose2Early))
W.sev <- with(ensemble.dat, cbind(Age, Sex, BMIDich,
                                  CalendarDateDose2Early))
complete <- apply(cbind(A, W), 1, function(x){all(!is.na(x))})

event.sev <- ensemble.dat$AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc
time.sev <- ensemble.dat$AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc

event.any.mod <- ensemble.dat$AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc
time.any.mod <- ensemble.dat$AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc

event.any <- ensemble.dat$AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc
time.any <- ensemble.dat$AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc

t0 <- 90 - 13

fit.sev <- CumulativeIncidenceFit(event = event.sev[complete], time = time.sev[complete],
                                  A = A[complete], W = W[complete,],
                                  t0 = t0, fit.times = seq(1, t0, 1), np = FALSE)
fit.any.mod <- CumulativeIncidenceFit(event = event.any.mod[complete], time = time.any.mod[complete],
                                      A = A[complete], W = W[complete,],
                                      t0 = t0, fit.times = seq(1, t0, 1), np = FALSE)
fit.any <- CumulativeIncidenceFit(event = event.any[complete], time = time.any[complete],
                                  A = A[complete], W = W[complete,],
                                  t0 = t0, fit.times = seq(1, t0, 1), np = FALSE)

out.sev <- CumulativeIncidenceSummary(time = time.sev[complete], event = event.sev[complete],
                                     W = W.sev[complete,],
                                     A = A[complete], t0 = t0, fit = fit.sev,
                                     t.start = 1, label = "Severe-critical (United States)",
                                     time.axis.label = "Time (Days after Second Ad26.COV2.S Dose)",
                                     path.prefix = paste0(spath, "/reports/"),
                                     path.suffix = paste0("severe-US_", Sys.Date()), y.axis.bound.risk = .25,
                                     time.axis.offset = 13)
out.any.mod <- CumulativeIncidenceSummary(time = time.any.mod[complete], event = event.any.mod[complete],
                                          W = W[complete,],
                                          A = A[complete], t0 = t0, fit = fit.any.mod,
                                          t.start = 1, label = "Moderate to Severe-critical (United States)",
                                          time.axis.label = "Time (Days after Second Ad26.COV2.S Dose)",
                                          path.prefix = paste0(spath, "/reports/"),
                                          path.suffix = paste0("moderate-to-severe-US_", Sys.Date()), y.axis.bound.risk = .25,
                                          time.axis.offset = 13)
out.any <- CumulativeIncidenceSummary(time = time.any[complete], event = event.any[complete],
                                      W = W[complete,],
                                      A = A[complete], t0 = t0, fit = fit.any,
                                      t.start = 1, label = "Infection, Any Symptomatology (United States)",
                                      time.axis.label = "Time (Days after Second Ad26.COV2.S Dose)",
                                      path.prefix = paste0(spath, "/reports/"),
                                      path.suffix = paste0("any-infection-US_", Sys.Date()), y.axis.bound.risk = .25,
                                      time.axis.offset = 13)


#adaptive y-axis
# out.sev.aaxis <- CumulativeIncidenceSummary(time = time.sev[complete], event = event.sev[complete],
#                                             W = W[complete,],
#                                             A = A[complete], t0 = t0, fit = fit.sev,
#                                             t.start = 1, label = "Severe-critical (United States)",
#                                             time.axis.label = "Time (Days after Second Ad26.COV2.S Dose)",
#                                             path.suffix = "severe-US-adap-axis-1-22-2026",
#                                             time.axis.offset = 13)
# out.any.mod.aaxis <- CumulativeIncidenceSummary(time = time.any.mod[complete], event = event.any.mod[complete],
#                                                 W = W[complete,],
#                                                 A = A[complete], t0 = t0, fit = fit.any.mod,
#                                                 t.start = 1, label = "Moderate to Severe-critical (United States)",
#                                                 time.axis.label = "Time (Days after Second Ad26.COV2.S Dose)",
#                                                 path.suffix = "moderate-to-severe-US-adap-axis-1-22-2026",
#                                                 time.axis.offset = 13)
# out.any.aaxis <- CumulativeIncidenceSummary(time = time.any[complete], event = event.any[complete],
#                                             W = W[complete,],
#                                             A = A[complete], t0 = t0, fit = fit.any,
#                                             t.start = 1, label = "Infection, Any Symptomatology (United States)",
#                                             time.axis.label = "Time (Days after Second Ad26.COV2.S Dose)",
#                                             path.suffix = "any-infection-US-adap-axis-1-22-2026",
#                                             time.axis.offset = 13)


tab <- matrix(NA, nrow = 3, ncol = 4)

rownames(tab) <- c(paste0("Severe-critical COVID through ",  out.sev$t0 + 13, " Days"),
                   paste0("Moderate to Severe-critical COVID through ",  out.any.mod$t0 + 13, " Days"),
                   paste0("Infection, Any Symptomatology through ", out.any$t0 + 13, " Days"))
colnames(tab) <- c("Vaccine Incidence", "Hybrid Incidence", "Incidence Ratio", "P-value")
tab[1,] <- c(paste0(signif(out.sev$incidence.ctrl[1], 2),
                    " (", signif(out.sev$incidence.ctrl[2], 2),
                    ", ", signif(out.sev$incidence.ctrl[3], 2), ")"),
             paste0(signif(out.sev$incidence.tx[1], 2),
                    " (", signif(out.sev$incidence.tx[2], 2),
                    ", ", signif(out.sev$incidence.tx[3], 2), ")"),
             paste0(signif(out.sev$rr.tab[1], 2),
                    " (", signif(out.sev$rr.tab[2], 2),
                    ", ", signif(out.sev$rr.tab[3], 2), ")"),
             signif(out.sev$rr.tab[4], 4))
tab[2,] <- c(paste0(signif(out.any.mod$incidence.ctrl[1], 2),
                    " (", signif(out.any.mod$incidence.ctrl[2], 2),
                    ", ", signif(out.any.mod$incidence.ctrl[3], 2), ")"),
             paste0(signif(out.any.mod$incidence.tx[1], 2),
                    " (", signif(out.any.mod$incidence.tx[2], 2),
                    ", ", signif(out.any.mod$incidence.tx[3], 2), ")"),
             paste0(signif(out.any.mod$rr.tab[1], 2),
                    " (", signif(out.any.mod$rr.tab[2], 2),
                    ", ", signif(out.any.mod$rr.tab[3], 2), ")"),
             signif(out.any.mod$rr.tab[4], 4))
tab[3,] <- c(paste0(signif(out.any$incidence.ctrl[1], 2),
                    " (", signif(out.any$incidence.ctrl[2], 2),
                    ", ", signif(out.any$incidence.ctrl[3], 2), ")"),
             paste0(signif(out.any$incidence.tx[1], 2),
                    " (", signif(out.any$incidence.tx[2], 2),
                    ", ", signif(out.any$incidence.tx[3], 2), ")"),
             paste0(signif(out.any$rr.tab[1], 2),
                    " (", signif(out.any$rr.tab[2], 2),
                    ", ", signif(out.any$rr.tab[3], 2), ")"),
             signif(out.any$rr.tab[4], 4))
# write.csv(tab, "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Efficacy-US-1-22-2026.csv")
write.csv(tab, paste0(spath, "/reports/Efficacy-US_", Sys.Date(), ".csv"))
