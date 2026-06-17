require(survminer)
require(survival)
require(lubridate)
set.seed(206)
# ensemble.dat <- read.csv("/Volumes/trials/covpn/p3003/analysis/mapping_immune_correlates/post unblinded part B stage 2/adata/COVID_ENSEMBLE_stage2_mapped_20250909.csv")
# ensemble.dat <- read.csv("/Volumes/trials/covpn/p3003/analysis/mapping_immune_correlates/post unblinded part B stage 2/adata/COVID_ENSEMBLE_stage2_mapped_20250925.csv")
# ensemble.dat <- read.csv("/Volumes/trials/covpn/p3003/analysis/mapping_immune_correlates/post unblinded part B stage 2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv")
if(!exists("spath")) {
  spath <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2"
}
if(!exists("data.path")) {
  data.path <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"
}
ensemble.dat <- read.csv(data.path) # read.csv(paste0(spath, "/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"))
ensemble.dat <- ensemble.dat[ensemble.dat$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat$History %in% c("History a", "History b") &
                             ensemble.dat$Region == 2,]
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

# Event variables for efficacy analyses
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

########################################
### Event counts and event rate by group
########################################

Tab <- matrix(nrow = 6, ncol = 5)
colnames(Tab) <- c("Crossover/Vaccine", "Crossover/Hybrid", "Original/Vaccine", "Original/Hybrid", "Overall")
rownames(Tab) <- c("N", "Moderate to Severe-critical COVID (Total, Annual Incidence)",
                   "Moderate COVID (Total, Annual Incidence)",
                   "Severe-critical COVID (Total, Annual Incidence)",
                   "Asymptomatic Infection (Total, Annual Incidence)",
                   "Infection, Any Symptomatology (Total, Annual Incidence)")
# N
Tab[1,] <- c(sum(ensemble.dat$ImmediateVaccine), sum(ensemble.dat$ImmediateHybrid),
             sum(ensemble.dat$DelayedVaccine), sum(ensemble.dat$DelayedHybrid), nrow(ensemble.dat))

# Moderate to severe COVID
event.total <- with(ensemble.dat, 
                    c(sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc)))
annual.incidence <- with(ensemble.dat,
                         c(sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                           sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                           sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                           sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                           sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc)/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
Tab[2,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")


# Moderate COVID

# mutally exclusive defintion of moderate covid
#set event indicator to 0 when moderate and severe covid occur concurrently
ensemble.dat$AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME <- 
  ifelse(ensemble.dat$AdjModEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
           ensemble.dat$AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc == 1,
         0, ensemble.dat$AdjModEventIndPrimaryBD14_CensorAtOutsideVacc)
#set censoring time as time of severe covid when moderate nad severe occur concurrently
ensemble.dat$AdjModEventTimePrimaryBD14_CensorAtOutsideVacc_ME <-
  ifelse(ensemble.dat$AdjModEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
           ensemble.dat$AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc == 1,
         ensemble.dat$AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc, 
         ensemble.dat$AdjModEventTimePrimaryBD14_CensorAtOutsideVacc)

event.total <- with(ensemble.dat, 
                    c(sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME[ImmediateVaccine == 1]),
                      sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME[ImmediateHybrid == 1]),
                      sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME[DelayedVaccine == 1]),
                      sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME[DelayedHybrid == 1]),
                      sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME)))
annual.incidence <- with(ensemble.dat,
                         c(sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME[ImmediateVaccine == 1])/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc_ME[ImmediateVaccine == 1]),
                           sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME[ImmediateHybrid == 1])/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc_ME[ImmediateHybrid == 1]),
                           sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME[DelayedVaccine == 1])/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc_ME[DelayedVaccine == 1]),
                           sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME[DelayedHybrid == 1])/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc_ME[DelayedHybrid == 1]),
                           sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc_ME)/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc_ME))) * 365
Tab[3,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# event.total <- with(ensemble.dat, 
#                     c(sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                       sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                       sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                       sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                       sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc)))
# annual.incidence <- with(ensemble.dat,
#                          c(sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc)/
#                              sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
# Tab[3,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Severe COVID
event.total <- with(ensemble.dat,
                    c(sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc)))
annual.incidence <- with(ensemble.dat,
                         c(sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                           sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                           sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                           sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                           sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc)/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
Tab[4,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Asymp (or other infection)
# mutally exclusive defintion of asymptomatic or other COVID
#set event indicator to 0 when asym and mod-severe covid occur concurrently
ensemble.dat$AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME <- 
  ifelse(ensemble.dat$AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
           ensemble.dat$AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc == 1,
         0, ensemble.dat$AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc)
#set censoring time as time of mod-severe covid when mod-severe and asym occur concurrently
ensemble.dat$AdjAsymOtherEventTimePrimaryBD14_CensorAtOutsideVacc_ME <-
  ifelse(ensemble.dat$AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
           ensemble.dat$AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc == 1,
         ensemble.dat$AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc, 
         ensemble.dat$AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc)

event.total <- with(ensemble.dat,
                    c(sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME[ImmediateVaccine == 1]),
                      sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME[ImmediateHybrid == 1]),
                      sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME[DelayedVaccine == 1]),
                      sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME[DelayedHybrid == 1]),
                      sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME)))
annual.incidence <- with(ensemble.dat,
                         c(sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME[ImmediateVaccine == 1])/
                             sum(AdjAsymOtherEventTimePrimaryBD14_CensorAtOutsideVacc_ME[ImmediateVaccine == 1]),
                           sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME[ImmediateHybrid == 1])/
                             sum(AdjAsymOtherEventTimePrimaryBD14_CensorAtOutsideVacc_ME[ImmediateHybrid == 1]),
                           sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME[DelayedVaccine == 1])/
                             sum(AdjAsymOtherEventTimePrimaryBD14_CensorAtOutsideVacc_ME[DelayedVaccine == 1]),
                           sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME[DelayedHybrid == 1])/
                             sum(AdjAsymOtherEventTimePrimaryBD14_CensorAtOutsideVacc_ME[DelayedHybrid == 1]),
                           sum(AdjAsymOtherEventIndPrimaryBD14_CensorAtOutsideVacc_ME)/
                             sum(AdjAsymOtherEventTimePrimaryBD14_CensorAtOutsideVacc_ME))) * 365
Tab[5,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Infection, Any Symptomatology
event.total <- with(ensemble.dat, 
                    c(sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc)))
annual.incidence <- with(ensemble.dat,
                         c(sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                           sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                           sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                           sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                           sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc)/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
Tab[6,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

write.csv(Tab, paste0(spath, "/reports/AdjEndpoints-SouthAfrica_", Sys.Date(), ".csv"))
# write.csv(Tab, "/Volumes/ahudson/Ensemble-stage2-correlates/SummaryTables/AdjEndpoints-SouthAfrica.csv")




########################################
### Event counts and event rate by group
########################################

Tab <- matrix(nrow = 6, ncol = 5)
colnames(Tab) <- c("Immediate Vacc", "Immedate Hybrid", "Delayed Vacc", "Delayed Hybrid", "Overall")
rownames(Tab) <- c("N", "Moderate to Severe-critical COVID (Total, Annual Incidence)",
                   "Moderate COVID (Total, Annual Incidence)",
                   "Severe-critical COVID (Total, Annual Incidence)",
                   "Asymptomatic Infection (Total, Annual Incidence)",
                   "Any Infection (Total, Annual Incidence)")
# N
Tab[1,] <- c(sum(ensemble.dat$ImmediateVaccine), sum(ensemble.dat$ImmediateHybrid),
             sum(ensemble.dat$DelayedVaccine), sum(ensemble.dat$DelayedHybrid), nrow(ensemble.dat))

# Moderate to Severe-critical COVID
event.total <- with(ensemble.dat[!(ensemble.dat$AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                   ensemble.dat$AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc)))
annual.incidence <- with(ensemble.dat[!(ensemble.dat$AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                        ensemble.dat$AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc > 90),],
                         c(sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                           sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                           sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                           sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                           sum(AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc)/
                             sum(AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
Tab[2,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")


# Moderate COVID
event.total <- with(ensemble.dat[!(ensemble.dat$AdjModEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                   ensemble.dat$AdjModEventTimePrimaryBD14_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc)))
annual.incidence <- with(ensemble.dat[!(ensemble.dat$AdjModEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                        ensemble.dat$AdjModEventTimePrimaryBD14_CensorAtOutsideVacc > 90),],
                         c(sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                           sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                           sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                           sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                           sum(AdjModEventIndPrimaryBD14_CensorAtOutsideVacc)/
                             sum(AdjModEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
Tab[3,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Severe-critical COVID
event.total <- with(ensemble.dat[!(ensemble.dat$AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                   ensemble.dat$AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc)))
annual.incidence <- with(ensemble.dat[!(ensemble.dat$AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                        ensemble.dat$AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc > 90),],
                         c(sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                           sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                           sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                           sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                           sum(AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc)/
                             sum(AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
Tab[4,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Asymp
event.total <- with(ensemble.dat[!(ensemble.dat$AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                   ensemble.dat$AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc)))
annual.incidence <- with(ensemble.dat[!(ensemble.dat$AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                        ensemble.dat$AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc > 90),],
                         c(sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
                             sum(AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                           sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
                             sum(AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                           sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
                             sum(AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                           sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
                             sum(AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                           sum(AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc)/
                             sum(AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
Tab[5,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Any infection
event.total <- with(ensemble.dat[!(ensemble.dat$AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                   ensemble.dat$AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc)))
annual.incidence <- with(ensemble.dat[!(ensemble.dat$AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc == 1 &
                                        ensemble.dat$AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc > 90),],
                         c(sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1])/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                           sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1])/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                           sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1])/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedVaccine == 1]),
                           sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1])/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc[DelayedHybrid == 1]),
                           sum(AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc)/
                             sum(AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc))) * 365
Tab[6,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# write.csv(Tab, "/Volumes/ahudson/Ensemble-stage2-correlates/SummaryTables/AdjEndpoints90Days-SouthAfrica.csv")
write.csv(Tab, paste0(spath, "/reports/AdjEndpoints90Days-SouthAfrica_", Sys.Date(), ".csv"))


###############################################
### Raw unadjusted Kaplan-Meier curves by group
###############################################

# pdf("/Volumes/ahudson/Ensemble-stage2-correlates/Summary Figures/AdjSurvivalCurves-SouthAfrica.pdf")
pdf(paste0(spath, "/reports/AdjSurvivalCurves-SouthAfrica_", Sys.Date(), ".pdf"))
fit <- survfit(Surv(event = AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc, time = AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Primary COVID", fun = "event",
           conf.int = TRUE, xlim = c(0, 210))

fit <- survfit(Surv(event = AdjModEventIndPrimaryBD14_CensorAtOutsideVacc, time = AdjModEventTimePrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Moderate COVID", fun = "event",
           conf.int = TRUE, xlim = c(0, 210))

fit <- survfit(Surv(event = AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc, time = AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Severe-critical COVID", fun = "event",
           conf.int = TRUE, xlim = c(0, 210))

fit <- survfit(Surv(event = AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc, time = AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Asymptomatic Infection", fun = "event",
           conf.int = TRUE, xlim = c(0, 210))

fit <- survfit(Surv(time = AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc, event = AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Any Infection", fun = "event",
           conf.int = TRUE, xlim = c(0, 210))
dev.off()

# pdf("/Volumes/ahudson/Ensemble-stage2-correlates/Summary Figures/AdjSurvivalCurves90Days-SouthAfrica.pdf")
pdf(paste0(spath, "/reports/AdjSurvivalCurves90Days-SouthAfrica_", Sys.Date(), ".pdf"))
fit <- survfit(Surv(event = AdjModtoSevEventIndPrimaryBD14_CensorAtOutsideVacc, time = AdjModtoSevEventTimePrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Primary COVID", fun = "event",
           conf.int = TRUE, xlim = c(0, 90))

fit <- survfit(Surv(event = AdjModEventIndPrimaryBD14_CensorAtOutsideVacc, time = AdjModEventTimePrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Moderate COVID", fun = "event",
           conf.int = TRUE, xlim = c(0, 90))

fit <- survfit(Surv(event = AdjSevEventIndPrimaryBD14_CensorAtOutsideVacc, time = AdjSevEventTimePrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Severe-critical COVID", fun = "event",
           conf.int = TRUE, xlim = c(0, 90))

fit <- survfit(Surv(event = AdjAsymEventIndPrimaryBD14_CensorAtOutsideVacc, time = AdjAsymEventTimePrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Asymptomatic Infection", fun = "event",
           conf.int = TRUE, xlim = c(0, 90))

fit <- survfit(Surv(time = AdjAnyEventTimePrimaryBD14_CensorAtOutsideVacc, event = AdjAnyEventIndPrimaryBD14_CensorAtOutsideVacc) ~ group,
               data = ensemble.dat)
ggsurvplot(fit, data = ensemble.dat, risk.table = TRUE, censor = FALSE, title = "Any Infection", fun = "event",
           conf.int = TRUE, xlim = c(0, 90))
dev.off()

##########################################
### Distribution of cases by calendar time
##########################################

# calendar time
ensemble.dat$CalendarDateDose2OLD <- Date(nrow(ensemble.dat)) 
for(i in 1:nrow(ensemble.dat)) {
  if(ensemble.dat$group[i] %in% c("Delayed Hybrid", "Delayed Vaccine")) {
    ensemble.dat$CalendarDateDose2OLD[i] <- as.Date(ensemble.dat$TRTSDT[i]) + ensemble.dat$NumberdaysD1toYear1booster[i]
  } else {
    ensemble.dat$CalendarDateDose2OLD[i] <- as.Date(ensemble.dat$CROSSDT[i]) + ensemble.dat$NumberdaysD1toYear1booster[i]
  }
}

for(i in 1:nrow(ensemble.dat)) {
  if(ensemble.dat$group[i] %in% c("Delayed Hybrid", "Delayed Vaccine")) {
    ensemble.dat$NumberdaysD1toBoosterNew[i] <- as.Date(ensemble.dat$BOOSTDT[i]) - as.Date(ensemble.dat$TRTSDT[i])
  } else {
    ensemble.dat$NumberdaysD1toBoosterNew[i] <- as.Date(ensemble.dat$BOOSTDT[i]) - as.Date(ensemble.dat$CROSSDT[i]) 
  }
}

ensemble.dat$CalendarDateDose2 <- as.Date(ensemble.dat$BOOSTDT)

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

# pdf("/Volumes/ahudson/Ensemble-stage2-correlates/Summary Figures/DateDose2-SouthAfrica.pdf")
pdf(paste0(spath, "/reports/DateDose2-SouthAfrica_", Sys.Date(), ".pdf"))
with(ensemble.dat, hist(CalendarDateDose2, breaks = "quarters", freq = TRUE, xlab = "Date of Booster Dose", 
                        main = "", col = "lightgrey"))
hist.hybrid <- with(ensemble.dat, hist(CalendarDateDose2[group %in% c("Delayed Hybrid", "Immediate Hybrid")], breaks = "quarters", plot = FALSE))
hist.vaccine <- with(ensemble.dat, hist(CalendarDateDose2[group %in% c("Delayed Vaccine", "Immediate Vaccine")], breaks = "quarters", plot = FALSE))
plot(hist.vaccine, add = TRUE, col = c1)
plot(hist.hybrid, add = TRUE, col = c2)

legend("topright", col = c("lightgrey", c1, c2), legend = c("overall", "vaccine", "hybrid"),
       pch = c(19, 19, 19))
dev.off()

