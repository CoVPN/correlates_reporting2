set.seed(206)
require(latex2exp)
require(vioplot)
require(lubridate)
if(!exists("spath")) {
  spath <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2"
}
source(paste0(spath, "/code/Methods-AverageTreatmentEffect.R"))
source(paste0(spath, "/code/Methods-NPCorrelates.R"))
if(!exists("data.path")) {
  data.path <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"
}
ensemble.dat <- read.csv(data.path)
ensemble.dat.immuno <- read.csv(data.path)

# ensemble.dat <- read.csv("/Volumes/trials/covpn/p3003/analysis/mapping_immune_correlates/post unblinded part B stage 2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv")
ensemble.dat <- ensemble.dat[ensemble.dat$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat$History %in% c("History a", "History b"),]
ensemble.dat <- ensemble.dat[ensemble.dat$BOOSTDT < as.Date("2022-05-22"),]
ensemble.dat$White <- with(ensemble.dat, ifelse(Black == 0 & Asian == 0 & NatAmer == 0 & PacIsl == 0 &
                                                  Multiracial == 0 & Notreported == 0 & Unknown == 0, 1, 0))
ensemble.dat$OtherRace <- with(ensemble.dat, ifelse(Notreported == 1 | Unknown == 1, 1, 0))

# ensemble.dat.immuno <- read.csv("/Volumes/trials/covpn/p3003/analysis/mapping_immune_correlates/post unblinded part B stage 2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv")
ensemble.dat.immuno <- ensemble.dat.immuno[ensemble.dat.immuno$BOOSTDT < as.Date("2022-05-22"),]
ensemble.dat.immuno <- ensemble.dat.immuno[ensemble.dat.immuno$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat.immuno$History %in% c("History a", "History b") &
                                             !(is.na(ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1) & is.na(ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.1)),]
ensemble.dat.immuno$White <- with(ensemble.dat.immuno, ifelse(Black == 0 & Asian == 0 & NatAmer == 0 & PacIsl == 0 &
                                                                Multiracial == 0 & Notreported == 0 & Unknown == 0, 1, 0))
ensemble.dat.immuno$OtherRace <- with(ensemble.dat.immuno, ifelse(Notreported == 1 | Unknown == 1, 1, 0))

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

ensemble.dat.immuno$DelayedVaccine <- ifelse(ensemble.dat.immuno$Cohort == "Cohort 1" & ensemble.dat.immuno$History == "History a", 1, 0)
ensemble.dat.immuno$ImmediateVaccine <- ifelse(ensemble.dat.immuno$Cohort == "Cohort 1" & ensemble.dat.immuno$History == "History b", 1, 0)
ensemble.dat.immuno$DelayedHybrid <- ifelse(ensemble.dat.immuno$Cohort == "Cohort 2" & ensemble.dat.immuno$History == "History a", 1, 0)
ensemble.dat.immuno$ImmediateHybrid <- ifelse(ensemble.dat.immuno$Cohort == "Cohort 2" & ensemble.dat.immuno$History == "History b", 1, 0)

ensemble.dat.immuno$group <- character(nrow(ensemble.dat.immuno))
ensemble.dat.immuno$group <- ifelse(ensemble.dat.immuno$ImmediateVaccine == 0,
                                    ensemble.dat.immuno$group, rep("Immediate Vaccine", nrow(ensemble.dat.immuno)))
ensemble.dat.immuno$group <- ifelse(ensemble.dat.immuno$DelayedVaccine == 0,
                                    ensemble.dat.immuno$group, rep("Delayed Vaccine", nrow(ensemble.dat.immuno)))
ensemble.dat.immuno$group <- ifelse(ensemble.dat.immuno$ImmediateHybrid == 0,
                                    ensemble.dat.immuno$group, rep("Immediate Hybrid", nrow(ensemble.dat.immuno)))
ensemble.dat.immuno$group <- ifelse(ensemble.dat.immuno$DelayedHybrid == 0,
                                    ensemble.dat.immuno$group, rep("Delayed Hybrid", nrow(ensemble.dat.immuno)))

### Event variables for correlates analyses

## day 28
D28.gap.immuno <- as.numeric(as.Date(ensemble.dat.immuno$SVSTDTC_YEAR.1_28.DAYS)) - 
  as.numeric(as.Date(ensemble.dat.immuno$BOOSTDT))
D28.gap <- as.numeric(as.Date(ensemble.dat$SVSTDTC_YEAR.1_28.DAYS)) - 
  as.numeric(as.Date(ensemble.dat$BOOSTDT))


# correlates cohort
ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc)


# per-protocol
ensemble.dat$AdjModtoSevEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)



## day 35

# correlates cohort
ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc)


# per-protocol
ensemble.dat$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjModEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)


# weights
ensemble.dat.immuno$wt.immuno <- numeric(nrow(ensemble.dat.immuno))

# Immediate Vaccine
for(i in which(ensemble.dat.immuno$group == "Immediate Vaccine")) {
  if(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Vaccine")
  } else if(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Vaccine")
  } else if(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Vaccine")
  } else {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
                                              ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc >= 1 &
                                              ensemble.dat$group == "Immediate Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
            ensemble.dat.immuno$group == "Immediate Vaccine")
  }
}
# Delayed Vaccine
for(i in which(ensemble.dat.immuno$group == "Delayed Vaccine")) {
  if(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Vaccine")
  } else if(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Vaccine")
  } else if(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Vaccine")
  } else {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
                                              ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc >= 1 &
                                              ensemble.dat$group == "Delayed Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
            ensemble.dat.immuno$group == "Delayed Vaccine")
  }
}
# Immediate Hybrid
for(i in which(ensemble.dat.immuno$group == "Immediate Hybrid")) {
  if(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Hybrid")
  } else if(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Hybrid")
  } else if(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Hybrid")
  } else {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
                                              ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc >= 1 &
                                              ensemble.dat$group == "Immediate Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
            ensemble.dat.immuno$group == "Immediate Hybrid")
  }
}
# Delayed Hybrid
for(i in which(ensemble.dat.immuno$group == "Delayed Hybrid")) {
  
  if(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Hybrid")
  } else if(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Hybrid")
  } else if(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Hybrid")
  } else {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
                                              ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc >= 1 &
                                              ensemble.dat$group == "Delayed Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
            ensemble.dat.immuno$group == "Delayed Hybrid")
  }
}

ensemble.dat.immuno$wt.immuno <- ensemble.dat.immuno$wt.immuno/mean(ensemble.dat.immuno$wt.immuno)

## covariates
ensemble.dat.immuno$BMIDich <- ifelse(ensemble.dat.immuno$BMI %in% c(1, 2), 0, 1)
ensemble.dat.immuno$BMIDich[is.na(ensemble.dat.immuno$BMI)] <- NA

ensemble.dat.immuno$RegionUS <- ifelse(ensemble.dat.immuno$Region == 0, 1, 0)
ensemble.dat.immuno$RegionLA <- ifelse(ensemble.dat.immuno$Region == 1, 1, 0)
ensemble.dat.immuno$RegionRSA <- ifelse(ensemble.dat.immuno$Region == 2, 1, 0)

ensemble.dat.immuno$HistoryA <- ifelse(ensemble.dat.immuno$History == "History a", 1, 0)
ensemble.dat.immuno$CalendarDateDose2 <- as.Date(ensemble.dat.immuno$BOOSTDT)
ensemble.dat.immuno$CalendarDateDose2Early <- ifelse(ensemble.dat.immuno$CalendarDateDose2 < as.Date("2021-12-01"), 1, 0)

################## NP Correlates BA-45

t0 <- c(39, 53)
t0.sev <- c(24, 38)
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

########### BA 4.5; Hybrid;
X <-  ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.4.5
X.eval <- quantile(X, seq(.01, .99, .01), na.rm = TRUE)
A <- with(ensemble.dat.immuno, ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
W <- with(ensemble.dat.immuno, cbind(Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex)))

# W <- with(ensemble.dat.immuno, cbind(Age, BMIDich * Sex, BMIDich * (1 - Sex),
#                                      RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))
# W <- with(ensemble.dat.immuno, cbind(RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))
# W <- with(ensemble.dat.immuno, cbind(HistoryA, CalendarDateDose2Early, Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex)))

sweights <- ensemble.dat.immuno$wt.immuno
complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))})

hyb.complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 1
# hyb.complete.sev <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 1 & 
#                     ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
hyb.complete.sev <- apply(cbind(A, X), 1, function(x){all(!is.na(x))}) & A == 1 & 
  ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
hyb.complete.modtosev <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 1 & 
                         ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
hyb.complete.any <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 1 & 
                    ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 0

hyb.time.modtosev <- ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[hyb.complete.modtosev]
hyb.event.modtosev <- ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[hyb.complete.modtosev]
hyb.time.sev <- ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[hyb.complete.sev]
hyb.event.sev <- ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[hyb.complete.sev]
hyb.time.any <- ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[hyb.complete.any]
hyb.event.any <- ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[hyb.complete.any]

# hyb.prop.score.fit <- EstimateConditionalDensity(a = X[hyb.complete], w = W[hyb.complete,], weights = sweights[hyb.complete])
hyb.prop.score.fit.modtosev <- EstimateConditionalDensity(a = X[hyb.complete.modtosev], w = W[hyb.complete.modtosev,], weights = sweights[hyb.complete.modtosev])
# hyb.prop.score.fit.sev <- EstimateConditionalDensity(a = X[hyb.complete.sev], w = W[hyb.complete.sev,], weights = sweights[hyb.complete.sev])
hyb.prop.score.fit.any <- EstimateConditionalDensity(a = X[hyb.complete.any], w = W[hyb.complete.any,], weights = sweights[hyb.complete.any])

hyb.modtosev.np.surv.fit <- EstimateConditionalSurvival(time = hyb.time.modtosev, event = hyb.event.modtosev, 
                                                        a = X[hyb.complete.modtosev], w = W[hyb.complete.modtosev,],
                                                        weights = sweights[hyb.complete.modtosev], t0 = t0)
hyb.any.np.surv.fit <- EstimateConditionalSurvival(time = hyb.time.any, event = hyb.event.any, 
                                                   a = X[hyb.complete.any], w = W[hyb.complete.any,],
                                                   weights = sweights[hyb.complete.any], t0 = t0)
hyb.sev.np.surv.fit <- EstimateConditionalSurvival(time = hyb.time.sev, event = hyb.event.sev, 
                                                   a = X[hyb.complete.sev], w = NULL, #w = W[hyb.complete.sev,],
                                                   weights = sweights[hyb.complete.sev], t0 = t0.sev,
                                                   nonparametric = FALSE)


hyb.modtosev.ksmooth.fit <-KernelSmoothFit(time = hyb.time.modtosev, event = hyb.event.modtosev, 
                                           a = X[hyb.complete.modtosev], a.eval = X.eval, w = W[hyb.complete.modtosev,], 
                                           weights = sweights[hyb.complete.modtosev],
                                           t0 = t0, 
                                           prop.score = hyb.prop.score.fit.modtosev$prop.score,
                                           fitted.surv = hyb.modtosev.np.surv.fit$fitted.surv,
                                           fitted.cens = hyb.modtosev.np.surv.fit$fitted.cens,
                                           cond.surv = hyb.modtosev.np.surv.fit$cond.surv,
                                           cond.cens = hyb.modtosev.np.surv.fit$cond.cens,
                                           fitted.haz.event = hyb.modtosev.np.surv.fit$fitted.haz.event,
                                           fitted.haz.cens = hyb.modtosev.np.surv.fit$fitted.haz.cens,
                                           alpha = .05)
hyb.modtosev.test <- MonotoneTest(time = hyb.time.modtosev, event = hyb.event.modtosev, 
                                  a = X[hyb.complete.modtosev], w = W[hyb.complete.modtosev,], 
                                  weights = sweights[hyb.complete.modtosev],
                                  t0 = t0, 
                                  prop.score = hyb.prop.score.fit.modtosev$prop.score,
                                  fitted.surv = hyb.modtosev.np.surv.fit$fitted.surv,
                                  fitted.cens = hyb.modtosev.np.surv.fit$fitted.cens,
                                  cond.surv = hyb.modtosev.np.surv.fit$cond.surv,
                                  cond.cens = hyb.modtosev.np.surv.fit$cond.cens,
                                  fitted.haz.event = hyb.modtosev.np.surv.fit$fitted.haz.event,
                                  fitted.haz.cens = hyb.modtosev.np.surv.fit$fitted.haz.cens)


hyb.any.ksmooth.fit <-KernelSmoothFit(time = hyb.time.any, event = hyb.event.any, 
                                      a = X[hyb.complete.any], a.eval = X.eval, w = W[hyb.complete.any,], 
                                      weights = sweights[hyb.complete.any],
                                      t0 = t0,
                                      prop.score = hyb.prop.score.fit.any$prop.score,
                                      fitted.surv = hyb.any.np.surv.fit$fitted.surv,
                                      fitted.cens = hyb.any.np.surv.fit$fitted.cens,
                                      cond.surv = hyb.any.np.surv.fit$cond.surv,
                                      cond.cens = hyb.any.np.surv.fit$cond.cens,
                                      fitted.haz.event = hyb.any.np.surv.fit$fitted.haz.event,
                                      fitted.haz.cens = hyb.any.np.surv.fit$fitted.haz.cens,
                                      alpha = .05)
hyb.any.test <- MonotoneTest(time = hyb.time.any, event = hyb.event.any, 
                             a = X[hyb.complete.any], w = W[hyb.complete.any,], 
                             weights = sweights[hyb.complete.any],
                             t0 = t0, 
                             prop.score = hyb.prop.score.fit.any$prop.score,
                             fitted.surv = hyb.any.np.surv.fit$fitted.surv,
                             fitted.cens = hyb.any.np.surv.fit$fitted.cens,
                             cond.surv = hyb.any.np.surv.fit$cond.surv,
                             cond.cens = hyb.any.np.surv.fit$cond.cens,
                             fitted.haz.event = hyb.any.np.surv.fit$fitted.haz.event,
                             fitted.haz.cens = hyb.any.np.surv.fit$fitted.haz.cens)

hyb.sev.ksmooth.fit <-KernelSmoothFit(time = hyb.time.sev, event = hyb.event.sev, 
                                      a = X[hyb.complete.sev], a.eval = X.eval, w = NULL, #W[hyb.complete.sev,], 
                                      weights = sweights[hyb.complete.sev],
                                      t0 = t0.sev,
                                      prop.score = rep(1, length(hyb.time.sev)), #hyb.prop.score.fit.sev$prop.score,
                                      fitted.surv = hyb.sev.np.surv.fit$fitted.surv,
                                      fitted.cens = hyb.sev.np.surv.fit$fitted.cens,
                                      cond.surv = hyb.sev.np.surv.fit$cond.surv,
                                      cond.cens = hyb.sev.np.surv.fit$cond.cens,
                                      fitted.haz.event = hyb.sev.np.surv.fit$fitted.haz.event,
                                      fitted.haz.cens = hyb.sev.np.surv.fit$fitted.haz.cens,
                                      alpha = .05)
hyb.sev.test <- MonotoneTest(time = hyb.time.sev, event = hyb.event.sev, 
                             a = X[hyb.complete.sev], w = NULL, #W[hyb.complete.sev,], 
                             weights = sweights[hyb.complete.sev],
                             t0 = t0.sev, 
                             prop.score = rep(1, length(hyb.time.sev)),#hyb.prop.score.fit.sev$prop.score,
                             fitted.surv = hyb.sev.np.surv.fit$fitted.surv,
                             fitted.cens = hyb.sev.np.surv.fit$fitted.cens,
                             cond.surv = hyb.sev.np.surv.fit$cond.surv,
                             cond.cens = hyb.sev.np.surv.fit$cond.cens,
                             fitted.haz.event = hyb.sev.np.surv.fit$fitted.haz.event,
                             fitted.haz.cens = hyb.sev.np.surv.fit$fitted.haz.cens)

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Corr-Survival-BA45Hybrid.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Corr-Survival-BA45Hybrid.pdf"), width = 9, height = 7)

par(mar = c(5,4,4,3) + .1)
X.dist <- hist(X[hyb.complete], plot = FALSE)
num.breaks <- length(X.dist$breaks)
dens.scale <- max(X.dist$density) * 10

# primary
plot(hyb.modtosev.ksmooth.fit$a.eval, hyb.modtosev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Hybrid; Moderate to Severe-critical COVID",
     xlim = c(min(X[hyb.complete]), max(X[hyb.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(hyb.modtosev.ksmooth.fit$a.eval, rev(hyb.modtosev.ksmooth.fit$a.eval)),
        y = c(1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2], rev(1-hyb.modtosev.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
legend("topright", paste0("p-value: ", signif(hyb.modtosev.test$p.vals[2])), cex = .75)

# any infection
plot(hyb.any.ksmooth.fit$a.eval, hyb.any.ksmooth.fit$drf.est[,2], ylim = c(0, max(1-hyb.any.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Infection, Any Symptomatology; Hybrid",
     xlim = c(min(X[hyb.complete]), max(X[hyb.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-hyb.any.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(hyb.any.ksmooth.fit$a.eval, rev(hyb.any.ksmooth.fit$a.eval)),
        y = c(1-hyb.any.ksmooth.fit$drf.cb.lower[,2], rev(1-hyb.any.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
legend("topright", paste0("p-value: ", signif(hyb.any.test$p.vals[2])), cex = .75)

# severe
plot(hyb.sev.ksmooth.fit$a.eval, hyb.sev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1-hyb.sev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 45 Days after BD-28",
     main = "Hybrid; Severe-critical COVID",
     xlim = c(min(X[hyb.complete]), max(X[hyb.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-hyb.sev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(hyb.sev.ksmooth.fit$a.eval, rev(hyb.sev.ksmooth.fit$a.eval)),
        y = c(1-hyb.sev.ksmooth.fit$drf.cb.lower[,2], rev(1-hyb.sev.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
legend("topright", paste0("p-value: ", signif(hyb.sev.test$p.vals[2])), cex = .75)




dev.off()



########### BA 4.5; Vaccine
X <-  ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.4.5
A <- with(ensemble.dat.immuno, ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
W <- with(ensemble.dat.immuno, cbind(Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex)))
# W <- with(ensemble.dat.immuno, cbind(RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))
sweights <- ensemble.dat.immuno$wt.immuno
complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))})

vacc.complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 0
# vacc.complete.sev <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 0 & 
#   ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
vacc.complete.sev <- apply(cbind(A, X), 1, function(x){all(!is.na(x))}) & A == 0 & 
  ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
vacc.complete.modtosev <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 0 & 
  ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
vacc.complete.any <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 0 & 
  ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 0

vacc.time.modtosev <- ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[vacc.complete.modtosev]
vacc.event.modtosev <- ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[vacc.complete.modtosev]
vacc.time.sev <- ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[vacc.complete.sev]
vacc.event.sev <- ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[vacc.complete.sev]
vacc.time.any <- ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[vacc.complete.any]
vacc.event.any <- ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[vacc.complete.any]

# vacc.prop.score.fit <- EstimateConditionalDensity(a = X[vacc.complete], w = W[vacc.complete,], weights = sweights[vacc.complete])
vacc.prop.score.fit.modtosev <- EstimateConditionalDensity(a = X[vacc.complete.modtosev], w = W[vacc.complete.modtosev,], weights = sweights[vacc.complete.modtosev])
# vacc.prop.score.fit.sev <- EstimateConditionalDensity(a = X[vacc.complete.sev], w = W[vacc.complete.sev,], weights = sweights[vacc.complete.sev])
vacc.prop.score.fit.any <- EstimateConditionalDensity(a = X[vacc.complete.any], w = W[vacc.complete.any,], weights = sweights[vacc.complete.any])

vacc.modtosev.np.surv.fit <- EstimateConditionalSurvival(time = vacc.time.modtosev, event = vacc.event.modtosev, 
                                                        a = X[vacc.complete.modtosev], w = W[vacc.complete.modtosev,],
                                                        weights = sweights[vacc.complete.modtosev], t0 = t0)
vacc.any.np.surv.fit <- EstimateConditionalSurvival(time = vacc.time.any, event = vacc.event.any, 
                                                   a = X[vacc.complete.any], w = W[vacc.complete.any,],
                                                   weights = sweights[vacc.complete.any], t0 = t0)
vacc.sev.np.surv.fit <- EstimateConditionalSurvival(time = vacc.time.sev, event = vacc.event.sev, 
                                                   a = X[vacc.complete.sev], w = NULL, #W[vacc.complete.sev,],
                                                   weights = sweights[vacc.complete.sev], t0 = t0.sev,
                                                   nonparametric = FALSE)


vacc.modtosev.ksmooth.fit <-KernelSmoothFit(time = vacc.time.modtosev, event = vacc.event.modtosev, 
                                            a = X[vacc.complete.modtosev], a.eval = X.eval, w = W[vacc.complete.modtosev,], 
                                            weights = sweights[vacc.complete.modtosev],
                                            t0 = t0, 
                                            prop.score = vacc.prop.score.fit.modtosev$prop.score,
                                            fitted.surv = vacc.modtosev.np.surv.fit$fitted.surv,
                                            fitted.cens = vacc.modtosev.np.surv.fit$fitted.cens,
                                            cond.surv = vacc.modtosev.np.surv.fit$cond.surv,
                                            cond.cens = vacc.modtosev.np.surv.fit$cond.cens,
                                            fitted.haz.event = vacc.modtosev.np.surv.fit$fitted.haz.event,
                                            fitted.haz.cens = vacc.modtosev.np.surv.fit$fitted.haz.cens,
                                            alpha = .05)
vacc.modtosev.test <- MonotoneTest(time = vacc.time.modtosev, event = vacc.event.modtosev, 
                                  a = X[vacc.complete.modtosev], w = W[vacc.complete.modtosev,], 
                                  weights = sweights[vacc.complete.modtosev],
                                  t0 = t0, 
                                  prop.score = vacc.prop.score.fit.modtosev$prop.score,
                                  fitted.surv = vacc.modtosev.np.surv.fit$fitted.surv,
                                  fitted.cens = vacc.modtosev.np.surv.fit$fitted.cens,
                                  cond.surv = vacc.modtosev.np.surv.fit$cond.surv,
                                  cond.cens = vacc.modtosev.np.surv.fit$cond.cens,
                                  fitted.haz.event = vacc.modtosev.np.surv.fit$fitted.haz.event,
                                  fitted.haz.cens = vacc.modtosev.np.surv.fit$fitted.haz.cens)


vacc.any.ksmooth.fit <-KernelSmoothFit(time = vacc.time.any, event = vacc.event.any, 
                                      a = X[vacc.complete.any], a.eval = X.eval, w = W[vacc.complete.any,], 
                                      weights = sweights[vacc.complete.any],
                                      t0 = t0,
                                      prop.score = vacc.prop.score.fit.any$prop.score,
                                      fitted.surv = vacc.any.np.surv.fit$fitted.surv,
                                      fitted.cens = vacc.any.np.surv.fit$fitted.cens,
                                      cond.surv = vacc.any.np.surv.fit$cond.surv,
                                      cond.cens = vacc.any.np.surv.fit$cond.cens,
                                      fitted.haz.event = vacc.any.np.surv.fit$fitted.haz.event,
                                      fitted.haz.cens = vacc.any.np.surv.fit$fitted.haz.cens,
                                      alpha = .05)
vacc.any.test <- MonotoneTest(time = vacc.time.any, event = vacc.event.any, 
                             a = X[vacc.complete.any], w = W[vacc.complete.any,], 
                             weights = sweights[vacc.complete.any],
                             t0 = t0, 
                             prop.score = vacc.prop.score.fit.any$prop.score,
                             fitted.surv = vacc.any.np.surv.fit$fitted.surv,
                             fitted.cens = vacc.any.np.surv.fit$fitted.cens,
                             cond.surv = vacc.any.np.surv.fit$cond.surv,
                             cond.cens = vacc.any.np.surv.fit$cond.cens,
                             fitted.haz.event = vacc.any.np.surv.fit$fitted.haz.event,
                             fitted.haz.cens = vacc.any.np.surv.fit$fitted.haz.cens)

vacc.sev.ksmooth.fit <-KernelSmoothFit(time = vacc.time.sev, event = vacc.event.sev, 
                                      a = X[vacc.complete.sev], a.eval = X.eval, w = NULL,#W[vacc.complete.sev,], 
                                      weights = sweights[vacc.complete.sev],
                                      t0 = t0.sev,
                                      prop.score = rep(1, length(vacc.time.sev)), #vacc.prop.score.fit.sev$prop.score,
                                      fitted.surv = vacc.sev.np.surv.fit$fitted.surv,
                                      fitted.cens = vacc.sev.np.surv.fit$fitted.cens,
                                      cond.surv = vacc.sev.np.surv.fit$cond.surv,
                                      cond.cens = vacc.sev.np.surv.fit$cond.cens,
                                      fitted.haz.event = vacc.sev.np.surv.fit$fitted.haz.event,
                                      fitted.haz.cens = vacc.sev.np.surv.fit$fitted.haz.cens,
                                      alpha = .05)
vacc.sev.test <- MonotoneTest(time = vacc.time.sev, event = vacc.event.sev, 
                             a = X[vacc.complete.sev], w = NULL,#W[vacc.complete.sev,], 
                             weights = sweights[vacc.complete.sev],
                             t0 = t0.sev, 
                             prop.score = rep(1, length(vacc.time.sev)),#vacc.prop.score.fit.sev$prop.score,
                             fitted.surv = vacc.sev.np.surv.fit$fitted.surv,
                             fitted.cens = vacc.sev.np.surv.fit$fitted.cens,
                             cond.surv = vacc.sev.np.surv.fit$cond.surv,
                             cond.cens = vacc.sev.np.surv.fit$cond.cens,
                             fitted.haz.event = vacc.sev.np.surv.fit$fitted.haz.event,
                             fitted.haz.cens = vacc.sev.np.surv.fit$fitted.haz.cens)

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Corr-Survival-BA45vaccine.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Corr-Survival-BA45vaccine.pdf"), width = 9, height = 7)


par(mar = c(5,4,4,3) + .1)

X.dist <- hist(X[vacc.complete], plot = FALSE)
num.breaks <- length(X.dist$breaks)
dens.scale <- max(X.dist$density) * 10

#primary
plot(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Moderate to Severe-critical COVID",
     xlim = c(min(X[complete]), max(X[complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-vacc.modtosev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(vacc.modtosev.ksmooth.fit$a.eval, rev(vacc.modtosev.ksmooth.fit$a.eval)),
        y = c(1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2], rev(1 - vacc.modtosev.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
legend("topright", paste0("p-value: ", signif(vacc.modtosev.test$p.vals[2])), cex = .75)

#any
plot(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.any.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Infection, Any Symptomatology",
     xlim = c(min(X[complete]), max(X[complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-vacc.any.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(vacc.any.ksmooth.fit$a.eval, rev(vacc.any.ksmooth.fit$a.eval)),
        y = c(1 - vacc.any.ksmooth.fit$drf.cb.lower[,2], rev(1 - vacc.any.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
# legend("bottomright", paste0("p-value: ", signif(test$p.vals[2])), cex = .75)
legend("topright", paste0("p-value: ", signif(vacc.any.test$p.vals[2])), cex = .75)

#severe
plot(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 45 Days after BD-28",
     main = "Severe-critical COVID",
     xlim = c(min(X[complete]), max(X[complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-vacc.sev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(vacc.sev.ksmooth.fit$a.eval, rev(vacc.sev.ksmooth.fit$a.eval)),
        y = c(1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2], rev(1 - vacc.sev.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
# legend("bottomright", paste0("p-value: ", signif(test$p.vals[2])), cex = .75)
legend("topright", paste0("p-value: ", signif(vacc.sev.test$p.vals[2])), cex = .75)

dev.off()



##### Simultaneous figures
# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Corr-Survival-BA45.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Corr-Survival-BA45.pdf"), width = 9, height = 7)
par(mar = c(5,4,4,3) + .1)
vacc.X.dist <- hist(X[vacc.complete], plot = FALSE)
vacc.num.breaks <- length(vacc.X.dist$breaks)
vacc.dens.scale <- max(vacc.X.dist$density) * 10
hyb.X.dist <- hist(X[hyb.complete], plot = FALSE)
hyb.num.breaks <- length(hyb.X.dist$breaks)
hyb.dens.scale <- max(hyb.X.dist$density) * 10

# modtosev
plot(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Moderate to Severe-critical COVID",
     xlim = c(min(X[vacc.complete]), max(X[vacc.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(c(vacc.X.dist$density, hyb.X.dist$density)) * 
  1/max(1-vacc.modtosev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(vacc.num.breaks-1)) {
  polygon(x = c(vacc.X.dist$breaks[i], vacc.X.dist$breaks[i], 
                vacc.X.dist$breaks[i+1], vacc.X.dist$breaks[i+1]), 
          y = c(-dens.scale, vacc.X.dist$density[i], vacc.X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
for(i in 1:(hyb.num.breaks-1)) {
  polygon(x = c(hyb.X.dist$breaks[i], hyb.X.dist$breaks[i], 
                hyb.X.dist$breaks[i+1], hyb.X.dist$breaks[i+1]), 
          y = c(-dens.scale, hyb.X.dist$density[i], hyb.X.dist$density[i], -dens.scale)/dens.scale,
          col = c1, border = c1) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "salmon")
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "salmon")
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "salmon")

lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")

legend("topright", 
       legend = c(paste0("Hybrid (p-value: ", signif(hyb.modtosev.test$p.vals[2]), ")"),
                  paste0("Vaccine (p-value: ", signif(vacc.modtosev.test$p.vals[2]), ")")),
       lty = c(1, 1), lwd = c(2,2),
       col = c("blue", "salmon"),
       cex = .75)

# any
plot(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.any.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Infection, Any Symptomatology",
     xlim = c(min(X[vacc.complete]), max(X[vacc.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(c(vacc.X.dist$density, hyb.X.dist$density)) * 
  1/max(1-vacc.any.ksmooth.fit$drf.cb.lower[,2])
vacc.dens.scale <- max(X.dist$density) * 1/max(1-vacc.any.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(vacc.num.breaks-1)) {
  polygon(x = c(vacc.X.dist$breaks[i], vacc.X.dist$breaks[i], 
                vacc.X.dist$breaks[i+1], vacc.X.dist$breaks[i+1]), 
          y = c(-dens.scale, vacc.X.dist$density[i], vacc.X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
for(i in 1:(hyb.num.breaks-1)) {
  polygon(x = c(hyb.X.dist$breaks[i], hyb.X.dist$breaks[i], 
                hyb.X.dist$breaks[i+1], hyb.X.dist$breaks[i+1]), 
          y = c(-dens.scale, hyb.X.dist$density[i], hyb.X.dist$density[i], -dens.scale)/dens.scale,
          col = c1, border = c1) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "salmon")
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "salmon")
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "salmon")

lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")

legend("topright", 
       legend = c(paste0("Hybrid (p-value: ", signif(hyb.any.test$p.vals[2]), ")"),
                  paste0("Vaccine (p-value: ", signif(vacc.any.test$p.vals[2]), ")")),
       lty = c(1, 1), lwd = c(2,2),
       col = c("blue", "salmon"),
       cex = .75)

# sev
plot(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 BA.4/5 titer (AU/mL)", ylab = "Cumulative Incidence 45 Days after BD-28",
     main = "Severe-critical COVID",
     xlim = c(min(X[vacc.complete]), max(X[vacc.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(c(vacc.X.dist$density, hyb.X.dist$density)) * 
  1/max(1-vacc.sev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(vacc.num.breaks-1)) {
  polygon(x = c(vacc.X.dist$breaks[i], vacc.X.dist$breaks[i], 
                vacc.X.dist$breaks[i+1], vacc.X.dist$breaks[i+1]), 
          y = c(-dens.scale, vacc.X.dist$density[i], vacc.X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
for(i in 1:(hyb.num.breaks-1)) {
  polygon(x = c(hyb.X.dist$breaks[i], hyb.X.dist$breaks[i], 
                hyb.X.dist$breaks[i+1], hyb.X.dist$breaks[i+1]), 
          y = c(-dens.scale, hyb.X.dist$density[i], hyb.X.dist$density[i], -dens.scale)/dens.scale,
          col = c1, border = c1) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "salmon")
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "salmon")
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "salmon")

lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")

legend("topright", 
       legend = c(paste0("Hybrid (p-value: ", signif(hyb.sev.test$p.vals[2]), ")"),
                  paste0("Vaccine (p-value: ", signif(vacc.sev.test$p.vals[2]), ")")),
       lty = c(1, 1), lwd = c(2,2),
       col = c("blue", "salmon"),
       cex = .75)

dev.off()











################## NP Correlates B1

t0 <- c(39, 53)
t0.sev <- c(24, 38)
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

########### B1; Hybrid;
X <-  ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1
X.eval <- quantile(X, seq(.01, .99, .01), na.rm = TRUE)
A <- with(ensemble.dat.immuno, ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
W <- with(ensemble.dat.immuno, cbind(Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex)))
W.sev <- with(ensemble.dat.immuno, cbind(Age, Sex))
W.sev <- matrix(ensemble.dat.immuno$Age, ncol = 1)
# W <- with(ensemble.dat.immuno, cbind(RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))
sweights <- ensemble.dat.immuno$wt.immuno
complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))})

hyb.complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 1
# hyb.complete.sev <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 1 & 
#   ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
hyb.complete.sev <- apply(cbind(A, X), 1, function(x){all(!is.na(x))}) & A == 1 & 
  ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
hyb.complete.modtosev <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 1 & 
  ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
hyb.complete.any <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 1 & 
  ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 0

hyb.time.modtosev <- ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[hyb.complete.modtosev]
hyb.event.modtosev <- ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[hyb.complete.modtosev]
hyb.time.sev <- ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[hyb.complete.sev]
hyb.event.sev <- ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[hyb.complete.sev]
hyb.time.any <- ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[hyb.complete.any]
hyb.event.any <- ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[hyb.complete.any]

# hyb.prop.score.fit <- EstimateConditionalDensity(a = X[hyb.complete], w = W[hyb.complete,], weights = sweights[hyb.complete])
hyb.prop.score.fit.modtosev <- EstimateConditionalDensity(a = X[hyb.complete.modtosev], w = W[hyb.complete.modtosev,], weights = sweights[hyb.complete.modtosev])
# hyb.prop.score.fit.sev <- EstimateConditionalDensity(a = X[hyb.complete.sev], w = W.sev[hyb.complete.sev,], weights = sweights[hyb.complete.sev])
hyb.prop.score.fit.any <- EstimateConditionalDensity(a = X[hyb.complete.any], w = W[hyb.complete.any,], weights = sweights[hyb.complete.any])

hyb.modtosev.np.surv.fit <- EstimateConditionalSurvival(time = hyb.time.modtosev, event = hyb.event.modtosev, 
                                                        a = X[hyb.complete.modtosev], w = W[hyb.complete.modtosev,],
                                                        weights = sweights[hyb.complete.modtosev], t0 = t0)
hyb.any.np.surv.fit <- EstimateConditionalSurvival(time = hyb.time.any, event = hyb.event.any, 
                                                   a = X[hyb.complete.any], w = W[hyb.complete.any,],
                                                   weights = sweights[hyb.complete.any], t0 = t0)
hyb.sev.np.surv.fit <- EstimateConditionalSurvival(time = hyb.time.sev, event = hyb.event.sev, 
                                                   a = X[hyb.complete.sev], w = NULL, #= W[hyb.complete.sev,1:2],
                                                   weights = sweights[hyb.complete.sev], t0 = t0.sev,
                                                   nonparametric = FALSE)


hyb.modtosev.ksmooth.fit <-KernelSmoothFit(time = hyb.time.modtosev, event = hyb.event.modtosev, 
                                           a = X[hyb.complete.modtosev], a.eval = X.eval, w = W[hyb.complete.modtosev,], 
                                           weights = sweights[hyb.complete.modtosev],
                                           t0 = t0, 
                                           prop.score = hyb.prop.score.fit.modtosev$prop.score,
                                           fitted.surv = hyb.modtosev.np.surv.fit$fitted.surv,
                                           fitted.cens = hyb.modtosev.np.surv.fit$fitted.cens,
                                           cond.surv = hyb.modtosev.np.surv.fit$cond.surv,
                                           cond.cens = hyb.modtosev.np.surv.fit$cond.cens,
                                           fitted.haz.event = hyb.modtosev.np.surv.fit$fitted.haz.event,
                                           fitted.haz.cens = hyb.modtosev.np.surv.fit$fitted.haz.cens,
                                           alpha = .05)
hyb.modtosev.test <- MonotoneTest(time = hyb.time.modtosev, event = hyb.event.modtosev, 
                                  a = X[hyb.complete.modtosev], w = W[hyb.complete.modtosev,], 
                                  weights = sweights[hyb.complete.modtosev],
                                  t0 = t0, 
                                  prop.score = hyb.prop.score.fit.modtosev$prop.score,
                                  fitted.surv = hyb.modtosev.np.surv.fit$fitted.surv,
                                  fitted.cens = hyb.modtosev.np.surv.fit$fitted.cens,
                                  cond.surv = hyb.modtosev.np.surv.fit$cond.surv,
                                  cond.cens = hyb.modtosev.np.surv.fit$cond.cens,
                                  fitted.haz.event = hyb.modtosev.np.surv.fit$fitted.haz.event,
                                  fitted.haz.cens = hyb.modtosev.np.surv.fit$fitted.haz.cens)


hyb.any.ksmooth.fit <-KernelSmoothFit(time = hyb.time.any, event = hyb.event.any, 
                                      a = X[hyb.complete.any], a.eval = X.eval, w = W[hyb.complete.any,], 
                                      weights = sweights[hyb.complete.any],
                                      t0 = t0,
                                      prop.score = hyb.prop.score.fit.any$prop.score,
                                      fitted.surv = hyb.any.np.surv.fit$fitted.surv,
                                      fitted.cens = hyb.any.np.surv.fit$fitted.cens,
                                      cond.surv = hyb.any.np.surv.fit$cond.surv,
                                      cond.cens = hyb.any.np.surv.fit$cond.cens,
                                      fitted.haz.event = hyb.any.np.surv.fit$fitted.haz.event,
                                      fitted.haz.cens = hyb.any.np.surv.fit$fitted.haz.cens,
                                      alpha = .05)
hyb.any.test <- MonotoneTest(time = hyb.time.any, event = hyb.event.any, 
                             a = X[hyb.complete.any], w = W[hyb.complete.any,], 
                             weights = sweights[hyb.complete.any],
                             t0 = t0, 
                             prop.score = hyb.prop.score.fit.any$prop.score,
                             fitted.surv = hyb.any.np.surv.fit$fitted.surv,
                             fitted.cens = hyb.any.np.surv.fit$fitted.cens,
                             cond.surv = hyb.any.np.surv.fit$cond.surv,
                             cond.cens = hyb.any.np.surv.fit$cond.cens,
                             fitted.haz.event = hyb.any.np.surv.fit$fitted.haz.event,
                             fitted.haz.cens = hyb.any.np.surv.fit$fitted.haz.cens)

hyb.sev.ksmooth.fit <-KernelSmoothFit(time = hyb.time.sev, event = hyb.event.sev, 
                                      a = X[hyb.complete.sev], a.eval = X.eval, w = NULL, #W.sev[hyb.complete.sev,], 
                                      weights = sweights[hyb.complete.sev],
                                      t0 = t0.sev,
                                      prop.score = rep(1, length(hyb.time.sev)), #hyb.prop.score.fit.sev$prop.score,
                                      fitted.surv = hyb.sev.np.surv.fit$fitted.surv,
                                      fitted.cens = hyb.sev.np.surv.fit$fitted.cens,
                                      cond.surv = hyb.sev.np.surv.fit$cond.surv,
                                      cond.cens = hyb.sev.np.surv.fit$cond.cens,
                                      fitted.haz.event = hyb.sev.np.surv.fit$fitted.haz.event,
                                      fitted.haz.cens = hyb.sev.np.surv.fit$fitted.haz.cens,
                                      alpha = .05)
hyb.sev.test <- MonotoneTest(time = hyb.time.sev, event = hyb.event.sev, 
                             a = X[hyb.complete.sev], w = NULL, #W[hyb.complete.sev,], 
                             weights = sweights[hyb.complete.sev],
                             t0 = t0.sev, 
                             prop.score = rep(1, length(hyb.time.sev)), #hyb.prop.score.fit.sev$prop.score,
                             fitted.surv = hyb.sev.np.surv.fit$fitted.surv,
                             fitted.cens = hyb.sev.np.surv.fit$fitted.cens,
                             cond.surv = hyb.sev.np.surv.fit$cond.surv,
                             cond.cens = hyb.sev.np.surv.fit$cond.cens,
                             fitted.haz.event = hyb.sev.np.surv.fit$fitted.haz.event,
                             fitted.haz.cens = hyb.sev.np.surv.fit$fitted.haz.cens)

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Corr-Survival-B1Hybrid.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Corr-Survival-B1Hybrid.pdf"), width = 9, height = 7)

par(mar = c(5,4,4,3) + .1)
X.dist <- hist(X[hyb.complete], plot = FALSE)
num.breaks <- length(X.dist$breaks)
dens.scale <- max(X.dist$density) * 10

# primary
plot(hyb.modtosev.ksmooth.fit$a.eval, hyb.modtosev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Hybrid; Moderate to Severe-critical COVID",
     xlim = c(min(X[hyb.complete]), max(X[hyb.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(hyb.modtosev.ksmooth.fit$a.eval, rev(hyb.modtosev.ksmooth.fit$a.eval)),
        y = c(1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2], rev(1-hyb.modtosev.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
legend("topright", paste0("p-value: ", signif(hyb.modtosev.test$p.vals[2])), cex = .75)

# any infection
plot(hyb.any.ksmooth.fit$a.eval, hyb.any.ksmooth.fit$drf.est[,2], ylim = c(0, max(1-hyb.any.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Hybrid; Infection, Any Symptomatology",
     xlim = c(min(X[hyb.complete]), max(X[hyb.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-hyb.any.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(hyb.any.ksmooth.fit$a.eval, rev(hyb.any.ksmooth.fit$a.eval)),
        y = c(1-hyb.any.ksmooth.fit$drf.cb.lower[,2], rev(1-hyb.any.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
legend("topright", paste0("p-value: ", signif(hyb.any.test$p.vals[2])), cex = .75)

# severe
plot(hyb.sev.ksmooth.fit$a.eval, hyb.sev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1-hyb.sev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 45 Days after BD-28",
     main = "Hybrid; Severe-critical COVID",
     xlim = c(min(X[hyb.complete]), max(X[hyb.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-hyb.sev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(hyb.sev.ksmooth.fit$a.eval, rev(hyb.sev.ksmooth.fit$a.eval)),
        y = c(1-hyb.sev.ksmooth.fit$drf.cb.lower[,2], rev(1-hyb.sev.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
legend("topright", paste0("p-value: ", signif(hyb.sev.test$p.vals[2])), cex = .75)




dev.off()



########### B1; Vaccine
X <-  ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1
A <- with(ensemble.dat.immuno, ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
W <- with(ensemble.dat.immuno, cbind(Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex)))
# W <- with(ensemble.dat.immuno, cbind(RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))
complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))})

vacc.complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 0
# vacc.complete.sev <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 0 & 
#   ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
vacc.complete.sev <- apply(cbind(A, X), 1, function(x){all(!is.na(x))}) & A == 0 & 
  ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
vacc.complete.modtosev <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 0 & 
  ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
vacc.complete.any <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & A == 0 & 
  ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 0

vacc.time.modtosev <- ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[vacc.complete.modtosev]
vacc.event.modtosev <- ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[vacc.complete.modtosev]
vacc.time.sev <- ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[vacc.complete.sev]
vacc.event.sev <- ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[vacc.complete.sev]
vacc.time.any <- ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[vacc.complete.any]
vacc.event.any <- ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[vacc.complete.any]

# vacc.prop.score.fit <- EstimateConditionalDensity(a = X[vacc.complete], w = W[vacc.complete,], weights = sweights[vacc.complete])
vacc.prop.score.fit.modtosev <- EstimateConditionalDensity(a = X[vacc.complete.modtosev], w = W[vacc.complete.modtosev,], weights = sweights[vacc.complete.modtosev])
# vacc.prop.score.fit.sev <- EstimateConditionalDensity(a = X[vacc.complete.sev], w = W[vacc.complete.sev,], weights = sweights[vacc.complete.sev])
vacc.prop.score.fit.any <- EstimateConditionalDensity(a = X[vacc.complete.any], w = W[vacc.complete.any,], weights = sweights[vacc.complete.any])

vacc.modtosev.np.surv.fit <- EstimateConditionalSurvival(time = vacc.time.modtosev, event = vacc.event.modtosev, 
                                                         a = X[vacc.complete.modtosev], w = W[vacc.complete.modtosev,],
                                                         weights = sweights[vacc.complete.modtosev], t0 = t0)
vacc.any.np.surv.fit <- EstimateConditionalSurvival(time = vacc.time.any, event = vacc.event.any, 
                                                    a = X[vacc.complete.any], w = W[vacc.complete.any,],
                                                    weights = sweights[vacc.complete.any], t0 = t0)
vacc.sev.np.surv.fit <- EstimateConditionalSurvival(time = vacc.time.sev, event = vacc.event.sev, 
                                                    a = X[vacc.complete.sev], w = NULL,#= W[vacc.complete.sev,],
                                                    weights = sweights[vacc.complete.sev], t0 = t0.sev,
                                                    nonparametric = FALSE)


vacc.modtosev.ksmooth.fit <-KernelSmoothFit(time = vacc.time.modtosev, event = vacc.event.modtosev, 
                                            a = X[vacc.complete.modtosev], a.eval = X.eval, w = W[vacc.complete.modtosev,], 
                                            weights = sweights[vacc.complete.modtosev],
                                            t0 = t0, 
                                            prop.score = vacc.prop.score.fit.modtosev$prop.score,
                                            fitted.surv = vacc.modtosev.np.surv.fit$fitted.surv,
                                            fitted.cens = vacc.modtosev.np.surv.fit$fitted.cens,
                                            cond.surv = vacc.modtosev.np.surv.fit$cond.surv,
                                            cond.cens = vacc.modtosev.np.surv.fit$cond.cens,
                                            fitted.haz.event = vacc.modtosev.np.surv.fit$fitted.haz.event,
                                            fitted.haz.cens = vacc.modtosev.np.surv.fit$fitted.haz.cens,
                                            alpha = .05)
vacc.modtosev.test <- MonotoneTest(time = vacc.time.modtosev, event = vacc.event.modtosev, 
                                   a = X[vacc.complete.modtosev], w = W[vacc.complete.modtosev,], 
                                   weights = sweights[vacc.complete.modtosev],
                                   t0 = t0, 
                                   prop.score = vacc.prop.score.fit.modtosev$prop.score,
                                   fitted.surv = vacc.modtosev.np.surv.fit$fitted.surv,
                                   fitted.cens = vacc.modtosev.np.surv.fit$fitted.cens,
                                   cond.surv = vacc.modtosev.np.surv.fit$cond.surv,
                                   cond.cens = vacc.modtosev.np.surv.fit$cond.cens,
                                   fitted.haz.event = vacc.modtosev.np.surv.fit$fitted.haz.event,
                                   fitted.haz.cens = vacc.modtosev.np.surv.fit$fitted.haz.cens)


vacc.any.ksmooth.fit <-KernelSmoothFit(time = vacc.time.any, event = vacc.event.any, 
                                       a = X[vacc.complete.any], a.eval = X.eval, w = W[vacc.complete.any,], 
                                       weights = sweights[vacc.complete.any],
                                       t0 = t0,
                                       prop.score = vacc.prop.score.fit.any$prop.score,
                                       fitted.surv = vacc.any.np.surv.fit$fitted.surv,
                                       fitted.cens = vacc.any.np.surv.fit$fitted.cens,
                                       cond.surv = vacc.any.np.surv.fit$cond.surv,
                                       cond.cens = vacc.any.np.surv.fit$cond.cens,
                                       fitted.haz.event = vacc.any.np.surv.fit$fitted.haz.event,
                                       fitted.haz.cens = vacc.any.np.surv.fit$fitted.haz.cens,
                                       alpha = .05)
vacc.any.test <- MonotoneTest(time = vacc.time.any, event = vacc.event.any, 
                              a = X[vacc.complete.any], w = W[vacc.complete.any,], 
                              weights = sweights[vacc.complete.any],
                              t0 = t0, 
                              prop.score = vacc.prop.score.fit.any$prop.score,
                              fitted.surv = vacc.any.np.surv.fit$fitted.surv,
                              fitted.cens = vacc.any.np.surv.fit$fitted.cens,
                              cond.surv = vacc.any.np.surv.fit$cond.surv,
                              cond.cens = vacc.any.np.surv.fit$cond.cens,
                              fitted.haz.event = vacc.any.np.surv.fit$fitted.haz.event,
                              fitted.haz.cens = vacc.any.np.surv.fit$fitted.haz.cens)

vacc.sev.ksmooth.fit <-KernelSmoothFit(time = vacc.time.sev, event = vacc.event.sev, 
                                       a = X[vacc.complete.sev], a.eval = X.eval, w = NULL,# W[vacc.complete.sev,], 
                                       weights = sweights[vacc.complete.sev],
                                       t0 = t0.sev,
                                       prop.score = rep(1, length(vacc.time.sev)),#vacc.prop.score.fit.sev$prop.score,
                                       fitted.surv = vacc.sev.np.surv.fit$fitted.surv,
                                       fitted.cens = vacc.sev.np.surv.fit$fitted.cens,
                                       cond.surv = vacc.sev.np.surv.fit$cond.surv,
                                       cond.cens = vacc.sev.np.surv.fit$cond.cens,
                                       fitted.haz.event = vacc.sev.np.surv.fit$fitted.haz.event,
                                       fitted.haz.cens = vacc.sev.np.surv.fit$fitted.haz.cens,
                                       alpha = .05)
vacc.sev.test <- MonotoneTest(time = vacc.time.sev, event = vacc.event.sev, 
                              a = X[vacc.complete.sev], w = NULL, #W[vacc.complete.sev,], 
                              weights = sweights[vacc.complete.sev],
                              t0 = t0.sev, 
                              prop.score = rep(1, length(vacc.time.sev)),#vacc.prop.score.fit.sev$prop.score,
                              fitted.surv = vacc.sev.np.surv.fit$fitted.surv,
                              fitted.cens = vacc.sev.np.surv.fit$fitted.cens,
                              cond.surv = vacc.sev.np.surv.fit$cond.surv,
                              cond.cens = vacc.sev.np.surv.fit$cond.cens,
                              fitted.haz.event = vacc.sev.np.surv.fit$fitted.haz.event,
                              fitted.haz.cens = vacc.sev.np.surv.fit$fitted.haz.cens)

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Corr-Survival-B1vaccine.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Corr-Survival-B1vaccine.pdf"), width = 9, height = 7)


par(mar = c(5,4,4,3) + .1)

X.dist <- hist(X[vacc.complete], plot = FALSE)
num.breaks <- length(X.dist$breaks)
dens.scale <- max(X.dist$density) * 10

#primary
plot(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Moderate to Severe-critical COVID",
     xlim = c(min(X[complete]), max(X[complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-vacc.modtosev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(vacc.modtosev.ksmooth.fit$a.eval, rev(vacc.modtosev.ksmooth.fit$a.eval)),
        y = c(1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2], rev(1 - vacc.modtosev.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
legend("topright", paste0("p-value: ", signif(vacc.modtosev.test$p.vals[2])), cex = .75)

#any
plot(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.any.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Infection, Any Symptomatology",
     xlim = c(min(X[complete]), max(X[complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-vacc.any.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(vacc.any.ksmooth.fit$a.eval, rev(vacc.any.ksmooth.fit$a.eval)),
        y = c(1 - vacc.any.ksmooth.fit$drf.cb.lower[,2], rev(1 - vacc.any.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
# legend("bottomright", paste0("p-value: ", signif(test$p.vals[2])), cex = .75)
legend("topright", paste0("p-value: ", signif(vacc.any.test$p.vals[2])), cex = .75)

#severe
plot(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.est[,2], ylim = c(0, max(1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2])),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 45 Days after BD-28",
     main = "Severe-critical COVID",
     xlim = c(min(X[complete]), max(X[complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(X.dist$density) * 1/max(1-vacc.sev.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(num.breaks-1)) {
  polygon(x = c(X.dist$breaks[i], X.dist$breaks[i], 
                X.dist$breaks[i+1], X.dist$breaks[i+1]), 
          y = c(-dens.scale, X.dist$density[i], X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

polygon(x = c(vacc.sev.ksmooth.fit$a.eval, rev(vacc.sev.ksmooth.fit$a.eval)),
        y = c(1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2], rev(1 - vacc.sev.ksmooth.fit$drf.cb.upper[,2])),
        col = c1, border = c1)
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")
# legend("bottomright", paste0("p-value: ", signif(test$p.vals[2])), cex = .75)
legend("topright", paste0("p-value: ", signif(vacc.sev.test$p.vals[2])), cex = .75)

dev.off()



##### Simultaneous figures
# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Corr-Survival-B1.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Corr-Survival-B1.pdf"), width = 9, height = 7)
par(mar = c(5,4,4,3) + .1)
vacc.X.dist <- hist(X[vacc.complete], plot = FALSE)
vacc.num.breaks <- length(vacc.X.dist$breaks)
vacc.dens.scale <- max(vacc.X.dist$density) * 10
hyb.X.dist <- hist(X[hyb.complete], plot = FALSE)
hyb.num.breaks <- length(hyb.X.dist$breaks)
hyb.dens.scale <- max(hyb.X.dist$density) * 10

# modtosev
plot(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.est[,2], 
     ylim = c(0, max(c(1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2], 1 - hyb.modtosev.ksmooth.fit$drf.cb.lower[,2]))),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Moderate to Severe-critical COVID",
     xlim = c(min(X[vacc.complete]), max(X[vacc.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(c(vacc.X.dist$density, hyb.X.dist$density)) * 
  1.2/max(c(1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2], 1 - hyb.modtosev.ksmooth.fit$drf.cb.lower[,2]))
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(vacc.num.breaks-1)) {
  polygon(x = c(vacc.X.dist$breaks[i], vacc.X.dist$breaks[i], 
                vacc.X.dist$breaks[i+1], vacc.X.dist$breaks[i+1]), 
          y = c(-dens.scale, vacc.X.dist$density[i], vacc.X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
for(i in 1:(hyb.num.breaks-1)) {
  polygon(x = c(hyb.X.dist$breaks[i], hyb.X.dist$breaks[i], 
                hyb.X.dist$breaks[i+1], hyb.X.dist$breaks[i+1]), 
          y = c(-dens.scale, hyb.X.dist$density[i], hyb.X.dist$density[i], -dens.scale)/dens.scale,
          col = c1, border = c1) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "salmon")
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "salmon")
lines(vacc.modtosev.ksmooth.fit$a.eval, 1 - vacc.modtosev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "salmon")

lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.modtosev.ksmooth.fit$a.eval, 1-hyb.modtosev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")

legend("topright", 
       legend = c(paste0("Hybrid (p-value: ", signif(hyb.modtosev.test$p.vals[2]), ")"),
                  paste0("Vaccine (p-value: ", signif(vacc.modtosev.test$p.vals[2]), ")")),
       lty = c(1, 1), lwd = c(2,2),
       col = c("blue", "salmon"),
       cex = .75)

# any
plot(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.est[,2], 
     ylim = c(0, max(c(1 - vacc.any.ksmooth.fit$drf.cb.lower[,2], 1 - hyb.any.ksmooth.fit$drf.cb.lower[,2]))),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 60 Days after BD-28",
     main = "Infection, Any Symptomatology",
     xlim = c(min(X[vacc.complete]), max(X[vacc.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(c(vacc.X.dist$density, hyb.X.dist$density)) * 
  1.2/max(c(1 - vacc.any.ksmooth.fit$drf.cb.lower[,2], 1 - hyb.any.ksmooth.fit$drf.cb.lower[,2]))
vacc.dens.scale <- max(X.dist$density) * 1/max(1-vacc.any.ksmooth.fit$drf.cb.lower[,2])
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(vacc.num.breaks-1)) {
  polygon(x = c(vacc.X.dist$breaks[i], vacc.X.dist$breaks[i], 
                vacc.X.dist$breaks[i+1], vacc.X.dist$breaks[i+1]), 
          y = c(-dens.scale, vacc.X.dist$density[i], vacc.X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
for(i in 1:(hyb.num.breaks-1)) {
  polygon(x = c(hyb.X.dist$breaks[i], hyb.X.dist$breaks[i], 
                hyb.X.dist$breaks[i+1], hyb.X.dist$breaks[i+1]), 
          y = c(-dens.scale, hyb.X.dist$density[i], hyb.X.dist$density[i], -dens.scale)/dens.scale,
          col = c1, border = c1) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "salmon")
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "salmon")
lines(vacc.any.ksmooth.fit$a.eval, 1 - vacc.any.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "salmon")

lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.any.ksmooth.fit$a.eval, 1-hyb.any.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")

legend("topright", 
       legend = c(paste0("Hybrid (p-value: ", signif(hyb.any.test$p.vals[2]), ")"),
                  paste0("Vaccine (p-value: ", signif(vacc.any.test$p.vals[2]), ")")),
       lty = c(1, 1), lwd = c(2,2),
       col = c("blue", "salmon"),
       cex = .75)

# sev
plot(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.est[,2],
     ylim = c(0, max(c(1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2], 1 - hyb.sev.ksmooth.fit$drf.cb.lower[,2]))),
     lwd = 2, col = "blue",
     xlab = "BD-28 nAb-ID50 B.1 titer (AU/mL)", ylab = "Cumulative Incidence 45 Days after BD-28",
     main = "Severe-critical COVID",
     xlim = c(min(X[vacc.complete]), max(X[vacc.complete])), type = "n",
     cex.main = .9, xaxt = "n")
dens.scale <- max(c(vacc.X.dist$density, hyb.X.dist$density)) * 
  1.2/max(c(1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2], 1 - hyb.sev.ksmooth.fit$drf.cb.lower[,2]))
axis(side = 1, at = 0:4, c(TeX(paste0("$10^", 0:4, "$"))))
for(i in 1:(vacc.num.breaks-1)) {
  polygon(x = c(vacc.X.dist$breaks[i], vacc.X.dist$breaks[i], 
                vacc.X.dist$breaks[i+1], vacc.X.dist$breaks[i+1]), 
          y = c(-dens.scale, vacc.X.dist$density[i], vacc.X.dist$density[i], -dens.scale)/dens.scale,
          col = c2, border = c2) 
}
for(i in 1:(hyb.num.breaks-1)) {
  polygon(x = c(hyb.X.dist$breaks[i], hyb.X.dist$breaks[i], 
                hyb.X.dist$breaks[i+1], hyb.X.dist$breaks[i+1]), 
          y = c(-dens.scale, hyb.X.dist$density[i], hyb.X.dist$density[i], -dens.scale)/dens.scale,
          col = c1, border = c1) 
}

axis(side = 4, at = seq(0, max(X.dist$density)/dens.scale, length.out = 6),
     labels = c(0, .5, 1, 1.5, 2, 2.5))

mtext("Density", side=4, line = 2)

lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "salmon")
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "salmon")
lines(vacc.sev.ksmooth.fit$a.eval, 1 - vacc.sev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "salmon")

lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.est[,2], lty = 1,
      lwd = 2, col = "blue")
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.cb.lower[,2], lty = 2, col = "blue")
lines(hyb.sev.ksmooth.fit$a.eval, 1-hyb.sev.ksmooth.fit$drf.cb.upper[,2], lty = 2, col = "blue")

legend("topright", 
       legend = c(paste0("Hybrid (p-value: ", signif(hyb.sev.test$p.vals[2]), ")"),
                  paste0("Vaccine (p-value: ", signif(vacc.sev.test$p.vals[2]), ")")),
       lty = c(1, 1), lwd = c(2,2),
       col = c("blue", "salmon"),
       cex = .75)

dev.off()
