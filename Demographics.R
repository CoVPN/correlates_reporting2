set.seed(206)
require(lubridate)
current.date <- format(Sys.time(), "%e%b%Y")

if(!exists("spath")) {
  spath <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2"
}
if(!exists("data.path")) {
  data.path <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"
}
ensemble.dat <- read.csv(data.path) #read.csv(paste0(spath, "/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"))
# ensemble.dat <- read.csv("/Volumes/trials/covpn/p3003/analysis/mapping_immune_correlates/post unblinded part B stage 2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv")
ensemble.dat <- ensemble.dat[ensemble.dat$BOOSTDT < as.Date("2022-05-22"),]
ensemble.dat <- ensemble.dat[ensemble.dat$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat$History %in% c("History a", "History b"),]
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

Tab1 <- matrix(nrow = 25, ncol = 5)
colnames(Tab1) <- c("Crossover/Vaccine", "Crossover/Hybrid", "Original/Vaccine", "Original/Hybrid", "Overall")
rownames(Tab1) <- c("N", "Age (Mean, SD)", "Female", "Male", "Undifferentiated Sex", 
                    "BMI Underweight", "BMI Normal weight", "BMI Overweight", "BMI Obese",
                    "Race Black", "Race Asian", "Race Native American", "Race Pacific Islander",
                    "Race White", "Race Multiracial", "Race Not Determined",
                    "USA", "Argentina", "Brazil", "Chile", "Colombia", "Mexico",
                    "Peru", "South Africa", "Comorbidities")

# N
Tab1[1,] <- c(sum(ensemble.dat$ImmediateVaccine), sum(ensemble.dat$ImmediateHybrid),
              sum(ensemble.dat$DelayedVaccine), sum(ensemble.dat$DelayedHybrid), nrow(ensemble.dat))

# Age
age.mean <- with(ensemble.dat, c(mean(Age[ImmediateVaccine == 1]), mean(Age[ImmediateHybrid == 1]),
                                 mean(Age[DelayedVaccine == 1]), mean(Age[DelayedHybrid == 1]), mean(Age)))
age.sd <- with(ensemble.dat, c(sd(Age[ImmediateVaccine == 1]), sd(Age[ImmediateHybrid == 1]),
                               sd(Age[DelayedVaccine == 1]), sd(Age[DelayedHybrid == 1]), sd(Age)))
Tab1[2,] <- paste0(round(age.mean, digits = 1), " (", round(age.sd, digits = 1), ")")

# female sex (female == 1; male == 0; undifferentiated == 2)
sex.total <- with(ensemble.dat, c(sum((Sex == 1)[ImmediateVaccine == 1]), sum((Sex == 1)[ImmediateHybrid == 1]),
                                 sum((Sex == 1)[DelayedVaccine == 1]), sum((Sex == 1)[DelayedHybrid == 1]), sum((Sex == 1))))
sex.perecent <- with(ensemble.dat, c(mean((Sex == 1)[ImmediateVaccine == 1]), mean((Sex == 1)[ImmediateHybrid == 1]),
                                     mean((Sex == 1)[DelayedVaccine == 1]), mean((Sex == 1)[DelayedHybrid == 1]), mean((Sex == 1)))) * 100
Tab1[3,] <- paste0(round(sex.total, digits = 1), " (", round(sex.perecent, digits = 1), "%)")
# male sex (female == 1; male == 0; undifferentiated == 2)
sex.total <- with(ensemble.dat, c(sum((Sex == 0)[ImmediateVaccine == 1]), sum((Sex == 0)[ImmediateHybrid == 1]),
                                  sum((Sex == 0)[DelayedVaccine == 1]), sum((Sex == 0)[DelayedHybrid == 1]), sum((Sex == 0))))
sex.perecent <- with(ensemble.dat, c(mean((Sex == 0)[ImmediateVaccine == 1]), mean((Sex == 0)[ImmediateHybrid == 1]),
                                     mean((Sex == 0)[DelayedVaccine == 1]), mean((Sex == 0)[DelayedHybrid == 1]), mean((Sex == 0)))) * 100
Tab1[4,] <- paste0(round(sex.total, digits = 1), " (", round(sex.perecent, digits = 1), "%)")
# undetermined sex (female == 1; male == 0; undifferentiated == 2)
sex.total <- with(ensemble.dat, c(sum((Sex == 2)[ImmediateVaccine == 1]), sum((Sex == 2)[ImmediateHybrid == 1]),
                                  sum((Sex == 2)[DelayedVaccine == 1]), sum((Sex == 2)[DelayedHybrid == 1]), sum((Sex == 2))))
sex.perecent <- with(ensemble.dat, c(mean((Sex == 2)[ImmediateVaccine == 1]), mean((Sex == 2)[ImmediateHybrid == 1]),
                                     mean((Sex == 2)[DelayedVaccine == 1]), mean((Sex == 2)[DelayedHybrid == 1]), mean((Sex ==2)))) * 100
Tab1[5,] <- paste0(round(sex.total, digits = 1), " (", round(sex.perecent, digits = 1), "%)")


# BMI (1 = Underweight; 2 = Normal; 3 = Overweight; 4 = Obese; Else = Missing)
bmi.total <- with(ensemble.dat, c(sum((BMI == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((BMI == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                  sum((BMI == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((BMI == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                  sum((BMI == 1), na.rm = TRUE)))
bmi.perecent <- with(ensemble.dat, c(mean((BMI == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((BMI == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                     mean((BMI == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((BMI == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                     mean((BMI == 1), na.rm = TRUE))) * 100
Tab1[6,] <- paste0(round(bmi.total, digits = 1), " (", round(bmi.perecent, digits = 1), "%)")
# BMI
bmi.total <- with(ensemble.dat, c(sum((BMI == 2)[ImmediateVaccine == 1], na.rm = TRUE), sum((BMI == 2)[ImmediateHybrid == 1], na.rm = TRUE),
                                  sum((BMI == 2)[DelayedVaccine == 1], na.rm = TRUE), sum((BMI == 2)[DelayedHybrid == 1], na.rm = TRUE),
                                  sum((BMI == 2), na.rm = TRUE)))
bmi.perecent <- with(ensemble.dat, c(mean((BMI == 2)[ImmediateVaccine == 1], na.rm = TRUE), mean((BMI == 2)[ImmediateHybrid == 1], na.rm = TRUE),
                                     mean((BMI == 2)[DelayedVaccine == 1], na.rm = TRUE), mean((BMI == 2)[DelayedHybrid == 1], na.rm = TRUE),
                                     mean((BMI == 2), na.rm = TRUE))) * 100
Tab1[7,] <- paste0(round(bmi.total, digits = 1), " (", round(bmi.perecent, digits = 1), "%)")

# BMI
bmi.total <- with(ensemble.dat, c(sum((BMI == 3)[ImmediateVaccine == 1], na.rm = TRUE), sum((BMI == 3)[ImmediateHybrid == 1], na.rm = TRUE),
                                  sum((BMI == 3)[DelayedVaccine == 1], na.rm = TRUE), sum((BMI == 3)[DelayedHybrid == 1], na.rm = TRUE),
                                  sum((BMI == 3), na.rm = TRUE)))
bmi.perecent <- with(ensemble.dat, c(mean((BMI == 3)[ImmediateVaccine == 1], na.rm = TRUE), mean((BMI == 3)[ImmediateHybrid == 1], na.rm = TRUE),
                                     mean((BMI == 3)[DelayedVaccine == 1], na.rm = TRUE), mean((BMI == 3)[DelayedHybrid == 1], na.rm = TRUE),
                                     mean((BMI == 3), na.rm = TRUE))) * 100
Tab1[8,] <- paste0(round(bmi.total, digits = 1), " (", round(bmi.perecent, digits = 1), "%)")

# BMI
bmi.total <- with(ensemble.dat, c(sum((BMI == 4)[ImmediateVaccine == 1], na.rm = TRUE), sum((BMI == 4)[ImmediateHybrid == 1], na.rm = TRUE),
                                  sum((BMI == 4)[DelayedVaccine == 1], na.rm = TRUE), sum((BMI == 4)[DelayedHybrid == 1], na.rm = TRUE),
                                  sum((BMI == 4), na.rm = TRUE)))
bmi.perecent <- with(ensemble.dat, c(mean((BMI == 4)[ImmediateVaccine == 1], na.rm = TRUE), mean((BMI == 4)[ImmediateHybrid == 1], na.rm = TRUE),
                                     mean((BMI == 4)[DelayedVaccine == 1], na.rm = TRUE), mean((BMI == 4)[DelayedHybrid == 1], na.rm = TRUE),
                                     mean((BMI == 4), na.rm = TRUE))) * 100
Tab1[9,] <- paste0(round(bmi.total, digits = 1), " (", round(bmi.perecent, digits = 1), "%)")

# Race (Black)
black.total <- with(ensemble.dat, c(sum((Black == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((Black == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                  sum((Black == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((Black == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                  sum((Black == 1), na.rm = TRUE)))
black.perecent <- with(ensemble.dat, c(mean((Black == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((Black == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                     mean((Black == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((Black == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                     mean((Black == 1), na.rm = TRUE))) * 100
Tab1[10,] <- paste0(round(black.total, digits = 1), " (", round(black.perecent, digits = 1), "%)")
# Race  (Asian)
asian.total <- with(ensemble.dat, c(sum((Asian == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((Asian == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                    sum((Asian == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((Asian == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                    sum((Asian == 1), na.rm = TRUE)))
asian.perecent <- with(ensemble.dat, c(mean((Asian == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((Asian == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                       mean((Asian == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((Asian == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                       mean((Asian == 1), na.rm = TRUE))) * 100
Tab1[11,] <- paste0(round(asian.total, digits = 1), " (", round(asian.perecent, digits = 1), "%)")
# Race  (NatAmer)
natamer.total <- with(ensemble.dat, c(sum((NatAmer == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((NatAmer == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                    sum((NatAmer == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((NatAmer == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                    sum((NatAmer == 1), na.rm = TRUE)))
natamer.perecent <- with(ensemble.dat, c(mean((NatAmer == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((NatAmer == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                       mean((NatAmer == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((NatAmer == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                       mean((NatAmer == 1), na.rm = TRUE))) * 100
Tab1[12,] <- paste0(round(natamer.total, digits = 1), " (", round(natamer.perecent, digits = 1), "%)")
# Race  (Pacific Islander)
pacisl.total <- with(ensemble.dat, c(sum((PacIsl == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((PacIsl == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                    sum((PacIsl == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((PacIsl == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                    sum((PacIsl == 1), na.rm = TRUE)))
pacisl.perecent <- with(ensemble.dat, c(mean((PacIsl == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((PacIsl == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                       mean((PacIsl == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((PacIsl == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                       mean((PacIsl == 1), na.rm = TRUE))) * 100
Tab1[13,] <- paste0(round(pacisl.total, digits = 1), " (", round(pacisl.perecent, digits = 1), "%)")
# Race  (White)
white.total <- with(ensemble.dat, c(sum((White == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((White == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                    sum((White == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((White == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                    sum((White == 1), na.rm = TRUE)))
white.perecent <- with(ensemble.dat, c(mean((White == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((White == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                       mean((White == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((White == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                       mean((White == 1), na.rm = TRUE))) * 100
Tab1[14,] <- paste0(round(white.total, digits = 1), " (", round(white.perecent, digits = 1), "%)")
# Race  (Multiracial)
multi.total <- with(ensemble.dat, c(sum((Multiracial == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((Multiracial == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                    sum((Multiracial == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((Multiracial == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                    sum((Multiracial == 1), na.rm = TRUE)))
multi.perecent <- with(ensemble.dat, c(mean((Multiracial == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((Multiracial == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                       mean((Multiracial == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((Multiracial == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                       mean((Multiracial == 1), na.rm = TRUE))) * 100
Tab1[15,] <- paste0(round(multi.total, digits = 1), " (", round(multi.perecent, digits = 1), "%)")
# Race  (Unknown or not reported)
otherrace.total <- with(ensemble.dat, c(sum((OtherRace == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((OtherRace == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                    sum((OtherRace == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((OtherRace == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                    sum((OtherRace == 1), na.rm = TRUE)))
otherrace.perecent <- with(ensemble.dat, c(mean((OtherRace == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((OtherRace == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                       mean((OtherRace == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((OtherRace == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                       mean((OtherRace == 1), na.rm = TRUE))) * 100
Tab1[16,] <- paste0(round(otherrace.total, digits = 1), " (", round(otherrace.perecent, digits = 1), "%)")

# Country (USA)
country.total <- with(ensemble.dat, c(sum((Country == 0)[ImmediateVaccine == 1], na.rm = TRUE), sum((Country == 0)[ImmediateHybrid == 1], na.rm = TRUE),
                                        sum((Country == 0)[DelayedVaccine == 1], na.rm = TRUE), sum((Country == 0)[DelayedHybrid == 1], na.rm = TRUE),
                                        sum((Country == 0), na.rm = TRUE)))
country.perecent <- with(ensemble.dat, c(mean((Country == 0)[ImmediateVaccine == 1], na.rm = TRUE), mean((Country == 0)[ImmediateHybrid == 1], na.rm = TRUE),
                                           mean((Country == 0)[DelayedVaccine == 1], na.rm = TRUE), mean((Country == 0)[DelayedHybrid == 1], na.rm = TRUE),
                                           mean((Country == 0), na.rm = TRUE))) * 100
Tab1[17,] <- paste0(round(country.total, digits = 1), " (", round(country.perecent, digits = 1), "%)")
# Country (Argentina)
country.total <- with(ensemble.dat, c(sum((Country == 1)[ImmediateVaccine == 1], na.rm = TRUE), sum((Country == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                      sum((Country == 1)[DelayedVaccine == 1], na.rm = TRUE), sum((Country == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                      sum((Country == 1), na.rm = TRUE)))
country.perecent <- with(ensemble.dat, c(mean((Country == 1)[ImmediateVaccine == 1], na.rm = TRUE), mean((Country == 1)[ImmediateHybrid == 1], na.rm = TRUE),
                                         mean((Country == 1)[DelayedVaccine == 1], na.rm = TRUE), mean((Country == 1)[DelayedHybrid == 1], na.rm = TRUE),
                                         mean((Country == 1), na.rm = TRUE))) * 100
Tab1[18,] <- paste0(round(country.total, digits = 1), " (", round(country.perecent, digits = 1), "%)")
# Country (Brazil)
country.total <- with(ensemble.dat, c(sum((Country == 2)[ImmediateVaccine == 1], na.rm = TRUE), sum((Country == 2)[ImmediateHybrid == 1], na.rm = TRUE),
                                      sum((Country == 2)[DelayedVaccine == 1], na.rm = TRUE), sum((Country == 2)[DelayedHybrid == 1], na.rm = TRUE),
                                      sum((Country == 2), na.rm = TRUE)))
country.perecent <- with(ensemble.dat, c(mean((Country == 2)[ImmediateVaccine == 1], na.rm = TRUE), mean((Country == 2)[ImmediateHybrid == 1], na.rm = TRUE),
                                         mean((Country == 2)[DelayedVaccine == 1], na.rm = TRUE), mean((Country == 2)[DelayedHybrid == 1], na.rm = TRUE),
                                         mean((Country == 2), na.rm = TRUE))) * 100
Tab1[19,] <- paste0(round(country.total, digits = 1), " (", round(country.perecent, digits = 1), "%)")
# Country (Chile)
country.total <- with(ensemble.dat, c(sum((Country == 3)[ImmediateVaccine == 1], na.rm = TRUE), sum((Country == 3)[ImmediateHybrid == 1], na.rm = TRUE),
                                      sum((Country == 3)[DelayedVaccine == 1], na.rm = TRUE), sum((Country == 3)[DelayedHybrid == 1], na.rm = TRUE),
                                      sum((Country == 3), na.rm = TRUE)))
country.perecent <- with(ensemble.dat, c(mean((Country == 3)[ImmediateVaccine == 1], na.rm = TRUE), mean((Country == 3)[ImmediateHybrid == 1], na.rm = TRUE),
                                         mean((Country == 3)[DelayedVaccine == 1], na.rm = TRUE), mean((Country == 3)[DelayedHybrid == 1], na.rm = TRUE),
                                         mean((Country == 3), na.rm = TRUE))) * 100
Tab1[20,] <- paste0(round(country.total, digits = 1), " (", round(country.perecent, digits = 1), "%)")
# Country (Colombia)
country.total <- with(ensemble.dat, c(sum((Country == 4)[ImmediateVaccine == 1], na.rm = TRUE), sum((Country == 4)[ImmediateHybrid == 1], na.rm = TRUE),
                                      sum((Country == 4)[DelayedVaccine == 1], na.rm = TRUE), sum((Country == 4)[DelayedHybrid == 1], na.rm = TRUE),
                                      sum((Country == 4), na.rm = TRUE)))
country.perecent <- with(ensemble.dat, c(mean((Country == 4)[ImmediateVaccine == 1], na.rm = TRUE), mean((Country == 4)[ImmediateHybrid == 1], na.rm = TRUE),
                                         mean((Country == 4)[DelayedVaccine == 1], na.rm = TRUE), mean((Country == 4)[DelayedHybrid == 1], na.rm = TRUE),
                                         mean((Country == 4), na.rm = TRUE))) * 100
Tab1[21,] <- paste0(round(country.total, digits = 1), " (", round(country.perecent, digits = 1), "%)")
# Country (Mexico)
country.total <- with(ensemble.dat, c(sum((Country == 5)[ImmediateVaccine == 1], na.rm = TRUE), sum((Country == 5)[ImmediateHybrid == 1], na.rm = TRUE),
                                      sum((Country == 5)[DelayedVaccine == 1], na.rm = TRUE), sum((Country == 5)[DelayedHybrid == 1], na.rm = TRUE),
                                      sum((Country == 5), na.rm = TRUE)))
country.perecent <- with(ensemble.dat, c(mean((Country == 5)[ImmediateVaccine == 1], na.rm = TRUE), mean((Country == 5)[ImmediateHybrid == 1], na.rm = TRUE),
                                         mean((Country == 5)[DelayedVaccine == 1], na.rm = TRUE), mean((Country == 5)[DelayedHybrid == 1], na.rm = TRUE),
                                         mean((Country == 5), na.rm = TRUE))) * 100
Tab1[22,] <- paste0(round(country.total, digits = 1), " (", round(country.perecent, digits = 1), "%)")
# Country (Peru)
country.total <- with(ensemble.dat, c(sum((Country == 6)[ImmediateVaccine == 1], na.rm = TRUE), sum((Country == 6)[ImmediateHybrid == 1], na.rm = TRUE),
                                      sum((Country == 6)[DelayedVaccine == 1], na.rm = TRUE), sum((Country == 6)[DelayedHybrid == 1], na.rm = TRUE),
                                      sum((Country == 6), na.rm = TRUE)))
country.perecent <- with(ensemble.dat, c(mean((Country == 6)[ImmediateVaccine == 1], na.rm = TRUE), mean((Country == 6)[ImmediateHybrid == 1], na.rm = TRUE),
                                         mean((Country == 6)[DelayedVaccine == 1], na.rm = TRUE), mean((Country == 6)[DelayedHybrid == 1], na.rm = TRUE),
                                         mean((Country == 6), na.rm = TRUE))) * 100
Tab1[23,] <- paste0(round(country.total, digits = 1), " (", round(country.perecent, digits = 1), "%)")
# Country (RSA)
country.total <- with(ensemble.dat, c(sum((Country == 7)[ImmediateVaccine == 1], na.rm = TRUE), sum((Country == 7)[ImmediateHybrid == 1], na.rm = TRUE),
                                      sum((Country == 7)[DelayedVaccine == 1], na.rm = TRUE), sum((Country == 7)[DelayedHybrid == 1], na.rm = TRUE),
                                      sum((Country == 7), na.rm = TRUE)))
country.perecent <- with(ensemble.dat, c(mean((Country == 7)[ImmediateVaccine == 1], na.rm = TRUE), mean((Country == 7)[ImmediateHybrid == 1], na.rm = TRUE),
                                         mean((Country == 7)[DelayedVaccine == 1], na.rm = TRUE), mean((Country == 7)[DelayedHybrid == 1], na.rm = TRUE),
                                         mean((Country == 7), na.rm = TRUE))) * 100
Tab1[24,] <- paste0(round(country.total, digits = 1), " (", round(country.perecent, digits = 1), "%)")

#comorbidities
comor.total <- with(ensemble.dat, c(sum((HighRiskInd == 1)[ImmediateVaccine == 1]), sum((HighRiskInd == 1)[ImmediateHybrid == 1]),
                                  sum((HighRiskInd == 1)[DelayedVaccine == 1]), sum((HighRiskInd == 1)[DelayedHybrid == 1]), sum((HighRiskInd == 1))))
comor.perecent <- with(ensemble.dat, c(mean((HighRiskInd == 1)[ImmediateVaccine == 1]), mean((HighRiskInd == 1)[ImmediateHybrid == 1]),
                                     mean((HighRiskInd == 1)[DelayedVaccine == 1]), mean((HighRiskInd == 1)[DelayedHybrid == 1]), mean((HighRiskInd == 1)))) * 100
Tab1[25,] <- paste0(round(comor.total, digits = 1), " (", round(comor.perecent, digits = 1), "%)")


# write.csv(Tab1, "/Volumes/ahudson/Ensemble-stage2-correlates/SummaryTables/Demographics.csv")
write.csv(Tab1, paste0(spath, "/reports/Demographics_", Sys.Date(), ".csv"))


##########################################
### Distribution of cases by calendar time
##########################################

ensemble.dat$CalendarDateDose2 <- as.Date(ensemble.dat$BOOSTDT)

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
c3 <- rgb(235,192,52,max = 255, alpha = 80, names = "lt.gold")
c4 <- rgb(91,235,52, max = 255, alpha = 80, names = "lt.green")

# pdf("/Volumes/ahudson/Ensemble-stage2-correlates/Summary Figures/DateDose2-FourGroups-Jan-14-2026.pdf")
pdf(paste0(spath, "/reports/DateDose2-FourGroups_", Sys.Date(), ".pdf"))

hist.del.hybrid <- with(ensemble.dat, hist(CalendarDateDose2[group == "Delayed Hybrid"], breaks = "months", plot = FALSE))
hist.imm.hybrid <- with(ensemble.dat, hist(CalendarDateDose2[group == "Immediate Hybrid"], breaks = "months", plot = FALSE))
hist.del.vaccine <- with(ensemble.dat, hist(CalendarDateDose2[group == "Delayed Vaccine"], breaks = "months", plot = FALSE))
hist.imm.vaccine <- with(ensemble.dat, hist(CalendarDateDose2[group == "Immediate Vaccine"], breaks = "months", plot = FALSE))

# with(ensemble.dat[ensemble.dat$group %in% c("Delayed Hybrid", "Immediate Hybrid"),],
#      hist(CalendarDateDose2, breaks = "quarters", freq = TRUE, xlab = "Date of Booster Dose", 
#           main = "", col = "lightgrey"))
# plot(hist.imm.hybrid, add = TRUE, col = c1)
# plot(hist.del.hybrid, add = TRUE, col = c2)
# legend("topright", col = c("lightgrey", c1, c2), legend = c("Hybrid Immunity", "Crossover/Hybrid", "Original/Hybrid"),
#        pch = c(19, 19, 19))
# 
# with(ensemble.dat[ensemble.dat$group %in% c("Delayed Vaccine", "Immediate Vaccine"),],
#      hist(CalendarDateDose2, breaks = "quarters", freq = TRUE, xlab = "Date of Booster Dose", 
#                         main = "", col = "lightgrey"))
# plot(hist.imm.vaccine, add = TRUE, col = c1)
# plot(hist.del.vaccine, add = TRUE, col = c2)
# legend("topright", col = c("lightgrey", c1, c2), legend = c("Vaccine Immunity", "Crossover/Vaccine", "Original/Vaccine"),
#        pch = c(19, 19, 19))

with(ensemble.dat[ensemble.dat$group == c("Delayed Vaccine"),],
     hist(CalendarDateDose2, breaks = "months", freq = TRUE, xlab = "Date of Booster Dose", 
          main = "Original/Vaccine", col = "salmon"))

with(ensemble.dat[ensemble.dat$group == c("Immediate Vaccine"),],
     hist(CalendarDateDose2, breaks = "months", freq = TRUE, xlab = "Date of Booster Dose", 
          main = "Crossover/Vaccine", col = "salmon"))

with(ensemble.dat[ensemble.dat$group == c("Delayed Hybrid"),],
     hist(CalendarDateDose2, breaks = "months", freq = TRUE, xlab = "Date of Booster Dose", 
          main = "Original/Hybrid", col = "salmon"))

with(ensemble.dat[ensemble.dat$group == c("Immediate Hybrid"),],
     hist(CalendarDateDose2, breaks = "months", freq = TRUE, xlab = "Date of Booster Dose", 
          main = "Crossover/Hybrid", col = "salmon"))

dev.off()

##########################################
### Distribution of number of days between 1st and second doses
##########################################

# pdf("/Volumes/ahudson/Ensemble-stage2-correlates/Summary Figures/DaysBetweenDoses.pdf")
pdf(paste0(spath, "/reports/DaysBetweenDoses_", Sys.Date(), ".pdf"))

ensemble.dat$NumberDaysD1toBoosterNew <- numeric(nrow(ensemble.dat))
for(i in 1:nrow(ensemble.dat)) {
  if(ensemble.dat$group[i] %in% c("Delayed Hybrid", "Delayed Vaccine")) {
    ensemble.dat$NumberDaysD1toBoosterNew[i] <- as.Date(ensemble.dat$BOOSTDT[i]) - as.Date(ensemble.dat$TRTSDT[i])
  } else {
    ensemble.dat$NumberDaysD1toBoosterNew[i] <- as.Date(ensemble.dat$BOOSTDT[i]) - as.Date(ensemble.dat$CROSSDT[i]) 
  }
}

with(subset(ensemble.dat, group == "Delayed Hybrid" | group == "Delayed Vaccine"), median(NumberDaysD1toBoosterNew))
with(subset(ensemble.dat, group == "Immediate Hybrid" | group == "Immediate Vaccine"), median(NumberDaysD1toBoosterNew))

hist.del.hybrid <- with(ensemble.dat, hist(NumberDaysD1toBoosterNew[group == "Delayed Hybrid"], plot = FALSE))
hist.imm.hybrid <- with(ensemble.dat, hist(NumberDaysD1toBoosterNew[group == "Immediate Hybrid"], plot = FALSE))
hist.del.vaccine <- with(ensemble.dat, hist(NumberDaysD1toBoosterNew[group == "Delayed Vaccine"], plot = FALSE))
hist.imm.vaccine <- with(ensemble.dat, hist(NumberDaysD1toBoosterNew[group == "Immediate Vaccine"], plot = FALSE))

with(ensemble.dat[ensemble.dat$group %in% c("Delayed Hybrid", "Immediate Hybrid"),],
     hist(NumberDaysD1toBoosterNew, xlab = "Days between First and Second Booster Dose",
          main = "Hybrid Immunity", col = NULL, border = NA, ylim = c(0, 700)))
plot(hist.imm.hybrid, add = TRUE, col = c1)
plot(hist.del.hybrid, add = TRUE, col = c2)
legend("topleft", col = c(c1, c2), legend = c("Crossover/Hybrid", "Original/Hybrid"),
       pch = c(19, 19))

with(ensemble.dat[ensemble.dat$group %in% c("Delayed Vaccine", "Immediate Vaccine"),],
     hist(NumberDaysD1toBoosterNew, freq = TRUE, xlab = "Days between First and Second Booster Dose",
          main = "Vaccine Immunity",  col = NULL, border = NA, ylim = c(0, 3000)))
plot(hist.imm.vaccine, add = TRUE, col = c1)
plot(hist.del.vaccine, add = TRUE, col = c2)
legend("topleft", col = c(c1, c2), legend = c("Crossover/Vaccine", "Original/Vaccine"),
       pch = c(19, 19))
dev.off()

### GISAID plot info
# ensemble.dat$CalendarDateDose2 <- Date(nrow(ensemble.dat))
# for(i in 1:nrow(ensemble.dat)) {
#   if(ensemble.dat$group[i] %in% c("Delayed Hybrid", "Delayed Vaccine")) {
#     ensemble.dat$CalendarDateDose2[i] <- as.Date(ensemble.dat$TRTSDT[i]) + ensemble.dat$NumberdaysD1toYear1booster[i]
#   } else {
#     ensemble.dat$CalendarDateDose2[i] <- as.Date(ensemble.dat$CROSSDT[i]) + ensemble.dat$NumberdaysD1toYear1booster[i]
#   }
# }
# 
# c(min(ensemble.dat$BOOSTDT[ensemble.dat$Country == 0]), 
#   max(ensemble.dat$BOOSTDT[ensemble.dat$Country == 0]))
# c(min(ensemble.dat$BOOSTDT[ensemble.dat$Country == 1]), 
#   max(ensemble.dat$BOOSTDT[ensemble.dat$Country == 1]))
# c(min(ensemble.dat$BOOSTDT[ensemble.dat$Country == 2]), 
#   max(ensemble.dat$BOOSTDT[ensemble.dat$Country == 2]))
# c(min(ensemble.dat$BOOSTDT[ensemble.dat$Country == 3]), 
#   max(ensemble.dat$BOOSTDT[ensemble.dat$Country == 3]))
# c(min(ensemble.dat$BOOSTDT[ensemble.dat$Country == 4]), 
#   max(ensemble.dat$BOOSTDT[ensemble.dat$Country == 4]))
# c(min(ensemble.dat$BOOSTDT[ensemble.dat$Country == 5]), 
#   max(ensemble.dat$BOOSTDT[ensemble.dat$Country == 5]))
# c(min(ensemble.dat$BOOSTDT[ensemble.dat$Country == 6]), 
#   max(ensemble.dat$BOOSTDT[ensemble.dat$Country == 6]))
# c(min(ensemble.dat$BOOSTDT[ensemble.dat$Country == 7]), 
#   max(ensemble.dat$BOOSTDT[ensemble.dat$Country == 7]))

