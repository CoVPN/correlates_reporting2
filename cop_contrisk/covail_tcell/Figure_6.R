dir.create("./Figure/Figure_6", showWarnings = FALSE, recursive=TRUE) 
 ###############################################################################
# Load packages
library(dplyr)
library(vaccine)
library(ggplot2)
library(patchwork)
library(Hmisc)
library(ggplot2)
###############################################################################
# Helper functions
###############################################################################
# Load the data
dt = read.csv('./covail_data_processed_20250612.csv')

# Define the PP correlates cohort
dt = dt %>%
  filter(ph1.D15.tcell == 1) %>% # n = 932
  mutate(risk_score = ifelse(is.na(risk_score), median(risk_score), risk_score)) %>%
  mutate(risk_score = scale(risk_score)) %>%
  mutate(casecontrol = TwophasesampIndD15.tcell) # all are included in the phase 2 sampling


dt = dt %>%
  mutate(booster_proximal_case = ifelse(COVIDIndD22toD91 == 1, 1, 0),
         non_case = ifelse(COVIDIndD22toD181 == 0, 1, 0))

## NN
dt_PP = dt %>%
  filter(TrtonedosemRNA == 1,
         ph2.D15.tcell == 1,
         naive == 0)

# The 3-month plot only includes booster-proximal and non-case
dt_PP = dt_PP %>%
  filter(booster_proximal_case == 1 | non_case == 1)

p = ggplot(dt_PP, aes(x = Day15pseudoneutid50_BA.4.BA.5, 
                  y = Bcd4_IFNg.IL2_BA.4.5.S)) +
  geom_point(aes(color = as.factor(COVIDIndD22toD91),
                 shape = as.factor(COVIDIndD22toD91),
                 size = as.factor(COVIDIndD22toD91))) +
  geom_hline(yintercept = median(dt_PP$Bcd4_IFNg.IL2_BA.4.5.S, na.rm = TRUE),
             color = 'black', alpha = 0.3, linetype = 'dashed', linewidth = 1.5) +
  geom_vline(xintercept = median(dt_PP$Day15pseudoneutid50_BA.4.BA.5, na.rm = TRUE),
             color = 'black', alpha = 0.3, linetype = 'dashed', linewidth = 1.5) +
  scale_color_manual('', values = c("steelblue", "#E69F00"), 
                        labels = c('Non-Case',
                                   'Booster-Proximal COVID-19 Case')) +
  scale_shape_manual('', values = c(16, 17), 
                     labels = c('Non-Case',
                                'Booster-Proximal COVID-19 Case')) +
  scale_size_manual('', values = c(4, 5.5), 
                    labels = c('Non-Case',
                               'Booster-Proximal COVID-19 Case')) +
  scale_x_continuous('D15 nAb-ID50 Omicron Spike BA.4/5', 
                     limits = c(1, 5),
                     breaks = c(1,2,3,4,5)) +
  scale_y_continuous(expression(paste('D1 CD4+ IFN-', gamma, ' and/or IL-2 Spike BA.4/5')),
                     limits = c(-3, 1),
                     breaks = c(-3, -2, -1, 0, 1),
                     labels = c('0.001%', '0.01%', '0.1%', '1%', '10%')) +
  theme_bw(base_size = 28) +
  theme(legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.position = "inside",
        legend.position.inside = c(0.3, 0.07))

ggsave('./Figure/Figure_6/two_marker_case_non_case_NN.pdf', p, width = 12, height = 12, units = 'in')

# NN all case + non-case
dt_PP = dt %>%
  filter(TrtonedosemRNA == 1,
         ph2.D15.tcell == 1,
         naive == 0)

p2 = ggplot(dt_PP, aes(x = Day15pseudoneutid50_BA.4.BA.5, 
                       y = Bcd4_IFNg.IL2_BA.4.5.S)) +
  geom_point(aes(color = as.factor(COVIDIndD22toD181),
                 shape = as.factor(COVIDIndD22toD181),
                 size = as.factor(COVIDIndD22toD181))) +
  geom_hline(yintercept = median(dt_PP$Bcd4_IFNg.IL2_BA.4.5.S, na.rm = TRUE),
             color = 'black', alpha = 0.3, linetype = 'dashed', linewidth = 1.5) +
  geom_vline(xintercept = median(dt_PP$Day15pseudoneutid50_BA.4.BA.5, na.rm = TRUE),
             color = 'black', alpha = 0.3, linetype = 'dashed', linewidth = 1.5) +
  scale_color_manual('', values = c("steelblue", "#E69F00"), 
                     labels = c('Non-Case',
                                'Overall (Booster-Proximal + Booster-Distal COVID-19 Case')) +
  scale_shape_manual('', values = c(16, 17), 
                     labels = c('Non-Case',
                                'Overall (Booster-Proximal + Booster-Distal COVID-19 Case')) +
  scale_size_manual('', values = c(4, 5.5), 
                    labels = c('Non-Case',
                               'Overall (Booster-Proximal + Booster-Distal COVID-19 Case')) +
  scale_x_continuous('D15 nAb-ID50 Omicron Spike BA.4/5', 
                     limits = c(1, 5),
                     breaks = c(1,2,3,4,5)) +
  scale_y_continuous(expression(paste('D1 CD4+ IFN-', gamma, ' and/or IL-2 Spike BA.4/5')),
                     limits = c(-3, 1),
                     breaks = c(-3, -2, -1, 0, 1),
                     labels = c('0.001%', '0.01%', '0.1%', '1%', '10%')) +
  theme_bw(base_size = 28) +
  theme(legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.position = "inside",
        legend.position.inside = c(0.473, 0.07))

ggsave('./Figure/Figure_6/two_marker_case_non_case_NN_all_case.pdf', 
       p2, width = 12, height = 12, units = 'in')



## Naive
dt_PP = dt %>%
  filter(TrtonedosemRNA == 1,
         ph2.D15.tcell == 1,
         naive == 1)

dt_PP = dt_PP %>%
  filter(booster_proximal_case == 1 | non_case == 1)

p = ggplot(dt_PP, aes(x = Day15pseudoneutid50_BA.4.BA.5, 
                      y = Bcd4_IFNg.IL2_BA.4.5.S)) +
  geom_point(aes(color = as.factor(COVIDIndD22toD91),
                 shape = as.factor(COVIDIndD22toD91),
                 size = as.factor(COVIDIndD22toD91))) +
  geom_hline(yintercept = median(dt_PP$Bcd4_IFNg.IL2_BA.4.5.S, na.rm = TRUE),
             color = 'black', alpha = 0.3, linetype = 'dashed', linewidth = 1.5) +
  geom_vline(xintercept = median(dt_PP$Day15pseudoneutid50_BA.4.BA.5, na.rm = TRUE),
             color = 'black', alpha = 0.3, linetype = 'dashed', linewidth = 1.5) +
  scale_color_manual('', values = c("steelblue", "#E69F00"), 
                     labels = c('Non-Case',
                                'Booster-Proximal COVID-19 Case')) +
  scale_shape_manual('', values = c(16, 17), 
                     labels = c('Non-Case',
                                'Booster-Proximal COVID-19 Case')) +
  scale_size_manual('', values = c(4, 5.5), 
                    labels = c('Non-Case',
                               'Booster-Proximal COVID-19 Case')) +
  scale_x_continuous('D15 nAb-ID50 Omicron Spike BA.4/5', 
                     limits = c(1, 5),
                     breaks = c(1,2,3,4,5)) +
  scale_y_continuous(expression(paste('D1 CD4+ IFN-', gamma, ' and/or IL-2 Spike BA.4/5')),
                     limits = c(-3, 1),
                     breaks = c(-3, -2, -1, 0, 1),
                     labels = c('0.001%', '0.01%', '0.1%', '1%', '10%')) +
  theme_bw(base_size = 28) +
  theme(legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.position = "inside",
        legend.position.inside = c(0.3, 0.07))

ggsave('./Figure/Figure_6/two_marker_case_non_case_N.pdf', p, width = 12, height = 12, units = 'in')


dt_PP = dt %>%
  filter(TrtonedosemRNA == 1,
         ph2.D15.tcell == 1,
         naive == 1)

p2 = ggplot(dt_PP, aes(x = Day15pseudoneutid50_BA.4.BA.5, 
                       y = Bcd4_IFNg.IL2_BA.4.5.S)) +
  geom_point(aes(color = as.factor(COVIDIndD22toD181),
                 shape = as.factor(COVIDIndD22toD181),
                 size = as.factor(COVIDIndD22toD181))) +
  geom_hline(yintercept = median(dt_PP$Bcd4_IFNg.IL2_BA.4.5.S, na.rm = TRUE),
             color = 'black', alpha = 0.3, linetype = 'dashed', linewidth = 1.5) +
  geom_vline(xintercept = median(dt_PP$Day15pseudoneutid50_BA.4.BA.5, na.rm = TRUE),
             color = 'black', alpha = 0.3, linetype = 'dashed', linewidth = 1.5) +
  scale_color_manual('', values = c("steelblue", "#E69F00"), 
                     labels = c('Non-Case',
                                'Overall (Booster-Proximal + Booster-Distal COVID-19 Case')) +
  scale_shape_manual('', values = c(16, 17), 
                     labels = c('Non-Case',
                                'Overall (Booster-Proximal + Booster-Distal COVID-19 Case')) +
  scale_size_manual('', values = c(4, 5.5), 
                    labels = c('Non-Case',
                               'Overall (Booster-Proximal + Booster-Distal COVID-19 Case')) +
  scale_x_continuous('D15 nAb-ID50 Omicron Spike BA.4/5', 
                     limits = c(1, 5),
                     breaks = c(1,2,3,4,5)) +
  scale_y_continuous(expression(paste('D1 CD4+ IFN-', gamma, ' and/or IL-2 Spike BA.4/5')),
                     limits = c(-3, 1),
                     breaks = c(-3, -2, -1, 0, 1),
                     labels = c('0.001%', '0.01%', '0.1%', '1%', '10%')) +
  theme_bw(base_size = 28) +
  theme(legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.position = "inside",
        legend.position.inside = c(0.473, 0.07))

ggsave('./Figure/Figure_6/two_marker_case_non_case_N_all_case.pdf', p2, width = 12, height = 12, units = 'in')

