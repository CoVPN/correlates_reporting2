dir.create("./Figure/Figure_5", showWarnings = FALSE, recursive=TRUE) 
###############################################################################
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
###############################################################################
# Load the CD4+ FS detailed data
dt_all_combn = read.csv('./CD4+T_dt-summary_COMPASS.csv')

# 19 different combinations
length(table(dt_all_combn$combn))


# These are coded as +/- for each of the following
# 8 cytokines:
# CD154
# IFNg
# IL2
# IL4
# IL5_OR_IL13
# IL17a
# IL21
# TNF

# Recode the expression combinations
code_combn <- function(combn_str){
  return(paste(unlist(str_extract_all(combn_str, '[+-]')), collapse = ''))
}

dt_all_combn$combn_b = unlist(lapply(str_extract_all(dt_all_combn$combn, '[+-]'), 
                                     function(x) paste0(x, collapse = '')))


dt2 = dt_all_combn %>%
  mutate(tp = ifelse(blinded_VISITNO == 3232324948, 1, 
                     ifelse(blinded_VISITNO == 3232325148, 15, 
                            ifelse(blinded_VISITNO == 3232325348, 91, 181))),
         pct = CYTNUM/NSUB) %>%
  pivot_wider(id_cols = c(PTID, tp, combn, combn_b, degree),
              names_from = 'STIM', values_from = 'pct') %>%
  mutate(BA_45_S1 = 100*(`BA.4-5 S1` - negctrl),
         BA_45_S2 = 100*(`BA.4-5 S2` - negctrl), 
         BA_45_S1_trunc = ifelse(BA_45_S1 >= 0.001, BA_45_S1, 0.001),
         BA_45_S2_trunc = ifelse(BA_45_S2 >= 0.001, BA_45_S2, 0.001),
         BA_Spike = BA_45_S1_trunc + BA_45_S2_trunc)

dt2 = dt2 %>%
  filter(tp == 1 | tp == 15)

# dim(dt2) = 21356
# This is approx 19 cytokine combinations for 569 ptids for 2 tps

# Beliw focus on BA.4/5 S1 and BA.4/5 S2
dt_BA45 = dt2 %>%
  select(PTID, tp, combn, combn_b, degree, BA_45_S1, BA_45_S2, BA_Spike) %>%
  pivot_longer(cols = c('BA_45_S1', 'BA_45_S2', 'BA_Spike'),
               names_to = 'STIM', values_to = 'pct_subtract')


###########################################################################
# To prepare for plotting, label different cytokine combinations
# CD154
# IFNg
# IL2
# IL4
# IL5_OR_IL13
# IL17a
# IL21
# TNF

# A total of 4 one-cyt
# 1 --> CD154
# 2 --> IFNg
# 3 --> IL2
# 4 --> IL21

# A total of 6 two-cyt
# 5 --> 154 + IFNg
# 6 --> 154 + IL2
# 7 --> 154 + TNF
# 8 --> IFNg + IL2
# 9 --> IFNg + TNF
# 10 --> IL2 + TNF

# A total of 5 three-cyt
# 11 --> 154 + IFNg + IL2
# 12 --> 154 + IFNg + TNF
# 13 --> 154 + IL2 + TNF
# 14 --> IFNg + IL2 + TNF
# 15 --> IL2 + IL21 + TNF

# A total of 3 four-cyt
# 16 --> 154 + IFNg + IL2 + TNF
# 17 --> 154 + IFNg + IL21 + TNF
# 18 --> 154 + IL2 + IL21 + TNF

# One five-cyt
# 19 --> 154 + IFNg + IL2 + IL21 + TNF

dt_BA45 = dt_BA45 %>%
mutate(combn_b_f =case_when(combn_b == '+-------' ~ 1,
                            combn_b == '-+------' ~ 2,
                            combn_b == '--+-----' ~ 3,
                            combn_b == '------+-' ~ 4,
                            combn_b == '++------' ~ 5,
                            combn_b == '+-+-----' ~ 6,
                            combn_b == '+------+' ~ 7,
                            combn_b == '-++-----' ~ 8,
                            combn_b == '-+-----+' ~ 9,
                            combn_b == '--+----+' ~ 10,
                            combn_b == '+++-----' ~ 11,
                            combn_b == '++-----+' ~ 12,
                            combn_b == '+-+----+' ~ 13,
                            combn_b == '-++----+' ~ 14,
                            combn_b == '--+---++' ~ 15,
                            combn_b == '+++----+' ~ 16,
                            combn_b == '++----++' ~ 17,
                            combn_b == '+-+---++' ~ 18,
                            combn_b == '+++---++' ~ 19,
                            TRUE ~ NA))


# Draw a heatmap corresponding to the cyt config
cyt = rep(c('CD154', 'IFNg', 'IL2','IL4', 'IL5 and/or IL13','IL17a', 'IL21', 'TNF'),
          19)

comb = rep(seq(1,19,1), each = 8)

col = c(1,0,0,0,0,0,0,0,
        0,1,0,0,0,0,0,0,
        0,0,1,0,0,0,0,0,
        0,0,0,0,0,0,1,0,
        2,2,0,0,0,0,0,0,
        2,0,2,0,0,0,0,0,
        2,0,0,0,0,0,0,2,
        0,2,2,0,0,0,0,0,
        0,2,0,0,0,0,0,2,
        0,0,2,0,0,0,0,2,
        3,3,3,0,0,0,0,0,
        3,3,0,0,0,0,0,3,
        3,0,3,0,0,0,0,3,
        0,3,3,0,0,0,0,3,
        0,0,3,0,0,0,3,3,
        4,4,4,0,0,0,0,4,
        4,4,0,0,0,0,4,4,
        4,0,4,0,0,0,4,4,
        5,5,5,0,0,0,5,5)

df = data.frame(cyt, comb, col)
df$cyt <- factor(df$cyt, levels = rev(unique(df$cyt)))

p_comb = ggplot(df, aes(x = comb, y = cyt, fill = as.factor(col))) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "0" = "white",
      "1" = "lightblue",
      "2" = "blue",
      "3" = "lightgreen",
      "4" = "green",
      "5" = "lightpink"
    ),
    na.value = "white"
  ) +
  scale_x_discrete(label = NULL) +
  theme_bw(base_size = 30) +
  theme(legend.position = 'none') +
  labs(x = NULL, y = NULL)

ggsave('./Figure/Figure_5/cyt_combination.pdf', p_comb, 
       width = 16, height = 6, units = 'in')


# Read cohort data 
dt_cohort = read.csv('./covail_data_processed_20250612.csv') 
dt_cohort_sub = dt_cohort %>%
  select(Ptid, TrtonedosemRNA, TrtSanofi, naive, COVIDIndD22toD91, COVIDIndD22toD181, 
         ph1.D15.tcell, ph2.D15.tcell, wt.D15.tcell)

dt_BA45_merge = merge(dt_BA45, dt_cohort_sub, by.x = 'PTID', by.y = 'Ptid')
dt_BA45_merge = dt_BA45_merge %>%
  mutate(pct_subtract = ifelse(pct_subtract > 0.001, pct_subtract, 0.001))

##########################################################################
# We do not plot distal cases
dt_BA45_merge = dt_BA45_merge %>%
  mutate(booster_proximal_case = ifelse(COVIDIndD22toD91 == 1, 1, 0),
         non_case = ifelse(COVIDIndD22toD181 == 0, 1, 0)) %>%
  filter(booster_proximal_case == 1 | non_case == 1)

# We only include one-dose mRNA
dt_BA45_merge = dt_BA45_merge %>%
  filter(TrtSanofi == 1)

##########################################################################
# Make some plots
# We are interested in the distribution of each marker in case vs non-case
# and also by naive/non-naive


dt_BA45_merge_naive = dt_BA45_merge %>%
  filter(naive == 1)

p_D1_N = dt_BA45_merge_naive %>%
  filter(ph1.D15.tcell == 1,
         tp == 1) %>%
  mutate(COVIDIndD22toD91_str = ifelse(COVIDIndD22toD91 == 1, 'Case', 'Non-Case')) %>%
  ggplot(aes(x = as.factor(combn_b_f), 
             y = pct_subtract,
             weight = wt.D15.tcell,
             fill = COVIDIndD22toD91_str)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual('Case status D22-D91', values = c('plum3', 'gold')) +
  xlab('') +
  scale_y_log10('D1 % CD4+ T cells expressing a cytokine \n combination against Spike BA.4/5',
                limits = c(0.001, 2),
                breaks = c(0.001, 0.01, 0.1, 1),
                labels = c('0.001%', '0.01%', '0.1%', '1%') ) +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 30) +
  theme(legend.position = 'none')

ggsave('./Figure/Figure_5/D1_N_BA45.pdf', p_D1_N, width = 16, height = 12, units = 'in')

p_D15_N = dt_BA45_merge_naive %>%
  filter(ph1.D15.tcell == 1,
         tp == 15) %>%
  mutate(COVIDIndD22toD91_str = ifelse(COVIDIndD22toD91 == 1, 'Case', 'Non-Case')) %>%
  ggplot(aes(x = as.factor(combn_b_f), 
             y = pct_subtract,
             weight = wt.D15.tcell,
             fill = COVIDIndD22toD91_str)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual('Case status D22-D91', values = c('plum3', 'gold')) +
  xlab('') +
  scale_y_log10('D15 % CD4+ T cells expressing a cytokine \n combination against Spike BA.4/5',
                limits = c(0.001, 2),
                breaks = c(0.001, 0.01, 0.1, 1),
                labels = c('0.001%', '0.01%', '0.1%', '1%') ) +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 30) +
  theme(legend.position = 'none')

ggsave('./Figure/Figure_5/D15_N_BA45.pdf', p_D15_N, width = 16, height = 12, units = 'in')

dt_BA45_merge_nonnaive = dt_BA45_merge %>%
  filter(naive == 0)

p_D1_NN = dt_BA45_merge_nonnaive %>%
  filter(ph1.D15.tcell == 1,
         tp == 1) %>%
  mutate(COVIDIndD22toD91_str = ifelse(COVIDIndD22toD91 == 1, 'Case', 'Non-Case')) %>%
  ggplot(aes(x = as.factor(combn_b_f), 
             y = pct_subtract,
             weight = wt.D15.tcell,
             fill = COVIDIndD22toD91_str)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual('Case status D22-D91', values = c('plum3', 'gold')) +
  xlab('') +
  scale_y_log10('D1 % CD4+ T cells expressing a cytokine \n combination against Spike BA.4/5',
                limits = c(0.001, 2),
                breaks = c(0.001, 0.01, 0.1, 1),
                labels = c('0.001%', '0.01%', '0.1%', '1%') ) +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 30) +
  theme(legend.position = 'none')

ggsave('./Figure/Figure_5/D1_NN_BA45.pdf', p_D1_NN, width = 16, height = 12, units = 'in')


p_D15_NN = dt_BA45_merge_nonnaive %>%
  filter(ph1.D15.tcell == 1,
         tp == 15) %>%
  mutate(COVIDIndD22toD91_str = ifelse(COVIDIndD22toD91 == 1, 'Case', 'Non-Case')) %>%
  ggplot(aes(x = as.factor(combn_b_f), 
             y = pct_subtract,
             weight = wt.D15.tcell,
             fill = COVIDIndD22toD91_str)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual('Case status D22-D91', values = c('plum3', 'gold')) +
  xlab('') +
  scale_y_log10('D15 % CD4+ T cells expressing a cytokine \n combination against Spike BA.4/5',
                limits = c(0.001, 2),
                breaks = c(0.001, 0.01, 0.1, 1),
                labels = c('0.001%', '0.01%', '0.1%', '1%') ) +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 30) +
  theme(legend.position = 'none')

ggsave('./Figure/Figure_5/D15_NN_BA45.pdf', p_D15_NN, width = 16, height = 12, units = 'in')





