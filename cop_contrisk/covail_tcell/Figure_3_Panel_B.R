dir.create("./Figure/Figure_3", showWarnings = FALSE, recursive=TRUE) 
################################################################################
# Load packages
library(dplyr)
library(haven)
library(forestplot)
library(Cairo)
################################################################################
################################################################################
################################################################################
# All mRNA; 6-Month
table_lables = c('\n Analysis \n Cohort \n',  
                 'Time\nPoint', 
                 'T cell \n Subset',
                 '\nPeptide \nPool\n',
                 'Cytokine \nCombination',
                 '\nHazard Ratio \n(95% CI)\n',
                 '\nP-value\n', 
                 '\nFWER\n')

# Naive

# D1 % CD4+ T-cells IFNg and/or IL-2 BA.4/5 Spike 
# D1 % CD8+ T-cells IFNg and/or IL-2 BA.4/5 Spike 1.03 (0.92, 1.15) 0.655


# D15 % CD4+ T-cells IFNg and/or IL-2 BA.4/5 Spike 
# D15 CD8+ IFNg and/or IL-2 BA.4/5 Spike 


# Non-Naive

# D1 % CD4+ T-cells IFNg and/or IL-2 BA.4/5 Spike 
# D1 % CD8+ T-cells IFNg and/or IL-2 BA.4/5 Spike BA.4/5

# D1 % CD4+ T-cells IFNg and/or IL-2 N Index 

# D15 % CD4+ T-cells IFNg and/or IL-2 BA.4/5 Spike 
# D15 CD8+ IFNg and/or IL-2 BA.4/5 Spike 


mean_vec = c(NA, 1.10, 1.03,
             1.03, 1.02,
             0.64, 0.92,
             0.75, 
             0.62, 0.92)

CI_low = c(NA, 0.91, 0.92,
           0.81, 0.91,
           0.44, 0.74,
           0.56, 
           0.39, 0.73)

CI_high = c(NA, 1.33, 1.15,
            1.32, 1.14,
            0.93, 1.14,
            1.01,  
            0.99, 1.15) 

p_val_vec = c(NA, 0.331, 0.655,
              0.803,0.732,  
              0.020, 0.426,
              0.055,  
              0.044, 0.459) 

p_adj_vec = c(NA, p.adjust(c(0.331, 0.559,
                             0.803,0.732), method ='holm'),
              p.adjust(c(0.020, 0.426,
                         0.055,  
                         0.044, 0.459), method = 'holm'))

# HR and CI
# Keep two digits, i.e., 1.40, 2.25, etc
CI = paste0(format(mean_vec, nsmall = 2), ' (', CI_low, ', ', CI_high, ')')
CI[is.na(mean_vec)] = NA
CI[1] = table_lables[6]

# P-value
p_val_vec[1] = table_lables[7]

# FWER
p_adj_vec[1] = table_lables[8]


# Make the table text
tabletext = cbind(c(table_lables[1], 'Naive', rep('', 3), 
                    'Non-Naive', rep('', 4)),
                  c(table_lables[2], 'D1','', 
                    'D15','', 
                    'D1','', '',  
                    'D15', ''),
                  c(table_lables[3], 'CD4+','CD8+', 'CD4+', 'CD8+',
                    'CD4+','CD8+', 'CD4+', 'CD4+', 'CD8+'),
                  c(table_lables[4], c('Spike BA.4/5', 'Spike BA.4/5',
                                       'Spike BA.4/5', 'Spike BA.4/5',
                                       'Spike BA.4/5', 'Spike BA.4/5', 
                                       'Nucleocapsid \nIndex',
                                       'Spike BA.4/5', 'Spike BA.4/5')),
                  c(table_lables[5], rep(c('IFN-\u03B3 and/or IL-2'), 9)),
                  CI,
                  p_val_vec,
                  p_adj_vec
)

# X ticks in the log-scale
xticks <- c(0.125, 0.25, 0.5, 1, 2, 4)
attr(xticks, "labels")  = c('0.125', '0.25','0.5', '1', '2', '4')


cairo_pdf('./Figure/Figure_3/ForestPlot_onedosemRNA_D22toD181.pdf', 
          width = 12, height = 6, pointsize = 8)
print(forestplot(
  tabletext,
  mean = mean_vec,
  low = CI_low,
  upper = CI_high,
  is.summary = FALSE,
  col = fpColors(
    box = "darkblue",
    line = "darkblue",
    summary = "royalblue",
    zero = 'white'
  ),
  graph.pos = 6,
  graphwidth = 'auto',
  hrzl_lines = list("2" = gpar(lty=1)),
  lwd.zero = 0.5,
  lwd.ci = 0.5,
  lwd.xaxis = 0.5,
  xticks = xticks,
  xlog = TRUE,
  boxsize = 0.1,
  lineheight = 'auto',
  grid = structure(c(1), 
                   gp = gpar(lty = 2, col = "#CCCCFF", size = 1.1)), 
  txt_gp = fpTxtGp(
    ticks = gpar(fontfamily = "Arial", cex = 1.5),
    label = gpar(fontfamily = "Arial", cex = 1.5),
    summary = gpar(fontfamily = "Arial", cex = 1)
  ),
  colgap = unit(4.5, "mm"),
  align = c("c", "c", "c", "c", "c"),
  mar = unit(c(3,4,5,5), "mm")
))
dev.off()


