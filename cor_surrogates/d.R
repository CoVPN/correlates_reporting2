d <- read_csv("d.csv")
# d = d %>% mutate(diff = Day210IgG3gp4140delta - Day1IgG3gp4140delta)
# mylogit <- glm(Delta.D210 ~ RSA + Age + BMI + Riskscore + diff, data = d, family = "binomial")
# summary(mylogit)
# 
# Call:
#   glm(formula = Delta.D210 ~ RSA + Age + BMI + Riskscore + Day1IgG3gp4140delta + 
#         Day210IgG3gp4140delta, family = "binomial", data = d)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.032  -0.657  -0.561  -0.451   2.255  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)           -1.61921    0.16238   -9.97   <2e-16 ***
#   RSA                    0.00257    0.18053    0.01    0.989    
#   Age                    0.08732    0.15725    0.56    0.579    
#   BMI                    0.01869    0.16985    0.11    0.912    
#   Riskscore              0.38564    0.15944    2.42    0.016 *  
#   Day1IgG3gp4140delta   -0.23019    0.17212   -1.34    0.181    
#   Day210IgG3gp4140delta -0.03507    0.17447   -0.20    0.841    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 270.89  on 292  degrees of freedom
# Residual deviance: 263.04  on 286  degrees of freedom
# AIC: 277
# 
# Number of Fisher Scoring iterations: 4
# 
# 
# For a one unit increase in Day210IgG3gp4140delta, the log odds of having COVID infection decreases by 0.03507.
# For a one unit increase in Day1IgG3gp4140delta, the log odds of having COVID infection decreases by 0.23.
# Being in RSA, versus other countries, the log odds of Covid infection increases by 0.00257.