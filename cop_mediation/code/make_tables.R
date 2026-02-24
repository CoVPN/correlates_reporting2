library(officer) # read_docx
library(flextable)
library(magrittr)
library(here)
library(glue)
library(kyotil)

Sys.setenv(TRIAL = "nextgen_mock")
source(here::here("..", "_common.R")) 

CORs=c("D31toM12_nextgen_mock_sera")

save.results.to = glue("output/{TRIAL}")

rr2ve = function(rr) formatDouble((1-rr)*100,1)%.%'%'

style_my_table <- function(data) {
  flextable(data) %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 8, part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%     # First col left
    align(j = -1, align = "center", part = "all") %>%  # Others center
    bold(part = "header") %>%
    autofit()
}


marker_sets = unique(assay_metadata$panel)
# hack, for now only work on one marker set
marker_sets=c('pseudoneutid50_sera')


doc <- read_docx("table_template.docx")
for (COR in CORs) {
for (marker_set in marker_sets) {
    
  res = readRDS(glue("{save.results.to}/{COR}/eff_{marker_set}.RDS"))
  
  tab = mysapply (res, function (eff) {
    c(
      # # Total
      # paste0(round(eff[1,1], 3), " (", round(eff[1,2], 3), ", ", round(eff[1,3], 3), ")"),
      
      # Direct, (1 - Direct RR)%
      paste0(rr2ve(eff[3,1]), " (", rr2ve(eff[3,3]), ", ", rr2ve(eff[3,2]), ")"),
      
      # Indirect, (1 - Indirect RR)%
      paste0(rr2ve(eff[2,1]), " (", rr2ve(eff[2,3]), ", ", rr2ve(eff[2,2]), ")"),
      
      # Prop_med
      paste0(formatDouble(eff[4,1], 2), " (", formatDouble(eff[4,2], 2), ", ", formatDouble(eff[4,3], 2), ")")
    )
  })
  colnames(tab) = c("Non-marker mediated rVE", "Marker mediated rVE", "Proportion of rVE mediated*")
  
  # turn to data frame for writing to doc
  tab = as.data.frame(tab)
  # add row names as a column
  tab = cbind("marker"=rownames(tab), tab)

  # Add table to doc
  doc <- body_add_flextable(doc, value = style_my_table(as.data.frame(tab)))
  doc <- body_add_par(doc, "", style = "Normal")
  
}}

print(doc, target = glue('cop_mediation_{TRIAL}_{format(Sys.Date(), "%Y%m%d")}.docx'))
