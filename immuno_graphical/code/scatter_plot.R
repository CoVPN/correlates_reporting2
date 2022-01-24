renv::activate(project = here::here(".."))
source(here::here("../_common.R"))

source(here::here("code", "ggally_cor_resample.R"))

# subset to baseline seronegs and vaccine
sub_dat <- dat.mock %>%
  dplyr::filter(Bserostatus == 0 & Trt == 1)

covid_corr_pairplots(
        plot_dat = sub_dat,
        time = "Day57",
        assays = assays,
        strata = "Bstratum",
        weight = "wt.subcohort",
        plot_title = paste0(
          c("D57"),
          " Ab markers: baseline negative", ", ",
          "vaccine arm"
        ),
        column_labels = labels.axis["Day57", seq_along(assays)] %>% unlist(),
        height = max(1.3 * length(assays) + 0.1, 5.5),
        width = max(1.3 * length(assays), 5.5),
        filename = here::here("figs", "scatter_baselineneg_vacc.pdf")
)
