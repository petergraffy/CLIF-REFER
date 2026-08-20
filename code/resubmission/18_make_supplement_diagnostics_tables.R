#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1 && nzchar(args[[1]])) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(file.path("output", "resubmission_pooled", "20260819_all_sites_with_michigan"), mustWork = TRUE)
}

site_dirs_path <- file.path(out_dir, "pooled_site_directories.csv")
stopifnot(file.exists(site_dirs_path))

site_dirs <- read_csv(site_dirs_path, show_col_types = FALSE)

read_sites_file <- function(filename) {
  bind_rows(lapply(seq_len(nrow(site_dirs)), function(i) {
    path <- file.path(site_dirs$site_dir[[i]], filename)
    if (!file.exists(path)) return(tibble())
    read_csv(path, show_col_types = FALSE) %>%
      mutate(
        site = site_dirs$site[[i]],
        site_label = site_dirs$site_label[[i]],
        site_dir = site_dirs$site_dir[[i]],
        source_file = filename
      )
  }))
}

fmt_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.0001 ~ "<0.0001",
    p < 0.001 ~ "<0.001",
    p < 0.01 ~ sprintf("%.3f", p),
    TRUE ~ sprintf("%.2f", p)
  )
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

fisher_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}

term_labels <- c(
  "pm25_per_5" = "PM2.5 per 5 ug/m3",
  "no2_per_10" = "NO2 per 10 ppb",
  "o3_per_10" = "O3 per 10 ppb",
  "age_10" = "Age per 10 years",
  "sex" = "Sex",
  "race_ethnicity" = "Race/ethnicity",
  "charlson_score" = "Charlson comorbidity index",
  "index_year_f" = "Index year",
  "acs_pct_poverty" = "ZCTA poverty",
  "acs_pct_unemployed" = "ZCTA unemployment",
  "acs_pct_no_vehicle" = "ZCTA households without vehicle",
  "acs_pct_nonwhite" = "ZCTA non-White population",
  "acs_median_household_income_10k" = "ZCTA median household income per $10,000",
  "acs_pct_bachelor_plus" = "ZCTA bachelor's degree or higher",
  "GLOBAL" = "Global test"
)

label_term <- function(x) coalesce(recode(x, !!!term_labels), x)

vif_site <- read_sites_file("resubmission_cause_specific_cox_vif_diagnostics.csv") %>%
  mutate(
    gvif = suppressWarnings(as.numeric(gvif)),
    gvif_adjusted = suppressWarnings(as.numeric(gvif_adjusted)),
    df = suppressWarnings(as.numeric(df)),
    n = suppressWarnings(as.numeric(n)),
    events = suppressWarnings(as.numeric(events)),
    term_label = label_term(term)
  )

vif_summary_by_term <- vif_site %>%
  filter(is.finite(gvif_adjusted)) %>%
  group_by(term, term_label, term_type, df) %>%
  summarise(
    n_sites = n_distinct(site),
    n_models = n(),
    median_adjusted_gvif = median(gvif_adjusted, na.rm = TRUE),
    q75_adjusted_gvif = quantile(gvif_adjusted, 0.75, na.rm = TRUE),
    max_adjusted_gvif = max(gvif_adjusted, na.rm = TRUE),
    n_adjusted_gvif_gt_2_5 = sum(gvif_adjusted > 2.5, na.rm = TRUE),
    n_adjusted_gvif_gt_5 = sum(gvif_adjusted > 5, na.rm = TRUE),
    n_adjusted_gvif_gt_10 = sum(gvif_adjusted > 10, na.rm = TRUE),
    sites_with_max = paste(sort(unique(site_label[gvif_adjusted == max(gvif_adjusted, na.rm = TRUE)])), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(max_adjusted_gvif), term_label)

vif_appendix_table <- vif_summary_by_term %>%
  transmute(
    Covariate = term_label,
    `Term type` = term_type,
    `Model df` = df,
    `Sites, n` = n_sites,
    `Model rows, n` = n_models,
    `Median adjusted GVIF` = fmt_num(median_adjusted_gvif),
    `75th percentile adjusted GVIF` = fmt_num(q75_adjusted_gvif),
    `Maximum adjusted GVIF` = fmt_num(max_adjusted_gvif),
    `Rows with adjusted GVIF >2.5` = n_adjusted_gvif_gt_2_5,
    `Rows with adjusted GVIF >5` = n_adjusted_gvif_gt_5,
    `Rows with adjusted GVIF >10` = n_adjusted_gvif_gt_10,
    `Site(s) with maximum` = sites_with_max
  )

fg_zph_site <- read_sites_file("fine_gray_ph_diagnostics_summary_for_pooling.csv") %>%
  mutate(
    p_value = suppressWarnings(as.numeric(p_value)),
    chisq = suppressWarnings(as.numeric(chisq)),
    df = suppressWarnings(as.numeric(df)),
    n = suppressWarnings(as.numeric(n)),
    events = suppressWarnings(as.numeric(events)),
    term_label = label_term(term)
  )

fg_zph_summary <- fg_zph_site %>%
  filter(is.finite(p_value)) %>%
  group_by(outcome, model, term, term_label, term_type, diagnostic, time_transform) %>%
  summarise(
    n_sites = n_distinct(site),
    total_n = sum(n, na.rm = TRUE),
    total_events = sum(events, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    min_p = min(p_value, na.rm = TRUE),
    fisher_p = fisher_p(p_value),
    n_sites_p_lt_0_05 = sum(p_value < 0.05, na.rm = TRUE),
    sites_p_lt_0_05 = paste(sort(unique(site_label[p_value < 0.05])), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(outcome, model, factor(term_type, levels = c("global", "exposure", "covariate")), term_label)

fg_zph_appendix_table <- fg_zph_summary %>%
  filter(term_type %in% c("global", "exposure")) %>%
  transmute(
    Outcome = outcome,
    Model = model,
    Term = term_label,
    `Term type` = term_type,
    `Sites, n` = n_sites,
    `Events, n` = total_events,
    `Median site p` = fmt_p(median_p),
    `Minimum site p` = fmt_p(min_p),
    `Fisher combined p` = fmt_p(fisher_p),
    `Sites with p<0.05, n` = n_sites_p_lt_0_05,
    `Sites with p<0.05` = sites_p_lt_0_05
  )

fg_time_site <- read_sites_file("fine_gray_ph_time_interaction_diagnostics.csv") %>%
  mutate(
    p_value = suppressWarnings(as.numeric(p_value)),
    beta = suppressWarnings(as.numeric(beta)),
    std_error = suppressWarnings(as.numeric(std_error)),
    n = suppressWarnings(as.numeric(n)),
    events = suppressWarnings(as.numeric(events)),
    exposure_label = label_term(exposure_term)
  )

fg_time_summary <- fg_time_site %>%
  filter(is.finite(p_value)) %>%
  group_by(outcome, model, exposure_term, exposure_label, diagnostic, time_transform) %>%
  summarise(
    n_sites = n_distinct(site),
    total_n = sum(n, na.rm = TRUE),
    total_events = sum(events, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    min_p = min(p_value, na.rm = TRUE),
    fisher_p = fisher_p(p_value),
    n_sites_p_lt_0_05 = sum(p_value < 0.05, na.rm = TRUE),
    sites_p_lt_0_05 = paste(sort(unique(site_label[p_value < 0.05])), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(outcome, model, exposure_label)

fg_time_appendix_table <- fg_time_summary %>%
  transmute(
    Outcome = outcome,
    Model = model,
    Exposure = exposure_label,
    `Sites, n` = n_sites,
    `Events, n` = total_events,
    `Median exposure x log(time) coefficient` = fmt_num(median_beta, 3),
    `Median site p` = fmt_p(median_p),
    `Minimum site p` = fmt_p(min_p),
    `Fisher combined p` = fmt_p(fisher_p),
    `Sites with p<0.05, n` = n_sites_p_lt_0_05,
    `Sites with p<0.05` = sites_p_lt_0_05
  )

standardize_model_diagnostics <- function(tbl) {
  if (!nrow(tbl)) return(tibble())
  if ("mean_vfd" %in% names(tbl)) {
    tbl <- tbl %>%
      rename(
        outcome_mean = mean_vfd,
        outcome_variance = variance_vfd,
        outcome_min = min_vfd,
        outcome_q25 = q25_vfd,
        outcome_median = median_vfd,
        outcome_q75 = q75_vfd,
        outcome_max = max_vfd,
        zero_outcome_n = zero_vfd_n,
        zero_outcome_percent = zero_vfd_percent,
        ceiling_outcome_n = ceiling_vfd_n,
        ceiling_outcome_percent = ceiling_vfd_percent,
        noninteger_outcome_n = noninteger_vfd_n,
        noninteger_outcome_percent = noninteger_vfd_percent
      )
  }
  tbl
}

vfd_diag_site <- standardize_model_diagnostics(read_sites_file("primary_vfd_model_diagnostics.csv"))
imv_diag_site <- standardize_model_diagnostics(read_sites_file("primary_imv_duration_model_diagnostics.csv"))
qp_diag_site <- bind_rows(vfd_diag_site, imv_diag_site) %>%
  mutate(
    across(
      any_of(c(
        "n", "outcome_mean", "outcome_variance", "variance_to_mean_ratio",
        "zero_outcome_percent", "ceiling_outcome_percent", "noninteger_outcome_percent",
        "pearson_dispersion", "deviance_dispersion", "quasipoisson_dispersion",
        "se_inflation_vs_poisson", "max_cooks_distance", "influential_cooks_percent",
        "poisson_overdispersion_p"
      )),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

qp_dispersion_summary <- qp_diag_site %>%
  filter(is.finite(quasipoisson_dispersion)) %>%
  group_by(outcome, model) %>%
  summarise(
    n_sites = n_distinct(site),
    total_n = sum(n, na.rm = TRUE),
    weighted_mean_outcome = weighted.mean(outcome_mean, n, na.rm = TRUE),
    weighted_mean_variance_to_mean_ratio = weighted.mean(variance_to_mean_ratio, n, na.rm = TRUE),
    median_quasipoisson_dispersion = median(quasipoisson_dispersion, na.rm = TRUE),
    min_quasipoisson_dispersion = min(quasipoisson_dispersion, na.rm = TRUE),
    max_quasipoisson_dispersion = max(quasipoisson_dispersion, na.rm = TRUE),
    median_se_inflation_vs_poisson = median(se_inflation_vs_poisson, na.rm = TRUE),
    n_sites_poisson_overdispersion_p_lt_0_05 = sum(poisson_overdispersion_p < 0.05, na.rm = TRUE),
    median_zero_percent = median(zero_outcome_percent, na.rm = TRUE),
    median_ceiling_percent = median(ceiling_outcome_percent, na.rm = TRUE),
    median_noninteger_percent = median(noninteger_outcome_percent, na.rm = TRUE),
    median_influential_cooks_percent = median(influential_cooks_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(outcome, model)

qp_dispersion_appendix_table <- qp_dispersion_summary %>%
  transmute(
    Outcome = outcome,
    Model = model,
    `Sites, n` = n_sites,
    `Participants, n` = total_n,
    `Weighted mean outcome` = fmt_num(weighted_mean_outcome),
    `Weighted variance-to-mean ratio` = fmt_num(weighted_mean_variance_to_mean_ratio),
    `Median quasi-Poisson dispersion` = fmt_num(median_quasipoisson_dispersion),
    `Dispersion range` = paste0(fmt_num(min_quasipoisson_dispersion), "-", fmt_num(max_quasipoisson_dispersion)),
    `Median SE inflation vs Poisson` = fmt_num(median_se_inflation_vs_poisson),
    `Sites with Poisson overdispersion p<0.05, n` = n_sites_poisson_overdispersion_p_lt_0_05,
    `Median zero outcome, %` = fmt_num(median_zero_percent),
    `Median ceiling outcome, %` = fmt_num(median_ceiling_percent),
    `Median non-integer outcome, %` = fmt_num(median_noninteger_percent),
    `Median influential Cook's distance, %` = fmt_num(median_influential_cooks_percent)
  )

qp_comparison_site <- bind_rows(
  read_sites_file("primary_vfd_poisson_vs_quasipoisson_diagnostics.csv"),
  read_sites_file("primary_imv_duration_poisson_vs_quasipoisson_diagnostics.csv")
) %>%
  mutate(
    quasipoisson_dispersion = suppressWarnings(as.numeric(quasipoisson_dispersion)),
    quasi_to_poisson_se_ratio = suppressWarnings(as.numeric(quasi_to_poisson_se_ratio)),
    quasi_p_value = suppressWarnings(as.numeric(quasi_p_value)),
    poisson_p_value = suppressWarnings(as.numeric(poisson_p_value)),
    quasi_ratio_of_means = suppressWarnings(as.numeric(quasi_ratio_of_means)),
    poisson_ratio_of_means = suppressWarnings(as.numeric(poisson_ratio_of_means)),
    term_label = label_term(term)
  )

qp_comparison_appendix_table <- qp_comparison_site %>%
  filter(is.finite(quasipoisson_dispersion)) %>%
  group_by(outcome, model, term, term_label) %>%
  summarise(
    n_sites = n_distinct(site),
    median_quasipoisson_dispersion = median(quasipoisson_dispersion, na.rm = TRUE),
    median_quasi_to_poisson_se_ratio = median(quasi_to_poisson_se_ratio, na.rm = TRUE),
    median_quasi_ratio_of_means = median(quasi_ratio_of_means, na.rm = TRUE),
    median_poisson_ratio_of_means = median(poisson_ratio_of_means, na.rm = TRUE),
    n_sites_quasi_p_lt_0_05 = sum(quasi_p_value < 0.05, na.rm = TRUE),
    n_sites_poisson_p_lt_0_05 = sum(poisson_p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(outcome, model, term_label) %>%
  transmute(
    Outcome = outcome,
    Model = model,
    Term = term_label,
    `Sites, n` = n_sites,
    `Median quasi-Poisson dispersion` = fmt_num(median_quasipoisson_dispersion),
    `Median quasi/Poisson SE ratio` = fmt_num(median_quasi_to_poisson_se_ratio),
    `Median quasi-Poisson ratio of means` = fmt_num(median_quasi_ratio_of_means),
    `Median Poisson ratio of means` = fmt_num(median_poisson_ratio_of_means),
    `Sites with quasi-Poisson p<0.05, n` = n_sites_quasi_p_lt_0_05,
    `Sites with Poisson p<0.05, n` = n_sites_poisson_p_lt_0_05
  )

write_csv(vif_site, file.path(out_dir, "site_vif_diagnostics.csv"))
write_csv(vif_summary_by_term, file.path(out_dir, "supplement_vif_diagnostics_summary.csv"))
write_csv(vif_appendix_table, file.path(out_dir, "supplement_vif_diagnostics_table.csv"))

write_csv(fg_zph_site, file.path(out_dir, "site_fine_gray_ph_diagnostics_cox_zph.csv"))
write_csv(fg_zph_summary, file.path(out_dir, "supplement_fine_gray_ph_diagnostics_summary.csv"))
write_csv(fg_zph_appendix_table, file.path(out_dir, "supplement_fine_gray_ph_diagnostics_table.csv"))
write_csv(fg_time_site, file.path(out_dir, "site_fine_gray_ph_time_interaction_diagnostics.csv"))
write_csv(fg_time_summary, file.path(out_dir, "supplement_fine_gray_time_interaction_diagnostics_summary.csv"))
write_csv(fg_time_appendix_table, file.path(out_dir, "supplement_fine_gray_time_interaction_diagnostics_table.csv"))

write_csv(qp_diag_site, file.path(out_dir, "site_quasipoisson_model_diagnostics.csv"))
write_csv(qp_dispersion_summary, file.path(out_dir, "supplement_quasipoisson_dispersion_summary.csv"))
write_csv(qp_dispersion_appendix_table, file.path(out_dir, "supplement_quasipoisson_dispersion_table.csv"))
write_csv(qp_comparison_site, file.path(out_dir, "site_poisson_vs_quasipoisson_diagnostics.csv"))
write_csv(qp_comparison_appendix_table, file.path(out_dir, "supplement_poisson_vs_quasipoisson_table.csv"))

notes <- tibble(
  Diagnostic = c(
    "VIF",
    "Fine-Gray proportional hazards",
    "Fine-Gray time interaction",
    "Quasi-Poisson dispersion"
  ),
  Note = c(
    "Adjusted GVIF is reported for multi-degree-of-freedom terms and is comparable with ordinary VIF for one-degree terms.",
    "The Fine-Gray proportionality table summarizes site-level cox.zph-style tests for the weighted Fine-Gray data, emphasizing global and exposure terms.",
    "The time-interaction table summarizes exposure-by-log(time) terms fit in Fine-Gray models; small p values indicate evidence that the exposure association varies over follow-up time.",
    "Quasi-Poisson dispersion >1 and quasi/Poisson SE ratios >1 indicate overdispersion relative to a Poisson variance assumption."
  )
)
write_csv(notes, file.path(out_dir, "supplement_diagnostics_notes.csv"))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  sheets <- list(
    "VIF" = vif_appendix_table,
    "FineGray PH" = fg_zph_appendix_table,
    "FineGray time interaction" = fg_time_appendix_table,
    "QuasiPoisson dispersion" = qp_dispersion_appendix_table,
    "Poisson vs QuasiPoisson" = qp_comparison_appendix_table,
    "Notes" = notes
  )
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]])
    openxlsx::setColWidths(wb, nm, cols = seq_along(sheets[[nm]]), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, file.path(out_dir, "supplement_model_diagnostics_tables.xlsx"), overwrite = TRUE)
}

message("Wrote supplement diagnostics tables to: ", out_dir)
