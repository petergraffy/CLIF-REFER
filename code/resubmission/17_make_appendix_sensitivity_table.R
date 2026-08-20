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

primary_path <- file.path(out_dir, "pooled_primary_effect_estimates.csv")
sensitivity_path <- file.path(out_dir, "pooled_sensitivity_effect_estimates.csv")
site_dirs_path <- file.path(out_dir, "pooled_site_directories.csv")

stopifnot(file.exists(primary_path), file.exists(sensitivity_path))

term_pollutant <- function(term) {
  case_when(
    str_detect(term, "pm25") ~ "PM2.5",
    str_detect(term, "no2") ~ "NO2",
    str_detect(term, "o3") ~ "O3",
    TRUE ~ term
  )
}

pool_log_effects <- function(tbl, group_cols) {
  if (!nrow(tbl)) return(tibble())
  tbl %>%
    filter(is.finite(estimate), estimate > 0, is.finite(conf_low), conf_low > 0, is.finite(conf_high), conf_high > 0) %>%
    mutate(
      beta_log = log(estimate),
      se_log = (log(conf_high) - log(conf_low)) / (2 * 1.959963984540054),
      weight = if_else(is.finite(se_log) & se_log > 0, 1 / se_log^2, NA_real_)
    ) %>%
    filter(is.finite(weight), weight > 0) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_sites = n_distinct(site),
      sites = paste(sort(unique(site)), collapse = "; "),
      n_total = sum(n, na.rm = TRUE),
      events_total = sum(events, na.rm = TRUE),
      log_estimate = sum(weight * beta_log) / sum(weight),
      se_log = sqrt(1 / sum(weight)),
      q_statistic = sum(weight * (beta_log - sum(weight * beta_log) / sum(weight))^2),
      q_df = n() - 1,
      i2 = if_else(q_df > 0 & q_statistic > q_df, 100 * (q_statistic - q_df) / q_statistic, 0),
      .groups = "drop"
    ) %>%
    mutate(
      estimate = exp(log_estimate),
      conf_low = exp(log_estimate - 1.959963984540054 * se_log),
      conf_high = exp(log_estimate + 1.959963984540054 * se_log),
      z = log_estimate / se_log,
      p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    )
}

read_covid_sensitivity <- function(site_dirs_path) {
  if (!file.exists(site_dirs_path)) return(tibble())

  site_dirs <- read_csv(site_dirs_path, show_col_types = FALSE)
  covid <- bind_rows(lapply(seq_len(nrow(site_dirs)), function(i) {
    site_dir <- site_dirs$site_dir[[i]]
    path <- file.path(site_dir, "primary_vs_exclude_peak_covid_12m_pooling_table.csv")
    if (!file.exists(path)) {
      path <- file.path(site_dir, "sensitivity_exclude_peak_covid_12m_pooling_table.csv")
    }
    if (!file.exists(path)) return(tibble())
    read_csv(path, show_col_types = FALSE) %>%
      mutate(site = site_dirs$site[[i]], site_label = site_dirs$site_label[[i]], source_file = path)
  }))

  if (!nrow(covid)) return(tibble())

  covid %>%
    filter(sensitivity == "exclude_peak_covid_12m") %>%
    mutate(
      analysis_family = "Sensitivity",
      estimate = as.numeric(estimate),
      conf_low = as.numeric(conf_low),
      conf_high = as.numeric(conf_high),
      p_value = suppressWarnings(as.numeric(p_value)),
      n = suppressWarnings(as.numeric(n)),
      events = suppressWarnings(as.numeric(events)),
      pollutant = term_pollutant(term),
      model = paste(sensitivity, analysis_model, pollutant_model, sep = " | "),
      comparison = paste(sensitivity, outcome, analysis_model, pollutant_model, term, sep = " | ")
    ) %>%
    select(site, analysis_family, sensitivity, outcome, analysis_model, model, term, pollutant, effect_measure,
           n, events, estimate, conf_low, conf_high, p_value, comparison, everything())
}

primary <- read_csv(primary_path, show_col_types = FALSE) %>%
  filter(analysis_family %in% c("Primary mortality", "Primary VFD", "Mortality Cox sensitivity")) %>%
  mutate(sensitivity = "Primary analysis")

sensitivity <- read_csv(sensitivity_path, show_col_types = FALSE)
covid_site <- read_covid_sensitivity(site_dirs_path)
covid_pooled <- pool_log_effects(
  covid_site,
  c("analysis_family", "sensitivity", "outcome", "analysis_model", "model", "term", "pollutant", "effect_measure")
)

if (nrow(covid_site)) {
  write_csv(covid_site, file.path(out_dir, "site_covid_sensitivity_effect_estimates.csv"))
  write_csv(covid_pooled, file.path(out_dir, "pooled_covid_sensitivity_effect_estimates.csv"))
}

all_results <- bind_rows(
  primary %>% mutate(source_table = "pooled_primary_effect_estimates.csv"),
  sensitivity %>% mutate(source_table = "pooled_sensitivity_effect_estimates.csv"),
  covid_pooled %>% mutate(source_table = "pooled_covid_sensitivity_effect_estimates.csv")
)

sensitivity_labels <- c(
  "Primary analysis" = "Primary analysis",
  "add_sofa_total_first_24h" = "Additional adjustment for first-24-hour SOFA total score",
  "remove_icu_los_24h_restriction" = "No restriction on ICU length of stay >=24 hours",
  "add_o3_copollutant" = "Co-pollutant model additionally adjusted for ozone",
  "exposure_window_3y" = "Alternative 3-year exposure window",
  "followup_window_14d" = "Alternative 14-day follow-up window",
  "exclude_peak_covid_12m" = "Excluded peak COVID-19 12-month period"
)

outcome_labels <- c(
  "Mortality by day 28 after ARF onset" = "Mortality by day 28",
  "Ventilator-free days through day 28" = "Ventilator-free days through day 28",
  "In-hospital mortality after ARF onset" = "In-hospital mortality after ARF onset",
  "mortality_day28_event" = "Mortality by day 28",
  "ventilator_free_days" = "Ventilator-free days through day 28",
  "mortality_time_to_event" = "In-hospital mortality after ARF onset",
  "mortality_day14_event" = "Mortality by day 14",
  "ventilator_free_days_day14" = "Ventilator-free days through day 14",
  "imv_days_through_vfd" = "IMV days through day 28",
  "imv_days_through_day14" = "IMV days through day 14",
  "IMV duration through day 28" = "IMV days through day 28"
)

effect_labels <- c(
  "odds_ratio" = "Odds ratio",
  "hazard_ratio" = "Hazard ratio",
  "ratio_of_means" = "Ratio of means"
)

pollutant_labels <- c(
  "NO2" = "NO2 per 10 ppb",
  "PM2.5" = "PM2.5 per 5 ug/m3",
  "O3" = "O3 per 10 ppb"
)

model_labels <- c(
  "NO2 single-pollutant" = "NO2 single-pollutant",
  "PM25 single-pollutant" = "PM2.5 single-pollutant",
  "PM25 + NO2" = "NO2 and PM2.5 co-pollutant",
  "NO2 3-year single-pollutant" = "NO2 3-year single-pollutant",
  "PM25 3-year single-pollutant" = "PM2.5 3-year single-pollutant",
  "PM25 + NO2 3-year" = "NO2 and PM2.5 3-year co-pollutant",
  "PM25 + NO2 + O3" = "NO2, PM2.5, and O3 co-pollutant"
)

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

fmt_ci <- function(estimate, low, high) {
  paste0(fmt_num(estimate), " (", fmt_num(low), "-", fmt_num(high), ")")
}

extract_pollutant_model <- function(model) {
  str_trim(str_replace(model, "^.*\\|", ""))
}

appendix_full <- all_results %>%
  mutate(
    analysis_model = case_when(
      !is.na(analysis_model) & nzchar(analysis_model) ~ analysis_model,
      effect_measure == "odds_ratio" ~ "Logistic regression",
      effect_measure == "hazard_ratio" ~ "Cox proportional hazards",
      effect_measure == "ratio_of_means" ~ "Quasi-Poisson regression",
      TRUE ~ analysis_model
    ),
    sensitivity_label = coalesce(recode(sensitivity, !!!sensitivity_labels), sensitivity),
    outcome_label = coalesce(recode(outcome, !!!outcome_labels), outcome),
    pollutant_label = coalesce(recode(pollutant, !!!pollutant_labels), pollutant),
    effect_label = coalesce(recode(effect_measure, !!!effect_labels), effect_measure),
    pollutant_model_raw = extract_pollutant_model(model),
    pollutant_model = coalesce(recode(pollutant_model_raw, !!!model_labels), pollutant_model_raw),
    estimate_ci = fmt_ci(estimate, conf_low, conf_high),
    p_value_formatted = fmt_p(p_value),
    i2_formatted = ifelse(is.na(i2), "", sprintf("%.1f", i2)),
    n_total_formatted = ifelse(is.na(n_total), "", format(n_total, big.mark = ",")),
    events_total_formatted = ifelse(is.na(events_total) | events_total == 0, "", format(events_total, big.mark = ","))
  ) %>%
  arrange(
    factor(
      sensitivity,
      levels = c(
        "Primary analysis",
        "add_sofa_total_first_24h",
        "remove_icu_los_24h_restriction",
        "add_o3_copollutant",
        "exposure_window_3y",
        "followup_window_14d",
        "exclude_peak_covid_12m"
      )
    ),
    factor(outcome_label, levels = c(
      "Mortality by day 28",
      "Ventilator-free days through day 28",
      "In-hospital mortality after ARF onset",
      "Mortality by day 14",
      "Ventilator-free days through day 14",
      "IMV days through day 28",
      "IMV days through day 14"
    )),
    analysis_model,
    pollutant_model,
    pollutant
  )

appendix_table <- appendix_full %>%
  transmute(
    `Sensitivity analysis` = sensitivity_label,
    Outcome = outcome_label,
    `Statistical model` = analysis_model,
    `Pollutant model` = pollutant_model,
    Contrast = pollutant_label,
    `Effect measure` = effect_label,
    `Estimate (95% CI)` = estimate_ci,
    `P value` = p_value_formatted,
    `I-squared, %` = i2_formatted,
    `Sites, n` = n_sites,
    `Participants, n` = n_total_formatted,
    `Events, n` = events_total_formatted
  )

appendix_table_primary_outcomes <- appendix_table %>%
  filter(!str_detect(Outcome, "^IMV days"))

comparison_table <- appendix_full %>%
  mutate(
    outcome_family = case_when(
      outcome_label %in% c("Mortality by day 28", "Mortality by day 14") ~ "Mortality",
      outcome_label %in% c("Ventilator-free days through day 28", "Ventilator-free days through day 14") ~ "Ventilator-free days",
      outcome_label == "In-hospital mortality after ARF onset" ~ "Mortality time-to-event",
      outcome_label %in% c("IMV days through day 28", "IMV days through day 14") ~ "IMV days",
      TRUE ~ outcome_label
    ),
    follow_up = case_when(
      outcome_label %in% c("Mortality by day 14", "Ventilator-free days through day 14", "IMV days through day 14") ~ "Day 14",
      outcome_label %in% c("Mortality by day 28", "Ventilator-free days through day 28", "IMV days through day 28") ~ "Day 28",
      outcome_label == "In-hospital mortality after ARF onset" ~ "In-hospital",
      TRUE ~ ""
    ),
    pollutant_model_comparable = case_when(
      str_detect(pollutant_model, "O3") ~ "NO2 + PM2.5 + O3 co-pollutant",
      str_detect(pollutant_model, "co-pollutant") ~ "NO2 + PM2.5 co-pollutant",
      str_detect(pollutant_model, "single-pollutant") ~ "Single-pollutant",
      TRUE ~ pollutant_model
    ),
    sensitivity_column = case_when(
      sensitivity == "Primary analysis" ~ "Primary",
      sensitivity == "add_sofa_total_first_24h" ~ "SOFA adjusted",
      sensitivity == "remove_icu_los_24h_restriction" ~ "No ICU LOS >=24h restriction",
      sensitivity == "add_o3_copollutant" ~ "O3 co-pollutant",
      sensitivity == "exposure_window_3y" ~ "3-year exposure",
      sensitivity == "followup_window_14d" ~ "14-day follow-up",
      sensitivity == "exclude_peak_covid_12m" ~ "Excluding peak COVID-19",
      TRUE ~ sensitivity_label
    ),
    comparison_value = paste0(
      estimate_ci,
      "; p=", p_value_formatted,
      "; I2=", i2_formatted, "%"
    )
  ) %>%
  transmute(
    `Outcome family` = outcome_family,
    `Statistical model` = analysis_model,
    `Effect measure` = effect_label,
    `Pollutant model` = pollutant_model_comparable,
    Contrast = pollutant_label,
    `Sensitivity column` = sensitivity_column,
    `Estimate (95% CI); p; I-squared` = comparison_value,
    `Sites, n` = n_sites,
    `Participants, n` = n_total_formatted,
    `Events, n` = events_total_formatted,
    follow_up
  ) %>%
  distinct()

comparison_wide <- comparison_table %>%
  select(
    `Outcome family`,
    `Statistical model`,
    `Effect measure`,
    `Pollutant model`,
    Contrast,
    `Sensitivity column`,
    `Estimate (95% CI); p; I-squared`
  ) %>%
  pivot_wider(
    names_from = `Sensitivity column`,
    values_from = `Estimate (95% CI); p; I-squared`
  ) %>%
  select(
    `Outcome family`,
    `Statistical model`,
    `Effect measure`,
    `Pollutant model`,
    Contrast,
    any_of(c(
      "Primary",
      "SOFA adjusted",
      "No ICU LOS >=24h restriction",
      "O3 co-pollutant",
      "3-year exposure",
      "14-day follow-up",
      "Excluding peak COVID-19"
    ))
  ) %>%
  arrange(
    factor(`Outcome family`, levels = c("Mortality", "Ventilator-free days", "Mortality time-to-event", "IMV days")),
    factor(`Statistical model`, levels = c("Logistic regression", "Quasi-Poisson regression", "Cox proportional hazards")),
    factor(`Pollutant model`, levels = c("Single-pollutant", "NO2 + PM2.5 co-pollutant", "NO2 + PM2.5 + O3 co-pollutant")),
    Contrast
  )

comparison_compact <- comparison_wide %>%
  filter(`Outcome family` %in% c("Mortality", "Ventilator-free days", "Mortality time-to-event"))

comparison_core <- comparison_compact %>%
  filter(`Pollutant model` != "NO2 + PM2.5 + O3 co-pollutant") %>%
  select(-any_of("O3 co-pollutant"))

comparison_o3 <- comparison_compact %>%
  filter(`Pollutant model` == "NO2 + PM2.5 + O3 co-pollutant")

write_csv(appendix_full, file.path(out_dir, "appendix_sensitivity_results_full.csv"))
write_csv(appendix_table, file.path(out_dir, "appendix_sensitivity_results_table.csv"))
write_csv(appendix_table_primary_outcomes, file.path(out_dir, "appendix_sensitivity_results_primary_outcomes_table.csv"))
write_csv(comparison_wide, file.path(out_dir, "appendix_sensitivity_results_comparison_wide.csv"))
write_csv(comparison_compact, file.path(out_dir, "appendix_sensitivity_results_comparison_wide_primary_outcomes.csv"))
write_csv(comparison_core, file.path(out_dir, "appendix_sensitivity_results_comparison_core.csv"))
write_csv(comparison_o3, file.path(out_dir, "appendix_sensitivity_results_o3_copollutant.csv"))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Appendix sensitivity")
  openxlsx::writeData(wb, "Appendix sensitivity", appendix_table)
  openxlsx::setColWidths(wb, "Appendix sensitivity", cols = seq_along(appendix_table), widths = "auto")

  openxlsx::addWorksheet(wb, "Primary outcomes only")
  openxlsx::writeData(wb, "Primary outcomes only", appendix_table_primary_outcomes)
  openxlsx::setColWidths(wb, "Primary outcomes only", cols = seq_along(appendix_table_primary_outcomes), widths = "auto")

  openxlsx::addWorksheet(wb, "Full numeric")
  openxlsx::writeData(wb, "Full numeric", appendix_full)
  openxlsx::setColWidths(wb, "Full numeric", cols = seq_along(appendix_full), widths = "auto")

  openxlsx::addWorksheet(wb, "Comparable wide")
  openxlsx::writeData(wb, "Comparable wide", comparison_wide)
  openxlsx::setColWidths(wb, "Comparable wide", cols = seq_along(comparison_wide), widths = "auto")

  openxlsx::addWorksheet(wb, "Comparable primary")
  openxlsx::writeData(wb, "Comparable primary", comparison_compact)
  openxlsx::setColWidths(wb, "Comparable primary", cols = seq_along(comparison_compact), widths = "auto")

  openxlsx::addWorksheet(wb, "Comparable core")
  openxlsx::writeData(wb, "Comparable core", comparison_core)
  openxlsx::setColWidths(wb, "Comparable core", cols = seq_along(comparison_core), widths = "auto")

  openxlsx::addWorksheet(wb, "O3 copollutant")
  openxlsx::writeData(wb, "O3 copollutant", comparison_o3)
  openxlsx::setColWidths(wb, "O3 copollutant", cols = seq_along(comparison_o3), widths = "auto")

  openxlsx::saveWorkbook(wb, file.path(out_dir, "appendix_sensitivity_results_table.xlsx"), overwrite = TRUE)
}

notes <- tibble(
  Note = c(
    "Estimates are pooled across CLIF sites using the same meta-analytic pooling approach as the primary pooled analysis.",
    "NO2 contrasts are scaled per 10 ppb, PM2.5 contrasts per 5 ug/m3, and O3 contrasts per 10 ppb.",
    "Mortality by day 28 uses logistic regression; ventilator-free days and IMV duration use quasi-Poisson regression; time-to-death sensitivity models use Cox proportional hazards regression.",
    "The COVID-19 sensitivity excludes ARF onsets from 2020-03-01 through 2021-02-28.",
    "The primary-outcomes-only CSV omits supplemental IMV duration rows; the full table preserves all sensitivity rows exported by sites."
  )
)

write_csv(notes, file.path(out_dir, "appendix_sensitivity_results_notes.csv"))

message("Wrote appendix sensitivity tables to: ", out_dir)
