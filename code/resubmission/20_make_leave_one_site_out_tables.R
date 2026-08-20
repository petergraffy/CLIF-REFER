#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1 && nzchar(args[[1]])) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(file.path("output", "resubmission_pooled", "20260819_all_sites_with_michigan"), mustWork = TRUE)
}

fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

site_primary_path <- file.path(out_dir, "site_primary_effect_estimates.csv")
site_fg_path <- file.path(out_dir, "site_fine_gray_effect_estimates.csv")
site_sens_path <- file.path(out_dir, "site_sensitivity_effect_estimates.csv")
stopifnot(file.exists(site_primary_path), file.exists(site_fg_path))

site_map <- read_csv(file.path(out_dir, "pooled_site_directories.csv"), show_col_types = FALSE) %>%
  arrange(site_label) %>%
  mutate(site_anon = paste0("Site ", row_number())) %>%
  select(site, site_label, site_anon)

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

fmt_ci <- function(est, lo, hi) {
  paste0(fmt_num(est), " (", fmt_num(lo), "-", fmt_num(hi), ")")
}

pollutant_label <- function(x) {
  recode(
    x,
    "NO2" = "NO2 per 10 ppb",
    "PM2.5" = "PM2.5 per 5 ug/m3",
    "O3" = "O3 per 10 ppb",
    .default = x
  )
}

model_label <- function(x) {
  recode(
    x,
    "NO2 single-pollutant" = "NO2 single-pollutant",
    "PM25 single-pollutant" = "PM2.5 single-pollutant",
    "PM25 + NO2" = "NO2 + PM2.5 co-pollutant",
    "PM25 + NO2 + O3" = "NO2 + PM2.5 + O3 co-pollutant",
    .default = x
  )
}

outcome_label <- function(x) {
  recode(
    x,
    "Mortality by day 28 after ARF onset" = "Mortality by day 28",
    "Ventilator-free days through day 28" = "Ventilator-free days through day 28",
    "IMV duration through day 28" = "IMV duration through day 28",
    "In-hospital mortality after ARF onset" = "In-hospital mortality after ARF onset",
    .default = x
  )
}

effect_label <- function(x) {
  recode(
    x,
    "odds_ratio" = "OR",
    "ratio_of_means" = "Ratio of means",
    "hazard_ratio" = "HR",
    "subdistribution_hazard_ratio" = "SHR",
    .default = x
  )
}

prep_site_effects <- function(df, source) {
  df %>%
    filter(is.finite(estimate), estimate > 0, is.finite(conf_low), conf_low > 0, is.finite(conf_high), conf_high > 0) %>%
    mutate(
      source = source,
      beta_log = log(estimate),
      se_log = (log(conf_high) - log(conf_low)) / (2 * 1.959963984540054),
      weight_fixed = if_else(is.finite(se_log) & se_log > 0, 1 / se_log^2, NA_real_),
      site_key = site
    ) %>%
    left_join(site_map, by = c("site" = "site")) %>%
    mutate(
      site_label = coalesce(site_label.y, site_label.x, site),
      site_anon = coalesce(site_anon, site_label, site)
    ) %>%
    select(-any_of(c("site_label.x", "site_label.y"))) %>%
    filter(is.finite(weight_fixed), weight_fixed > 0)
}

pool_one <- function(df) {
  if (!nrow(df)) return(tibble())
  k <- nrow(df)
  beta_fixed <- sum(df$weight_fixed * df$beta_log) / sum(df$weight_fixed)
  se_fixed <- sqrt(1 / sum(df$weight_fixed))
  q <- sum(df$weight_fixed * (df$beta_log - beta_fixed)^2)
  q_df <- k - 1
  c_val <- sum(df$weight_fixed) - sum(df$weight_fixed^2) / sum(df$weight_fixed)
  tau2_dl <- if (q_df > 0 && c_val > 0) max(0, (q - q_df) / c_val) else 0
  w_re <- 1 / (df$se_log^2 + tau2_dl)
  beta_re <- sum(w_re * df$beta_log) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))
  tibble(
    n_sites = n_distinct(df$site),
    sites = paste(sort(unique(df$site)), collapse = "; "),
    n_total = sum(df$n, na.rm = TRUE),
    events_total = sum(df$events, na.rm = TRUE),
    log_estimate = beta_fixed,
    se_log = se_fixed,
    q_statistic = q,
    q_df = q_df,
    i2 = if_else(q_df > 0 && q > q_df, 100 * (q - q_df) / q, 0),
    tau2_dl = tau2_dl,
    estimate = exp(beta_fixed),
    conf_low = exp(beta_fixed - 1.959963984540054 * se_fixed),
    conf_high = exp(beta_fixed + 1.959963984540054 * se_fixed),
    z = beta_fixed / se_fixed,
    p_value = 2 * stats::pnorm(abs(beta_fixed / se_fixed), lower.tail = FALSE),
    estimate_random_dl = exp(beta_re),
    conf_low_random_dl = exp(beta_re - 1.959963984540054 * se_re),
    conf_high_random_dl = exp(beta_re + 1.959963984540054 * se_re),
    p_value_random_dl = 2 * stats::pnorm(abs(beta_re / se_re), lower.tail = FALSE)
  )
}

loo_by_group <- function(df, group_cols) {
  groups <- df %>% distinct(across(all_of(group_cols)))
  bind_rows(lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, , drop = FALSE]
    one <- df
    for (col in group_cols) {
      one <- one %>% filter(.data[[col]] == g[[col]][[1]])
    }
    full <- pool_one(one) %>% mutate(omitted_site = "None", omitted_site_label = "None", omitted_site_anon = "None")
    loo <- bind_rows(lapply(sort(unique(one$site)), function(s) {
      left_out <- one %>% filter(site != s)
      omitted <- one %>% filter(site == s) %>% slice(1)
      pool_one(left_out) %>%
        mutate(
          omitted_site = s,
          omitted_site_label = omitted$site_label,
          omitted_site_anon = omitted$site_anon
        )
    }))
    bind_cols(g[rep(1, 1 + nrow(loo)), , drop = FALSE], bind_rows(full, loo))
  }))
}

primary_site <- read_csv(site_primary_path, show_col_types = FALSE) %>%
  filter(analysis_family %in% c("Primary mortality", "Primary VFD", "Mortality Cox sensitivity")) %>%
  prep_site_effects("Primary/Cox")

fg_site <- read_csv(site_fg_path, show_col_types = FALSE) %>%
  prep_site_effects("Fine-Gray")

sens_site <- if (file.exists(site_sens_path)) {
  read_csv(site_sens_path, show_col_types = FALSE) %>%
    filter(sensitivity %in% c("add_sofa_total_first_24h", "remove_icu_los_24h_restriction", "exposure_window_3y", "followup_window_14d")) %>%
    prep_site_effects("Sensitivity")
} else {
  tibble()
}

primary_loo <- loo_by_group(primary_site, c("source", "analysis_family", "outcome", "model", "term", "pollutant", "effect_measure"))
fg_loo <- loo_by_group(fg_site, c("source", "analysis_family", "outcome", "model", "term", "pollutant", "effect_measure"))
sens_loo <- if (nrow(sens_site)) {
  loo_by_group(sens_site, c("source", "analysis_family", "sensitivity", "outcome", "analysis_model", "model", "term", "pollutant", "effect_measure"))
} else {
  tibble()
}

all_loo <- bind_rows(primary_loo, fg_loo, sens_loo) %>%
  group_by(source, analysis_family, outcome, model, term, pollutant, effect_measure, across(any_of(c("sensitivity", "analysis_model")))) %>%
  mutate(
    full_estimate = estimate[omitted_site == "None"][1],
    full_conf_low = conf_low[omitted_site == "None"][1],
    full_conf_high = conf_high[omitted_site == "None"][1],
    full_i2 = i2[omitted_site == "None"][1],
    percent_change_vs_full = 100 * (estimate - full_estimate) / full_estimate,
    crosses_null = case_when(
      effect_measure %in% c("odds_ratio", "ratio_of_means", "hazard_ratio", "subdistribution_hazard_ratio") ~ conf_low <= 1 & conf_high >= 1,
      TRUE ~ NA
    )
  ) %>%
  ungroup() %>%
  mutate(
    outcome_display = outcome_label(outcome),
    model_display = model_label(model),
    contrast_display = pollutant_label(pollutant),
    effect_display = effect_label(effect_measure),
    estimate_ci = fmt_ci(estimate, conf_low, conf_high),
    random_dl_estimate_ci = fmt_ci(estimate_random_dl, conf_low_random_dl, conf_high_random_dl),
    p_display = fmt_p(p_value),
    i2_display = sprintf("%.1f", i2),
    full_estimate_ci = fmt_ci(full_estimate, full_conf_low, full_conf_high)
  )

loo_only <- all_loo %>% filter(omitted_site != "None")

loo_summary <- loo_only %>%
  group_by(source, analysis_family, outcome, model, term, pollutant, effect_measure, across(any_of(c("sensitivity", "analysis_model")))) %>%
  summarise(
    outcome_display = first(outcome_label(outcome)),
    model_display = first(model_label(model)),
    contrast_display = first(pollutant_label(pollutant)),
    effect_display = first(effect_label(effect_measure)),
    full_estimate = first(full_estimate),
    full_conf_low = first(full_conf_low),
    full_conf_high = first(full_conf_high),
    full_i2 = first(full_i2),
    min_leave_one_estimate = min(estimate, na.rm = TRUE),
    max_leave_one_estimate = max(estimate, na.rm = TRUE),
    min_leave_one_i2 = min(i2, na.rm = TRUE),
    max_leave_one_i2 = max(i2, na.rm = TRUE),
    max_abs_percent_change = max(abs(percent_change_vs_full), na.rm = TRUE),
    site_max_change = omitted_site_anon[which.max(abs(percent_change_vs_full))],
    site_label_max_change = omitted_site_label[which.max(abs(percent_change_vs_full))],
    n_leave_one_cross_null = sum(crosses_null, na.rm = TRUE),
    n_leave_one_significant = sum(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    full_estimate_ci = fmt_ci(full_estimate, full_conf_low, full_conf_high),
    leave_one_estimate_range = paste0(fmt_num(min_leave_one_estimate), "-", fmt_num(max_leave_one_estimate)),
    leave_one_i2_range = paste0(sprintf("%.1f", min_leave_one_i2), "-", sprintf("%.1f", max_leave_one_i2)),
    max_abs_percent_change_display = paste0(sprintf("%.1f", max_abs_percent_change), "%"),
    full_i2_display = sprintf("%.1f", full_i2)
  ) %>%
  arrange(source, outcome_display, model_display, contrast_display)

primary_summary_table <- loo_summary %>%
  filter(analysis_family %in% c("Primary mortality", "Primary VFD", "Mortality Cox sensitivity")) %>%
  transmute(
    Analysis = recode(
      analysis_family,
      "Primary mortality" = "Primary mortality",
      "Primary VFD" = "Primary ventilator-free days",
      "Mortality Cox sensitivity" = "Cox mortality sensitivity"
    ),
    Outcome = outcome_display,
    Model = model_display,
    Contrast = contrast_display,
    Estimand = effect_display,
    `All-site estimate (95% CI)` = full_estimate_ci,
    `Leave-one-site-out estimate range` = leave_one_estimate_range,
    `All-site I2, %` = full_i2_display,
    `Leave-one-site-out I2 range, %` = leave_one_i2_range,
    `Maximum absolute change` = max_abs_percent_change_display,
    `Omitted site with largest change` = site_max_change,
    `Leave-one estimates crossing null, n` = n_leave_one_cross_null,
    `Leave-one estimates p<0.05, n` = n_leave_one_significant
  )

fg_summary_table <- loo_summary %>%
  filter(analysis_family == "Adjusted Fine-Gray") %>%
  transmute(
    Outcome = outcome_display,
    Model = model_display,
    Contrast = contrast_display,
    Estimand = effect_display,
    `All-site estimate (95% CI)` = full_estimate_ci,
    `Leave-one-site-out estimate range` = leave_one_estimate_range,
    `All-site I2, %` = full_i2_display,
    `Leave-one-site-out I2 range, %` = leave_one_i2_range,
    `Maximum absolute change` = max_abs_percent_change_display,
    `Omitted site with largest change` = site_max_change,
    `Leave-one estimates crossing null, n` = n_leave_one_cross_null,
    `Leave-one estimates p<0.05, n` = n_leave_one_significant
  )

write_csv(all_loo, file.path(out_dir, "leave_one_site_out_all_estimates.csv"))
write_csv(loo_summary, file.path(out_dir, "leave_one_site_out_summary.csv"))
write_csv(primary_summary_table, file.path(out_dir, "leave_one_site_out_primary_summary_table.csv"))
write_csv(fg_summary_table, file.path(out_dir, "leave_one_site_out_fine_gray_summary_table.csv"))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  sheets <- list(
    "Primary summary" = primary_summary_table,
    "FineGray summary" = fg_summary_table,
    "All LOO summary" = loo_summary,
    "All LOO estimates" = all_loo
  )
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]])
    openxlsx::setColWidths(wb, nm, cols = seq_along(sheets[[nm]]), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, file.path(out_dir, "leave_one_site_out_tables.xlsx"), overwrite = TRUE)
}

plot_df <- all_loo %>%
  filter(
    omitted_site != "None",
    analysis_family %in% c("Primary mortality", "Primary VFD"),
    model %in% c("NO2 single-pollutant", "PM25 single-pollutant")
  ) %>%
  mutate(
    panel = paste(outcome_display, contrast_display, sep = "\n"),
    omitted_site_anon = factor(omitted_site_anon, levels = paste0("Site ", seq_len(n_distinct(site_map$site_anon))))
  )

if (nrow(plot_df)) {
  p <- ggplot(plot_df, aes(x = omitted_site_anon, y = estimate)) +
    geom_hline(aes(yintercept = full_estimate), color = "grey35", linewidth = 0.35) +
    geom_point(color = "#0072B2", size = 2.1) +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high), color = "#0072B2", width = 0, linewidth = 0.45) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    labs(x = "Omitted site", y = "Leave-one-site-out estimate") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1.1, "lines"),
      plot.margin = margin(10, 12, 10, 10)
    )
  ggsave(file.path(fig_dir, "leave_one_site_out_primary_influence.png"), p, width = 10, height = 8, dpi = 300)
  ggsave(file.path(fig_dir, "leave_one_site_out_primary_influence.pdf"), p, width = 10, height = 8)
}

message("Wrote leave-one-site-out outputs to: ", out_dir)
