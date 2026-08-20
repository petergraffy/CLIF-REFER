#!/usr/bin/env Rscript

# Pool aggregate site outputs from the REFER resubmission pipeline.
#
# Usage:
#   Rscript --vanilla code/resubmission/08_pool_site_outputs.R \
#     <ucmc_output_dir> <box_project_dir> [pooled_output_dir]
#
# The script intentionally uses aggregate result tables and prediction grids
# rather than row-level analysis datasets.

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(purrr)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

args <- commandArgs(trailingOnly = TRUE)
arg_or <- function(i, default = NA_character_) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
repo <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)

ucmc_dir <- normalizePath(arg_or(1, file.path(repo, "output", "resubmission", "20260807_075932")), mustWork = TRUE)
box_project_dir <- normalizePath(arg_or(2, "/Users/saborpete/Library/CloudStorage/Box-Box/CLIF/Projects/CLIF-REFER_v1"), mustWork = TRUE)
out_dir <- arg_or(3, file.path(repo, "output", "resubmission_pooled", format(Sys.time(), "%Y%m%d_%H%M%S")))
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("UCMC output: ", ucmc_dir)
message("Box project: ", box_project_dir)
message("Pooled output: ", out_dir)

read_optional <- function(site_dir, filename) {
  path <- file.path(site_dir, filename)
  if (!file.exists(path)) return(tibble())
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      across(matches("(formatted|display|label|notes|source|adjustment|terms)$"), as.character),
      source_dir = site_dir,
      source_file = filename
    )
}

latest_site_run <- function(site_root) {
  dirs <- list.dirs(site_root, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "resubmission_pipeline_manifest.csv")) |
                 file.exists(file.path(dirs, "primary_mortality_day28_logistic_results.csv"))]
  if (!length(dirs)) return(NA_character_)
  dirs[order(basename(dirs), decreasing = TRUE)][[1]]
}

discover_box_site_runs <- function(box_dir) {
  roots <- list.dirs(box_dir, recursive = FALSE, full.names = TRUE)
  roots <- roots[!basename(roots) %in% c("archive")]
  roots <- roots[!str_detect(basename(roots), regex("^\\.|ucmc", ignore_case = TRUE))]
  tibble(
    site_label = basename(roots),
    site_dir = map_chr(roots, latest_site_run)
  ) %>%
    filter(!is.na(site_dir), dir.exists(site_dir))
}

infer_site <- function(site_dir) {
  candidates <- c(
    "site_inclusion_flow_counts.csv",
    "primary_mortality_day28_logistic_results.csv",
    "cohort_summary_reviewer_optimized.csv",
    "resubmission_cohort_summary.csv"
  )
  for (filename in candidates) {
    x <- read_optional(site_dir, filename)
    if (nrow(x) && "site" %in% names(x)) {
      site <- x$site[!is.na(x$site)][[1]] %||% NA_character_
      if (!is.na(site)) return(as.character(site))
    }
  }
  basename(dirname(site_dir))
}

site_dirs <- bind_rows(
  tibble(site_label = "UCMC", site_dir = ucmc_dir),
  discover_box_site_runs(box_project_dir)
) %>%
  filter(!is.na(site_dir), dir.exists(site_dir)) %>%
  mutate(site = map_chr(site_dir, infer_site)) %>%
  distinct(site, .keep_all = TRUE)

if (nrow(site_dirs) < 2) {
  stop("Need at least two site output directories; found: ", nrow(site_dirs))
}

read_sites <- function(filename) {
  site_dirs %>%
    mutate(data = map(site_dir, read_optional, filename = filename)) %>%
    select(site_label, site_expected = site, site_dir, data) %>%
    tidyr::unnest(data) %>%
    mutate(site = if ("site" %in% names(.)) coalesce(as.character(.data$site), site_expected) else site_expected)
}

clean_pollutant <- function(x) {
  case_when(
    str_detect(x, regex("pm", ignore_case = TRUE)) ~ "PM2.5",
    str_detect(x, regex("no2|nitrogen", ignore_case = TRUE)) ~ "NO2",
    str_detect(x, regex("o3|ozone", ignore_case = TRUE)) ~ "O3",
    TRUE ~ as.character(x)
  )
}

term_pollutant <- function(term) {
  case_when(
    str_detect(term, "pm25") ~ "PM2.5",
    str_detect(term, "no2") ~ "NO2",
    str_detect(term, "o3") ~ "O3",
    TRUE ~ term
  )
}

normalise_effect_table <- function(tbl, analysis_family, estimate_col, effect_measure) {
  if (!nrow(tbl)) return(tibble())
  tbl %>%
    mutate(
      analysis_family = analysis_family,
      estimate = as.numeric(.data[[estimate_col]]),
      conf_low = as.numeric(conf_low),
      conf_high = as.numeric(conf_high),
      p_value = suppressWarnings(as.numeric(p_value)),
      effect_measure = effect_measure,
      pollutant = term_pollutant(term),
      comparison = paste(outcome, model, term, sep = " | "),
      n = suppressWarnings(as.numeric(n)),
      events = if ("events" %in% names(.)) suppressWarnings(as.numeric(events)) else NA_real_
    ) %>%
    select(site, analysis_family, outcome, model, term, pollutant, effect_measure, n, events,
           estimate, conf_low, conf_high, p_value, comparison, everything())
}

primary_effects <- bind_rows(
  normalise_effect_table(read_sites("primary_mortality_day28_logistic_results.csv"), "Primary mortality", "odds_ratio", "odds_ratio"),
  normalise_effect_table(read_sites("primary_vfd_quasipoisson_results.csv"), "Primary VFD", "ratio_of_means", "ratio_of_means"),
  normalise_effect_table(read_sites("primary_imv_duration_quasipoisson_results.csv"), "Supplemental IMV duration", "ratio_of_means", "ratio_of_means"),
  normalise_effect_table(read_sites("primary_mortality_cox_results.csv"), "Mortality Cox sensitivity", "hazard_ratio", "hazard_ratio")
)

fine_gray_effects <- normalise_effect_table(
  read_sites("fine_gray_results_same_covariates_as_cox.csv"),
  "Adjusted Fine-Gray",
  "subdistribution_hazard_ratio",
  "subdistribution_hazard_ratio"
)

sensitivity_effects <- read_sites("primary_sensitivity_models_pooling_table.csv") %>%
  mutate(
    analysis_family = "Sensitivity",
    estimate = as.numeric(estimate),
    conf_low = as.numeric(conf_low),
    conf_high = as.numeric(conf_high),
    p_value = suppressWarnings(as.numeric(p_value)),
    pollutant = term_pollutant(term),
    model = paste(sensitivity, analysis_model, pollutant_model, sep = " | "),
    comparison = paste(sensitivity, outcome, analysis_model, pollutant_model, term, sep = " | ")
  ) %>%
  select(site, analysis_family, sensitivity, outcome, analysis_model, model, term, pollutant, effect_measure,
         n, events, estimate, conf_low, conf_high, p_value, comparison, everything())

subgroup_effects <- read_sites("primary_subgroup_specific_estimates.csv") %>%
  filter(!str_detect(outcome, regex("imv|invasive|ventilation duration", ignore_case = TRUE))) %>%
  mutate(
    analysis_family = "Subgroup",
    model = paste(subgroup, subgroup_level, sep = ": "),
    comparison = paste(outcome, pollutant, subgroup, subgroup_level, sep = " | "),
    estimate = as.numeric(estimate),
    conf_low = as.numeric(conf_low),
    conf_high = as.numeric(conf_high),
    p_value = suppressWarnings(as.numeric(p_value)),
    effect_measure = estimand
  ) %>%
  select(site, analysis_family, outcome, model, pollutant, subgroup, subgroup_level, exposure_term,
         effect_measure, n, events, estimate, conf_low, conf_high, p_value, comparison, everything())

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

primary_pooled <- pool_log_effects(primary_effects, c("analysis_family", "outcome", "model", "term", "pollutant", "effect_measure"))
fine_gray_pooled <- pool_log_effects(fine_gray_effects, c("analysis_family", "outcome", "model", "term", "pollutant", "effect_measure"))
sensitivity_pooled <- pool_log_effects(sensitivity_effects, c("analysis_family", "sensitivity", "outcome", "analysis_model", "model", "term", "pollutant", "effect_measure"))
subgroup_pooled <- pool_log_effects(subgroup_effects, c("analysis_family", "outcome", "pollutant", "subgroup", "subgroup_level", "exposure_term", "effect_measure"))

readr::write_csv(site_dirs, file.path(out_dir, "pooled_site_directories.csv"))
readr::write_csv(primary_effects, file.path(out_dir, "site_primary_effect_estimates.csv"))
readr::write_csv(primary_pooled, file.path(out_dir, "pooled_primary_effect_estimates.csv"))
readr::write_csv(fine_gray_effects, file.path(out_dir, "site_fine_gray_effect_estimates.csv"))
readr::write_csv(fine_gray_pooled, file.path(out_dir, "pooled_fine_gray_effect_estimates.csv"))
readr::write_csv(sensitivity_effects, file.path(out_dir, "site_sensitivity_effect_estimates.csv"))
readr::write_csv(sensitivity_pooled, file.path(out_dir, "pooled_sensitivity_effect_estimates.csv"))
readr::write_csv(subgroup_effects, file.path(out_dir, "site_subgroup_effect_estimates.csv"))
readr::write_csv(subgroup_pooled, file.path(out_dir, "pooled_subgroup_effect_estimates.csv"))

format_n_pct <- function(n, denom, digits = 1) {
  ifelse(
    is.finite(n) & is.finite(denom) & denom > 0,
    paste0(format(round(n), big.mark = ","), " (", sprintf(paste0("%.", digits, "f"), 100 * n / denom), "%)"),
    NA_character_
  )
}

format_mean_sd <- function(mean, sd, digits = 1) {
  ifelse(
    is.finite(mean) & is.finite(sd),
    paste0(sprintf(paste0("%.", digits, "f"), mean), " +/- ", sprintf(paste0("%.", digits, "f"), sd)),
    NA_character_
  )
}

pool_continuous_table1 <- function(tbl) {
  if (!nrow(tbl)) return(tibble())
  tbl %>%
    mutate(
      site_n = as.numeric(n),
      mean_value = as.numeric(.data$mean),
      sd_value = as.numeric(.data$sd),
      q25_value = as.numeric(.data$q25),
      median_value = as.numeric(.data$median),
      q75_value = as.numeric(.data$q75)
    ) %>%
    filter(is.finite(site_n), site_n > 0) %>%
    group_by(table, variable, label) %>%
    summarise(
      n_sites = n_distinct(site),
      n = sum(site_n, na.rm = TRUE),
      mean = sum(site_n * mean_value, na.rm = TRUE) / sum(site_n, na.rm = TRUE),
      sd = {
        pooled_n <- sum(site_n, na.rm = TRUE)
        pooled_mean <- sum(site_n * mean_value, na.rm = TRUE) / pooled_n
        sqrt(sum((site_n - 1) * sd_value^2 + site_n * (mean_value - pooled_mean)^2, na.rm = TRUE) / pmax(pooled_n - 1, 1))
      },
      q25_site_weighted = sum(site_n * q25_value, na.rm = TRUE) / sum(site_n, na.rm = TRUE),
      median_site_weighted = sum(site_n * median_value, na.rm = TRUE) / sum(site_n, na.rm = TRUE),
      q75_site_weighted = sum(site_n * q75_value, na.rm = TRUE) / sum(site_n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      variable_type = "continuous",
      level = NA_character_,
      denominator = n,
      pct = NA_real_,
      display = paste0(
        format_mean_sd(mean, sd),
        " (median approx. ",
        sprintf("%.1f", median_site_weighted),
        "; IQR approx. ",
        sprintf("%.1f", q25_site_weighted),
        ", ",
        sprintf("%.1f", q75_site_weighted),
        ")"
      )
    ) %>%
    select(table, variable_type, variable, label, level, n_sites, n, denominator, pct,
           mean, sd, q25_site_weighted, median_site_weighted, q75_site_weighted, display)
}

pool_categorical_table1 <- function(tbl) {
  if (!nrow(tbl)) return(tibble())
  variable_denominators <- tbl %>%
    mutate(n = as.numeric(n), total_n = as.numeric(total_n)) %>%
    group_by(table, variable, label, site) %>%
    summarise(site_total_n = max(total_n, na.rm = TRUE), .groups = "drop") %>%
    group_by(table, variable, label) %>%
    summarise(denominator = sum(site_total_n, na.rm = TRUE), .groups = "drop")

  tbl %>%
    mutate(n = as.numeric(n), total_n = as.numeric(total_n)) %>%
    group_by(table, variable, label, level) %>%
    summarise(
      n_sites = n_distinct(site),
      n = sum(n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(variable_denominators, by = c("table", "variable", "label")) %>%
    mutate(
      variable_type = "categorical",
      pct = if_else(denominator > 0, 100 * n / denominator, NA_real_),
      mean = NA_real_,
      sd = NA_real_,
      q25_site_weighted = NA_real_,
      median_site_weighted = NA_real_,
      q75_site_weighted = NA_real_,
      display = format_n_pct(n, denominator)
    ) %>%
    select(table, variable_type, variable, label, level, n_sites, n, denominator, pct,
           mean, sd, q25_site_weighted, median_site_weighted, q75_site_weighted, display)
}

table1_continuous_site <- read_sites("table1_continuous_site_resubmission.csv")
table1_categorical_site <- read_sites("table1_categorical_site_resubmission.csv")
table1_chronic_site <- read_sites("table1_baseline_chronic_disease_prevalence.csv")

table1_continuous_pooled <- pool_continuous_table1(table1_continuous_site)
table1_categorical_pooled <- pool_categorical_table1(table1_categorical_site)
table1_chronic_pooled <- table1_chronic_site %>%
  transmute(table, domain, variable, label, n = as.numeric(n), denominator = as.numeric(denominator), site) %>%
  group_by(table, domain, variable, label) %>%
  summarise(
    n_sites = n_distinct(site),
    n = sum(n, na.rm = TRUE),
    denominator = sum(denominator, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prevalence = if_else(denominator > 0, n / denominator, NA_real_),
    prevalence_percent = 100 * prevalence,
    display = format_n_pct(n, denominator)
  )
table1_pooled <- bind_rows(table1_continuous_pooled, table1_categorical_pooled)

readr::write_csv(table1_continuous_site, file.path(out_dir, "site_table1_continuous_resubmission.csv"))
readr::write_csv(table1_categorical_site, file.path(out_dir, "site_table1_categorical_resubmission.csv"))
readr::write_csv(table1_chronic_site, file.path(out_dir, "site_table1_baseline_chronic_disease_prevalence.csv"))
readr::write_csv(table1_continuous_pooled, file.path(out_dir, "pooled_table1_continuous_resubmission.csv"))
readr::write_csv(table1_categorical_pooled, file.path(out_dir, "pooled_table1_categorical_resubmission.csv"))
readr::write_csv(table1_chronic_pooled, file.path(out_dir, "pooled_table1_baseline_chronic_disease_prevalence.csv"))
readr::write_csv(table1_pooled, file.path(out_dir, "pooled_table1_resubmission_long.csv"))

inclusion_flow <- read_sites("site_inclusion_flow_counts.csv")
readr::write_csv(inclusion_flow, file.path(out_dir, "pooled_site_inclusion_flow_counts_long.csv"))
inclusion_summary <- inclusion_flow %>%
  group_by(step_order, step) %>%
  summarise(
    n_sites = n_distinct(site),
    n_hospitalizations_or_encounters = sum(n_hospitalizations_or_encounters, na.rm = TRUE),
    n_patients = if (all(is.na(n_patients))) NA_integer_ else sum(n_patients, na.rm = TRUE),
    n_excluded_from_prior = if (all(is.na(n_excluded_from_prior))) NA_integer_ else sum(n_excluded_from_prior, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(inclusion_summary, file.path(out_dir, "pooled_inclusion_flow_counts.csv"))

subgroup_interactions <- read_sites("primary_subgroup_interaction_tests.csv") %>%
  filter(!str_detect(outcome, regex("imv|invasive|ventilation duration", ignore_case = TRUE))) %>%
  mutate(p_interaction = suppressWarnings(as.numeric(p_interaction))) %>%
  group_by(outcome, pollutant, exposure_term, exposure_scale, subgroup, test) %>%
  summarise(
    n_sites = n_distinct(site),
    sites = paste(sort(unique(site)), collapse = "; "),
    fisher_statistic = -2 * sum(log(p_interaction), na.rm = TRUE),
    fisher_df = 2 * sum(!is.na(p_interaction)),
    p_fisher = stats::pchisq(fisher_statistic, df = fisher_df, lower.tail = FALSE),
    .groups = "drop"
  )
readr::write_csv(
  read_sites("primary_subgroup_interaction_tests.csv") %>%
    filter(!str_detect(outcome, regex("imv|invasive|ventilation duration", ignore_case = TRUE))),
  file.path(out_dir, "site_subgroup_interaction_tests.csv")
)
readr::write_csv(subgroup_interactions, file.path(out_dir, "pooled_subgroup_interaction_tests_fisher.csv"))

pool_prediction_curve <- function(tbl, id_cols, grid_n = 100) {
  if (!nrow(tbl)) return(tibble())
  approx_safe <- function(x, y, xout) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 2) return(rep(NA_real_, length(xout)))
    approx(x[ok], y[ok], xout = xout, rule = 2, ties = "ordered")$y
  }
  derive_link <- function(endpoint, estimate, conf_low, conf_high) {
    endpoint <- as.character(endpoint)
    out_est <- rep(NA_real_, length(estimate))
    out_se <- rep(NA_real_, length(estimate))

    mort_ok <- endpoint == "mortality" &
      is.finite(estimate) & estimate > 0 & estimate < 1 &
      is.finite(conf_low) & conf_low > 0 & conf_low < 1 &
      is.finite(conf_high) & conf_high > 0 & conf_high < 1 &
      conf_high > conf_low
    out_est[mort_ok] <- stats::qlogis(estimate[mort_ok])
    out_se[mort_ok] <- (stats::qlogis(conf_high[mort_ok]) - stats::qlogis(conf_low[mort_ok])) / (2 * 1.959963984540054)

    vfd_ok <- endpoint == "vfd" &
      is.finite(estimate) & estimate > 0 &
      is.finite(conf_low) & conf_low > 0 &
      is.finite(conf_high) & conf_high > 0 &
      conf_high > conf_low
    out_est[vfd_ok] <- log(estimate[vfd_ok])
    out_se[vfd_ok] <- (log(conf_high[vfd_ok]) - log(conf_low[vfd_ok])) / (2 * 1.959963984540054)

    tibble(link_estimate_derived = out_est, link_se_derived = out_se)
  }
  if (!"link_estimate" %in% names(tbl)) tbl$link_estimate <- NA_real_
  if (!"link_se" %in% names(tbl)) tbl$link_se <- NA_real_
  tbl <- tbl %>%
    mutate(
      exposure_raw = as.numeric(exposure_raw),
      estimate = as.numeric(estimate),
      conf_low = as.numeric(conf_low),
      conf_high = as.numeric(conf_high),
      link_estimate = suppressWarnings(as.numeric(link_estimate)),
      link_se = suppressWarnings(as.numeric(link_se)),
      se = pmax((conf_high - conf_low) / (2 * 1.959963984540054), 1e-8)
    ) %>%
    bind_cols(derive_link(.$endpoint, .$estimate, .$conf_low, .$conf_high)) %>%
    mutate(
      link_estimate = coalesce(link_estimate, link_estimate_derived),
      link_se = coalesce(link_se, link_se_derived)
    ) %>%
    select(-link_estimate_derived, -link_se_derived) %>%
    filter(is.finite(exposure_raw), is.finite(estimate), is.finite(se), se > 0)

  tbl %>%
    group_by(across(all_of(id_cols))) %>%
    group_modify(function(.x, .y) {
      ranges <- .x %>%
        group_by(site) %>%
        summarise(x_min = min(exposure_raw), x_max = max(exposure_raw), .groups = "drop")
      x_min <- min(ranges$x_min)
      x_max <- max(ranges$x_max)
      if (!is.finite(x_min) || !is.finite(x_max) || x_max <= x_min || nrow(ranges) < 1) return(tibble())
      grid <- seq(x_min, x_max, length.out = grid_n)
      interp <- .x %>%
        arrange(site, exposure_raw) %>%
        group_by(site) %>%
        group_modify(function(site_df, site_key) {
          tibble(
            exposure_raw = grid,
            estimate = approx_safe(site_df$exposure_raw, site_df$estimate, grid),
            conf_low = approx_safe(site_df$exposure_raw, site_df$conf_low, grid),
            conf_high = approx_safe(site_df$exposure_raw, site_df$conf_high, grid),
            se = approx_safe(site_df$exposure_raw, site_df$se, grid),
            link_estimate = approx_safe(site_df$exposure_raw, site_df$link_estimate, grid),
            link_se = approx_safe(site_df$exposure_raw, site_df$link_se, grid),
            n = max(site_df$n, na.rm = TRUE),
            events = max(site_df$events, na.rm = TRUE)
          )
        }) %>%
        ungroup()
      endpoint_value <- as.character(.y$endpoint[[1]] %||% NA_character_)
      link_site_summary <- interp %>%
        group_by(site) %>%
        summarise(has_link = sum(is.finite(link_estimate) & is.finite(link_se) & link_se > 0) >= 2, .groups = "drop")
      has_link <- nrow(link_site_summary) > 0 && all(link_site_summary$has_link)
      if (has_link && endpoint_value %in% c("mortality", "vfd")) {
        interp %>%
          filter(is.finite(link_estimate), is.finite(link_se), link_se > 0) %>%
          mutate(weight = 1 / link_se^2) %>%
          group_by(exposure_raw) %>%
          summarise(
            n_sites = n_distinct(site),
            n_total = sum(n, na.rm = TRUE),
            events_total = sum(events, na.rm = TRUE),
            link_estimate = sum(weight * link_estimate) / sum(weight),
            link_se = sqrt(1 / sum(weight)),
            .groups = "drop"
          ) %>%
          mutate(
            estimate = if (identical(endpoint_value, "mortality")) plogis(link_estimate) else exp(link_estimate),
            conf_low = if (identical(endpoint_value, "mortality")) plogis(link_estimate - 1.959963984540054 * link_se) else exp(link_estimate - 1.959963984540054 * link_se),
            conf_high = if (identical(endpoint_value, "mortality")) plogis(link_estimate + 1.959963984540054 * link_se) else exp(link_estimate + 1.959963984540054 * link_se),
            se = NA_real_
          )
      } else {
        interp %>%
          filter(is.finite(estimate), is.finite(se), se > 0) %>%
          mutate(weight = 1 / pmax(se, 1e-8)^2) %>%
          group_by(exposure_raw) %>%
          summarise(
            n_sites = n_distinct(site),
            n_total = sum(n, na.rm = TRUE),
            events_total = sum(events, na.rm = TRUE),
            estimate = sum(weight * estimate) / sum(weight),
            se = sqrt(1 / sum(weight)),
            conf_low = estimate - 1.959963984540054 * se,
            conf_high = estimate + 1.959963984540054 * se,
            link_estimate = NA_real_,
            link_se = NA_real_,
            .groups = "drop"
          )
      }
    }) %>%
    ungroup()
}

primary_curve_sites <- bind_rows(
  read_sites("primary_exposure_response_predictions.csv") %>% mutate(prediction_family = "Whole cohort and ARF subtype"),
  read_sites("primary_exposure_response_by_sex_race_predictions.csv") %>% mutate(prediction_family = "Sex and race/ethnicity")
) %>%
  filter(endpoint %in% c("mortality", "vfd")) %>%
  mutate(
    group_var = if ("group_var" %in% names(.)) coalesce(group_var, "arf_subtype") else "arf_subtype",
    group_level = if ("group_level" %in% names(.)) coalesce(as.character(group_level), as.character(arf_subtype)) else as.character(arf_subtype),
    curve_label = case_when(
      !is.na(group_var) & !is.na(group_level) ~ paste(group_var, group_level, sep = ": "),
      TRUE ~ as.character(curve_type)
    ),
    pollutant = clean_pollutant(pollutant),
    y_label = case_when(
      endpoint == "mortality" ~ "Predicted mortality probability by day 28",
      endpoint == "vfd" ~ "Predicted ventilator-free days",
      TRUE ~ y_label
    )
  ) %>%
  filter(!(group_var == "sex" & !group_level %in% c("Female", "Male")))

primary_curve_ids <- c("prediction_family", "curve_type", "outcome", "endpoint", "pollutant", "pollutant_label", "group_var", "group_level", "curve_label", "y_label")
primary_curve_pooled <- pool_prediction_curve(primary_curve_sites, primary_curve_ids)
pollutant_label_levels_primary <- c(
  "Nitrogen~dioxide~(NO[2])",
  "Fine~particulate~matter~(PM[2.5])"
)
primary_curve_pooled <- primary_curve_pooled %>%
  mutate(
    pollutant = factor(pollutant, levels = c("NO2", "PM2.5")),
    pollutant_label = factor(
      case_when(
        as.character(pollutant) == "NO2" ~ pollutant_label_levels_primary[[1]],
        as.character(pollutant) == "PM2.5" ~ pollutant_label_levels_primary[[2]],
        TRUE ~ as.character(pollutant_label)
      ),
      levels = pollutant_label_levels_primary
    )
  )
readr::write_csv(primary_curve_sites, file.path(out_dir, "site_primary_exposure_response_predictions.csv"))
readr::write_csv(primary_curve_pooled, file.path(out_dir, "pooled_primary_exposure_response_predictions.csv"))

build_primary_exposure_rug <- function() {
  empty_rug <- tibble(
    site = character(),
    group_var = character(),
    group_level = character(),
    arf_subtype = character(),
    pollutant = factor(character(), levels = c("NO2", "PM2.5")),
    pollutant_label = factor(character(), levels = pollutant_label_levels_primary),
    exposure_raw = numeric(),
    n = numeric(),
    endpoint = character()
  )

  binned <- read_sites("primary_exposure_distribution_bins.csv")
  if (nrow(binned) && all(c("group_var", "group_level", "pollutant", "bin_mid", "n") %in% names(binned))) {
    return(
      binned %>%
        transmute(
          site = as.character(site),
          group_var = as.character(group_var),
          group_level = as.character(group_level),
          arf_subtype = if_else(group_var == "arf_subtype", group_level, NA_character_),
          pollutant = factor(clean_pollutant(pollutant), levels = c("NO2", "PM2.5")),
          pollutant_label = factor(
            case_when(
              clean_pollutant(pollutant) == "NO2" ~ pollutant_label_levels_primary[[1]],
              clean_pollutant(pollutant) == "PM2.5" ~ pollutant_label_levels_primary[[2]]
            ),
            levels = pollutant_label_levels_primary
          ),
          exposure_raw = as.numeric(bin_mid),
          n = as.numeric(n)
        ) %>%
        filter(is.finite(exposure_raw), is.finite(n), n > 0) %>%
        tidyr::crossing(endpoint = c("mortality", "vfd"))
    )
  }

  rug_data <- site_dirs %>%
    mutate(data = map2(site_dir, site, function(site_dir, site_expected) {
      path <- file.path(site_dir, "resubmission_analysis_dataset.csv")
      if (!file.exists(path)) return(tibble())
      readr::read_csv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        col_select = any_of(c("site", "arf_subtype", "no2_12m_zcta", "pm25_12m_zcta"))
      ) %>%
        mutate(site = if ("site" %in% names(.)) coalesce(as.character(site), site_expected) else site_expected)
    })) %>%
    pull(data) %>%
    bind_rows()

  if (!nrow(rug_data) || !all(c("site", "arf_subtype", "no2_12m_zcta", "pm25_12m_zcta") %in% names(rug_data))) {
    return(empty_rug)
  }

  rug_data %>%
    transmute(
      site,
      group_var = "arf_subtype",
      group_level = as.character(arf_subtype),
      arf_subtype = as.character(arf_subtype),
      NO2 = as.numeric(no2_12m_zcta),
      PM2.5 = as.numeric(pm25_12m_zcta)
    ) %>%
    pivot_longer(c(NO2, PM2.5), names_to = "pollutant", values_to = "exposure_raw") %>%
    filter(is.finite(exposure_raw)) %>%
    mutate(
      pollutant = factor(pollutant, levels = c("NO2", "PM2.5")),
      pollutant_label = factor(
        case_when(
          as.character(pollutant) == "NO2" ~ pollutant_label_levels_primary[[1]],
          as.character(pollutant) == "PM2.5" ~ pollutant_label_levels_primary[[2]]
        ),
        levels = pollutant_label_levels_primary
      ),
      n = 1
    ) %>%
    tidyr::crossing(endpoint = c("mortality", "vfd"))
}

primary_rug <- build_primary_exposure_rug()

exposure_plot_lower_quantile <- as.numeric(Sys.getenv("REFER_EXPOSURE_RESPONSE_X_Q_LOW", "0.025"))
exposure_plot_upper_quantile <- as.numeric(Sys.getenv("REFER_EXPOSURE_RESPONSE_X_Q_HIGH", "0.975"))

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) < 2 || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  o <- order(x)
  x <- x[o]
  w <- w[o]
  cw <- cumsum(w) / sum(w)
  stats::approx(c(0, cw), c(x[[1]], x), xout = probs, rule = 2, ties = "ordered")$y
}

exposure_plot_windows <- primary_rug %>%
  filter(group_var == "overall", is.finite(exposure_raw), is.finite(n), n > 0) %>%
  group_by(pollutant) %>%
  summarise(
    x_low = weighted_quantile(exposure_raw, n, exposure_plot_lower_quantile)[[1]],
    x_high = weighted_quantile(exposure_raw, n, exposure_plot_upper_quantile)[[1]],
    n_total = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(x_low), is.finite(x_high), x_high > x_low)

readr::write_csv(exposure_plot_windows, file.path(out_dir, "pooled_exposure_response_plot_windows.csv"))

pool_aj_curve <- function(tbl, id_cols, grid = seq(0, 28, by = 0.25)) {
  if (!nrow(tbl)) return(tibble())
  tbl <- tbl %>%
    mutate(
      time = as.numeric(time),
      cif = as.numeric(cif),
      var = suppressWarnings(as.numeric(var)),
      se = sqrt(pmax(var, 0))
    ) %>%
    filter(is.finite(time), time <= 28, is.finite(cif))

  tbl %>%
    group_by(across(all_of(id_cols))) %>%
    group_modify(function(.x, .y) {
      interp <- .x %>%
        arrange(site, time) %>%
        group_by(site) %>%
        group_modify(function(site_df, site_key) {
          tibble(
            time = grid,
            cif = approx(site_df$time, site_df$cif, xout = grid, rule = 2, ties = "ordered")$y,
            se = approx(site_df$time, site_df$se, xout = grid, rule = 2, ties = "ordered")$y
          )
        }) %>%
        ungroup()
      interp %>%
        mutate(weight = if_else(is.finite(se) & se > 0, 1 / se^2, NA_real_)) %>%
        group_by(time) %>%
        summarise(
          n_sites = n_distinct(site),
          estimate = if (all(is.na(weight))) mean(cif, na.rm = TRUE) else sum(weight * cif, na.rm = TRUE) / sum(weight, na.rm = TRUE),
          se = if (all(is.na(weight))) NA_real_ else sqrt(1 / sum(weight, na.rm = TRUE)),
          .groups = "drop"
        ) %>%
        mutate(conf_low = pmax(0, estimate - 1.959963984540054 * se), conf_high = pmin(1, estimate + 1.959963984540054 * se))
    }) %>%
    ungroup()
}

aj_sites <- read_sites("unadjusted_aalen_johansen_cif_by_exposure_quartile.csv") %>%
  mutate(
    site = if_else(is.na(site), site_label, site),
    pollutant = clean_pollutant(pollutant),
    outcome = recode(outcome, "Successful extubation" = "Extubation")
  )
aj_pooled <- pool_aj_curve(aj_sites, c("pollutant", "exposure_quartile", "event_code", "outcome"))
readr::write_csv(aj_sites, file.path(out_dir, "site_unadjusted_aalen_johansen_cif_by_exposure_quartile.csv"))
readr::write_csv(aj_pooled, file.path(out_dir, "pooled_unadjusted_aalen_johansen_cif_by_exposure_quartile.csv"))

pool_weighted_cif <- function(tbl, id_cols, grid = seq(0, 28, by = 0.25)) {
  if (!nrow(tbl)) return(tibble())
  tbl <- tbl %>%
    mutate(time = as.numeric(time), cif = as.numeric(cif), pollutant = clean_pollutant(pollutant)) %>%
    filter(is.finite(time), time <= 28, is.finite(cif))
  weights <- fine_gray_effects %>%
    distinct(site, outcome, model, n) %>%
    group_by(site, outcome, model) %>%
    summarise(n = max(n, na.rm = TRUE), .groups = "drop")
  tbl <- tbl %>% left_join(weights, by = c("site", "outcome", "model")) %>% mutate(n = if_else(is.finite(n), n, 1))

  tbl %>%
    group_by(across(all_of(id_cols))) %>%
    group_modify(function(.x, .y) {
      interp <- .x %>%
        arrange(site, time) %>%
        group_by(site) %>%
        group_modify(function(site_df, site_key) {
          tibble(
            time = grid,
            cif = approx(site_df$time, site_df$cif, xout = grid, rule = 2, ties = "ordered")$y,
            n = max(site_df$n, na.rm = TRUE)
          )
        }) %>%
        ungroup()
      interp %>%
        group_by(time) %>%
        summarise(n_sites = n_distinct(site), n_total = sum(n, na.rm = TRUE), cif = weighted.mean(cif, w = n, na.rm = TRUE), .groups = "drop")
    }) %>%
    ungroup()
}

fg_curve_sites <- read_sites("fine_gray_adjusted_cif_predictions_by_exposure_quartile.csv") %>%
  mutate(pollutant = clean_pollutant(pollutant))
fg_curve_pooled <- pool_weighted_cif(fg_curve_sites, c("pollutant", "exposure_quartile", "outcome", "event_code", "model"))
readr::write_csv(fg_curve_sites, file.path(out_dir, "site_fine_gray_adjusted_cif_predictions_by_exposure_quartile.csv"))
readr::write_csv(fg_curve_pooled, file.path(out_dir, "pooled_fine_gray_adjusted_cif_predictions_by_exposure_quartile.csv"))

format_effect <- function(x) {
  paste0(number(x$estimate, accuracy = 0.01), " (", number(x$conf_low, accuracy = 0.01), ", ", number(x$conf_high, accuracy = 0.01), ")")
}

write_plot <- function(plot, filename, width = 12, height = 8) {
  ggsave(file.path(fig_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 320)
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot, width = width, height = height)
}

theme_refer <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "white", colour = "grey55", linewidth = 0.4),
      strip.text = element_text(face = "bold"),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      legend.position = "bottom",
      plot.title = element_blank()
    )
}

forest_primary <- bind_rows(
  primary_effects %>% mutate(row_type = "Site-specific"),
  primary_pooled %>% mutate(site = "Pooled", row_type = "Pooled")
) %>%
  filter(model %in% c("PM25 single-pollutant", "NO2 single-pollutant", "PM25 + NO2")) %>%
  filter(!str_detect(outcome, regex("imv|invasive|ventilation duration", ignore_case = TRUE))) %>%
  mutate(
    site = fct_relevel(factor(site), "Pooled", after = Inf),
    label = paste(site, term, sep = " | "),
    outcome = factor(outcome)
  )

p_forest_primary <- ggplot(forest_primary, aes(x = estimate, y = fct_rev(label), xmin = conf_low, xmax = conf_high, colour = site)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(height = 0.16, linewidth = 0.6) +
  geom_point(size = 2.2) +
  scale_x_log10() +
  facet_grid(outcome ~ model, scales = "free_y", space = "free_y") +
  labs(x = "Effect estimate with 95% CI", y = NULL, colour = NULL) +
  theme_refer(10)
write_plot(p_forest_primary, "pooled_primary_effect_forest", width = 15, height = 9)

curve_colours <- c("Mortality by day 28" = "#C82536", "Ventilator-free days" = "#2C77BF")
subtype_colours <- c("Hypoxemic" = "#0072B2", "Hypercapnic" = "#009E73", "Mixed" = "#D55E00")
sex_colours <- c("Female" = "#009E73", "Male" = "#D55E00")
race_levels <- primary_curve_pooled %>%
  filter(prediction_family == "Sex and race/ethnicity", group_var == "race_ethnicity") %>%
  pull(group_level) %>%
  unique() %>%
  na.omit() %>%
  as.character()
race_levels <- intersect(
  c("Asian", "Other/Unknown", "Hispanic White", "Non-Hispanic Black", "Non-Hispanic White", race_levels),
  race_levels
)
race_colours <- setNames(
  c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")[seq_along(race_levels)],
  race_levels
)

pollutant_x_labels <- list(
  "NO2" = expression(NO[2]~"(ppb)"),
  "PM2.5" = expression(PM[2.5]~"("*mu*"g/m"^3*")")
)

pollutant_title <- function(pollutant_value) {
  if (as.character(pollutant_value) == "NO2") {
    expression("Nitrogen dioxide (NO"[2]*")")
  } else {
    expression("Fine particulate matter (PM"[2.5]*")")
  }
}

theme_response <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(face = "bold", size = 8.5, margin = margin(t = 8, r = 8, b = 8, l = 8)),
      strip.placement = "outside",
      legend.position = "bottom",
      plot.margin = margin(5.5, 5.5, 5.5, 14)
    )
}

label_outcome_axis <- function(x) {
  ifelse(
    abs(x) < 1e-8,
    "0",
    if_else(
      abs(x) >= 1,
      label_number(accuracy = 1, trim = TRUE)(x),
      label_number(accuracy = 0.01, trim = TRUE)(x)
    )
  )
}

apply_outcome_y_windows <- function(df) {
  df %>%
    mutate(
      y_window_low = case_when(endpoint == "mortality" ~ 0, TRUE ~ 15),
      y_window_high = case_when(endpoint == "mortality" ~ 0.40, TRUE ~ 28),
      estimate_plot = pmin(pmax(estimate, y_window_low), y_window_high),
      conf_low_plot = pmin(pmax(conf_low, y_window_low), y_window_high),
      conf_high_plot = pmin(pmax(conf_high, y_window_low), y_window_high),
      y_label = factor(
        y_label,
        levels = c("Predicted mortality probability by day 28", "Predicted ventilator-free days")
      )
    )
}

make_y_window_df <- function(plot_df) {
  plot_df %>%
    distinct(y_label, endpoint) %>%
    mutate(
      exposure_raw = min(plot_df$exposure_raw, na.rm = TRUE),
      y_window_low = case_when(endpoint == "mortality" ~ 0, TRUE ~ 15),
      y_window_high = case_when(endpoint == "mortality" ~ 0.40, TRUE ~ 28)
    ) %>%
    pivot_longer(c(y_window_low, y_window_high), names_to = "window_bound", values_to = "y_value")
}

wrap_pollutant_panels <- function(plot_fun) {
  plots <- lapply(c("NO2", "PM2.5"), plot_fun)
  patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom")
}

make_separate_rug_panel <- function(rugs, group_col, palette_values, pollutant_value) {
  group_col <- rlang::ensym(group_col)
  rugs <- rugs %>%
    mutate(rug_level = factor(as.character(!!group_col), levels = names(palette_values))) %>%
    filter(!is.na(rug_level))
  if (!"n" %in% names(rugs)) {
    rugs <- rugs %>% mutate(n = 1)
  }
  rugs <- rugs %>%
    group_by(rug_level, exposure_raw) %>%
    summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      density_alpha = if (n_distinct(n) > 1) {
        scales::rescale(log1p(n), to = c(0.08, 1.00))
      } else {
        0.70
      }
    )

  ggplot(rugs, aes(x = exposure_raw, y = rug_level, color = rug_level, alpha = density_alpha)) +
    geom_point(shape = "|", size = 3.3, show.legend = FALSE) +
    scale_color_manual(values = palette_values, drop = TRUE) +
    scale_alpha_identity(guide = "none") +
    scale_x_continuous(limits = exposure_plot_windows %>%
                         filter(as.character(pollutant) == pollutant_value) %>%
                         select(x_low, x_high) %>%
                         unlist(use.names = FALSE),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = pollutant_x_labels[[as.character(pollutant_value)]], y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10, color = "grey20"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.25),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 13),
      plot.margin = margin(0, 5.5, 5.5, 5.5),
      legend.position = "none"
    )
}

make_combined_curve_plot <- function(df, rugs, palette_values, filename_stub, plot_width = 12.5, plot_height = 11,
                                     ribbon_alpha = 0.14) {
  df <- df %>%
    filter(as.character(group_level) %in% names(palette_values)) %>%
    mutate(group_level = factor(as.character(group_level), levels = names(palette_values)))
  rugs <- rugs %>%
    filter(as.character(group_level) %in% names(palette_values)) %>%
    mutate(group_level = factor(as.character(group_level), levels = names(palette_values)))

  p <- wrap_pollutant_panels(function(pollutant_value) {
    x_window <- exposure_plot_windows %>%
      filter(as.character(pollutant) == pollutant_value) %>%
      select(x_low, x_high) %>%
      unlist(use.names = FALSE)
    if (length(x_window) != 2 || any(!is.finite(x_window)) || x_window[[2]] <= x_window[[1]]) {
      plot_df_window <- df %>% filter(as.character(pollutant) == pollutant_value)
      x_window <- range(plot_df_window$exposure_raw, na.rm = TRUE)
    }
    plot_df <- df %>%
      filter(
        as.character(pollutant) == pollutant_value,
        exposure_raw >= x_window[[1]],
        exposure_raw <= x_window[[2]]
      )
    plot_rugs <- rugs %>%
      filter(
        as.character(pollutant) == pollutant_value,
        exposure_raw >= x_window[[1]],
        exposure_raw <= x_window[[2]]
      )

    curve_plot <- ggplot(plot_df, aes(exposure_raw, estimate_plot, color = group_level, fill = group_level)) +
      geom_blank(data = make_y_window_df(plot_df), aes(x = exposure_raw, y = y_value), inherit.aes = FALSE) +
      geom_ribbon(aes(ymin = conf_low_plot, ymax = conf_high_plot), alpha = ribbon_alpha, color = NA, show.legend = FALSE) +
      geom_line(linewidth = 1.05) +
      facet_grid(y_label ~ ., scales = "free_y", switch = "y") +
      scale_color_manual(values = palette_values, drop = TRUE) +
      scale_fill_manual(values = palette_values, drop = TRUE) +
      scale_x_continuous(limits = x_window, expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(labels = label_outcome_axis, expand = expansion(mult = c(0, 0))) +
      labs(title = pollutant_title(pollutant_value), x = NULL, y = NULL, color = NULL, fill = NULL) +
      theme_response()

    curve_plot / make_separate_rug_panel(plot_rugs, group_level, palette_values, pollutant_value) +
      plot_layout(heights = c(1, 0.18))
  })

  write_plot(p, filename_stub, width = plot_width, height = plot_height)
  p
}

primary_rug_by_subtype <- primary_rug %>%
  filter(group_var == "arf_subtype") %>%
  filter(!is.na(arf_subtype)) %>%
  mutate(group_level = factor(as.character(arf_subtype), levels = names(subtype_colours)))

build_primary_exposure_rug_by_sex_race <- function() {
  empty_rug <- tibble(
    site = character(),
    sex = character(),
    race_ethnicity = character(),
    pollutant = factor(character(), levels = c("NO2", "PM2.5")),
    exposure_raw = numeric(),
    n = numeric()
  )

  binned <- primary_rug %>%
    filter(group_var %in% c("sex", "race_ethnicity")) %>%
    transmute(
      site,
      sex = if_else(group_var == "sex", group_level, NA_character_),
      race_ethnicity = if_else(group_var == "race_ethnicity", group_level, NA_character_),
      pollutant,
      exposure_raw,
      n
    )
  if (nrow(binned)) {
    return(binned)
  }

  rug_data <- site_dirs %>%
    mutate(data = map2(site_dir, site, function(site_dir, site_expected) {
    path <- file.path(site_dir, "resubmission_analysis_dataset.csv")
    if (!file.exists(path)) return(tibble())
    readr::read_csv(
      path,
      show_col_types = FALSE,
      progress = FALSE,
      col_select = any_of(c("site", "sex", "race_ethnicity", "no2_12m_zcta", "pm25_12m_zcta"))
    ) %>%
      mutate(site = if ("site" %in% names(.)) coalesce(as.character(site), site_expected) else site_expected)
    })) %>%
    pull(data) %>%
    bind_rows()

  if (!nrow(rug_data) || !all(c("site", "sex", "race_ethnicity", "no2_12m_zcta", "pm25_12m_zcta") %in% names(rug_data))) {
    return(empty_rug)
  }

  rug_data %>%
    mutate(
      race_ethnicity = forcats::fct_collapse(
        factor(race_ethnicity),
        "Other/Unknown" = c("Hispanic Other/Unknown", "Other/Unknown")
      )
    ) %>%
    transmute(
      site,
      sex = as.character(sex),
      race_ethnicity = as.character(race_ethnicity),
      NO2 = as.numeric(no2_12m_zcta),
      PM2.5 = as.numeric(pm25_12m_zcta)
    ) %>%
    pivot_longer(c(NO2, PM2.5), names_to = "pollutant", values_to = "exposure_raw") %>%
    filter(is.finite(exposure_raw)) %>%
    mutate(
      pollutant = factor(pollutant, levels = c("NO2", "PM2.5")),
      n = 1
    )
}

primary_rug_by_sex <- build_primary_exposure_rug_by_sex_race()

p_curve_whole <- primary_curve_pooled %>%
  filter(prediction_family == "Whole cohort and ARF subtype", curve_type == "Whole cohort") %>%
  apply_outcome_y_windows() %>%
  mutate(group_level = factor(as.character(outcome), levels = names(curve_colours)))

primary_rug_whole <- primary_rug %>%
  filter(group_var == "overall") %>%
  mutate(group_level = factor("Mortality by day 28", levels = names(curve_colours)))

p_curve_whole <- make_combined_curve_plot(
  p_curve_whole,
  primary_rug_whole,
  curve_colours,
  "pooled_primary_exposure_response_whole_cohort",
  plot_width = 12.5,
  plot_height = 11
)

p_curve_subtype <- primary_curve_pooled %>%
  filter(prediction_family == "Whole cohort and ARF subtype", curve_type != "Whole cohort") %>%
  apply_outcome_y_windows() %>%
  mutate(group_level = factor(group_level, levels = names(subtype_colours)))

p_curve_subtype <- make_combined_curve_plot(
  p_curve_subtype,
  primary_rug_by_subtype,
  subtype_colours,
  "pooled_primary_exposure_response_by_arf_subtype",
  plot_width = 12.5,
  plot_height = 11
)
write_plot(p_curve_subtype, "figure1_primary_outcomes_with_arf_subtype", width = 12.5, height = 11)

p_curve_sex <- primary_curve_pooled %>%
  filter(prediction_family == "Sex and race/ethnicity", group_var == "sex") %>%
  apply_outcome_y_windows() %>%
  mutate(group_level = factor(group_level, levels = names(sex_colours)))

p_curve_sex <- make_combined_curve_plot(
  p_curve_sex,
  primary_rug_by_sex %>% transmute(pollutant, exposure_raw, group_level = factor(sex, levels = names(sex_colours))),
  sex_colours,
  "pooled_primary_exposure_response_by_sex",
  plot_width = 12.5,
  plot_height = 11
)

p_curve_race <- primary_curve_pooled %>%
  filter(prediction_family == "Sex and race/ethnicity", group_var == "race_ethnicity") %>%
  apply_outcome_y_windows() %>%
  mutate(group_level = factor(group_level, levels = names(race_colours)))

p_curve_race <- make_combined_curve_plot(
  p_curve_race,
  primary_rug_by_sex %>% transmute(pollutant, exposure_raw, group_level = factor(race_ethnicity, levels = names(race_colours))),
  race_colours,
  "pooled_primary_exposure_response_by_race_ethnicity",
  plot_width = 16,
  plot_height = 11,
  ribbon_alpha = 0.06
)

quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
quartile_cols <- c(
  "Q1 lowest" = "#2166AC",
  "Q2" = "#67A9CF",
  "Q3" = "#F4A582",
  "Q4 highest" = "#B2182B"
)
tick_times <- c(0, 7, 14, 21, 28)
x_limits_cif <- c(-2.5, 31.0)
x_limits_cif_curves_only <- c(-0.35, 28.35)
outcome_levels_cif <- c("Extubation", "Death", "Persistent RF")
pollutant_levels_cif <- c("NO2", "PM2.5")

label_quartile <- function(x) {
  recode(as.character(x), "Q1" = "Q1 lowest", "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4 highest")
}

theme_cif_transplant <- function() {
  theme_minimal(base_size = 18) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      legend.key.width = unit(30, "pt"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16, color = "grey20"),
      panel.spacing.x = unit(28, "pt"),
      panel.spacing.y = unit(22, "pt"),
      plot.margin = margin(4, 8, 4, 8)
    )
}

theme_table_box <- function() {
  theme_minimal(base_size = 16) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.title = element_text(size = 17, color = "grey10", margin = margin(t = 7)),
      axis.text.x = element_text(size = 16, color = "grey20", margin = margin(t = 2)),
      axis.text.y = element_text(size = 17, face = "bold"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.25),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 17, face = "bold", color = "grey20", margin = margin(b = 3)),
      plot.margin = margin(3, 8, 4, 8)
    )
}

build_events_at_risk_from_analysis <- function() {
  analysis_sites <- site_dirs %>%
    mutate(data = map2(site_dir, site, function(site_dir, site_expected) {
      path <- file.path(site_dir, "resubmission_analysis_dataset.csv")
      if (!file.exists(path)) return(tibble())
      readr::read_csv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        col_select = any_of(c("site", "has_imv_after_arf", "ftime_days", "event_code", "no2_q", "pm25_q"))
      ) %>%
        mutate(site = if ("site" %in% names(.)) coalesce(as.character(site), site_expected) else site_expected)
    })) %>%
    select(site_expected = site, site_dir, data) %>%
    tidyr::unnest(data)
  if (!nrow(analysis_sites)) return(tibble())
  if ("has_imv_after_arf" %in% names(analysis_sites)) {
    analysis_sites <- analysis_sites %>% filter(has_imv_after_arf %in% c(TRUE, 1, "TRUE", "true"))
  }
  needed <- c("ftime_days", "event_code", "no2_q", "pm25_q")
  if (!all(needed %in% names(analysis_sites))) return(tibble())

  make_one <- function(data, exposure_var, pollutant_label) {
    base <- data %>%
      filter(!is.na(.data[[exposure_var]])) %>%
      mutate(
        exposure_quartile = factor(label_quartile(.data[[exposure_var]]), levels = quartile_labels),
        ftime_days = as.numeric(ftime_days),
        event_code = as.integer(event_code)
      ) %>%
      filter(!is.na(exposure_quartile), is.finite(ftime_days))

    tidyr::expand_grid(
      pollutant = pollutant_label,
      outcome = outcome_levels_cif,
      exposure_quartile = factor(quartile_labels, levels = quartile_labels),
      time = tick_times
    ) %>%
      mutate(event_code = match(outcome, outcome_levels_cif)) %>%
      rowwise() %>%
      mutate(
        n_risk = sum(base$exposure_quartile == exposure_quartile & base$ftime_days >= time, na.rm = TRUE),
        n_events = sum(base$exposure_quartile == exposure_quartile & base$event_code == event_code & base$ftime_days <= time, na.rm = TRUE),
        label = paste0(n_events, "/", n_risk)
      ) %>%
      ungroup()
  }

  bind_rows(
    make_one(analysis_sites, "no2_q", "NO2"),
    make_one(analysis_sites, "pm25_q", "PM2.5")
  )
}

risk_event_table <- read_sites("unadjusted_aalen_johansen_events_at_risk_by_tick.csv")
if (nrow(risk_event_table)) {
  risk_event_table <- risk_event_table %>%
    mutate(
      time = as.numeric(time),
      n_events = as.numeric(n_events),
      n_risk = as.numeric(n_risk),
      pollutant = clean_pollutant(pollutant),
      outcome = recode(outcome, "Successful extubation" = "Extubation", "Persistent respiratory failure" = "Persistent RF"),
      exposure_quartile = factor(label_quartile(exposure_quartile), levels = quartile_labels)
    ) %>%
    group_by(pollutant, outcome, exposure_quartile, time) %>%
    summarise(n_events = sum(n_events, na.rm = TRUE), n_risk = sum(n_risk, na.rm = TRUE), .groups = "drop") %>%
    mutate(label = paste0(n_events, "/", n_risk))
} else {
  risk_event_table <- build_events_at_risk_from_analysis()
}
readr::write_csv(risk_event_table, file.path(out_dir, "pooled_unadjusted_aalen_johansen_events_at_risk_by_tick.csv"))

aj_plot_data <- aj_pooled %>%
  mutate(
    pollutant = factor(clean_pollutant(pollutant), levels = pollutant_levels_cif),
    outcome = recode(outcome, "Successful extubation" = "Extubation", "Persistent respiratory failure" = "Persistent RF"),
    outcome = factor(outcome, levels = outcome_levels_cif),
    exposure_quartile = factor(label_quartile(exposure_quartile), levels = quartile_labels)
  ) %>%
  filter(!is.na(pollutant), !is.na(outcome), !is.na(exposure_quartile), time <= 28)

aj_panel_scales <- aj_plot_data %>%
  group_by(pollutant, outcome) %>%
  summarise(y_max = max(conf_high, estimate, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    y_limit = case_when(
      outcome == "Persistent RF" ~ pmax(ceiling(y_max / 0.005) * 0.005, 0.025),
      outcome == "Death" ~ pmax(ceiling(y_max / 0.05) * 0.05, 0.25),
      TRUE ~ pmax(ceiling(y_max / 0.10) * 0.10, 0.70)
    )
  )

make_column_header <- function(label_expr) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "grey45", linewidth = 0.55) +
    annotate("text", x = 0.5, y = 0.42, label = label_expr, parse = TRUE, hjust = 0.5, vjust = 0.5, fontface = "bold", size = 6.4) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(12, 8, -18, 8))
}

make_row_header <- function(label_text) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "grey45", linewidth = 0.55) +
    annotate("text", x = 0.5, y = 0.5, label = label_text, angle = -90, hjust = 0.5, vjust = 0.5, fontface = "bold", size = 6.8) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(4, 8, 2, -2))
}

make_aj_curve_panel <- function(pollutant_value, outcome_value, show_y_axis = TRUE, show_x_axis = FALSE, x_limits_use = x_limits_cif) {
  dat <- aj_plot_data %>% filter(pollutant == pollutant_value, outcome == outcome_value)
  y_limit <- aj_panel_scales %>% filter(pollutant == pollutant_value, outcome == outcome_value) %>% pull(y_limit)
  y_breaks <- if (identical(as.character(outcome_value), "Persistent RF")) {
    seq(0, y_limit, by = 0.005)
  } else if (identical(as.character(outcome_value), "Death")) {
    seq(0, y_limit, by = 0.05)
  } else {
    seq(0, y_limit, by = 0.20)
  }
  y_labels <- if (identical(as.character(outcome_value), "Persistent RF")) {
    label_percent(accuracy = 0.5)
  } else {
    label_percent(accuracy = 1)
  }

  ggplot(dat, aes(x = time, y = estimate, color = exposure_quartile)) +
    geom_step(linewidth = 1.05) +
    scale_color_manual(values = quartile_cols, drop = FALSE) +
    scale_x_continuous(breaks = tick_times, labels = label_number(accuracy = 1, trim = TRUE), limits = x_limits_use, expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = c(0, y_limit), expand = expansion(mult = c(0.015, 0.035))) +
    labs(x = if (show_x_axis) "Days since ARF onset" else NULL, y = if (show_y_axis) "Cumulative incidence" else NULL, color = "Quartile") +
    theme_cif_transplant() +
    theme(
      axis.title.y = if (show_y_axis) element_text(size = 19) else element_blank(),
      axis.text.y = if (show_y_axis) element_text(size = 17, color = "grey20") else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line(color = "black", linewidth = 0.35) else element_blank(),
      axis.title.x = if (show_x_axis) element_text(size = 18, margin = margin(t = 7)) else element_blank(),
      legend.position = "bottom",
      plot.margin = margin(4, 8, if (show_x_axis) 8 else 1, 8)
    )
}

make_aj_table_panel <- function(pollutant_value, outcome_value, show_x_axis = FALSE) {
  dat <- risk_event_table %>%
    filter(pollutant == pollutant_value, outcome == outcome_value) %>%
    mutate(exposure_quartile = factor(exposure_quartile, levels = rev(quartile_labels)))

  ggplot(dat, aes(x = time, y = exposure_quartile, label = label, color = exposure_quartile)) +
    geom_text(size = 4.75, fontface = "bold") +
    scale_color_manual(values = quartile_cols, drop = FALSE, guide = "none") +
    scale_x_continuous(breaks = tick_times, labels = label_number(accuracy = 1, trim = TRUE), limits = x_limits_cif, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_discrete(labels = function(x) str_replace(x, " lowest| highest", "")) +
    coord_cartesian(clip = "off") +
    labs(title = "Events / at risk", x = if (show_x_axis) "Days since ARF onset" else NULL, y = NULL) +
    theme_table_box() +
    theme(
      axis.text.x = if (show_x_axis) element_text(size = 16, color = "grey20", margin = margin(t = 2)) else element_blank(),
      axis.ticks.x = if (show_x_axis) element_line(color = "black", linewidth = 0.25) else element_blank(),
      legend.position = "none",
      plot.margin = margin(1, 38, if (show_x_axis) 6 else 2, 8)
    )
}

make_aj_cell <- function(pollutant_value, outcome_value, show_y_axis = TRUE, show_x_axis = FALSE) {
  make_aj_curve_panel(pollutant_value, outcome_value, show_y_axis = show_y_axis) /
    make_aj_table_panel(pollutant_value, outcome_value, show_x_axis = show_x_axis) +
    plot_layout(heights = c(4.2, 1.28))
}

make_cif_layout <- function(cell_fun, header_margin = c(12, 8, -18, 8)) {
  header_row <- wrap_plots(
    make_column_header("'Nitrogen dioxide'~(NO[2])"),
    make_column_header("'Fine particulate matter'~(PM[2.5])"),
    plot_spacer(),
    nrow = 1,
    widths = c(1, 1, 0.075)
  )
  row_plots <- lapply(seq_along(outcome_levels_cif), function(i) {
    outcome_value <- outcome_levels_cif[[i]]
    wrap_plots(
      cell_fun(pollutant_levels_cif[[1]], outcome_value, show_y_axis = TRUE, show_x_axis = i == length(outcome_levels_cif)),
      cell_fun(pollutant_levels_cif[[2]], outcome_value, show_y_axis = FALSE, show_x_axis = i == length(outcome_levels_cif)),
      make_row_header(outcome_value),
      nrow = 1,
      widths = c(1, 1, 0.075)
    )
  })
  (header_row / wrap_plots(row_plots, ncol = 1)) +
    plot_layout(heights = c(0.075, 3), guides = "collect") +
    plot_annotation(theme = theme(legend.position = "bottom", legend.text = element_text(size = 19), plot.margin = margin(8, 14, 8, 14))) &
    theme(legend.position = "bottom")
}

p_aj <- make_cif_layout(make_aj_cell)
write_plot(p_aj, "pooled_unadjusted_aalen_johansen_cif", width = 21, height = 19)

p_aj_curves_only <- make_cif_layout(function(pollutant_value, outcome_value, show_y_axis = TRUE, show_x_axis = FALSE) {
  make_aj_curve_panel(
    pollutant_value,
    outcome_value,
    show_y_axis = show_y_axis,
    show_x_axis = show_x_axis,
    x_limits_use = x_limits_cif_curves_only
  )
})
write_plot(p_aj_curves_only, "pooled_unadjusted_aalen_johansen_cif_curves_only", width = 21, height = 16)

fg_plot_data <- fg_curve_pooled %>%
  mutate(
    pollutant = factor(clean_pollutant(pollutant), levels = pollutant_levels_cif),
    outcome = recode(outcome, "Successful extubation" = "Extubation", "Persistent respiratory failure" = "Persistent RF"),
    outcome = factor(outcome, levels = outcome_levels_cif),
    exposure_quartile = factor(label_quartile(exposure_quartile), levels = quartile_labels)
  ) %>%
  filter(!is.na(pollutant), !is.na(outcome), !is.na(exposure_quartile), time <= 28)

fg_panel_scales <- fg_plot_data %>%
  group_by(pollutant, outcome) %>%
  summarise(y_max = max(cif, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_limit = pmax(y_max * 1.04, 0.01))

make_fg_curve_panel <- function(pollutant_value, outcome_value, show_y_axis = TRUE, show_x_axis = FALSE) {
  dat <- fg_plot_data %>% filter(pollutant == pollutant_value, outcome == outcome_value)
  y_limit <- fg_panel_scales %>% filter(pollutant == pollutant_value, outcome == outcome_value) %>% pull(y_limit)
  ggplot(dat, aes(time, cif, color = exposure_quartile)) +
    geom_line(linewidth = 1.05) +
    scale_color_manual(values = quartile_cols, drop = FALSE) +
    scale_x_continuous(limits = c(0, 28), breaks = tick_times, expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, y_limit), expand = expansion(mult = c(0.015, 0.035))) +
    labs(
      x = if (show_x_axis) "Days since ARF onset" else NULL,
      y = if (show_y_axis) "Adjusted cumulative incidence" else NULL,
      color = "Quartile"
    ) +
    theme_cif_transplant() +
    theme(
      axis.title.y = if (show_y_axis) element_text(size = 19) else element_blank(),
      axis.text.y = if (show_y_axis) element_text(size = 17, color = "grey20") else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line(color = "black", linewidth = 0.35) else element_blank(),
      axis.title.x = if (show_x_axis) element_text(size = 18, margin = margin(t = 7)) else element_blank(),
      legend.position = "bottom"
    )
}

p_fg <- make_cif_layout(function(pollutant_value, outcome_value, show_y_axis = TRUE, show_x_axis = FALSE) {
  make_fg_curve_panel(pollutant_value, outcome_value, show_y_axis = show_y_axis, show_x_axis = show_x_axis)
})
write_plot(p_fg, "pooled_adjusted_fine_gray_cif", width = 21, height = 19)

forest_secondary <- bind_rows(
  fine_gray_effects %>% mutate(row_type = "Site-specific"),
  fine_gray_pooled %>% mutate(site = "Pooled", row_type = "Pooled")
) %>%
  mutate(label = paste(site, term, sep = " | "))

p_forest_fg <- ggplot(forest_secondary, aes(x = estimate, y = fct_rev(label), xmin = conf_low, xmax = conf_high, colour = site)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(height = 0.16, linewidth = 0.6) +
  geom_point(size = 2.1) +
  scale_x_log10() +
  facet_grid(outcome ~ model, scales = "free_y", space = "free_y") +
  labs(x = "Subdistribution hazard ratio with 95% CI", y = NULL, colour = NULL) +
  theme_refer(9)
write_plot(p_forest_fg, "pooled_fine_gray_effect_forest", width = 15, height = 10)

p_forest_sens <- sensitivity_pooled %>%
  filter(term %in% c("pm25_per_5", "no2_per_10")) %>%
  filter(!str_detect(outcome, regex("imv|invasive|ventilation duration", ignore_case = TRUE))) %>%
  mutate(label = paste(sensitivity, outcome, analysis_model, term, sep = " | ")) %>%
  ggplot(aes(x = estimate, y = fct_rev(label), xmin = conf_low, xmax = conf_high, colour = pollutant)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(height = 0.18, linewidth = 0.6) +
  geom_point(size = 2.0) +
  scale_x_log10() +
  labs(x = "Pooled effect estimate with 95% CI", y = NULL, colour = NULL) +
  theme_refer(9)
write_plot(p_forest_sens, "pooled_sensitivity_effect_forest", width = 14, height = 11)

p_forest_subgroup <- subgroup_pooled %>%
  filter(pollutant %in% c("PM2.5", "NO2")) %>%
  mutate(label = paste(outcome, pollutant, subgroup, subgroup_level, sep = " | ")) %>%
  ggplot(aes(x = estimate, y = fct_rev(label), xmin = conf_low, xmax = conf_high, colour = pollutant)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(height = 0.15, linewidth = 0.5) +
  geom_point(size = 1.8) +
  scale_x_log10() +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(x = "Pooled subgroup-specific effect with 95% CI", y = NULL, colour = NULL) +
  theme_refer(8)
write_plot(p_forest_subgroup, "pooled_subgroup_effect_forest", width = 16, height = 12)

manifest <- tibble(
  pooled_output_dir = normalizePath(out_dir, mustWork = FALSE),
  generated_at = as.character(Sys.time()),
  n_sites = nrow(site_dirs),
  sites = paste(site_dirs$site, collapse = "; "),
  ucmc_dir = ucmc_dir,
  box_project_dir = box_project_dir
)
readr::write_csv(manifest, file.path(out_dir, "pooled_manifest.csv"))

message("Pooled outputs complete.")
message("Output directory: ", normalizePath(out_dir, mustWork = FALSE))
message("Figures: ", normalizePath(fig_dir, mustWork = FALSE))
