#!/usr/bin/env Rscript

# Subgroup-specific estimates and formal interaction tests for reviewer response.
#
# Primary outcomes:
#   - Mortality by day 28 after ARF onset: logistic regression, reported as OR
#   - Ventilator-free days through day 28: quasi-Poisson regression, reported as
#     mean ratio
#   - IMV duration through day 28: quasi-Poisson regression, reported as mean ratio
#
# Subgroups:
#   - Sex
#   - Race/ethnicity
#   - ARF subtype

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(forcats)
  library(glue)
  library(janitor)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !file.exists(args[[1]])) {
  stop("Usage: Rscript code/resubmission/04_subgroup_interaction_estimates.R <analysis_dataset.csv> [output_dir]")
}
arg_or <- function(i, default = NA_character_) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
repo <- normalizePath(file.path(dirname(input_path), "..", ".."), mustWork = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- arg_or(2, file.path(repo, "output", "resubmission", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading analysis dataset: ", input_path)
message("Output: ", out_dir)

drop_uninformative_covars <- function(model_df, covars) {
  purrr::keep(covars, function(covar) {
    x <- model_df[[covar]]
    if (is.factor(x) || is.character(x)) return(dplyr::n_distinct(x, na.rm = TRUE) >= 2)
    isTRUE(stats::var(as.numeric(x), na.rm = TRUE) > 0)
  }) %>% unname()
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

format_estimate <- function(estimate, conf_low, conf_high) {
  sprintf("%.2f (%.2f, %.2f)", estimate, conf_low, conf_high)
}

make_formula <- function(outcome, rhs_terms) {
  stats::as.formula(paste(outcome, "~", paste(rhs_terms, collapse = " + ")))
}

fit_glm_safe <- function(formula, data, family) {
  tryCatch(
    withCallingHandlers(
      stats::glm(formula, data = data, family = family),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    warning = function(w) {
      warning(conditionMessage(w))
      NULL
    },
    error = function(e) {
      warning(conditionMessage(e))
      NULL
    }
  )
}

extract_term <- function(fit, term, outcome_label) {
  if (is.null(fit)) return(tibble())

  broom::tidy(fit) %>%
    filter(.data$term == !!term) %>%
    mutate(
      beta = .data$estimate,
      std_error = .data$std.error
    ) %>%
    transmute(
      beta = .data$beta,
      std_error = .data$std_error,
      statistic = .data$statistic,
      p_value = .data$p.value,
      estimate = exp(.data$beta),
      conf_low = exp(.data$beta - 1.96 * .data$std_error),
      conf_high = exp(.data$beta + 1.96 * .data$std_error),
      estimand = if_else(outcome_label == "Mortality by day 28", "odds ratio", "mean ratio")
    )
}

model_df_for <- function(data, outcome, terms) {
  data %>%
    dplyr::select(all_of(c(outcome, terms))) %>%
    tidyr::drop_na() %>%
    mutate(
      across(any_of(c("sex", "race_ethnicity", "arf_subtype", "index_year_f")), droplevels)
    )
}

subgroup_estimate <- function(data, outcome_label, outcome, family, pollutant, term,
                              other_term, subgroup_var, subgroup_level, base_covars) {
  subgroup_df <- data %>%
    filter(.data[[subgroup_var]] == subgroup_level)

  candidate_covars <- setdiff(c(other_term, base_covars), subgroup_var)
  model_df <- model_df_for(subgroup_df, outcome, c(term, candidate_covars))
  covars <- drop_uninformative_covars(model_df, candidate_covars)
  rhs_terms <- c(term, covars)

  if (nrow(model_df) < 50 || dplyr::n_distinct(model_df[[outcome]], na.rm = TRUE) < 2) {
    return(tibble(
      outcome = outcome_label,
      pollutant = pollutant,
      exposure_term = term,
      exposure_scale = if_else(term == "no2_per_10", "per 10 ppb", "per 5 ug/m3"),
      subgroup = subgroup_var,
      subgroup_level = as.character(subgroup_level),
      n = nrow(model_df),
      events = if_else(.env$outcome == "mortality_day28_event", sum(model_df[[.env$outcome]] == 1L), NA_integer_),
      mean_vfd = if_else(.env$outcome == "ventilator_free_days", mean(model_df[[.env$outcome]], na.rm = TRUE), NA_real_),
      mean_imv_days = if_else(.env$outcome == "imv_days_through_vfd", mean(model_df[[.env$outcome]], na.rm = TRUE), NA_real_),
      beta = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      estimate = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      estimand = if_else(outcome_label == "Mortality by day 28", "odds ratio", "mean ratio"),
      estimate_ci = NA_character_,
      p_value_formatted = NA_character_,
      adjustment_covariates = paste(covars, collapse = "; ")
    ))
  }

  fit <- fit_glm_safe(make_formula(outcome, rhs_terms), model_df, family)
  extract_term(fit, term, outcome_label) %>%
    mutate(
      outcome = outcome_label,
      pollutant = pollutant,
      exposure_term = term,
      exposure_scale = if_else(term == "no2_per_10", "per 10 ppb", "per 5 ug/m3"),
      subgroup = subgroup_var,
      subgroup_level = as.character(subgroup_level),
      n = nrow(model_df),
      events = if_else(.env$outcome == "mortality_day28_event", sum(model_df[[.env$outcome]] == 1L), NA_integer_),
      mean_vfd = if_else(.env$outcome == "ventilator_free_days", mean(model_df[[.env$outcome]], na.rm = TRUE), NA_real_),
      mean_imv_days = if_else(.env$outcome == "imv_days_through_vfd", mean(model_df[[.env$outcome]], na.rm = TRUE), NA_real_),
      estimate_ci = format_estimate(.data$estimate, .data$conf_low, .data$conf_high),
      p_value_formatted = format_p(.data$p_value),
      adjustment_covariates = paste(covars, collapse = "; "),
      .before = 1
    )
}

interaction_test <- function(data, outcome_label, outcome, family, pollutant, term,
                             other_term, subgroup_var, base_covars) {
  candidate_covars <- setdiff(c(other_term, subgroup_var, base_covars), subgroup_var)
  candidate_terms <- c(term, subgroup_var, candidate_covars)
  model_df <- model_df_for(data, outcome, candidate_terms)
  covars <- drop_uninformative_covars(model_df, candidate_covars)
  rhs_base <- c(term, subgroup_var, covars)
  rhs_interaction <- c(paste0(term, " * ", subgroup_var), covars)

  if (
    nrow(model_df) < 50 ||
      dplyr::n_distinct(model_df[[outcome]], na.rm = TRUE) < 2 ||
      dplyr::n_distinct(model_df[[subgroup_var]], na.rm = TRUE) < 2
  ) {
    return(tibble(
      outcome = outcome_label,
      pollutant = pollutant,
      exposure_term = term,
      exposure_scale = if_else(term == "no2_per_10", "per 10 ppb", "per 5 ug/m3"),
      subgroup = subgroup_var,
      n = nrow(model_df),
      events = if_else(.env$outcome == "mortality_day28_event", sum(model_df[[.env$outcome]] == 1L), NA_integer_),
      test = NA_character_,
      df = NA_real_,
      statistic = NA_real_,
      p_interaction = NA_real_,
      p_interaction_formatted = NA_character_,
      base_model_terms = paste(rhs_base, collapse = " + "),
      interaction_model_terms = paste(rhs_interaction, collapse = " + ")
    ))
  }

  fit_base <- fit_glm_safe(make_formula(outcome, rhs_base), model_df, family)
  fit_interaction <- fit_glm_safe(make_formula(outcome, rhs_interaction), model_df, family)
  if (is.null(fit_base) || is.null(fit_interaction)) {
    p_int <- df <- stat <- NA_real_
    test_name <- NA_character_
  } else if (outcome == "mortality_day28_event") {
    test_tbl <- stats::anova(fit_base, fit_interaction, test = "Chisq")
    p_int <- test_tbl$`Pr(>Chi)`[[2]]
    df <- test_tbl$Df[[2]]
    stat <- test_tbl$Deviance[[2]]
    test_name <- "Likelihood ratio test"
  } else {
    test_tbl <- stats::anova(fit_base, fit_interaction, test = "F")
    p_int <- test_tbl$`Pr(>F)`[[2]]
    df <- test_tbl$Df[[2]]
    stat <- test_tbl$F[[2]]
    test_name <- "Partial F test"
  }

  tibble(
    outcome = outcome_label,
    pollutant = pollutant,
    exposure_term = term,
    exposure_scale = if_else(term == "no2_per_10", "per 10 ppb", "per 5 ug/m3"),
    subgroup = subgroup_var,
    n = nrow(model_df),
    events = if_else(.env$outcome == "mortality_day28_event", sum(model_df[[.env$outcome]] == 1L), NA_integer_),
    test = test_name,
    df = df,
    statistic = stat,
    p_interaction = p_int,
    p_interaction_formatted = format_p(p_int),
    base_model_terms = paste(rhs_base, collapse = " + "),
    interaction_model_terms = paste(rhs_interaction, collapse = " + ")
  )
}

interaction_coefficients <- function(data, outcome_label, outcome, family, pollutant, term,
                                     other_term, subgroup_var, base_covars) {
  candidate_covars <- setdiff(c(other_term, subgroup_var, base_covars), subgroup_var)
  candidate_terms <- c(term, subgroup_var, candidate_covars)
  model_df <- model_df_for(data, outcome, candidate_terms)
  covars <- drop_uninformative_covars(model_df, candidate_covars)
  rhs_interaction <- c(paste0(term, " * ", subgroup_var), covars)

  fit <- fit_glm_safe(make_formula(outcome, rhs_interaction), model_df, family)
  if (is.null(fit)) return(tibble())

  broom::tidy(fit) %>%
    filter(str_detect(.data$term, fixed(":")), str_detect(.data$term, fixed(term))) %>%
    mutate(
      beta = .data$estimate,
      std_error = .data$std.error
    ) %>%
    transmute(
      outcome = outcome_label,
      pollutant = pollutant,
      exposure_term = .env$term,
      exposure_scale = if_else(.env$term == "no2_per_10", "per 10 ppb", "per 5 ug/m3"),
      subgroup = subgroup_var,
      interaction_term = .data$term,
      n = nrow(model_df),
      events = if_else(.env$outcome == "mortality_day28_event", sum(model_df[[.env$outcome]] == 1L), NA_integer_),
      beta = .data$beta,
      std_error = .data$std_error,
      statistic = .data$statistic,
      p_value = .data$p.value,
      estimate = exp(.data$beta),
      conf_low = exp(.data$beta - 1.96 * .data$std_error),
      conf_high = exp(.data$beta + 1.96 * .data$std_error),
      estimand = if_else(outcome_label == "Mortality by day 28", "interaction odds ratio ratio", "interaction mean ratio ratio"),
      interaction_model_terms = paste(rhs_interaction, collapse = " + ")
    )
}

analysis_df <- readr::read_csv(input_path, show_col_types = FALSE, progress = FALSE) %>%
  janitor::clean_names() %>%
  mutate(
    mortality_ftime_days = as.numeric(mortality_ftime_days),
    mortality_event = as.integer(mortality_event),
    mortality_day28_event = as.integer(
      if ("mortality_day28_event" %in% names(.)) {
        mortality_day28_event
      } else {
        mortality_event == 1L & mortality_ftime_days <= 28
      }
    ),
    ventilator_free_days = as.numeric(ventilator_free_days),
    imv_days_through_vfd = as.numeric(imv_days_through_vfd),
    age_10 = age_10 %||% (age / 10),
    pm25_per_5 = pm25_per_5 %||% (pm25_12m_zcta / 5),
    no2_per_10 = no2_per_10 %||% (no2_12m_zcta / 10),
    acs_median_household_income_10k = acs_median_household_income_10k %||% (acs_median_household_income / 10000),
    sex = factor(sex),
    race_ethnicity = forcats::fct_collapse(
      factor(race_ethnicity),
      "Other/Unknown" = c("Hispanic Other/Unknown", "Other/Unknown")
    ),
    race_ethnicity = forcats::fct_lump_min(race_ethnicity, min = 100, other_level = "Other/Unknown"),
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed")),
    index_year_f = factor(index_year_f %||% index_year)
  ) %>%
  filter(
    is.finite(mortality_ftime_days),
    mortality_ftime_days > 0,
    !is.na(mortality_day28_event),
    !is.na(ventilator_free_days),
    !is.na(imv_days_through_vfd)
  )

site_name <- dplyr::first(stats::na.omit(analysis_df$site)) %||% "site"

social_covars <- c(
  "acs_pct_poverty",
  "acs_pct_unemployed",
  "acs_pct_no_vehicle",
  "acs_pct_nonwhite",
  "acs_median_household_income_10k",
  "acs_pct_bachelor_plus"
)

base_covars <- c(
  "age_10",
  "sex",
  "race_ethnicity",
  "charlson_score",
  "index_year_f",
  intersect(social_covars, names(analysis_df))
)

pollutant_specs <- tibble::tribble(
  ~pollutant, ~term, ~other_term,
  "NO2", "no2_per_10", "pm25_per_5",
  "PM2.5", "pm25_per_5", "no2_per_10"
)

outcome_specs <- tibble::tribble(
  ~outcome_label, ~outcome, ~family,
  "Mortality by day 28", "mortality_day28_event", list(stats::binomial(link = "logit")),
  "Ventilator-free days", "ventilator_free_days", list(stats::quasipoisson(link = "log")),
  "IMV duration", "imv_days_through_vfd", list(stats::quasipoisson(link = "log"))
)

subgroup_vars <- c("sex", "race_ethnicity", "arf_subtype")

subgroup_grid <- tidyr::crossing(
  outcome_specs,
  pollutant_specs,
  subgroup = subgroup_vars
) %>%
  mutate(
    subgroup_levels = purrr::map(.data$subgroup, ~ levels(droplevels(analysis_df[[.x]])))
  ) %>%
  tidyr::unnest(subgroup_levels)

subgroup_estimates <- pmap_dfr(
  subgroup_grid,
  function(outcome_label, outcome, family, pollutant, term, other_term, subgroup, subgroup_levels) {
    subgroup_estimate(
      data = analysis_df,
      outcome_label = outcome_label,
      outcome = outcome,
      family = family[[1]],
      pollutant = pollutant,
      term = term,
      other_term = other_term,
      subgroup_var = subgroup,
      subgroup_level = subgroup_levels,
      base_covars = base_covars
    )
  }
)

interaction_grid <- tidyr::crossing(
  outcome_specs,
  pollutant_specs,
  subgroup = subgroup_vars
)

interaction_tests <- pmap_dfr(
  interaction_grid,
  function(outcome_label, outcome, family, pollutant, term, other_term, subgroup) {
    interaction_test(
      data = analysis_df,
      outcome_label = outcome_label,
      outcome = outcome,
      family = family[[1]],
      pollutant = pollutant,
      term = term,
      other_term = other_term,
      subgroup_var = subgroup,
      base_covars = base_covars
    )
  }
)

interaction_coef <- pmap_dfr(
  interaction_grid,
  function(outcome_label, outcome, family, pollutant, term, other_term, subgroup) {
    interaction_coefficients(
      data = analysis_df,
      outcome_label = outcome_label,
      outcome = outcome,
      family = family[[1]],
      pollutant = pollutant,
      term = term,
      other_term = other_term,
      subgroup_var = subgroup,
      base_covars = base_covars
    )
  }
)

subgroup_estimates <- subgroup_estimates %>%
  mutate(site = site_name, .before = 1)

interaction_tests <- interaction_tests %>%
  mutate(site = site_name, .before = 1)

interaction_coef <- interaction_coef %>%
  mutate(site = site_name, .before = 1)

manuscript_summary <- subgroup_estimates %>%
  left_join(
    interaction_tests %>%
      select(site, outcome, pollutant, subgroup, p_interaction, p_interaction_formatted),
    by = c("site", "outcome", "pollutant", "subgroup")
  ) %>%
  transmute(
    site,
    outcome,
    pollutant,
    exposure_scale,
    subgroup,
    subgroup_level,
    n,
    events,
    mean_vfd,
    mean_imv_days,
    estimand,
    estimate_ci,
    p_value = p_value_formatted,
    p_interaction = p_interaction_formatted
  )

readr::write_csv(
  subgroup_estimates,
  file.path(out_dir, "primary_subgroup_specific_estimates.csv")
)
readr::write_csv(
  interaction_tests,
  file.path(out_dir, "primary_subgroup_interaction_tests.csv")
)
readr::write_csv(
  interaction_coef,
  file.path(out_dir, "primary_subgroup_interaction_coefficients.csv")
)
readr::write_csv(
  manuscript_summary,
  file.path(out_dir, "primary_subgroup_estimates_for_manuscript.csv")
)

message("Wrote:")
message("  ", file.path(out_dir, "primary_subgroup_specific_estimates.csv"))
message("  ", file.path(out_dir, "primary_subgroup_interaction_tests.csv"))
message("  ", file.path(out_dir, "primary_subgroup_interaction_coefficients.csv"))
message("  ", file.path(out_dir, "primary_subgroup_estimates_for_manuscript.csv"))

interaction_tests %>%
  arrange(.data$p_interaction) %>%
  select(outcome, pollutant, subgroup, n, events, test, statistic, df, p_interaction_formatted) %>%
  print(n = Inf)
