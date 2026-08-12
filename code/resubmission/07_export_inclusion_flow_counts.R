#!/usr/bin/env Rscript

# Consolidate site-level inclusion-flow counts from the resubmission pipeline.
#
# Usage:
#   Rscript --vanilla code/resubmission/07_export_inclusion_flow_counts.R <output_dir>
#
# This is intentionally an aggregate-only export. It reads the site-specific
# CSVs already produced by the pipeline and writes standardized long and wide
# flow tables that can be pooled across CLIF sites without row-level data.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else {
  stop("Usage: Rscript --vanilla code/resubmission/07_export_inclusion_flow_counts.R <output_dir>")
}
out_dir <- normalizePath(out_dir, mustWork = TRUE)
generated_at <- as.character(Sys.time())

read_optional_csv <- function(filename) {
  path <- file.path(out_dir, filename)
  if (!file.exists(path)) return(tibble())
  readr::read_csv(path, show_col_types = FALSE)
}

first_present_value <- function(df, col, default = NA) {
  if (!nrow(df) || !(col %in% names(df))) return(default)
  x <- df[[col]]
  x <- x[!is.na(x)]
  if (!length(x)) default else x[[1]]
}

count_complete_rows <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (!length(vars) || !nrow(df)) return(NA_integer_)
  sum(stats::complete.cases(df[, vars, drop = FALSE]))
}

add_flow_row <- function(step_order, step, n_hospitalizations = NA_integer_,
                         n_patients = NA_integer_, n_excluded_from_prior = NA_integer_,
                         source_file = NA_character_, notes = NA_character_) {
  tibble(
    site = site_name,
    step_order = step_order,
    step = step,
    n_hospitalizations_or_encounters = as.integer(round(n_hospitalizations)),
    n_patients = as.integer(round(n_patients)),
    n_excluded_from_prior = as.integer(round(n_excluded_from_prior)),
    source_file = source_file,
    notes = notes
  )
}

analysis_df <- read_optional_csv("analysis_dataset_reviewer_optimized.csv")
primary_dataset <- read_optional_csv("resubmission_analysis_dataset.csv")
no_icu_dataset <- read_optional_csv("resubmission_analysis_dataset_no_icu_los_restriction.csv")
all_clif_overall <- read_optional_csv("supplement_all_clif_raw_mortality_overall.csv")
raw_cohort_summary <- read_optional_csv("resubmission_cohort_summary.csv")
raw_no_icu_summary <- read_optional_csv("resubmission_cohort_summary_no_icu_los_restriction.csv")
model_cohort_summary <- read_optional_csv("cohort_summary_reviewer_optimized.csv")
imv_assignment <- read_optional_csv("competing_risk_exhaustive_imv_assignment_diagnostics.csv")
fine_gray_results <- read_optional_csv("fine_gray_results_same_covariates_as_cox.csv")
primary_mortality_results <- read_optional_csv("primary_mortality_day28_logistic_results.csv")
primary_vfd_results <- read_optional_csv("primary_vfd_quasipoisson_results.csv")
sensitivity_cohorts <- read_optional_csv("primary_sensitivity_models_cohort_summaries.csv")

site_name <- first_present_value(analysis_df, "site") %||%
  first_present_value(primary_mortality_results, "site") %||%
  first_present_value(model_cohort_summary, "site") %||%
  first_present_value(raw_cohort_summary, "site") %||%
  first_present_value(all_clif_overall, "site") %||%
  "site"

default_adjustment_covars <- c(
  "age_10",
  "sex",
  "race_ethnicity",
  "charlson_score",
  "index_year_f",
  "acs_pct_poverty",
  "acs_pct_unemployed",
  "acs_pct_no_vehicle",
  "acs_pct_nonwhite",
  "acs_median_household_income_10k",
  "acs_pct_bachelor_plus"
)

adjustment_covars <- first_present_value(model_cohort_summary, "adjustment_covariates", NA_character_)
adjustment_covars <- if (!is.na(adjustment_covars)) {
  trimws(strsplit(adjustment_covars, "\\s*\\+\\s*")[[1]])
} else {
  default_adjustment_covars
}
if (nrow(analysis_df)) {
  adjustment_covars <- intersect(adjustment_covars, names(analysis_df))
}

primary_complete_vars <- c(
  "mortality_day28_event",
  "ventilator_free_days",
  "imv_days_through_vfd",
  "pm25_per_5",
  "no2_per_10",
  adjustment_covars
)
primary_complete <- if (nrow(analysis_df)) {
  stats::complete.cases(analysis_df[, intersect(primary_complete_vars, names(analysis_df)), drop = FALSE])
} else {
  logical()
}

competing_complete_vars <- c(
  "ftime_days",
  "event_code",
  "pm25_per_5",
  "no2_per_10",
  adjustment_covars
)
competing_pool <- if (nrow(analysis_df) && "has_imv_after_arf" %in% names(analysis_df)) {
  analysis_df %>% filter(.data$has_imv_after_arf %in% TRUE)
} else {
  tibble()
}
competing_complete <- if (nrow(competing_pool)) {
  stats::complete.cases(competing_pool[, intersect(competing_complete_vars, names(competing_pool)), drop = FALSE])
} else {
  logical()
}

n_all_clif_hosp <- first_present_value(all_clif_overall, "n_hospitalizations")
n_all_clif_patients <- first_present_value(all_clif_overall, "n_patients")
n_arf_no_icu <- first_present_value(raw_no_icu_summary, "n_arf", nrow(no_icu_dataset))
n_patients_arf_no_icu <- first_present_value(
  raw_no_icu_summary, "n_patients",
  if (nrow(no_icu_dataset) && "patient_id" %in% names(no_icu_dataset)) n_distinct(no_icu_dataset$patient_id) else NA_integer_
)
n_arf_primary <- first_present_value(raw_cohort_summary, "n_arf", nrow(analysis_df))
n_patients_arf_primary <- first_present_value(
  raw_cohort_summary, "n_patients",
  if (nrow(analysis_df) && "patient_id" %in% names(analysis_df)) n_distinct(analysis_df$patient_id) else NA_integer_
)
n_excluded_icu_los <- if (!is.na(n_arf_no_icu) && !is.na(n_arf_primary)) n_arf_no_icu - n_arf_primary else NA_integer_
n_with_pm25 <- first_present_value(raw_cohort_summary, "n_with_pm25", count_complete_rows(analysis_df, "pm25_per_5"))
n_with_no2 <- first_present_value(raw_cohort_summary, "n_with_no2", count_complete_rows(analysis_df, "no2_per_10"))
n_with_pm25_no2 <- first_present_value(
  model_cohort_summary,
  "n_with_pm25_no2",
  count_complete_rows(analysis_df, c("pm25_per_5", "no2_per_10"))
)
n_with_primary_covars <- first_present_value(
  model_cohort_summary,
  "n_with_primary_adjustment_covariates",
  count_complete_rows(analysis_df, adjustment_covars)
)
n_primary_complete <- if (nrow(primary_mortality_results) && "n" %in% names(primary_mortality_results)) {
  max(primary_mortality_results$n, na.rm = TRUE)
} else if (!is.na(first_present_value(model_cohort_summary, "n_primary_complete_cases"))) {
  first_present_value(model_cohort_summary, "n_primary_complete_cases")
} else {
  sum(primary_complete)
}
n_primary_day28_deaths <- if (nrow(primary_mortality_results) && "events" %in% names(primary_mortality_results)) {
  max(primary_mortality_results$events, na.rm = TRUE)
} else if (n_primary_complete && "mortality_day28_event" %in% names(analysis_df)) {
  sum(analysis_df$mortality_day28_event[primary_complete] %in% TRUE, na.rm = TRUE)
} else {
  NA_integer_
}
n_vfd_complete <- if (nrow(primary_vfd_results) && "n" %in% names(primary_vfd_results)) {
  max(primary_vfd_results$n, na.rm = TRUE)
} else {
  count_complete_rows(analysis_df, c("ventilator_free_days", "imv_days_through_vfd", "pm25_per_5", "no2_per_10", adjustment_covars))
}
n_imv_after_arf <- first_present_value(raw_cohort_summary, "n_with_imv_after_arf", first_present_value(imv_assignment, "n_with_imv_after_arf"))
n_competing_complete <- if (nrow(fine_gray_results) && "n" %in% names(fine_gray_results)) {
  max(fine_gray_results$n, na.rm = TRUE)
} else if (!is.na(first_present_value(model_cohort_summary, "n_imv_competing_risk_complete_cases"))) {
  first_present_value(model_cohort_summary, "n_imv_competing_risk_complete_cases")
} else {
  sum(competing_complete)
}
n_primary_complete_patients <- first_present_value(model_cohort_summary, "n_primary_complete_case_patients")
n_competing_complete_patients <- first_present_value(model_cohort_summary, "n_imv_competing_risk_complete_case_patients")
n_fg_events <- if (nrow(fine_gray_results)) {
  fine_gray_results %>%
    group_by(outcome) %>%
    summarise(events = max(events, na.rm = TRUE), .groups = "drop") %>%
    summarise(total = sum(events, na.rm = TRUE), .groups = "drop") %>%
    pull(total)
} else {
  NA_integer_
}

flow_long <- bind_rows(
  add_flow_row(
    1,
    "All CLIF hospitalizations available at site",
    n_all_clif_hosp,
    n_all_clif_patients,
    source_file = "supplement_all_clif_raw_mortality_overall.csv",
    notes = "Aggregate site-level CLIF denominator; not restricted to ARF."
  ),
  add_flow_row(
    2,
    "ARF onset within 24 hours of admission, before ICU LOS restriction",
    n_arf_no_icu,
    n_patients_arf_no_icu,
    source_file = "resubmission_cohort_summary_no_icu_los_restriction.csv",
    notes = "Sensitivity cohort used to quantify the 24-hour ICU LOS restriction."
  ),
  add_flow_row(
    3,
    "Excluded by ICU length of stay less than 24 hours",
    n_excluded_icu_los,
    NA_integer_,
    n_excluded_icu_los,
    source_file = "resubmission_cohort_summary_no_icu_los_restriction.csv; resubmission_cohort_summary.csv",
    notes = "Difference between the no-ICU-LOS-restriction ARF cohort and the primary ARF cohort."
  ),
  add_flow_row(
    4,
    "Primary ARF cohort with ICU length of stay at least 24 hours",
    n_arf_primary,
    n_patients_arf_primary,
    source_file = "resubmission_cohort_summary.csv",
    notes = "Primary analytic cohort before model complete-case restrictions."
  ),
  add_flow_row(
    5,
    "Primary ARF cohort with PM2.5 exposure",
    n_with_pm25,
    NA_integer_,
    if (!is.na(n_arf_primary) && !is.na(n_with_pm25)) n_arf_primary - n_with_pm25 else NA_integer_,
    source_file = "resubmission_cohort_summary.csv",
    notes = "ZCTA-level 12-month PM2.5 exposure available."
  ),
  add_flow_row(
    6,
    "Primary ARF cohort with NO2 exposure",
    n_with_no2,
    NA_integer_,
    if (!is.na(n_arf_primary) && !is.na(n_with_no2)) n_arf_primary - n_with_no2 else NA_integer_,
    source_file = "resubmission_cohort_summary.csv",
    notes = "ZCTA-level 12-month NO2 exposure available; annual NO2 used before monthly data are available."
  ),
  add_flow_row(
    7,
    "Primary ARF cohort with both PM2.5 and NO2 exposure",
    n_with_pm25_no2,
    NA_integer_,
    if (!is.na(n_arf_primary) && !is.na(n_with_pm25_no2)) n_arf_primary - n_with_pm25_no2 else NA_integer_,
    source_file = "cohort_summary_reviewer_optimized.csv",
    notes = "Exposure-complete denominator exported as an aggregate count; no row-level dataset required."
  ),
  add_flow_row(
    8,
    "Primary ARF cohort with complete adjustment covariates",
    n_with_primary_covars,
    NA_integer_,
    if (!is.na(n_arf_primary) && !is.na(n_with_primary_covars)) n_arf_primary - n_with_primary_covars else NA_integer_,
    source_file = "cohort_summary_reviewer_optimized.csv",
    notes = paste(adjustment_covars, collapse = " + ")
  ),
  add_flow_row(
    9,
    "Final primary complete-case mortality and VFD model cohort",
    n_primary_complete,
    if (!is.na(n_primary_complete_patients)) {
      n_primary_complete_patients
    } else if (nrow(analysis_df) && "patient_id" %in% names(analysis_df)) {
      n_distinct(analysis_df$patient_id[primary_complete])
    } else {
      NA_integer_
    },
    if (!is.na(n_arf_primary) && !is.na(n_primary_complete)) n_arf_primary - n_primary_complete else NA_integer_,
    source_file = "primary_mortality_day28_logistic_results.csv; primary_vfd_quasipoisson_results.csv",
    notes = "Complete mortality day 28, VFD, IMV duration, PM2.5, NO2, and primary adjustment covariates."
  ),
  add_flow_row(
    10,
    "Day-28 deaths in final primary model cohort",
    n_primary_day28_deaths,
    NA_integer_,
    source_file = "primary_mortality_day28_logistic_results.csv",
    notes = "Primary mortality outcome count among complete cases."
  ),
  add_flow_row(
    11,
    "Patients with invasive mechanical ventilation after ARF onset",
    n_imv_after_arf,
    NA_integer_,
    if (!is.na(n_arf_primary) && !is.na(n_imv_after_arf)) n_arf_primary - n_imv_after_arf else NA_integer_,
    source_file = "resubmission_cohort_summary.csv; competing_risk_exhaustive_imv_assignment_diagnostics.csv",
    notes = "Denominator for secondary extubation/death/persistent RF competing-risk analyses."
  ),
  add_flow_row(
    12,
    "Final complete-case adjusted Fine-Gray model cohort",
    n_competing_complete,
    if (!is.na(n_competing_complete_patients)) {
      n_competing_complete_patients
    } else if (nrow(competing_pool) && "patient_id" %in% names(competing_pool)) {
      n_distinct(competing_pool$patient_id[competing_complete])
    } else {
      NA_integer_
    },
    if (!is.na(n_imv_after_arf) && !is.na(n_competing_complete)) n_imv_after_arf - n_competing_complete else NA_integer_,
    source_file = "fine_gray_results_same_covariates_as_cox.csv",
    notes = "Complete IMV competing-risk time/event, PM2.5, NO2, and primary adjustment covariates."
  ),
  add_flow_row(
    13,
    "Fine-Gray competing-risk events represented across extubation, death, and persistent RF",
    n_fg_events,
    NA_integer_,
    source_file = "fine_gray_results_same_covariates_as_cox.csv",
    notes = "Sum of outcome-specific event counts from adjusted Fine-Gray output."
  ),
  add_flow_row(
    14,
    "IMV competing-risk assignment: successful extubation",
    first_present_value(imv_assignment, "n_extubation"),
    NA_integer_,
    source_file = "competing_risk_exhaustive_imv_assignment_diagnostics.csv",
    notes = "Exhaustive IMV outcome assignment before Fine-Gray complete-case restriction."
  ),
  add_flow_row(
    15,
    "IMV competing-risk assignment: death or hospice",
    first_present_value(imv_assignment, "n_death"),
    NA_integer_,
    source_file = "competing_risk_exhaustive_imv_assignment_diagnostics.csv",
    notes = "Exhaustive IMV outcome assignment before Fine-Gray complete-case restriction."
  ),
  add_flow_row(
    16,
    "IMV competing-risk assignment: persistent respiratory failure",
    first_present_value(imv_assignment, "n_persistent_rf"),
    NA_integer_,
    source_file = "competing_risk_exhaustive_imv_assignment_diagnostics.csv",
    notes = "Exhaustive IMV outcome assignment before Fine-Gray complete-case restriction."
  ),
  add_flow_row(
    17,
    "IMV competing-risk assignment: unassigned/censored after exhaustive recode",
    first_present_value(imv_assignment, "n_event_code_0_after_exhaustive_recode"),
    NA_integer_,
    source_file = "competing_risk_exhaustive_imv_assignment_diagnostics.csv",
    notes = "Should be zero under the updated exhaustive IMV outcome definition."
  )
) %>%
  mutate(generated_at = generated_at, .after = site)

sensitivity_flow <- if (nrow(sensitivity_cohorts)) {
  sensitivity_cohorts %>%
    transmute(
      site = site_name,
      generated_at = generated_at,
      step_order = 100 + row_number(),
      step = paste0("Sensitivity cohort: ", sensitivity),
      n_hospitalizations_or_encounters = as.integer(n_arf),
      n_patients = as.integer(n_patients),
      n_excluded_from_prior = NA_integer_,
      source_file = "primary_sensitivity_models_cohort_summaries.csv",
      notes = paste0(
        "Day-28 deaths=", n_day28_deaths,
        "; day-14 deaths=", n_day14_deaths,
        "; with SOFA=", n_with_sofa_total,
        "; with O3=", n_with_o3
      )
    )
} else {
  tibble()
}

flow_long <- bind_rows(flow_long, sensitivity_flow) %>%
  arrange(step_order)

flow_wide <- flow_long %>%
  select(site, generated_at, step_order, step, n_hospitalizations_or_encounters) %>%
  mutate(step_key = paste0("step_", sprintf("%03d", step_order), "_n")) %>%
  select(site, generated_at, step_key, n_hospitalizations_or_encounters) %>%
  pivot_wider(names_from = step_key, values_from = n_hospitalizations_or_encounters)

flow_dictionary <- flow_long %>%
  select(step_order, step, source_file, notes) %>%
  distinct() %>%
  arrange(step_order)

readr::write_csv(flow_long, file.path(out_dir, "site_inclusion_flow_counts.csv"))
readr::write_csv(flow_wide, file.path(out_dir, "site_inclusion_flow_counts_wide.csv"))
readr::write_csv(flow_dictionary, file.path(out_dir, "site_inclusion_flow_counts_dictionary.csv"))

message("Wrote site inclusion-flow exports:")
message(" - ", file.path(out_dir, "site_inclusion_flow_counts.csv"))
message(" - ", file.path(out_dir, "site_inclusion_flow_counts_wide.csv"))
message(" - ", file.path(out_dir, "site_inclusion_flow_counts_dictionary.csv"))
