#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1 && nzchar(args[[1]])) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(file.path("output", "resubmission_pooled", "20260819_all_sites_with_michigan"), mustWork = TRUE)
}

pooled_long <- read_csv(file.path(out_dir, "pooled_table1_resubmission_long.csv"), show_col_types = FALSE)
chronic <- read_csv(file.path(out_dir, "pooled_table1_baseline_chronic_disease_prevalence.csv"), show_col_types = FALSE)
cohort_summary <- read_csv(file.path(out_dir, "site_primary_effect_estimates.csv"), show_col_types = FALSE)
site_dirs <- read_csv(file.path(out_dir, "pooled_site_directories.csv"), show_col_types = FALSE)

cohort_n <- pooled_long %>%
  filter(variable == "age") %>%
  summarise(n = max(n, na.rm = TRUE)) %>%
  pull(n)

fmt_mean_sd <- function(mean, sd, digits = 1) {
  if (!is.finite(mean) || !is.finite(sd)) return("")
  paste0(sprintf(paste0("%.", digits, "f"), mean), " +/- ", sprintf(paste0("%.", digits, "f"), sd))
}

fmt_median_iqr <- function(median, q25, q75, digits = 1) {
  if (!is.finite(median) || !is.finite(q25) || !is.finite(q75)) return("")
  paste0(sprintf(paste0("%.", digits, "f"), median), " [", sprintf(paste0("%.", digits, "f"), q25), ", ", sprintf(paste0("%.", digits, "f"), q75), "]")
}

fmt_n_pct <- function(n, denom, digits = 1) {
  if (!is.finite(n) || !is.finite(denom) || denom <= 0) return("")
  paste0(format(round(n), big.mark = ","), " (", sprintf(paste0("%.", digits, "f"), 100 * n / denom), "%)")
}

fmt_currency_mean_sd <- function(mean, sd) {
  if (!is.finite(mean) || !is.finite(sd)) return("")
  paste0("$", format(round(mean), big.mark = ","), " +/- $", format(round(sd), big.mark = ","))
}

continuous_row <- function(variable_name, label, digits = 1, style = c("mean_sd", "percent", "currency")) {
  style <- match.arg(style)
  x <- pooled_long %>% filter(variable_type == "continuous", .data$variable == variable_name) %>% slice(1)
  if (!nrow(x)) return(tibble(Characteristic = label, `Overall` = "Not exported"))
  value <- switch(
    style,
    mean_sd = fmt_mean_sd(x$mean, x$sd, digits),
    percent = fmt_mean_sd(100 * x$mean, 100 * x$sd, digits),
    currency = fmt_currency_mean_sd(x$mean, x$sd)
  )
  tibble(Characteristic = label, `Overall` = value)
}

continuous_median_row <- function(variable_name, label, digits = 1, style = c("raw", "percent")) {
  style <- match.arg(style)
  x <- pooled_long %>% filter(variable_type == "continuous", .data$variable == variable_name) %>% slice(1)
  if (!nrow(x)) return(tibble(Characteristic = label, `Overall` = "Not exported"))
  mult <- if (style == "percent") 100 else 1
  tibble(
    Characteristic = label,
    `Overall` = fmt_median_iqr(
      mult * x$median_site_weighted,
      mult * x$q25_site_weighted,
      mult * x$q75_site_weighted,
      digits
    )
  )
}

categorical_level_row <- function(variable_name, level_name, label = level_name) {
  x <- pooled_long %>%
    filter(
      variable_type == "categorical",
      .data$variable == variable_name,
      as.character(.data$level) == as.character(level_name)
    ) %>%
    slice(1)
  if (!nrow(x)) return(tibble(Characteristic = paste0("  ", label), `Overall` = "Not exported"))
  tibble(Characteristic = paste0("  ", label), `Overall` = fmt_n_pct(x$n, x$denominator))
}

binary_true_row <- function(variable_name, label) {
  categorical_level_row(variable_name, "TRUE", label)
}

race_other_unknown <- pooled_long %>%
  filter(
    variable_type == "categorical",
    variable == "race_ethnicity_simple",
    level %in% c("Other/Unknown", "Hispanic Other/Unknown", "(Missing)")
  ) %>%
  summarise(n = sum(n, na.rm = TRUE), denominator = cohort_n, .groups = "drop")

sex_other_unknown <- pooled_long %>%
  filter(
    variable_type == "categorical",
    variable == "sex_category",
    level %in% c("Other/Unknown", "Other", "(Missing)")
  ) %>%
  summarise(n = sum(n, na.rm = TRUE), denominator = cohort_n, .groups = "drop")

section <- function(label) tibble(Characteristic = label, `Overall` = "")

charlson_summary <- read_csv(file.path(out_dir, "pooled_site_directories.csv"), show_col_types = FALSE) %>%
  transmute(site) %>%
  left_join(
    bind_rows(lapply(unique(site_dirs$site_dir), function(site_dir) {
      path <- file.path(site_dir, "cohort_summary_reviewer_optimized.csv")
      if (!file.exists(path)) return(tibble())
      read_csv(path, show_col_types = FALSE)
    })),
    by = "site"
  )

charlson_row <- charlson_summary %>%
  filter(is.finite(n_arf), is.finite(mean_charlson), is.finite(median_charlson)) %>%
  summarise(
    mean_charlson = weighted.mean(mean_charlson, n_arf, na.rm = TRUE),
    median_charlson = weighted.mean(median_charlson, n_arf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  transmute(
    Characteristic = "Charlson comorbidity index",
    `Overall` = paste0(sprintf("%.1f", mean_charlson), " mean; ", sprintf("%.1f", median_charlson), " median")
  )

table1 <- bind_rows(
  section("Demographics"),
  continuous_row("age", "Age, years"),
  section("Sex, n (%)"),
  categorical_level_row("sex_category", "Female", "Female"),
  categorical_level_row("sex_category", "Male", "Male"),
  tibble(Characteristic = "  Other/Unknown", `Overall` = fmt_n_pct(sex_other_unknown$n, sex_other_unknown$denominator)),
  section("Race and ethnicity, n (%)"),
  categorical_level_row("race_ethnicity_simple", "Asian", "Asian"),
  categorical_level_row("race_ethnicity_simple", "Hispanic White", "Hispanic White"),
  categorical_level_row("race_ethnicity_simple", "Non-Hispanic Black", "Non-Hispanic Black"),
  categorical_level_row("race_ethnicity_simple", "Non-Hispanic White", "Non-Hispanic White"),
  tibble(Characteristic = "  Other/Unknown", `Overall` = fmt_n_pct(race_other_unknown$n, race_other_unknown$denominator)),
  section("Air pollution exposure"),
  continuous_row("pm25_mean", "PM\u2082\u2024\u2085 annual mean, ug/m3"),
  continuous_row("no2_mean", "NO2 annual mean, ppb"),
  section("Neighborhood social vulnerability and ZCTA characteristics"),
  continuous_row("svi_overall", "CDC SVI overall"),
  continuous_row("acs_pct_poverty", "Poverty, %", style = "percent"),
  continuous_row("acs_pct_unemployed", "Unemployed, %", style = "percent"),
  continuous_row("acs_pct_no_vehicle", "Households without vehicle, %", style = "percent"),
  continuous_row("acs_pct_nonwhite", "Non-White population, %", style = "percent"),
  continuous_row("acs_pct_black", "Black population, %", style = "percent"),
  continuous_row("acs_pct_asian", "Asian population, %", style = "percent"),
  continuous_row("acs_pct_hispanic", "Hispanic/Latino population, %", style = "percent"),
  continuous_row("acs_pct_bachelor_plus", "Bachelor's degree or higher, %", style = "percent"),
  continuous_row("acs_median_household_income", "Median household income", style = "currency"),
  section("Clinical severity, utilization, and outcomes"),
  charlson_row,
  continuous_row("sofa_total", "SOFA total score"),
  continuous_row("ventilator_free_days", "Ventilator-free days through day 28"),
  continuous_median_row("icu_los_days", "ICU length of stay, days, median [IQR]"),
  continuous_median_row("hosp_los_days", "Hospital length of stay, days, median [IQR]"),
  categorical_level_row("death_28d", "1", "Death by day 28"),
  categorical_level_row("in_hosp_death", "1", "In-hospital death"),
  section("Baseline chronic cardiopulmonary and comorbidity burden, n (%)"),
  binary_true_row("chronic_pulmonary_disease", "Chronic pulmonary disease"),
  binary_true_row("congestive_heart_failure", "Congestive heart failure"),
  binary_true_row("myocardial_infarction", "Myocardial infarction"),
  binary_true_row("peripheral_vascular_disease", "Peripheral vascular disease"),
  binary_true_row("cerebrovascular_disease", "Cerebrovascular disease"),
  binary_true_row("renal_disease", "Renal disease"),
  binary_true_row("diabetes_without_complication", "Diabetes without chronic complication"),
  binary_true_row("diabetes_with_complication", "Diabetes with chronic complication"),
  binary_true_row("any_malignancy", "Any malignancy"),
  binary_true_row("metastatic_solid_tumor", "Metastatic solid tumor")
) %>%
  mutate(`Overall` = str_replace_all(`Overall`, " \\+/- ", " +/- "))

names(table1)[2] <- paste0("Overall (N=", format(cohort_n, big.mark = ","), ")")

notes <- tibble(
  Note = c(
    "Values are pooled across 10 CLIF sites using aggregate site-level exports.",
    "Continuous variables are shown as pooled mean +/- SD unless otherwise noted.",
    "Exact pooled quantiles cannot be reconstructed from aggregate site exports; rows labeled median [IQR] use site-size-weighted medians and quartiles.",
    "Charlson comorbidity index is shown from cohort-summary exports as site-size-weighted mean and median because Charlson SD was not exported in the 10-site run. Future site Table 1 exports now include Charlson directly.",
    "Hispanic Other/Unknown, Other/Unknown, and missing race/ethnicity are combined as Other/Unknown for the manuscript-facing table."
  )
)

write_csv(table1, file.path(out_dir, "table1_manuscript_overall.csv"))
write_csv(notes, file.path(out_dir, "table1_manuscript_notes.csv"))

message("Wrote: ", file.path(out_dir, "table1_manuscript_overall.csv"))
message("Wrote: ", file.path(out_dir, "table1_manuscript_notes.csv"))
