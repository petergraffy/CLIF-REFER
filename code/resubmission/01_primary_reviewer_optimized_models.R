#!/usr/bin/env Rscript

# Reviewer-optimized REFER primary models
#
# Intended use:
#   Run the original cohort/linkage code through construction of `arf_exp`, then source
#   this script to replace the first-run model block with the reviewer-focused models.
#
# Core changes from the first-run analysis:
#   1. Use binary mortality by day 28 after ARF onset as the primary mortality
#      endpoint.
#   2. Keep Cox proportional hazards models for mortality after ARF onset as a
#      sensitivity analysis.
#   3. Make ventilator-free days through day 28 the primary ventilation endpoint,
#      while exporting supplemental raw IMV-duration models.
#   4. Use the expanded adjustment set: demographics, calendar year, Charlson, and
#      available social vulnerability covariates.

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(forcats)
  library(glue)
  library(janitor)
  library(lubridate)
  library(MASS)
  library(purrr)
  library(readr)
  library(rlang)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
site_name <- get0("site_name", inherits = TRUE) %||%
  (get0("config", inherits = TRUE)$site_name %||% "site")

repo <- get0("repo", inherits = TRUE) %||% normalizePath(getwd(), mustWork = FALSE)
cfg <- get0("cfg", inherits = TRUE)
base_output <- cfg$output_dir %||% file.path(repo, "output")
args <- commandArgs(trailingOnly = TRUE)
arg_or <- function(i, default = NA_character_) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}
out_dir <- arg_or(2, file.path(base_output, "resubmission", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vfd_days <- as.integer(get0("vfd_days", inherits = TRUE) %||% 28)
imv_gap_hours <- as.numeric(get0("imv_gap_hours", inherits = TRUE) %||% 6)
charlson_lookback_days <- as.integer(get0("charlson_lookback_days", inherits = TRUE) %||% 365)
charlson_include_index <- isTRUE(get0("charlson_include_index", inherits = TRUE) %||% TRUE)

message("Site: ", site_name)
message("Reviewer-optimized output: ", out_dir)

if (!exists("arf_exp", inherits = TRUE) && length(args) >= 1 && file.exists(args[[1]])) {
  message("Reading analytic frame from: ", normalizePath(args[[1]], mustWork = TRUE))
  arf_exp <- readr::read_csv(args[[1]], show_col_types = FALSE, progress = FALSE)
} else {
  stopifnot(exists("arf_exp", inherits = TRUE))
  arf_exp <- get("arf_exp", inherits = TRUE)
}
arf_exp <- arf_exp %>% janitor::clean_names()
if (identical(site_name, "site") && "site" %in% names(arf_exp)) {
  site_name <- dplyr::first(stats::na.omit(arf_exp$site)) %||% site_name
}

safe_ts <- get0("safe_ts", inherits = TRUE) %||% function(x, tz = Sys.timezone()) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = tz))
  if (is.numeric(x)) {
    x2 <- ifelse(x > 1e12, x / 1000, x)
    return(as.POSIXct(x2, origin = "1970-01-01", tz = tz))
  }
  suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c("ymd_HMS", "ymd_HM", "ymd", "ymdTz", "ymdT", "mdy_HMS", "mdy_HM", "mdy", "dmy_HMS", "dmy_HM", "dmy"),
    tz = tz,
    quiet = TRUE
  ))
}

first_present <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

icd_clean <- function(x) str_replace_all(str_to_upper(as.character(x)), "[^A-Z0-9]", "")

charlson_from_codes <- function(dx) {
  dx <- dx %>% mutate(code = icd_clean(diagnosis_code))
  has <- function(pattern) str_detect(dx$code, pattern)

  flags <- dx %>%
    transmute(
      patient_id,
      hospitalization_id,
      myocardial_infarction = has("^(I21|I22|I252)"),
      congestive_heart_failure = has("^(I50|I110|I130|I132|I255|I420|I425|I426|I427|I428|I429|P290)"),
      peripheral_vascular_disease = has("^(I70|I71|I731|I738|I739|I771|I790|I792|K551|K558|K559|Z958|Z959)"),
      cerebrovascular_disease = has("^(G45|G46|H340|I6[0-9])"),
      dementia = has("^(F00|F01|F02|F03|F051|G30|G311)"),
      chronic_pulmonary_disease = has("^(I278|I279|J4[0-7]|J6[0-7]|J684|J701|J703)"),
      rheumatic_disease = has("^(M05|M06|M315|M3[2-4]|M351|M353|M360)"),
      peptic_ulcer_disease = has("^(K25|K26|K27|K28)"),
      mild_liver_disease = has("^(B18|K700|K701|K702|K703|K709|K713|K714|K715|K717|K73|K74|K760|K762|K763|K764|K768|K769|Z944)"),
      diabetes_without_complication = has("^(E10[01689]|E11[01689]|E12[01689]|E13[01689]|E14[01689])"),
      diabetes_with_complication = has("^(E10[2-5]|E11[2-5]|E12[2-5]|E13[2-5]|E14[2-5])"),
      hemiplegia_paraplegia = has("^(G81|G82|G041|G114|G801|G802|G830|G831|G832|G833|G834|G839)"),
      renal_disease = has("^(I120|I131|N03|N052|N053|N054|N055|N056|N057|N18|N19|N250|Z490|Z491|Z492|Z940|Z992)"),
      any_malignancy = has("^(C0[0-9]|C1[0-9]|C2[0-6]|C3[0-4]|C37|C38|C39|C40|C41|C43|C4[5-9]|C5[0-8]|C6[0-9]|C7[0-6]|C81|C82|C83|C84|C85|C88|C90|C91|C92|C93|C94|C95|C96|C97)"),
      moderate_severe_liver_disease = has("^(I85|I864|I982|K704|K711|K721|K729|K765|K766|K767)"),
      metastatic_solid_tumor = has("^(C77|C78|C79|C80)"),
      aids_hiv = has("^(B20|B21|B22|B24)")
    ) %>%
    group_by(patient_id, hospitalization_id) %>%
    summarise(across(where(is.logical), ~ any(.x, na.rm = TRUE)), .groups = "drop")

  flags %>%
    mutate(
      diabetes_without_complication = diabetes_without_complication & !diabetes_with_complication,
      any_malignancy = any_malignancy & !metastatic_solid_tumor,
      mild_liver_disease = mild_liver_disease & !moderate_severe_liver_disease,
      charlson_score =
        1 * myocardial_infarction +
        1 * congestive_heart_failure +
        1 * peripheral_vascular_disease +
        1 * cerebrovascular_disease +
        1 * dementia +
        1 * chronic_pulmonary_disease +
        1 * rheumatic_disease +
        1 * peptic_ulcer_disease +
        1 * mild_liver_disease +
        1 * diabetes_without_complication +
        2 * diabetes_with_complication +
        2 * hemiplegia_paraplegia +
        2 * renal_disease +
        2 * any_malignancy +
        3 * moderate_severe_liver_disease +
        6 * metastatic_solid_tumor +
        6 * aids_hiv
    )
}

add_charlson_if_needed <- function(df) {
  if ("charlson_score" %in% names(df)) {
    return(df %>% mutate(charlson_score = coalesce(as.numeric(charlson_score), 0)))
  }

  if (!exists("diagnosis", inherits = TRUE)) {
    warning("No `charlson_score` column and no `diagnosis` object found; setting Charlson to 0.")
    return(df %>% mutate(charlson_score = 0))
  }

  diagnosis <- get("diagnosis", inherits = TRUE) %>% janitor::clean_names()
  diagnosis_code_col <- first_present(diagnosis, c("diagnosis_code", "dx_code", "icd_code", "code"))
  if (is.na(diagnosis_code_col)) {
    warning("Diagnosis table has no diagnosis code column; setting Charlson to 0.")
    return(df %>% mutate(charlson_score = 0))
  }

  diagnosis <- diagnosis %>%
    rename(diagnosis_code = all_of(diagnosis_code_col)) %>%
    filter(!is.na(patient_id), !is.na(hospitalization_id), !is.na(diagnosis_code))

  hosp_lookup <- df %>%
    dplyr::select(patient_id, hospitalization_id, t0) %>%
    distinct()

  dx_for_charlson <- diagnosis %>%
    inner_join(hosp_lookup, by = "patient_id", relationship = "many-to-many")

  if (exists("hospitalization", inherits = TRUE)) {
    hospitalization_raw <- get("hospitalization", inherits = TRUE) %>% janitor::clean_names()
    admission_col <- first_present(hospitalization_raw, c("admission_dttm", "hospital_admission_dttm"))
    if (!is.na(admission_col)) {
      hospitalization <- hospitalization_raw %>%
        transmute(
        hospitalization_id,
        dx_admission_dttm = safe_ts(.data[[admission_col]])
      )
      dx_for_charlson <- dx_for_charlson %>%
        left_join(hospitalization, by = c("hospitalization_id.x" = "hospitalization_id")) %>%
        filter(
          (hospitalization_id.x == hospitalization_id.y & charlson_include_index) |
            (!is.na(dx_admission_dttm) &
               dx_admission_dttm < t0 &
               dx_admission_dttm >= t0 - days(charlson_lookback_days))
        )
    }
  }

  charlson_tbl <- charlson_from_codes(
    dx_for_charlson %>%
      transmute(
        patient_id,
        hospitalization_id = hospitalization_id.y,
        diagnosis_code
      )
  ) %>%
    group_by(patient_id, hospitalization_id) %>%
    summarise(charlson_score = max(charlson_score, na.rm = TRUE), .groups = "drop")

  df %>%
    left_join(charlson_tbl, by = c("patient_id", "hospitalization_id")) %>%
    mutate(charlson_score = coalesce(charlson_score, 0))
}

coerce_analysis_frame <- function(df) {
  t0_col <- first_present(df, c("arf_onset", "index_admit", "admission_dttm"))
  discharge_col <- first_present(df, c("index_discharge", "discharge_dttm"))
  if (is.na(t0_col) || is.na(discharge_col)) {
    stop("Need ARF onset/index admission and discharge time columns in `arf_exp`.")
  }

  pm25_col <- first_present(df, c("pm25_12m_zcta", "pm25_12m_mean", "pm25_mean"))
  no2_col <- first_present(df, c("no2_12m_zcta", "no2_12m_mean", "no2_mean"))
  if (is.na(pm25_col) || is.na(no2_col)) {
    stop("Need PM2.5 and NO2 exposure columns in `arf_exp`.")
  }

  sex_col <- first_present(df, c("sex", "sex_category"))
  race_col <- first_present(df, c("race_ethnicity", "race_ethnicity_simple"))
  if (is.na(sex_col) || is.na(race_col)) {
    stop("Need sex and race/ethnicity columns in `arf_exp`.")
  }

  out <- df %>%
    mutate(
      t0 = safe_ts(.data[[t0_col]]),
      discharge_time = safe_ts(.data[[discharge_col]]),
      pm25_exposure = as.numeric(.data[[pm25_col]]),
      no2_exposure = as.numeric(.data[[no2_col]]),
      pm25_per_5 = pm25_exposure / 5,
      no2_per_10 = no2_exposure / 10,
      age_10 = as.numeric(age) / 10,
      sex = factor(.data[[sex_col]]),
      race_ethnicity = forcats::fct_lump_min(factor(.data[[race_col]]), min = 100, other_level = "Other/Unknown"),
      index_year = if ("index_year" %in% names(.)) as.integer(index_year) else lubridate::year(t0),
      index_year_f = factor(index_year)
    )

  out %>%
    add_charlson_if_needed()
}

add_mortality_time <- function(df) {
  if (all(c("mortality_ftime_days", "mortality_event") %in% names(df))) {
    return(df %>%
      mutate(
        mortality_ftime_days = pmax(as.numeric(mortality_ftime_days), 1 / 24),
        mortality_event = as.integer(mortality_event)
      ))
  }

  death_col <- first_present(df, c("death_time", "death_ts", "death_dttm_final", "death_dttm"))

  out <- df
  if (!is.na(death_col)) {
    out <- out %>% mutate(death_time = safe_ts(.data[[death_col]]))
  } else if (exists("patient", inherits = TRUE)) {
    patient <- get("patient", inherits = TRUE) %>%
      janitor::clean_names() %>%
      dplyr::select(patient_id, death_dttm)
    out <- out %>%
      left_join(patient, by = "patient_id") %>%
      mutate(death_time = safe_ts(death_dttm))
  } else {
    out <- out %>%
      mutate(death_time = as.POSIXct(NA_real_, origin = "1970-01-01", tz = Sys.timezone()))
  }

  if (!"in_hosp_death" %in% names(out)) {
    out <- out %>% mutate(in_hosp_death = as.integer(!is.na(death_time) & death_time <= discharge_time))
  }

  out %>%
    mutate(
      death_time = if_else(in_hosp_death == 1L & is.na(death_time), discharge_time, death_time),
      mortality_event = as.integer(!is.na(death_time) & death_time >= t0 & death_time <= discharge_time),
      mortality_end_time = if_else(mortality_event == 1L, death_time, discharge_time),
      mortality_ftime_days = as.numeric(difftime(mortality_end_time, t0, units = "days")),
      mortality_ftime_days = pmax(mortality_ftime_days, 1 / 24)
    )
}

add_day28_mortality <- function(df) {
  death_col <- first_present(df, c("death_time", "death_ts", "death_dttm_final", "death_dttm"))
  out <- df
  if (!is.na(death_col)) {
    out <- out %>% mutate(death_time_day28 = safe_ts(.data[[death_col]]))
  } else {
    out <- out %>% mutate(death_time_day28 = as.POSIXct(NA_real_, origin = "1970-01-01", tz = Sys.timezone()))
  }

  out %>%
    mutate(
      mortality_day28_event = as.integer(
        (
          !is.na(death_time_day28) &
            death_time_day28 >= t0 &
            death_time_day28 <= t0 + lubridate::days(28)
        ) |
          (
            mortality_event == 1L &
              !is.na(mortality_ftime_days) &
              mortality_ftime_days <= 28
          )
      ),
      alive_or_censored_before_day28 = as.integer(
        mortality_day28_event == 0L &
          !is.na(mortality_ftime_days) &
          mortality_ftime_days < 28
      )
    )
}

fallback_vfd_from_vent_hours <- function(df) {
  vent_hours <- if ("vent_hours" %in% names(df)) as.numeric(df$vent_hours) else rep(0, nrow(df))
  df %>%
    mutate(
      imv_days_through_vfd = pmin(coalesce(vent_hours, 0) / 24, vfd_days),
      imv_days_through_day14 = pmin(coalesce(vent_hours, 0) / 24, 14),
      died_before_vfd_horizon = mortality_event == 1L & mortality_ftime_days <= vfd_days,
      died_before_vfd_day14 = mortality_event == 1L & mortality_ftime_days <= 14,
      ventilator_free_days = if_else(
        died_before_vfd_horizon,
        0,
        pmax(vfd_days - imv_days_through_vfd, 0)
      ),
      ventilator_free_days_day14 = if_else(
        died_before_vfd_day14,
        0,
        pmax(14 - imv_days_through_day14, 0)
      )
    )
}

add_vfd <- function(df) {
  if ("ventilator_free_days" %in% names(df)) {
    return(df %>%
      mutate(
        ventilator_free_days = pmin(pmax(as.numeric(ventilator_free_days), 0), vfd_days),
        imv_days_through_vfd = if ("imv_days_through_vfd" %in% names(.)) as.numeric(imv_days_through_vfd) else NA_real_,
        died_before_vfd_horizon = if ("died_before_vfd_horizon" %in% names(.)) died_before_vfd_horizon else NA,
        imv_days_through_day14 = if ("imv_days_through_day14" %in% names(.)) as.numeric(imv_days_through_day14) else NA_real_,
        died_before_vfd_day14 = if ("died_before_vfd_day14" %in% names(.)) died_before_vfd_day14 else mortality_event == 1L & mortality_ftime_days <= 14,
        ventilator_free_days_day14 = if ("ventilator_free_days_day14" %in% names(.)) {
          pmin(pmax(as.numeric(ventilator_free_days_day14), 0), 14)
        } else if ("imv_days_through_day14" %in% names(.)) {
          if_else(died_before_vfd_day14, 0, pmax(14 - imv_days_through_day14, 0))
        } else {
          NA_real_
        }
      ))
  }

  if (!exists("support_class", inherits = TRUE)) {
    warning("No `support_class` object found; approximating VFD from `vent_hours`.")
    return(fallback_vfd_from_vent_hours(df))
  }

  support_class <- get("support_class", inherits = TRUE) %>%
    janitor::clean_names()
  rec_time_col <- first_present(support_class, c("rec_time", "recorded_dttm"))
  imv_col <- first_present(support_class, c("is_invasive_vent", "is_imv"))
  if (is.na(rec_time_col) || is.na(imv_col)) {
    warning("`support_class` lacks IMV/time columns; approximating VFD from `vent_hours`.")
    return(fallback_vfd_from_vent_hours(df))
  }

  imv_runs <- support_class %>%
    transmute(
      hospitalization_id,
      rec_time = safe_ts(.data[[rec_time_col]]),
      is_imv = as.logical(.data[[imv_col]])
    ) %>%
    inner_join(df %>% dplyr::select(hospitalization_id, t0), by = "hospitalization_id") %>%
    filter(is_imv, !is.na(rec_time), rec_time >= t0) %>%
    arrange(hospitalization_id, rec_time) %>%
    group_by(hospitalization_id) %>%
    mutate(
      gap_h = as.numeric(difftime(rec_time, lag(rec_time), units = "hours")),
      run_id = cumsum(if_else(is.na(gap_h) | gap_h > imv_gap_hours, 1L, 0L))
    ) %>%
    group_by(hospitalization_id, run_id) %>%
    summarise(imv_start = min(rec_time), imv_end = max(rec_time), t0 = first(t0), .groups = "drop")

  summarise_imv_days <- function(horizon_days, output_col) {
    imv_runs %>%
    mutate(
      vfd_end = t0 + days(horizon_days),
      run_start = pmax(imv_start, t0, na.rm = TRUE),
      run_end = pmin(imv_end + hours(imv_gap_hours), vfd_end, na.rm = TRUE),
      imv_days = pmax(as.numeric(difftime(run_end, run_start, units = "days")), 0)
    ) %>%
    filter(is.finite(imv_days), imv_days > 0) %>%
      group_by(hospitalization_id) %>%
      summarise("{output_col}" := sum(imv_days, na.rm = TRUE), .groups = "drop")
  }

  imv_days <- summarise_imv_days(vfd_days, "imv_days_through_vfd")
  imv_days_14 <- summarise_imv_days(14, "imv_days_through_day14")

  df %>%
    left_join(imv_days, by = "hospitalization_id") %>%
    left_join(imv_days_14, by = "hospitalization_id") %>%
    mutate(
      imv_days_through_vfd = coalesce(imv_days_through_vfd, 0),
      imv_days_through_day14 = coalesce(imv_days_through_day14, 0),
      died_before_vfd_horizon = mortality_event == 1L & mortality_ftime_days <= vfd_days,
      died_before_vfd_day14 = mortality_event == 1L & mortality_ftime_days <= 14,
      ventilator_free_days = if_else(
        died_before_vfd_horizon,
        0,
        pmax(vfd_days - imv_days_through_vfd, 0)
      ),
      ventilator_free_days = pmin(ventilator_free_days, vfd_days),
      ventilator_free_days_day14 = if_else(
        died_before_vfd_day14,
        0,
        pmax(14 - imv_days_through_day14, 0)
      ),
      ventilator_free_days_day14 = pmin(ventilator_free_days_day14, 14)
    )
}

social_covars <- function(df) {
  zcta_covars <- c(
    "acs_pct_poverty",
    "acs_pct_unemployed",
    "acs_pct_no_vehicle",
    "acs_pct_nonwhite",
    "acs_median_household_income_10k",
    "acs_pct_bachelor_plus"
  )
  tract_covars <- c(
    "svi_overall",
    "acs_median_income_10k",
    "acs_pct_lt_hs",
    "acs_unemp_rate_pct",
    "acs_pct_insured"
  )
  if (all(zcta_covars %in% names(df))) zcta_covars else intersect(tract_covars, names(df))
}

ensure_scaled_social_covars <- function(df) {
  df <- if ("acs_median_household_income_10k" %in% names(df)) {
    df
  } else if ("acs_median_household_income" %in% names(df)) {
    mutate(df, acs_median_household_income_10k = acs_median_household_income / 10000)
  } else {
    df
  }

  if ("acs_median_income_10k" %in% names(df)) {
    df
  } else if ("acs_median_income" %in% names(df)) {
    mutate(df, acs_median_income_10k = acs_median_income / 10000)
  } else {
    df
  }
}

drop_uninformative_covars <- function(model_df, covars, exposure_terms) {
  purrr::keep(covars, function(covar) {
    if (covar %in% exposure_terms) return(TRUE)
    x <- model_df[[covar]]
    if (is.factor(x) || is.character(x)) return(dplyr::n_distinct(x, na.rm = TRUE) >= 2)
    isTRUE(stats::var(as.numeric(x), na.rm = TRUE) > 0)
  }) %>% unname()
}

tidy_cox <- function(fit, model_df, exposure_terms, model_label, covars) {
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      site = site_name,
      outcome = "In-hospital mortality after ARF onset",
      model = model_label,
      term,
      n = nrow(model_df),
      events = sum(model_df$mortality_event == 1L),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      adjustment_set = paste(setdiff(covars, exposure_terms), collapse = " + ")
    )
}

tidy_logistic <- function(fit, model_df, exposure_terms, model_label, covars) {
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      site = site_name,
      outcome = "Mortality by day 28 after ARF onset",
      model = model_label,
      term,
      n = nrow(model_df),
      events = sum(model_df$mortality_day28_event == 1L),
      odds_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      adjustment_set = paste(setdiff(covars, exposure_terms), collapse = " + ")
    )
}

tidy_vfd <- function(fit, model_df, exposure_terms, model_label, covars) {
  broom::tidy(fit) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      site = site_name,
      outcome = glue("Ventilator-free days through day {vfd_days}"),
      model = model_label,
      term,
      n = nrow(model_df),
      mean_vfd = mean(model_df$ventilator_free_days, na.rm = TRUE),
      ratio_of_means = exp(estimate),
      conf_low = exp(estimate - 1.96 * std.error),
      conf_high = exp(estimate + 1.96 * std.error),
      p_value = p.value,
      dispersion = summary(fit)$dispersion,
      adjustment_set = paste(setdiff(covars, exposure_terms), collapse = " + ")
    )
}

continuous_outcome_specs <- list(
  vfd = list(
    outcome_col = "ventilator_free_days",
    label = function() glue("Ventilator-free days through day {vfd_days}"),
    mean_col = "mean_vfd",
    prefix = "vfd",
    ceiling = function() vfd_days
  ),
  imv_duration = list(
    outcome_col = "imv_days_through_vfd",
    label = function() glue("IMV duration through day {vfd_days}"),
    mean_col = "mean_imv_days",
    prefix = "imv_days",
    ceiling = function() vfd_days
  )
)

tidy_continuous_quasipoisson <- function(fit, model_df, exposure_terms, model_label, covars, spec) {
  y <- model_df[[spec$outcome_col]]
  out <- broom::tidy(fit) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      site = site_name,
      outcome = spec$label(),
      model = model_label,
      term,
      n = nrow(model_df),
      ratio_of_means = exp(estimate),
      conf_low = exp(estimate - 1.96 * std.error),
      conf_high = exp(estimate + 1.96 * std.error),
      p_value = p.value,
      dispersion = summary(fit)$dispersion,
      adjustment_set = paste(setdiff(covars, exposure_terms), collapse = " + ")
    )
  out[[spec$mean_col]] <- mean(y, na.rm = TRUE)
  out
}

diagnose_continuous_quasipoisson_model <- function(obj, spec) {
  fit <- obj$fit
  model_df <- obj$model_df
  y <- model_df[[spec$outcome_col]]
  fitted_values <- stats::fitted(fit)
  pearson_resid <- stats::residuals(fit, type = "pearson")
  deviance_resid <- stats::residuals(fit, type = "deviance")
  residual_df <- stats::df.residual(fit)
  pearson_chisq <- sum(pearson_resid^2, na.rm = TRUE)
  cooks <- tryCatch(stats::cooks.distance(fit), error = function(e) rep(NA_real_, nrow(model_df)))
  poisson_fit <- tryCatch(
    stats::glm(stats::formula(fit), data = model_df, family = poisson(link = "log")),
    warning = function(w) {
      suppressWarnings(stats::glm(stats::formula(fit), data = model_df, family = poisson(link = "log")))
    },
    error = function(e) NULL
  )

  tibble(
    site = site_name,
    outcome = spec$label(),
    model = obj$model_label,
    n = nrow(model_df),
    outcome_mean = mean(y, na.rm = TRUE),
    outcome_variance = stats::var(y, na.rm = TRUE),
    variance_to_mean_ratio = outcome_variance / outcome_mean,
    outcome_min = min(y, na.rm = TRUE),
    outcome_q25 = as.numeric(stats::quantile(y, 0.25, na.rm = TRUE)),
    outcome_median = stats::median(y, na.rm = TRUE),
    outcome_q75 = as.numeric(stats::quantile(y, 0.75, na.rm = TRUE)),
    outcome_max = max(y, na.rm = TRUE),
    zero_outcome_n = sum(y == 0, na.rm = TRUE),
    zero_outcome_percent = 100 * zero_outcome_n / n,
    ceiling_outcome_n = sum(y >= spec$ceiling(), na.rm = TRUE),
    ceiling_outcome_percent = 100 * ceiling_outcome_n / n,
    noninteger_outcome_n = sum(abs(y - round(y)) > 1e-8, na.rm = TRUE),
    noninteger_outcome_percent = 100 * noninteger_outcome_n / n,
    residual_df = residual_df,
    pearson_chisq = pearson_chisq,
    pearson_dispersion = pearson_chisq / residual_df,
    deviance_dispersion = stats::deviance(fit) / residual_df,
    quasipoisson_dispersion = summary(fit)$dispersion,
    poisson_overdispersion_p = stats::pchisq(pearson_chisq, df = residual_df, lower.tail = FALSE),
    se_inflation_vs_poisson = sqrt(summary(fit)$dispersion),
    poisson_aic = if (!is.null(poisson_fit)) stats::AIC(poisson_fit) else NA_real_,
    fitted_min = min(fitted_values, na.rm = TRUE),
    fitted_q25 = as.numeric(stats::quantile(fitted_values, 0.25, na.rm = TRUE)),
    fitted_median = stats::median(fitted_values, na.rm = TRUE),
    fitted_q75 = as.numeric(stats::quantile(fitted_values, 0.75, na.rm = TRUE)),
    fitted_max = max(fitted_values, na.rm = TRUE),
    pearson_resid_q05 = as.numeric(stats::quantile(pearson_resid, 0.05, na.rm = TRUE)),
    pearson_resid_median = stats::median(pearson_resid, na.rm = TRUE),
    pearson_resid_q95 = as.numeric(stats::quantile(pearson_resid, 0.95, na.rm = TRUE)),
    deviance_resid_q05 = as.numeric(stats::quantile(deviance_resid, 0.05, na.rm = TRUE)),
    deviance_resid_median = stats::median(deviance_resid, na.rm = TRUE),
    deviance_resid_q95 = as.numeric(stats::quantile(deviance_resid, 0.95, na.rm = TRUE)),
    max_cooks_distance = max(cooks, na.rm = TRUE),
    influential_cooks_n = sum(cooks > (4 / nrow(model_df)), na.rm = TRUE),
    influential_cooks_percent = 100 * influential_cooks_n / nrow(model_df),
    adjustment_set = paste(setdiff(obj$covars, obj$exposure_terms), collapse = " + ")
  )
}

diagnose_vfd_model <- function(obj) {
  fit <- obj$fit
  model_df <- obj$model_df
  y <- model_df$ventilator_free_days
  fitted_values <- stats::fitted(fit)
  pearson_resid <- stats::residuals(fit, type = "pearson")
  deviance_resid <- stats::residuals(fit, type = "deviance")
  residual_df <- stats::df.residual(fit)
  pearson_chisq <- sum(pearson_resid^2, na.rm = TRUE)
  cooks <- tryCatch(stats::cooks.distance(fit), error = function(e) rep(NA_real_, nrow(model_df)))
  poisson_fit <- tryCatch(
    stats::glm(stats::formula(fit), data = model_df, family = poisson(link = "log")),
    warning = function(w) {
      suppressWarnings(stats::glm(stats::formula(fit), data = model_df, family = poisson(link = "log")))
    },
    error = function(e) NULL
  )

  tibble(
    site = site_name,
    outcome = glue("Ventilator-free days through day {vfd_days}"),
    model = obj$model_label,
    n = nrow(model_df),
    mean_vfd = mean(y, na.rm = TRUE),
    variance_vfd = stats::var(y, na.rm = TRUE),
    variance_to_mean_ratio = variance_vfd / mean_vfd,
    min_vfd = min(y, na.rm = TRUE),
    q25_vfd = as.numeric(stats::quantile(y, 0.25, na.rm = TRUE)),
    median_vfd = stats::median(y, na.rm = TRUE),
    q75_vfd = as.numeric(stats::quantile(y, 0.75, na.rm = TRUE)),
    max_vfd = max(y, na.rm = TRUE),
    zero_vfd_n = sum(y == 0, na.rm = TRUE),
    zero_vfd_percent = 100 * zero_vfd_n / n,
    ceiling_vfd_n = sum(y >= vfd_days, na.rm = TRUE),
    ceiling_vfd_percent = 100 * ceiling_vfd_n / n,
    noninteger_vfd_n = sum(abs(y - round(y)) > 1e-8, na.rm = TRUE),
    noninteger_vfd_percent = 100 * noninteger_vfd_n / n,
    residual_df = residual_df,
    pearson_chisq = pearson_chisq,
    pearson_dispersion = pearson_chisq / residual_df,
    deviance_dispersion = stats::deviance(fit) / residual_df,
    quasipoisson_dispersion = summary(fit)$dispersion,
    poisson_overdispersion_p = stats::pchisq(pearson_chisq, df = residual_df, lower.tail = FALSE),
    se_inflation_vs_poisson = sqrt(summary(fit)$dispersion),
    poisson_aic = if (!is.null(poisson_fit)) stats::AIC(poisson_fit) else NA_real_,
    fitted_min = min(fitted_values, na.rm = TRUE),
    fitted_q25 = as.numeric(stats::quantile(fitted_values, 0.25, na.rm = TRUE)),
    fitted_median = stats::median(fitted_values, na.rm = TRUE),
    fitted_q75 = as.numeric(stats::quantile(fitted_values, 0.75, na.rm = TRUE)),
    fitted_max = max(fitted_values, na.rm = TRUE),
    pearson_resid_q05 = as.numeric(stats::quantile(pearson_resid, 0.05, na.rm = TRUE)),
    pearson_resid_median = stats::median(pearson_resid, na.rm = TRUE),
    pearson_resid_q95 = as.numeric(stats::quantile(pearson_resid, 0.95, na.rm = TRUE)),
    deviance_resid_q05 = as.numeric(stats::quantile(deviance_resid, 0.05, na.rm = TRUE)),
    deviance_resid_median = stats::median(deviance_resid, na.rm = TRUE),
    deviance_resid_q95 = as.numeric(stats::quantile(deviance_resid, 0.95, na.rm = TRUE)),
    max_cooks_distance = max(cooks, na.rm = TRUE),
    influential_cooks_n = sum(cooks > (4 / nrow(model_df)), na.rm = TRUE),
    influential_cooks_percent = 100 * influential_cooks_n / nrow(model_df),
    adjustment_set = paste(setdiff(obj$covars, obj$exposure_terms), collapse = " + ")
  )
}

compare_vfd_poisson_quasi <- function(obj) {
  poisson_fit <- tryCatch(
    stats::glm(stats::formula(obj$fit), data = obj$model_df, family = poisson(link = "log")),
    warning = function(w) {
      suppressWarnings(stats::glm(stats::formula(obj$fit), data = obj$model_df, family = poisson(link = "log")))
    },
    error = function(e) NULL
  )
  if (is.null(poisson_fit)) return(tibble())

  quasi_terms <- broom::tidy(obj$fit) %>%
    filter(term %in% obj$exposure_terms) %>%
    transmute(
      term,
      quasi_beta = estimate,
      quasi_std_error = std.error,
      quasi_p_value = p.value,
      quasi_ratio_of_means = exp(estimate),
      quasi_conf_low = exp(estimate - 1.96 * std.error),
      quasi_conf_high = exp(estimate + 1.96 * std.error)
    )

  poisson_terms <- broom::tidy(poisson_fit) %>%
    filter(term %in% obj$exposure_terms) %>%
    transmute(
      term,
      poisson_beta = estimate,
      poisson_std_error = std.error,
      poisson_p_value = p.value,
      poisson_ratio_of_means = exp(estimate),
      poisson_conf_low = exp(estimate - 1.96 * std.error),
      poisson_conf_high = exp(estimate + 1.96 * std.error)
    )

  quasi_terms %>%
    left_join(poisson_terms, by = "term") %>%
    mutate(
      site = site_name,
      outcome = glue("Ventilator-free days through day {vfd_days}"),
      model = obj$model_label,
      n = nrow(obj$model_df),
      quasipoisson_dispersion = summary(obj$fit)$dispersion,
      quasi_to_poisson_se_ratio = quasi_std_error / poisson_std_error,
      adjustment_set = paste(setdiff(obj$covars, obj$exposure_terms), collapse = " + "),
      .before = 1
  )
}

compare_continuous_poisson_quasi <- function(obj, spec) {
  poisson_fit <- tryCatch(
    stats::glm(stats::formula(obj$fit), data = obj$model_df, family = poisson(link = "log")),
    warning = function(w) {
      suppressWarnings(stats::glm(stats::formula(obj$fit), data = obj$model_df, family = poisson(link = "log")))
    },
    error = function(e) NULL
  )
  if (is.null(poisson_fit)) return(tibble())

  quasi_terms <- broom::tidy(obj$fit) %>%
    filter(term %in% obj$exposure_terms) %>%
    transmute(
      term,
      quasi_beta = estimate,
      quasi_std_error = std.error,
      quasi_p_value = p.value,
      quasi_ratio_of_means = exp(estimate),
      quasi_conf_low = exp(estimate - 1.96 * std.error),
      quasi_conf_high = exp(estimate + 1.96 * std.error)
    )

  poisson_terms <- broom::tidy(poisson_fit) %>%
    filter(term %in% obj$exposure_terms) %>%
    transmute(
      term,
      poisson_beta = estimate,
      poisson_std_error = std.error,
      poisson_p_value = p.value,
      poisson_ratio_of_means = exp(estimate),
      poisson_conf_low = exp(estimate - 1.96 * std.error),
      poisson_conf_high = exp(estimate + 1.96 * std.error)
    )

  quasi_terms %>%
    left_join(poisson_terms, by = "term") %>%
    mutate(
      site = site_name,
      outcome = spec$label(),
      model = obj$model_label,
      n = nrow(obj$model_df),
      quasipoisson_dispersion = summary(obj$fit)$dispersion,
      quasi_to_poisson_se_ratio = quasi_std_error / poisson_std_error,
      adjustment_set = paste(setdiff(obj$covars, obj$exposure_terms), collapse = " + "),
      .before = 1
    )
}

calibrate_vfd_model <- function(obj, n_bins = 10) {
  obj$model_df %>%
    mutate(
      fitted_vfd = stats::fitted(obj$fit),
      fitted_decile = dplyr::ntile(fitted_vfd, n_bins)
    ) %>%
    group_by(fitted_decile) %>%
    summarise(
      n = n(),
      fitted_mean = mean(fitted_vfd, na.rm = TRUE),
      observed_mean = mean(ventilator_free_days, na.rm = TRUE),
      observed_sd = stats::sd(ventilator_free_days, na.rm = TRUE),
      observed_median = stats::median(ventilator_free_days, na.rm = TRUE),
      observed_zero_percent = 100 * mean(ventilator_free_days == 0, na.rm = TRUE),
      observed_ceiling_percent = 100 * mean(ventilator_free_days >= vfd_days, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      site = site_name,
      outcome = glue("Ventilator-free days through day {vfd_days}"),
      model = obj$model_label,
      .before = 1
  )
}

calibrate_continuous_model <- function(obj, spec, n_bins = 10) {
  outcome_col <- spec$outcome_col
  obj$model_df %>%
    mutate(
      fitted_outcome = stats::fitted(obj$fit),
      fitted_decile = dplyr::ntile(fitted_outcome, n_bins)
    ) %>%
    group_by(fitted_decile) %>%
    summarise(
      n = n(),
      fitted_mean = mean(fitted_outcome, na.rm = TRUE),
      observed_mean = mean(.data[[outcome_col]], na.rm = TRUE),
      observed_sd = stats::sd(.data[[outcome_col]], na.rm = TRUE),
      observed_median = stats::median(.data[[outcome_col]], na.rm = TRUE),
      observed_zero_percent = 100 * mean(.data[[outcome_col]] == 0, na.rm = TRUE),
      observed_ceiling_percent = 100 * mean(.data[[outcome_col]] >= spec$ceiling(), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      site = site_name,
      outcome = spec$label(),
      model = obj$model_label,
      .before = 1
    )
}

summarise_binary_prevalence <- function(df, variables, labels) {
  purrr::map2_dfr(variables, labels, function(variable, label) {
    if (!variable %in% names(df)) return(tibble())
    x <- df[[variable]]
    x <- if (is.logical(x)) {
      x
    } else if (is.numeric(x) || is.integer(x)) {
      x == 1
    } else {
      str_to_lower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
    }

    denom <- sum(!is.na(x))
    count <- sum(x, na.rm = TRUE)
    tibble(
      site = site_name,
      table = "Table 1",
      domain = "Baseline chronic cardiopulmonary and comorbidity burden",
      variable = variable,
      label = label,
      n = count,
      denominator = denom,
      prevalence = count / denom,
      prevalence_percent = 100 * count / denom,
      display = sprintf("%s (%.1f%%)", format(count, big.mark = ","), 100 * count / denom)
    )
  })
}

prepare_table1_frame <- function(df) {
  df %>%
    mutate(
      sex_category = if ("sex_category" %in% names(.)) sex_category else sex,
      race_ethnicity_simple = if ("race_ethnicity_simple" %in% names(.)) race_ethnicity_simple else race_ethnicity,
      pm25_mean = if ("pm25_mean" %in% names(.)) pm25_mean else pm25_exposure,
      no2_mean = if ("no2_mean" %in% names(.)) no2_mean else no2_exposure,
      svi_overall = if ("svi_overall" %in% names(.)) svi_overall else acs_zcta_svi_proxy,
      icu_los_days = if ("icu_los_days" %in% names(.)) icu_los_days else icu_los_hours / 24,
      hosp_los_days = if ("hosp_los_days" %in% names(.)) {
        hosp_los_days
      } else {
        as.numeric(difftime(discharge_time, t0, units = "days"))
      },
      imv_days_through_vfd = as.numeric(imv_days_through_vfd),
      in_hosp_death = if ("in_hosp_death" %in% names(.)) in_hosp_death else mortality_event,
      death_28d = if ("death_28d" %in% names(.)) {
        death_28d
      } else {
        mortality_day28_event
      }
    ) %>%
    mutate(
      sex_category = forcats::fct_na_value_to_level(factor(sex_category), level = "(Missing)"),
      race_ethnicity_simple = forcats::fct_na_value_to_level(factor(race_ethnicity_simple), level = "(Missing)")
    )
}

summarise_table1_continuous <- function(df, variables, labels) {
  purrr::map_dfr(variables, function(variable) {
    if (!variable %in% names(df)) return(tibble())
    x <- suppressWarnings(as.numeric(df[[variable]]))
    x <- x[is.finite(x)]
    tibble(
      site = site_name,
      table = "Table 1",
      variable = variable,
      label = unname(labels[[.env$variable]]),
      n = length(x),
      mean = if (length(x)) mean(x, na.rm = TRUE) else NA_real_,
      sd = if (length(x)) stats::sd(x, na.rm = TRUE) else NA_real_,
      q25 = if (length(x)) as.numeric(stats::quantile(x, 0.25, na.rm = TRUE)) else NA_real_,
      median = if (length(x)) as.numeric(stats::quantile(x, 0.50, na.rm = TRUE)) else NA_real_,
      q75 = if (length(x)) as.numeric(stats::quantile(x, 0.75, na.rm = TRUE)) else NA_real_,
      display = if (length(x)) {
        sprintf("%.1f +/- %.1f (%.1f; %.1f, %.1f)", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE), stats::median(x, na.rm = TRUE), stats::quantile(x, 0.25, na.rm = TRUE), stats::quantile(x, 0.75, na.rm = TRUE))
      } else {
        NA_character_
      }
    )
  })
}

summarise_table1_categorical <- function(df, variables, labels) {
  purrr::map_dfr(variables, function(variable) {
    if (!variable %in% names(df)) return(tibble())
    x <- df[[variable]]
    x_fac <- if (is.factor(x)) x else factor(x)
    cnt <- as.data.frame(table(x_fac, useNA = "ifany"), stringsAsFactors = FALSE)
    names(cnt) <- c("level", "n")
    cnt %>%
      mutate(
        site = site_name,
        table = "Table 1",
        variable = variable,
        label = unname(labels[[.env$variable]]),
        total_n = sum(n),
        pct = if_else(total_n > 0, 100 * n / total_n, NA_real_),
        display = sprintf("%s (%.1f%%)", format(n, big.mark = ","), pct)
      ) %>%
      dplyr::select(site, table, variable, label, level, n, pct, total_n, display)
  })
}

fit_cox_model <- function(df, exposure_terms, model_label, adjustment_covars) {
  covars <- c(exposure_terms, adjustment_covars)
  model_df <- df %>%
    dplyr::select(mortality_ftime_days, mortality_event, all_of(covars)) %>%
    drop_na()
  covars <- drop_uninformative_covars(model_df, covars, exposure_terms)
  model_df <- model_df %>%
    dplyr::select(mortality_ftime_days, mortality_event, all_of(covars))

  fit <- survival::coxph(
    as.formula(paste0("Surv(mortality_ftime_days, mortality_event) ~ ", paste(covars, collapse = " + "))),
    data = model_df,
    ties = "efron",
    x = TRUE
  )

  list(fit = fit, model_df = model_df, exposure_terms = exposure_terms, model_label = model_label, covars = covars)
}

fit_logistic_model <- function(df, exposure_terms, model_label, adjustment_covars) {
  covars <- c(exposure_terms, adjustment_covars)
  model_df <- df %>%
    dplyr::select(mortality_day28_event, all_of(covars)) %>%
    drop_na()
  covars <- drop_uninformative_covars(model_df, covars, exposure_terms)
  model_df <- model_df %>%
    dplyr::select(mortality_day28_event, all_of(covars))

  fit <- stats::glm(
    as.formula(paste0("mortality_day28_event ~ ", paste(covars, collapse = " + "))),
    data = model_df,
    family = binomial(link = "logit")
  )

  list(fit = fit, model_df = model_df, exposure_terms = exposure_terms, model_label = model_label, covars = covars)
}

fit_vfd_model <- function(df, exposure_terms, model_label, adjustment_covars) {
  covars <- c(exposure_terms, adjustment_covars)
  model_df <- df %>%
    dplyr::select(ventilator_free_days, all_of(covars)) %>%
    drop_na()
  covars <- drop_uninformative_covars(model_df, covars, exposure_terms)
  model_df <- model_df %>%
    dplyr::select(ventilator_free_days, all_of(covars))

  fit <- stats::glm(
    as.formula(paste0("ventilator_free_days ~ ", paste(covars, collapse = " + "))),
    data = model_df,
    family = quasipoisson(link = "log")
  )

  list(fit = fit, model_df = model_df, exposure_terms = exposure_terms, model_label = model_label, covars = covars)
}

fit_continuous_quasipoisson_model <- function(df, exposure_terms, model_label, adjustment_covars, outcome_col) {
  covars <- c(exposure_terms, adjustment_covars)
  model_df <- df %>%
    dplyr::select(all_of(outcome_col), all_of(covars)) %>%
    drop_na()
  covars <- drop_uninformative_covars(model_df, covars, exposure_terms)
  model_df <- model_df %>%
    dplyr::select(all_of(outcome_col), all_of(covars))

  fit <- stats::glm(
    as.formula(paste0(outcome_col, " ~ ", paste(covars, collapse = " + "))),
    data = model_df,
    family = quasipoisson(link = "log")
  )

  list(fit = fit, model_df = model_df, exposure_terms = exposure_terms, model_label = model_label, covars = covars)
}

fit_primary_model_set <- function(df) {
  logistic_fits <- purrr::map(
    model_specs,
    ~ fit_logistic_model(df, .x$exposure_terms, .x$model_label, adjustment_covars)
  )

  cox_fits <- purrr::map(
    model_specs,
    ~ fit_cox_model(df, .x$exposure_terms, .x$model_label, adjustment_covars)
  )

  vfd_fits <- purrr::map(
    model_specs,
    ~ fit_vfd_model(df, .x$exposure_terms, .x$model_label, adjustment_covars)
  )

  imv_duration_fits <- purrr::map(
    model_specs,
    ~ fit_continuous_quasipoisson_model(
      df,
      .x$exposure_terms,
      .x$model_label,
      adjustment_covars,
      continuous_outcome_specs$imv_duration$outcome_col
    )
  )

  list(
    logistic_fits = logistic_fits,
    cox_fits = cox_fits,
    vfd_fits = vfd_fits,
    imv_duration_fits = imv_duration_fits,
    logistic_results = purrr::map_dfr(logistic_fits, ~ tidy_logistic(.x$fit, .x$model_df, .x$exposure_terms, .x$model_label, .x$covars)),
    cox_results = purrr::map_dfr(cox_fits, ~ tidy_cox(.x$fit, .x$model_df, .x$exposure_terms, .x$model_label, .x$covars)),
    vfd_results = purrr::map_dfr(vfd_fits, ~ tidy_vfd(.x$fit, .x$model_df, .x$exposure_terms, .x$model_label, .x$covars)),
    imv_duration_results = purrr::map_dfr(
      imv_duration_fits,
      ~ tidy_continuous_quasipoisson(.x$fit, .x$model_df, .x$exposure_terms, .x$model_label, .x$covars, continuous_outcome_specs$imv_duration)
    )
  )
}

summarise_modeling_cohort <- function(df, sensitivity_label, primary_n = NA_integer_) {
  both_exposure_complete <- stats::complete.cases(df[, intersect(c("pm25_per_5", "no2_per_10"), names(df)), drop = FALSE])
  adjustment_complete <- stats::complete.cases(df[, intersect(adjustment_covars, names(df)), drop = FALSE])
  primary_complete_vars <- c(
    "mortality_day28_event",
    "ventilator_free_days",
    "imv_days_through_vfd",
    "pm25_per_5",
    "no2_per_10",
    adjustment_covars
  )
  primary_complete <- stats::complete.cases(df[, intersect(primary_complete_vars, names(df)), drop = FALSE])
  competing_complete_vars <- c("ftime_days", "event_code", "pm25_per_5", "no2_per_10", adjustment_covars)
  competing_pool <- df %>% filter(.data$has_imv_after_arf %in% TRUE)
  competing_complete <- if (nrow(competing_pool)) {
    stats::complete.cases(competing_pool[, intersect(competing_complete_vars, names(competing_pool)), drop = FALSE])
  } else {
    logical()
  }

  tibble(
    site = site_name,
    sensitivity = sensitivity_label,
    n_arf = nrow(df),
    n_patients = n_distinct(df$patient_id),
    n_excluded_vs_primary = if (is.na(primary_n)) NA_integer_ else primary_n - nrow(df),
    mortality_events = sum(df$mortality_event == 1L, na.rm = TRUE),
    mortality_day28_events = sum(df$mortality_day28_event == 1L, na.rm = TRUE),
    alive_or_censored_before_day28 = sum(df$alive_or_censored_before_day28 == 1L, na.rm = TRUE),
    median_mortality_followup_days = median(df$mortality_ftime_days, na.rm = TRUE),
    mean_vfd = mean(df$ventilator_free_days, na.rm = TRUE),
    median_vfd = median(df$ventilator_free_days, na.rm = TRUE),
    mean_imv_days = mean(df$imv_days_through_vfd, na.rm = TRUE),
    median_imv_days = median(df$imv_days_through_vfd, na.rm = TRUE),
    vfd_horizon_days = vfd_days,
    mean_charlson = mean(df$charlson_score, na.rm = TRUE),
    median_charlson = median(df$charlson_score, na.rm = TRUE),
    n_with_pm25_no2 = sum(both_exposure_complete, na.rm = TRUE),
    n_with_primary_adjustment_covariates = sum(adjustment_complete, na.rm = TRUE),
    n_primary_complete_cases = sum(primary_complete, na.rm = TRUE),
    n_primary_complete_case_patients = n_distinct(df$patient_id[primary_complete]),
    n_imv_competing_risk_complete_cases = sum(competing_complete, na.rm = TRUE),
    n_imv_competing_risk_complete_case_patients = if (nrow(competing_pool)) {
      n_distinct(competing_pool$patient_id[competing_complete])
    } else {
      0L
    },
    pm25_exposure_column = first_present(arf_exp, c("pm25_12m_zcta", "pm25_12m_mean", "pm25_mean")),
    no2_exposure_column = first_present(arf_exp, c("no2_12m_zcta", "no2_12m_mean", "no2_mean")),
    adjustment_covariates = paste(adjustment_covars, collapse = " + ")
  )
}

make_primary_pooling_table <- function(logistic_tbl, cox_tbl, vfd_tbl, imv_duration_tbl, sensitivity_label,
                                       covid_exclude_start = as.Date(NA),
                                       covid_exclude_end = as.Date(NA)) {
  bind_rows(
    logistic_tbl %>%
      transmute(
        site,
        sensitivity = sensitivity_label,
        covid_exclude_start,
        covid_exclude_end,
        outcome,
        analysis_model = "Logistic regression",
        pollutant_model = model,
        term,
        effect_measure = "odds_ratio",
        n,
        events,
        mean_vfd = NA_real_,
        mean_imv_days = NA_real_,
        estimate = odds_ratio,
        conf_low,
        conf_high,
        p_value,
        dispersion = NA_real_,
        adjustment_set
      ),
    cox_tbl %>%
      transmute(
        site,
        sensitivity = sensitivity_label,
        covid_exclude_start,
        covid_exclude_end,
        outcome,
        analysis_model = "Cox proportional hazards",
        pollutant_model = model,
        term,
        effect_measure = "hazard_ratio",
        n,
        events,
        mean_vfd = NA_real_,
        mean_imv_days = NA_real_,
        estimate = hazard_ratio,
        conf_low,
        conf_high,
        p_value,
        dispersion = NA_real_,
        adjustment_set
      ),
    vfd_tbl %>%
      transmute(
        site,
        sensitivity = sensitivity_label,
        covid_exclude_start,
        covid_exclude_end,
        outcome,
        analysis_model = "Quasi-Poisson regression",
        pollutant_model = model,
        term,
        effect_measure = "ratio_of_means",
        n,
        events = NA_integer_,
        mean_vfd,
        mean_imv_days = NA_real_,
        estimate = ratio_of_means,
        conf_low,
        conf_high,
        p_value,
        dispersion,
        adjustment_set
      ),
    imv_duration_tbl %>%
      transmute(
        site,
        sensitivity = sensitivity_label,
        covid_exclude_start,
        covid_exclude_end,
        outcome,
        analysis_model = "Quasi-Poisson regression",
        pollutant_model = model,
        term,
        effect_measure = "ratio_of_means",
        n,
        events = NA_integer_,
        mean_vfd = NA_real_,
        mean_imv_days,
        estimate = ratio_of_means,
        conf_low,
        conf_high,
        p_value,
        dispersion,
        adjustment_set
      )
  )
}

analysis_df <- arf_exp %>%
  coerce_analysis_frame() %>%
  add_mortality_time() %>%
  add_day28_mortality() %>%
  ensure_scaled_social_covars() %>%
  add_vfd()

adjustment_covars <- c(
  "age_10",
  "sex",
  "race_ethnicity",
  "charlson_score",
  "index_year_f",
  social_covars(analysis_df)
)

model_specs <- list(
  list(exposure_terms = "pm25_per_5", model_label = "PM25 single-pollutant"),
  list(exposure_terms = "no2_per_10", model_label = "NO2 single-pollutant"),
  list(exposure_terms = c("pm25_per_5", "no2_per_10"), model_label = "PM25 + NO2")
)

primary_model_set <- fit_primary_model_set(analysis_df)
logistic_fits <- primary_model_set$logistic_fits
cox_fits <- primary_model_set$cox_fits
vfd_fits <- primary_model_set$vfd_fits
imv_duration_fits <- primary_model_set$imv_duration_fits
logistic_results <- primary_model_set$logistic_results
cox_results <- primary_model_set$cox_results
vfd_results <- primary_model_set$vfd_results
imv_duration_results <- primary_model_set$imv_duration_results
vfd_model_diagnostics <- purrr::map_dfr(vfd_fits, diagnose_vfd_model)
vfd_poisson_quasi_comparison <- purrr::map_dfr(vfd_fits, compare_vfd_poisson_quasi)
vfd_calibration_by_decile <- purrr::map_dfr(vfd_fits, calibrate_vfd_model)
imv_duration_model_diagnostics <- purrr::map_dfr(
  imv_duration_fits,
  ~ diagnose_continuous_quasipoisson_model(.x, continuous_outcome_specs$imv_duration)
)
imv_duration_poisson_quasi_comparison <- purrr::map_dfr(
  imv_duration_fits,
  ~ compare_continuous_poisson_quasi(.x, continuous_outcome_specs$imv_duration)
)
imv_duration_calibration_by_decile <- purrr::map_dfr(
  imv_duration_fits,
  ~ calibrate_continuous_model(.x, continuous_outcome_specs$imv_duration)
)

covid_exclude_start <- as.Date(Sys.getenv("REFER_COVID_EXCLUDE_START", "2020-03-01"))
covid_exclude_end <- as.Date(Sys.getenv("REFER_COVID_EXCLUDE_END", "2021-02-28"))
analysis_df_no_peak_covid <- analysis_df %>%
  mutate(covid_sensitivity_index_date = as.Date(t0)) %>%
  filter(
    !is.na(covid_sensitivity_index_date),
    covid_sensitivity_index_date < covid_exclude_start |
      covid_sensitivity_index_date > covid_exclude_end
  )

no_peak_covid_model_set <- fit_primary_model_set(analysis_df_no_peak_covid)
no_peak_covid_logistic_results <- no_peak_covid_model_set$logistic_results
no_peak_covid_cox_results <- no_peak_covid_model_set$cox_results
no_peak_covid_vfd_results <- no_peak_covid_model_set$vfd_results
no_peak_covid_imv_duration_results <- no_peak_covid_model_set$imv_duration_results

no_peak_covid_cohort_summary <- summarise_modeling_cohort(
  analysis_df_no_peak_covid,
  sensitivity_label = "exclude_peak_covid_12m",
  primary_n = nrow(analysis_df)
) %>%
  mutate(
    covid_exclude_start = covid_exclude_start,
    covid_exclude_end = covid_exclude_end,
    covid_index_date = "ARF onset",
    .after = sensitivity
  )

primary_pooling_table <- make_primary_pooling_table(
  logistic_results,
  cox_results,
  vfd_results,
  imv_duration_results,
  sensitivity_label = "primary"
)

no_peak_covid_pooling_table <- make_primary_pooling_table(
  no_peak_covid_logistic_results,
  no_peak_covid_cox_results,
  no_peak_covid_vfd_results,
  no_peak_covid_imv_duration_results,
  sensitivity_label = "exclude_peak_covid_12m",
  covid_exclude_start = covid_exclude_start,
  covid_exclude_end = covid_exclude_end
)

primary_vs_no_peak_covid_pooling_table <- bind_rows(
  primary_pooling_table,
  no_peak_covid_pooling_table
)

ph_results <- purrr::map_dfr(cox_fits, function(obj) {
  zph <- tryCatch(survival::cox.zph(obj$fit, terms = TRUE), error = function(e) NULL)
  if (is.null(zph)) return(tibble())
  as_tibble(zph$table, rownames = "term") %>%
    transmute(
      site = site_name,
      outcome = "In-hospital mortality after ARF onset",
      model = obj$model_label,
      term,
      chisq,
      p_value = p,
      term_type = case_when(
        term %in% obj$exposure_terms ~ "exposure",
        term == "GLOBAL" ~ "global",
        TRUE ~ "covariate"
      )
    )
})

baseline_burden_vars <- c(
  "chronic_pulmonary_disease",
  "congestive_heart_failure",
  "myocardial_infarction",
  "peripheral_vascular_disease",
  "cerebrovascular_disease",
  "renal_disease",
  "diabetes_without_complication",
  "diabetes_with_complication",
  "any_malignancy",
  "metastatic_solid_tumor"
)

baseline_burden_labels <- c(
  "Chronic pulmonary disease",
  "Congestive heart failure",
  "Myocardial infarction",
  "Peripheral vascular disease",
  "Cerebrovascular disease",
  "Renal disease",
  "Diabetes without chronic complication",
  "Diabetes with chronic complication",
  "Any malignancy",
  "Metastatic solid tumor"
)

baseline_burden_prevalence <- summarise_binary_prevalence(
  analysis_df,
  baseline_burden_vars,
  baseline_burden_labels
)

table1_df <- prepare_table1_frame(analysis_df)

table1_continuous_vars <- c(
  "age",
  "pm25_mean",
  "no2_mean",
  "svi_overall",
  "acs_pct_poverty",
  "acs_pct_unemployed",
  "acs_pct_no_vehicle",
  "acs_pct_nonwhite",
  "acs_pct_black",
  "acs_pct_asian",
  "acs_pct_hispanic",
  "acs_median_household_income",
  "acs_pct_bachelor_plus",
  "icu_los_days",
  "hosp_los_days",
  "vital_coverage_days",
  "ventilator_free_days",
  "sofa_total",
  "prior_icu_stays"
)

table1_continuous_labels <- c(
  age = "Age (years)",
  pm25_mean = "PM2.5 (annual mean, ug/m3)",
  no2_mean = "NO2 (annual mean, ppb)",
  svi_overall = "SVI (overall)",
  acs_pct_poverty = "ZCTA poverty (%)",
  acs_pct_unemployed = "ZCTA unemployment (%)",
  acs_pct_no_vehicle = "ZCTA households without vehicle (%)",
  acs_pct_nonwhite = "ZCTA non-White population (%)",
  acs_pct_black = "ZCTA Black population (%)",
  acs_pct_asian = "ZCTA Asian population (%)",
  acs_pct_hispanic = "ZCTA Hispanic/Latino population (%)",
  acs_median_household_income = "ZCTA median household income",
  acs_pct_bachelor_plus = "ZCTA bachelor's degree or higher (%)",
  icu_los_days = "ICU length of stay (days)",
  hosp_los_days = "Hospital length of stay (days)",
  vital_coverage_days = "Vital-sign coverage span (days)",
  ventilator_free_days = "Ventilator-free days through day 28",
  sofa_total = "SOFA total score",
  prior_icu_stays = "Prior ICU stays (count)"
)

table1_categorical_vars <- c(
  "sex_category",
  "race_ethnicity_simple",
  "in_hosp_death",
  "death_28d",
  "vaso_flag",
  "aki_flag",
  "any_acei_arb",
  "any_diuretic",
  "any_bb",
  baseline_burden_vars
)

table1_categorical_labels <- c(
  sex_category = "Sex",
  race_ethnicity_simple = "Race/Ethnicity",
  in_hosp_death = "In-hospital death",
  death_28d = "Death by day 28",
  vaso_flag = "Vasopressor use",
  aki_flag = "AKI",
  any_acei_arb = "ACEi/ARB during lookback",
  any_diuretic = "Diuretic during lookback",
  any_bb = "Beta-blocker during lookback",
  stats::setNames(baseline_burden_labels, baseline_burden_vars)
)

table1_continuous_stats <- summarise_table1_continuous(
  table1_df,
  table1_continuous_vars,
  table1_continuous_labels
)

table1_categorical_stats <- summarise_table1_categorical(
  table1_df,
  table1_categorical_vars,
  table1_categorical_labels
)

table1_missing_variables <- tibble(
  site = site_name,
  table = "Table 1",
  variable = c(table1_continuous_vars, table1_categorical_vars),
  variable_type = c(
    rep("continuous", length(table1_continuous_vars)),
    rep("categorical", length(table1_categorical_vars))
  ),
  present = variable %in% names(table1_df)
) %>%
  filter(!present)

table1_long <- bind_rows(
  table1_continuous_stats %>%
    transmute(
      site,
      table,
      variable_type = "continuous",
      variable,
      label,
      level = NA_character_,
      n,
      denominator = n,
      pct = NA_real_,
      mean,
      sd,
      q25,
      median,
      q75,
      display
    ),
  table1_categorical_stats %>%
    transmute(
      site,
      table,
      variable_type = "categorical",
      variable,
      label,
      level = as.character(level),
      n,
      denominator = total_n,
      pct,
      mean = NA_real_,
      sd = NA_real_,
      q25 = NA_real_,
      median = NA_real_,
      q75 = NA_real_,
      display
    )
)

baseline_burden_wide <- baseline_burden_prevalence %>%
  transmute(
    name = paste0(variable, "_prevalence_percent"),
    value = prevalence_percent
  ) %>%
  tidyr::pivot_wider(names_from = name, values_from = value)

cohort_summary <- summarise_modeling_cohort(
  analysis_df,
  sensitivity_label = "primary",
  primary_n = NA_integer_
) %>%
  dplyr::select(-sensitivity, -n_excluded_vs_primary) %>%
  bind_cols(baseline_burden_wide)

internal_analysis_dataset <- arg_or(
  3,
  Sys.getenv(
    "REFER_INTERNAL_ANALYSIS_DATASET",
    unset = file.path(out_dir, "analysis_dataset_reviewer_optimized.csv")
  )
)
export_row_level_datasets <- tolower(Sys.getenv("REFER_EXPORT_ROW_LEVEL_DATASETS", unset = "false")) %in%
  c("1", "true", "yes", "y")

if (nzchar(internal_analysis_dataset)) {
  dir.create(dirname(internal_analysis_dataset), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(analysis_df, internal_analysis_dataset)
}

if (export_row_level_datasets && !identical(normalizePath(internal_analysis_dataset, mustWork = FALSE), normalizePath(file.path(out_dir, "analysis_dataset_reviewer_optimized.csv"), mustWork = FALSE))) {
  readr::write_csv(analysis_df, file.path(out_dir, "analysis_dataset_reviewer_optimized.csv"))
}
readr::write_csv(cohort_summary, file.path(out_dir, "cohort_summary_reviewer_optimized.csv"))
readr::write_csv(baseline_burden_prevalence, file.path(out_dir, "table1_baseline_chronic_disease_prevalence.csv"))
readr::write_csv(table1_continuous_stats, file.path(out_dir, "table1_continuous_site_resubmission.csv"))
readr::write_csv(table1_categorical_stats, file.path(out_dir, "table1_categorical_site_resubmission.csv"))
readr::write_csv(table1_long, file.path(out_dir, "table1_resubmission_long.csv"))
readr::write_csv(table1_missing_variables, file.path(out_dir, "table1_resubmission_missing_variables.csv"))
readr::write_csv(logistic_results, file.path(out_dir, "primary_mortality_day28_logistic_results.csv"))
readr::write_csv(cox_results, file.path(out_dir, "primary_mortality_cox_results.csv"))
readr::write_csv(ph_results, file.path(out_dir, "primary_mortality_cox_ph_diagnostics.csv"))
readr::write_csv(vfd_results, file.path(out_dir, "primary_vfd_quasipoisson_results.csv"))
readr::write_csv(vfd_model_diagnostics, file.path(out_dir, "primary_vfd_model_diagnostics.csv"))
readr::write_csv(vfd_poisson_quasi_comparison, file.path(out_dir, "primary_vfd_poisson_vs_quasipoisson_diagnostics.csv"))
readr::write_csv(vfd_calibration_by_decile, file.path(out_dir, "primary_vfd_calibration_by_fitted_decile.csv"))
readr::write_csv(imv_duration_results, file.path(out_dir, "primary_imv_duration_quasipoisson_results.csv"))
readr::write_csv(imv_duration_model_diagnostics, file.path(out_dir, "primary_imv_duration_model_diagnostics.csv"))
readr::write_csv(imv_duration_poisson_quasi_comparison, file.path(out_dir, "primary_imv_duration_poisson_vs_quasipoisson_diagnostics.csv"))
readr::write_csv(imv_duration_calibration_by_decile, file.path(out_dir, "primary_imv_duration_calibration_by_fitted_decile.csv"))
readr::write_csv(no_peak_covid_cohort_summary, file.path(out_dir, "sensitivity_exclude_peak_covid_12m_cohort_summary.csv"))
readr::write_csv(no_peak_covid_logistic_results, file.path(out_dir, "sensitivity_exclude_peak_covid_12m_mortality_day28_logistic_results.csv"))
readr::write_csv(no_peak_covid_cox_results, file.path(out_dir, "sensitivity_exclude_peak_covid_12m_mortality_cox_results.csv"))
readr::write_csv(no_peak_covid_vfd_results, file.path(out_dir, "sensitivity_exclude_peak_covid_12m_vfd_quasipoisson_results.csv"))
readr::write_csv(no_peak_covid_imv_duration_results, file.path(out_dir, "sensitivity_exclude_peak_covid_12m_imv_duration_quasipoisson_results.csv"))
readr::write_csv(no_peak_covid_pooling_table, file.path(out_dir, "sensitivity_exclude_peak_covid_12m_pooling_table.csv"))
readr::write_csv(primary_vs_no_peak_covid_pooling_table, file.path(out_dir, "primary_vs_exclude_peak_covid_12m_pooling_table.csv"))

message("Wrote reviewer-optimized primary outputs:")
message(" - ", file.path(out_dir, "cohort_summary_reviewer_optimized.csv"))
message(" - ", file.path(out_dir, "table1_baseline_chronic_disease_prevalence.csv"))
message(" - ", file.path(out_dir, "table1_continuous_site_resubmission.csv"))
message(" - ", file.path(out_dir, "table1_categorical_site_resubmission.csv"))
message(" - ", file.path(out_dir, "table1_resubmission_long.csv"))
message(" - ", file.path(out_dir, "table1_resubmission_missing_variables.csv"))
message(" - ", file.path(out_dir, "primary_mortality_day28_logistic_results.csv"))
message(" - ", file.path(out_dir, "primary_mortality_cox_results.csv"))
message(" - ", file.path(out_dir, "primary_mortality_cox_ph_diagnostics.csv"))
message(" - ", file.path(out_dir, "primary_vfd_quasipoisson_results.csv"))
message(" - ", file.path(out_dir, "primary_vfd_model_diagnostics.csv"))
message(" - ", file.path(out_dir, "primary_vfd_poisson_vs_quasipoisson_diagnostics.csv"))
message(" - ", file.path(out_dir, "primary_vfd_calibration_by_fitted_decile.csv"))
message(" - ", file.path(out_dir, "primary_imv_duration_quasipoisson_results.csv"))
message(" - ", file.path(out_dir, "primary_imv_duration_model_diagnostics.csv"))
message(" - ", file.path(out_dir, "primary_imv_duration_poisson_vs_quasipoisson_diagnostics.csv"))
message(" - ", file.path(out_dir, "primary_imv_duration_calibration_by_fitted_decile.csv"))
message(" - ", file.path(out_dir, "sensitivity_exclude_peak_covid_12m_cohort_summary.csv"))
message(" - ", file.path(out_dir, "sensitivity_exclude_peak_covid_12m_mortality_day28_logistic_results.csv"))
message(" - ", file.path(out_dir, "sensitivity_exclude_peak_covid_12m_mortality_cox_results.csv"))
message(" - ", file.path(out_dir, "sensitivity_exclude_peak_covid_12m_vfd_quasipoisson_results.csv"))
message(" - ", file.path(out_dir, "sensitivity_exclude_peak_covid_12m_imv_duration_quasipoisson_results.csv"))
message(" - ", file.path(out_dir, "sensitivity_exclude_peak_covid_12m_pooling_table.csv"))
message(" - ", file.path(out_dir, "primary_vs_exclude_peak_covid_12m_pooling_table.csv"))
