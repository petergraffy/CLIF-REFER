#!/usr/bin/env Rscript

# Reviewer-optimized REFER primary models
#
# Intended use:
#   Run the original cohort/linkage code through construction of `arf_exp`, then source
#   this script to replace the first-run model block with the reviewer-focused models.
#
# Core changes from the first-run analysis:
#   1. Replace logistic mortality models with Cox proportional hazards models.
#   2. Use mortality after ARF onset when available; otherwise fall back to index admission.
#   3. Replace invasive ventilation duration with ventilator-free days through day 28.
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
out_dir <- file.path(base_output, "code_resubmit", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vfd_days <- as.integer(get0("vfd_days", inherits = TRUE) %||% 28)
imv_gap_hours <- as.numeric(get0("imv_gap_hours", inherits = TRUE) %||% 6)
charlson_lookback_days <- as.integer(get0("charlson_lookback_days", inherits = TRUE) %||% 365)
charlson_include_index <- isTRUE(get0("charlson_include_index", inherits = TRUE) %||% TRUE)

message("Site: ", site_name)
message("Reviewer-optimized output: ", out_dir)

args <- commandArgs(trailingOnly = TRUE)
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

fallback_vfd_from_vent_hours <- function(df) {
  vent_hours <- if ("vent_hours" %in% names(df)) as.numeric(df$vent_hours) else rep(0, nrow(df))
  df %>%
    mutate(
      imv_days_through_vfd = pmin(coalesce(vent_hours, 0) / 24, vfd_days),
      died_before_vfd_horizon = mortality_event == 1L & mortality_ftime_days <= vfd_days,
      ventilator_free_days = if_else(
        died_before_vfd_horizon,
        0,
        pmax(vfd_days - imv_days_through_vfd, 0)
      )
    )
}

add_vfd <- function(df) {
  if ("ventilator_free_days" %in% names(df)) {
    return(df %>%
      mutate(
        ventilator_free_days = pmin(pmax(as.numeric(ventilator_free_days), 0), vfd_days),
        imv_days_through_vfd = if ("imv_days_through_vfd" %in% names(.)) as.numeric(imv_days_through_vfd) else NA_real_,
        died_before_vfd_horizon = if ("died_before_vfd_horizon" %in% names(.)) died_before_vfd_horizon else NA
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

  imv_days <- imv_runs %>%
    mutate(
      vfd_end = t0 + days(vfd_days),
      run_start = pmax(imv_start, t0, na.rm = TRUE),
      run_end = pmin(imv_end + hours(imv_gap_hours), vfd_end, na.rm = TRUE),
      imv_days = pmax(as.numeric(difftime(run_end, run_start, units = "days")), 0)
    ) %>%
    filter(is.finite(imv_days), imv_days > 0) %>%
    group_by(hospitalization_id) %>%
    summarise(imv_days_through_vfd = sum(imv_days, na.rm = TRUE), .groups = "drop")

  df %>%
    left_join(imv_days, by = "hospitalization_id") %>%
    mutate(
      imv_days_through_vfd = coalesce(imv_days_through_vfd, 0),
      died_before_vfd_horizon = mortality_event == 1L & mortality_ftime_days <= vfd_days,
      ventilator_free_days = if_else(
        died_before_vfd_horizon,
        0,
        pmax(vfd_days - imv_days_through_vfd, 0)
      ),
      ventilator_free_days = pmin(ventilator_free_days, vfd_days)
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

analysis_df <- arf_exp %>%
  coerce_analysis_frame() %>%
  add_mortality_time() %>%
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

cox_fits <- purrr::map(
  model_specs,
  ~ fit_cox_model(analysis_df, .x$exposure_terms, .x$model_label, adjustment_covars)
)

vfd_fits <- purrr::map(
  model_specs,
  ~ fit_vfd_model(analysis_df, .x$exposure_terms, .x$model_label, adjustment_covars)
)

cox_results <- purrr::map_dfr(cox_fits, ~ tidy_cox(.x$fit, .x$model_df, .x$exposure_terms, .x$model_label, .x$covars))
vfd_results <- purrr::map_dfr(vfd_fits, ~ tidy_vfd(.x$fit, .x$model_df, .x$exposure_terms, .x$model_label, .x$covars))

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

cohort_summary <- tibble(
  site = site_name,
  n_arf = nrow(analysis_df),
  n_patients = n_distinct(analysis_df$patient_id),
  mortality_events = sum(analysis_df$mortality_event == 1L, na.rm = TRUE),
  median_mortality_followup_days = median(analysis_df$mortality_ftime_days, na.rm = TRUE),
  mean_vfd = mean(analysis_df$ventilator_free_days, na.rm = TRUE),
  median_vfd = median(analysis_df$ventilator_free_days, na.rm = TRUE),
  vfd_horizon_days = vfd_days,
  mean_charlson = mean(analysis_df$charlson_score, na.rm = TRUE),
  median_charlson = median(analysis_df$charlson_score, na.rm = TRUE),
  pm25_exposure_column = first_present(arf_exp, c("pm25_12m_zcta", "pm25_12m_mean", "pm25_mean")),
  no2_exposure_column = first_present(arf_exp, c("no2_12m_zcta", "no2_12m_mean", "no2_mean")),
  adjustment_covariates = paste(adjustment_covars, collapse = " + ")
)

readr::write_csv(analysis_df, file.path(out_dir, "analysis_dataset_reviewer_optimized.csv"))
readr::write_csv(cohort_summary, file.path(out_dir, "cohort_summary_reviewer_optimized.csv"))
readr::write_csv(cox_results, file.path(out_dir, "primary_mortality_cox_results.csv"))
readr::write_csv(ph_results, file.path(out_dir, "primary_mortality_cox_ph_diagnostics.csv"))
readr::write_csv(vfd_results, file.path(out_dir, "primary_vfd_quasipoisson_results.csv"))

message("Wrote reviewer-optimized primary outputs:")
message(" - ", file.path(out_dir, "cohort_summary_reviewer_optimized.csv"))
message(" - ", file.path(out_dir, "primary_mortality_cox_results.csv"))
message(" - ", file.path(out_dir, "primary_mortality_cox_ph_diagnostics.csv"))
message(" - ", file.path(out_dir, "primary_vfd_quasipoisson_results.csv"))
