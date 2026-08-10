#!/usr/bin/env Rscript

# Poolable sensitivity models for the reviewer-optimized REFER primary analysis.

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(forcats)
  library(glue)
  library(janitor)
  library(lubridate)
  library(purrr)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !file.exists(args[[1]])) {
  stop("Usage: Rscript code/resubmission/06_primary_sensitivity_models.R <primary_analysis_dataset.csv> [output_dir] [no_icu_los_dataset.csv]")
}
arg_or <- function(i, default = NA_character_) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
repo_guess <- normalizePath(file.path(dirname(input_path), "..", ".."), mustWork = FALSE)
repo <- if (dir.exists(file.path(getwd(), "code", "resubmission"))) normalizePath(getwd(), mustWork = TRUE) else repo_guess
resolve_repo_path <- function(path, default) {
  path <- path %||% default
  if (is.na(path) || !nzchar(path)) return(default)
  if (grepl("^(/|[A-Za-z]:[/\\\\]|~)", path)) path else file.path(repo, path)
}
out_dir <- arg_or(2, file.path(repo, "output", "resubmission", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_no_icu <- c(
  arg_or(3),
  Sys.getenv("REFER_NO_ICU_LOS_DATASET", NA_character_),
  file.path(dirname(input_path), "resubmission_analysis_dataset_no_icu_los_restriction.csv"),
  file.path(dirname(input_path), "analysis_dataset_no_icu_los_restriction.csv")
)
no_icu_hits <- candidate_no_icu[file.exists(candidate_no_icu)]
no_icu_path <- if (length(no_icu_hits)) no_icu_hits[[1]] else NA_character_

site_name <- "site"
primary_followup_days <- 28L
short_followup_days <- as.integer(Sys.getenv("REFER_SHORT_FOLLOWUP_DAYS", "14"))

message("Reading primary analysis dataset: ", input_path)
message("Output: ", out_dir)

safe_ts <- function(x, tz = Sys.timezone()) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = tz))
  if (is.numeric(x)) {
    x2 <- ifelse(x > 1e12, x / 1000, x)
    return(as.POSIXct(x2, origin = "1970-01-01", tz = tz))
  }
  suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c("ymd_HMS", "ymd_HM", "ymd", "ymdTz", "ymdT", "mdy_HMS", "mdy_HM", "mdy", "dmy_HMS", "dmy_HM", "dmy"),
    quiet = TRUE
  ))
}

normalize_zip <- function(x) {
  x <- str_replace_all(as.character(x), "[^0-9]", "")
  x <- ifelse(nchar(x) >= 5, substr(x, 1, 5), x)
  ifelse(nchar(x) == 5, x, NA_character_)
}

first_present <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

col_or <- function(df, candidates, default = NA) {
  hit <- first_present(df, candidates)
  if (is.na(hit)) rep(default, nrow(df)) else df[[hit]]
}

read_analysis_dataset <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    janitor::clean_names()

  t0_col <- first_present(df, c("t0", "arf_onset", "index_admit", "admission_dttm"))
  discharge_col <- first_present(df, c("discharge_time", "index_discharge", "discharge_dttm"))
  pm25_col <- first_present(df, c("pm25_12m_zcta", "pm25_12m_mean", "pm25_mean", "pm25_exposure"))
  no2_col <- first_present(df, c("no2_12m_zcta", "no2_12m_mean", "no2_mean", "no2_exposure"))
  o3_col <- first_present(df, c("o3_12m_zcta", "o3_12m_mean", "o3_mean", "o3_exposure"))
  if (is.na(t0_col) || is.na(discharge_col) || is.na(pm25_col) || is.na(no2_col)) {
    stop("Analysis dataset is missing required date or PM2.5/NO2 columns: ", path)
  }

  df %>%
    mutate(
      site = col_or(df, "site", "site"),
      t0 = safe_ts(.data[[t0_col]]),
      discharge_time = safe_ts(.data[[discharge_col]]),
      zipcode_five_digit = normalize_zip(col_or(df, c("zipcode_five_digit", "zip", "postal_code"))),
      mortality_ftime_days = as.numeric(mortality_ftime_days),
      mortality_event = as.integer(mortality_event),
      mortality_day28_event = as.integer(
        col_or(df, "mortality_day28_event") %||%
          (mortality_event == 1L & mortality_ftime_days <= primary_followup_days)
      ),
      mortality_day14_event = as.integer(mortality_event == 1L & mortality_ftime_days <= short_followup_days),
      mortality_ftime_days_14 = pmax(pmin(mortality_ftime_days, short_followup_days), 1 / 24),
      mortality_event_14 = as.integer(mortality_event == 1L & mortality_ftime_days <= short_followup_days),
      ventilator_free_days = as.numeric(ventilator_free_days),
      imv_days_through_vfd = if ("imv_days_through_vfd" %in% names(.)) as.numeric(imv_days_through_vfd) else NA_real_,
      imv_days_through_day14 = if ("imv_days_through_day14" %in% names(.)) as.numeric(imv_days_through_day14) else NA_real_,
      ventilator_free_days_day14 = if ("ventilator_free_days_day14" %in% names(.)) {
        as.numeric(ventilator_free_days_day14)
      } else if ("imv_days_through_day14" %in% names(.)) {
        if_else(
          mortality_day14_event == 1L,
          0,
          pmax(short_followup_days - as.numeric(imv_days_through_day14), 0)
        )
      } else {
        NA_real_
      },
      pm25_exposure = as.numeric(.data[[pm25_col]]),
      no2_exposure = as.numeric(.data[[no2_col]]),
      o3_exposure = if (!is.na(o3_col)) as.numeric(.data[[o3_col]]) else NA_real_,
      pm25_per_5 = pm25_exposure / 5,
      no2_per_10 = no2_exposure / 10,
      o3_per_10 = o3_exposure / 10,
      pm25_36m_per_5 = if ("pm25_36m_zcta" %in% names(.)) as.numeric(pm25_36m_zcta) / 5 else NA_real_,
      no2_36m_per_10 = if ("no2_36m_zcta" %in% names(.)) as.numeric(no2_36m_zcta) / 10 else NA_real_,
      age_10 = as.numeric(col_or(df, "age_10") %||% (col_or(df, "age") / 10)),
      sex = factor(sex),
      race_ethnicity = forcats::fct_lump_min(factor(race_ethnicity), min = 100, other_level = "Other/Unknown"),
      index_year = as.integer(col_or(df, "index_year") %||% lubridate::year(t0)),
      index_year_f = factor(col_or(df, "index_year_f") %||% index_year),
      acs_median_household_income_10k = as.numeric(
        col_or(df, "acs_median_household_income_10k") %||%
          (col_or(df, "acs_median_household_income") / 10000)
      ),
      sofa_total = if ("sofa_total" %in% names(.)) as.numeric(sofa_total) else NA_real_
    )
}

primary_df <- read_analysis_dataset(input_path)
site_name <- dplyr::first(stats::na.omit(primary_df$site)) %||% site_name

no_icu_df <- if (!is.na(no_icu_path)) {
  message("Reading no-ICU-LOS-restriction dataset: ", no_icu_path)
  read_analysis_dataset(no_icu_path)
} else {
  tibble()
}

social_covars <- function(df) {
  covars <- c(
    "acs_pct_poverty",
    "acs_pct_unemployed",
    "acs_pct_no_vehicle",
    "acs_pct_nonwhite",
    "acs_median_household_income_10k",
    "acs_pct_bachelor_plus"
  )
  intersect(covars, names(df))
}

base_adjustment_covars <- function(df, include_sofa = FALSE) {
  covars <- c(
    "age_10",
    "sex",
    "race_ethnicity",
    "charlson_score",
    "index_year_f",
    social_covars(df)
  )
  if (include_sofa) covars <- append(covars, "sofa_total", after = match("charlson_score", covars))
  covars
}

drop_uninformative_covars <- function(model_df, covars, exposure_terms) {
  purrr::keep(covars, function(covar) {
    if (covar %in% exposure_terms) return(TRUE)
    x <- model_df[[covar]]
    if (is.factor(x) || is.character(x)) return(dplyr::n_distinct(x, na.rm = TRUE) >= 2)
    isTRUE(stats::var(as.numeric(x), na.rm = TRUE) > 0)
  }) %>%
    unname()
}

fit_model_set <- function(df, sensitivity, exposure_terms, pollutant_model, include_sofa = FALSE,
                          mortality_outcome = "mortality_day28_event",
                          vfd_outcome = "ventilator_free_days",
                          imv_duration_outcome = "imv_days_through_vfd",
                          cox_time = "mortality_ftime_days",
                          cox_event = "mortality_event") {
  covars <- c(exposure_terms, base_adjustment_covars(df, include_sofa))
  needed <- c(mortality_outcome, vfd_outcome, imv_duration_outcome, cox_time, cox_event, covars)
  missing_needed <- setdiff(needed, names(df))
  if (length(missing_needed)) {
    return(list(results = tibble(), status = tibble(
      site = site_name, sensitivity, pollutant_model,
      status = "skipped", reason = paste("Missing columns:", paste(missing_needed, collapse = ", "))
    )))
  }

  run_one <- function(outcome_type) {
    response_cols <- switch(
      outcome_type,
      logistic = mortality_outcome,
      cox = c(cox_time, cox_event),
      vfd = vfd_outcome,
      imv_duration = imv_duration_outcome
    )
    model_df <- df %>%
      dplyr::select(all_of(c(response_cols, covars))) %>%
      tidyr::drop_na()
    if (!nrow(model_df)) return(tibble())
    if (outcome_type %in% c("logistic", "cox") && sum(model_df[[if (outcome_type == "cox") cox_event else mortality_outcome]] == 1L) < 10) {
      return(tibble())
    }
    covars_i <- drop_uninformative_covars(model_df, covars, exposure_terms)
    model_df <- model_df %>% dplyr::select(all_of(c(response_cols, covars_i)))
    rhs <- paste(covars_i, collapse = " + ")

    if (outcome_type == "logistic") {
      fit <- stats::glm(as.formula(paste0(mortality_outcome, " ~ ", rhs)), data = model_df, family = binomial("logit"))
      return(broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term %in% exposure_terms) %>%
        transmute(
          site = site_name, sensitivity, outcome = mortality_outcome,
          analysis_model = "Logistic regression", pollutant_model, term,
          effect_measure = "odds_ratio", n = nrow(model_df),
          events = sum(model_df[[mortality_outcome]] == 1L), mean_vfd = NA_real_,
          estimate, conf_low = conf.low, conf_high = conf.high, p_value = p.value,
          dispersion = NA_real_, adjustment_set = paste(setdiff(covars_i, exposure_terms), collapse = " + ")
        ))
    }

    if (outcome_type == "cox") {
      fit <- survival::coxph(
        as.formula(paste0("Surv(", cox_time, ", ", cox_event, ") ~ ", rhs)),
        data = model_df,
        ties = "efron",
        x = TRUE
      )
      return(broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term %in% exposure_terms) %>%
        transmute(
          site = site_name, sensitivity, outcome = "mortality_time_to_event",
          analysis_model = "Cox proportional hazards", pollutant_model, term,
          effect_measure = "hazard_ratio", n = nrow(model_df),
          events = sum(model_df[[cox_event]] == 1L), mean_vfd = NA_real_,
          estimate, conf_low = conf.low, conf_high = conf.high, p_value = p.value,
          dispersion = NA_real_, adjustment_set = paste(setdiff(covars_i, exposure_terms), collapse = " + ")
        ))
    }

    continuous_outcome <- if (outcome_type == "imv_duration") imv_duration_outcome else vfd_outcome
    fit <- stats::glm(as.formula(paste0(continuous_outcome, " ~ ", rhs)), data = model_df, family = quasipoisson("log"))
    broom::tidy(fit) %>%
      filter(term %in% exposure_terms) %>%
      mutate(beta = estimate) %>%
      transmute(
        site = site_name, sensitivity, outcome = continuous_outcome,
        analysis_model = "Quasi-Poisson regression", pollutant_model, term,
        effect_measure = "ratio_of_means", n = nrow(model_df), events = NA_integer_,
        mean_vfd = if (outcome_type == "vfd") mean(model_df[[vfd_outcome]], na.rm = TRUE) else NA_real_,
        mean_imv_days = if (outcome_type == "imv_duration") mean(model_df[[imv_duration_outcome]], na.rm = TRUE) else NA_real_,
        estimate = exp(beta),
        conf_low = exp(beta - 1.96 * std.error),
        conf_high = exp(beta + 1.96 * std.error),
        p_value = p.value,
        dispersion = summary(fit)$dispersion,
        adjustment_set = paste(setdiff(covars_i, exposure_terms), collapse = " + ")
      )
  }

  logistic_results <- run_one("logistic")
  cox_results <- run_one("cox")
  vfd_results <- run_one("vfd")
  imv_duration_results <- run_one("imv_duration")
  results <- bind_rows(logistic_results, cox_results, vfd_results, imv_duration_results)
  completed_models <- c(
    if (nrow(logistic_results)) "logistic" else character(),
    if (nrow(cox_results)) "cox" else character(),
    if (nrow(vfd_results)) "quasi_poisson_vfd" else character(),
    if (nrow(imv_duration_results)) "quasi_poisson_imv_duration" else character()
  )
  status_value <- if (length(completed_models) == 4L) {
    "completed"
  } else if (length(completed_models) > 0L) {
    "completed_partial"
  } else {
    "skipped"
  }
  reason_value <- if (identical(status_value, "completed")) {
    NA_character_
  } else if (identical(status_value, "completed_partial")) {
    "At least one model type had no complete cases or too few events; check completed_models."
  } else {
    "No complete cases or too few mortality events"
  }
  status <- tibble(
    site = site_name, sensitivity, pollutant_model,
    status = status_value,
    completed_models = paste(completed_models, collapse = ";"),
    reason = reason_value
  )
  list(results = results, status = status)
}

model_specs <- tibble::tribble(
  ~pollutant_model, ~exposure_terms,
  "PM25 single-pollutant", list("pm25_per_5"),
  "NO2 single-pollutant", list("no2_per_10"),
  "PM25 + NO2", list(c("pm25_per_5", "no2_per_10"))
)

derive_36m_exposures <- function(df) {
  if (all(c("pm25_36m_per_5", "no2_36m_per_10") %in% names(df)) &&
      any(!is.na(df$pm25_36m_per_5)) &&
      any(!is.na(df$no2_36m_per_10))) {
    return(df)
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    warning("Package `arrow` is not available; skipping 3-year exposure derivation.")
    return(df)
  }

  pm25_path <- resolve_repo_path(Sys.getenv(
    "REFER_ZCTA_PM25_MONTHLY_PATH",
    file.path(repo, "exposome", "zcta", "air_pollution_zcta_pm25_monthly_2005_2023.parquet")
  ), file.path(repo, "exposome", "zcta", "air_pollution_zcta_pm25_monthly_2005_2023.parquet"))
  no2_monthly_dir <- resolve_repo_path(Sys.getenv(
    "REFER_ZCTA_NO2_MONTHLY_DIR",
    file.path(repo, "exposome", "zcta", "no2_monthly")
  ), file.path(repo, "exposome", "zcta", "no2_monthly"))
  no2_annual_path <- resolve_repo_path(Sys.getenv(
    "REFER_ZCTA_NO2_ANNUAL_PATH",
    file.path(repo, "exposome", "zcta", "air_pollution_zcta_no2_annual_2005_2025.parquet")
  ), file.path(repo, "exposome", "zcta", "air_pollution_zcta_no2_annual_2005_2025.parquet"))

  if (!file.exists(pm25_path) || !file.exists(no2_annual_path)) {
    warning("Required PM2.5/NO2 ZCTA files are absent; skipping 3-year exposure derivation.")
    return(df)
  }

  index_tbl <- df %>%
    transmute(
      hospitalization_id,
      zipcode_five_digit,
      exposure_start = floor_date(as.Date(t0), "month") %m-% months(36),
      exposure_end = floor_date(as.Date(t0), "month") - days(1),
      index_year = lubridate::year(t0)
    ) %>%
    filter(!is.na(zipcode_five_digit), !is.na(exposure_start), !is.na(exposure_end))

  pm25_monthly <- arrow::read_parquet(pm25_path) %>%
    transmute(
      zipcode_five_digit = normalize_zip(zip),
      exposure_month = as.Date(sprintf("%04d-%02d-01", as.integer(year), as.integer(month))),
      pm25_ug_m3 = as.numeric(pm25_ug_m3)
    ) %>%
    filter(!is.na(zipcode_five_digit), !is.na(exposure_month), !is.na(pm25_ug_m3))

  pm25_36m <- index_tbl %>%
    left_join(pm25_monthly, by = "zipcode_five_digit", relationship = "many-to-many") %>%
    filter(exposure_month >= exposure_start, exposure_month <= exposure_end) %>%
    group_by(hospitalization_id) %>%
    summarise(
      pm25_36m_zcta = mean(pm25_ug_m3, na.rm = TRUE),
      pm25_36m_months_observed = n_distinct(exposure_month),
      .groups = "drop"
    )

  no2_monthly_files <- if (dir.exists(no2_monthly_dir)) {
    list.files(no2_monthly_dir, pattern = "^no2_zcta_monthly_[0-9]{4}\\.parquet$", full.names = TRUE)
  } else {
    character()
  }

  no2_36m_monthly <- tibble(hospitalization_id = character(), no2_36m_monthly_zcta = numeric(), no2_36m_months_observed = integer())
  if (length(no2_monthly_files)) {
    no2_monthly_raw <- purrr::map_dfr(no2_monthly_files, arrow::read_parquet) %>%
      janitor::clean_names()
    no2_monthly_value_col <- first_present(no2_monthly_raw, c("no2_ppbv", "no2_ppb", "no2"))
    if (is.na(no2_monthly_value_col)) {
      warning("Monthly NO2 files do not contain a recognized NO2 column; using annual fallback only.")
    } else {
      no2_monthly <- no2_monthly_raw %>%
      transmute(
        zipcode_five_digit = normalize_zip(zip),
        exposure_month = as.Date(sprintf("%04d-%02d-01", as.integer(year), as.integer(month))),
        no2_ppb = as.numeric(.data[[no2_monthly_value_col]])
      ) %>%
      filter(!is.na(zipcode_five_digit), !is.na(exposure_month), !is.na(no2_ppb)) %>%
      distinct(zipcode_five_digit, exposure_month, .keep_all = TRUE)

      no2_36m_monthly <- index_tbl %>%
        left_join(no2_monthly, by = "zipcode_five_digit", relationship = "many-to-many") %>%
        filter(exposure_month >= exposure_start, exposure_month <= exposure_end) %>%
        group_by(hospitalization_id) %>%
        summarise(
          no2_36m_monthly_zcta = mean(no2_ppb, na.rm = TRUE),
          no2_36m_months_observed = n_distinct(exposure_month),
          .groups = "drop"
        ) %>%
        mutate(no2_36m_monthly_zcta = if_else(no2_36m_months_observed == 36L, no2_36m_monthly_zcta, NA_real_))
    }
  }

  no2_annual_raw <- arrow::read_parquet(no2_annual_path) %>%
    janitor::clean_names()
  no2_annual_value_col <- first_present(no2_annual_raw, c("no2", "no2_ppb", "no2_ppbv"))
  if (is.na(no2_annual_value_col)) {
    warning("Annual NO2 file does not contain a recognized NO2 column; skipping 3-year NO2.")
    return(df)
  }
  no2_annual <- no2_annual_raw %>%
    transmute(
      zipcode_five_digit = normalize_zip(zip),
      no2_exposure_year = as.integer(year),
      no2_ppb = as.numeric(.data[[no2_annual_value_col]])
    ) %>%
    filter(!is.na(zipcode_five_digit), !is.na(no2_exposure_year), !is.na(no2_ppb)) %>%
    distinct(zipcode_five_digit, no2_exposure_year, .keep_all = TRUE)

  no2_36m_annual <- index_tbl %>%
    transmute(
      hospitalization_id,
      zipcode_five_digit,
      no2_exposure_year = purrr::map(index_year, ~ (.x - 3L):(.x - 1L))
    ) %>%
    tidyr::unnest(no2_exposure_year) %>%
    left_join(no2_annual, by = c("zipcode_five_digit", "no2_exposure_year")) %>%
    group_by(hospitalization_id) %>%
    summarise(
      no2_36m_annual_zcta = mean(no2_ppb, na.rm = TRUE),
      no2_36m_annual_years_observed = sum(!is.na(no2_ppb)),
      .groups = "drop"
    ) %>%
    mutate(no2_36m_annual_zcta = if_else(no2_36m_annual_years_observed > 0L, no2_36m_annual_zcta, NA_real_))

  no2_36m <- index_tbl %>%
    select(hospitalization_id) %>%
    left_join(no2_36m_monthly, by = "hospitalization_id") %>%
    left_join(no2_36m_annual, by = "hospitalization_id") %>%
    mutate(
      no2_36m_zcta = coalesce(no2_36m_monthly_zcta, no2_36m_annual_zcta),
      no2_36m_source = case_when(
        !is.na(no2_36m_monthly_zcta) ~ "monthly_complete_36m",
        !is.na(no2_36m_annual_zcta) ~ "annual_prior_3y_fallback",
        TRUE ~ NA_character_
      )
    ) %>%
    select(hospitalization_id, no2_36m_zcta, no2_36m_source, everything())

  df %>%
    select(-any_of(c(
      "pm25_36m_zcta", "pm25_36m_per_5", "pm25_36m_months_observed",
      "no2_36m_zcta", "no2_36m_per_10", "no2_36m_source",
      "no2_36m_monthly_zcta", "no2_36m_months_observed",
      "no2_36m_annual_zcta", "no2_36m_annual_years_observed"
    ))) %>%
    left_join(pm25_36m, by = "hospitalization_id") %>%
    left_join(no2_36m, by = "hospitalization_id") %>%
    mutate(
      pm25_36m_per_5 = pm25_36m_zcta / 5,
      no2_36m_per_10 = no2_36m_zcta / 10
    )
}

run_specs <- function(df, sensitivity, specs = model_specs, include_sofa = FALSE,
                      mortality_outcome = "mortality_day28_event",
                      vfd_outcome = "ventilator_free_days",
                      imv_duration_outcome = "imv_days_through_vfd",
                      cox_time = "mortality_ftime_days",
                      cox_event = "mortality_event") {
  out <- purrr::pmap(specs, function(pollutant_model, exposure_terms) {
    fit_model_set(
      df = df,
      sensitivity = sensitivity,
      exposure_terms = unlist(exposure_terms),
      pollutant_model = pollutant_model,
      include_sofa = include_sofa,
      mortality_outcome = mortality_outcome,
      vfd_outcome = vfd_outcome,
      imv_duration_outcome = imv_duration_outcome,
      cox_time = cox_time,
      cox_event = cox_event
    )
  })
  list(
    results = bind_rows(purrr::map(out, "results")),
    status = bind_rows(purrr::map(out, "status"))
  )
}

o3_specs <- tibble::tribble(
  ~pollutant_model, ~exposure_terms,
  "PM25 + NO2 + O3", list(c("pm25_per_5", "no2_per_10", "o3_per_10"))
)

exposure_36m_specs <- tibble::tribble(
  ~pollutant_model, ~exposure_terms,
  "PM25 3-year single-pollutant", list("pm25_36m_per_5"),
  "NO2 3-year single-pollutant", list("no2_36m_per_10"),
  "PM25 + NO2 3-year", list(c("pm25_36m_per_5", "no2_36m_per_10"))
)

primary_df_36m <- derive_36m_exposures(primary_df)

runs <- list(
  sofa_total_24h = run_specs(primary_df, "add_sofa_total_first_24h", include_sofa = TRUE),
  o3_copollutant = run_specs(primary_df, "add_o3_copollutant", specs = o3_specs),
  exposure_window_3y = run_specs(primary_df_36m, "exposure_window_3y", specs = exposure_36m_specs),
  followup_window_14d = run_specs(
    primary_df,
    "followup_window_14d",
    mortality_outcome = "mortality_day14_event",
    vfd_outcome = "ventilator_free_days_day14",
    imv_duration_outcome = "imv_days_through_day14",
    cox_time = "mortality_ftime_days_14",
    cox_event = "mortality_event_14"
  )
)

if (nrow(no_icu_df)) {
  runs$remove_icu_los_24h_restriction <- run_specs(no_icu_df, "remove_icu_los_24h_restriction")
} else {
  runs$remove_icu_los_24h_restriction <- list(
    results = tibble(),
    status = tibble(
      site = site_name,
      sensitivity = "remove_icu_los_24h_restriction",
      pollutant_model = NA_character_,
      status = "skipped",
      reason = "No companion no-ICU-LOS-restriction dataset found"
    )
  )
}

sensitivity_results <- bind_rows(purrr::map(runs, "results")) %>%
  mutate(
    exposure_window = case_when(
      str_detect(sensitivity, "3y") ~ "36 months before ARF onset",
      TRUE ~ "12 months before ARF onset"
    ),
    followup_window_days = case_when(
      sensitivity == "followup_window_14d" ~ short_followup_days,
      TRUE ~ primary_followup_days
    ),
    .after = sensitivity
  )

sensitivity_status <- bind_rows(purrr::map(runs, "status"))

cohort_summaries <- bind_rows(
  primary_df %>%
    summarise(
      site = site_name,
      sensitivity = "primary_dataset",
      n_arf = n(),
      n_patients = n_distinct(patient_id),
      n_day28_deaths = sum(mortality_day28_event == 1L, na.rm = TRUE),
      n_day14_deaths = sum(mortality_day14_event == 1L, na.rm = TRUE),
      n_with_sofa_total = sum(!is.na(sofa_total)),
      n_with_o3 = sum(!is.na(o3_per_10)),
      n_with_pm25_36m = sum(!is.na(pm25_36m_per_5)),
      n_with_no2_36m = sum(!is.na(no2_36m_per_10)),
      n_with_exact_vfd14 = sum(!is.na(ventilator_free_days_day14)),
      mean_vfd = mean(ventilator_free_days, na.rm = TRUE),
      mean_imv_days = mean(imv_days_through_vfd, na.rm = TRUE),
      mean_imv_days_day14 = mean(imv_days_through_day14, na.rm = TRUE)
    ),
  if (nrow(no_icu_df)) {
    no_icu_df %>%
      summarise(
        site = site_name,
        sensitivity = "remove_icu_los_24h_restriction",
        n_arf = n(),
        n_patients = n_distinct(patient_id),
        n_day28_deaths = sum(mortality_day28_event == 1L, na.rm = TRUE),
        n_day14_deaths = sum(mortality_day14_event == 1L, na.rm = TRUE),
        n_with_sofa_total = sum(!is.na(sofa_total)),
        n_with_o3 = sum(!is.na(o3_per_10)),
        n_with_pm25_36m = sum(!is.na(pm25_36m_per_5)),
        n_with_no2_36m = sum(!is.na(no2_36m_per_10)),
        n_with_exact_vfd14 = sum(!is.na(ventilator_free_days_day14)),
        mean_vfd = mean(ventilator_free_days, na.rm = TRUE),
        mean_imv_days = mean(imv_days_through_vfd, na.rm = TRUE),
        mean_imv_days_day14 = mean(imv_days_through_day14, na.rm = TRUE)
      )
  } else {
    tibble()
  },
  primary_df_36m %>%
    summarise(
      site = site_name,
      sensitivity = "exposure_window_3y",
      n_arf = n(),
      n_patients = n_distinct(patient_id),
      n_day28_deaths = sum(mortality_day28_event == 1L, na.rm = TRUE),
      n_day14_deaths = sum(mortality_day14_event == 1L, na.rm = TRUE),
      n_with_sofa_total = sum(!is.na(sofa_total)),
      n_with_o3 = sum(!is.na(o3_per_10)),
      n_with_pm25_36m = sum(!is.na(pm25_36m_per_5)),
      n_with_no2_36m = sum(!is.na(no2_36m_per_10)),
      n_with_exact_vfd14 = sum(!is.na(ventilator_free_days_day14)),
      mean_vfd = mean(ventilator_free_days, na.rm = TRUE),
      mean_imv_days = mean(imv_days_through_vfd, na.rm = TRUE),
      mean_imv_days_day14 = mean(imv_days_through_day14, na.rm = TRUE)
    )
)

readr::write_csv(sensitivity_results, file.path(out_dir, "primary_sensitivity_models_pooling_table.csv"))
readr::write_csv(sensitivity_status, file.path(out_dir, "primary_sensitivity_models_run_status.csv"))
readr::write_csv(cohort_summaries, file.path(out_dir, "primary_sensitivity_models_cohort_summaries.csv"))
readr::write_csv(
  primary_df_36m %>%
    select(
      any_of(c(
        "site", "patient_id", "hospitalization_id", "zipcode_five_digit", "t0",
        "pm25_36m_zcta", "pm25_36m_months_observed",
        "no2_36m_zcta", "no2_36m_source", "no2_36m_months_observed",
        "no2_36m_annual_years_observed"
      ))
    ),
  file.path(out_dir, "primary_sensitivity_3y_exposure_diagnostics.csv")
)

message("Wrote sensitivity outputs:")
message(" - ", file.path(out_dir, "primary_sensitivity_models_pooling_table.csv"))
message(" - ", file.path(out_dir, "primary_sensitivity_models_run_status.csv"))
message(" - ", file.path(out_dir, "primary_sensitivity_models_cohort_summaries.csv"))
message(" - ", file.path(out_dir, "primary_sensitivity_3y_exposure_diagnostics.csv"))
