#!/usr/bin/env Rscript

# Federated supplementary export: all-CLIF patient counts and raw mortality rates.
#
# This deliberately does not restrict to the ARF cohort. Each site can run this
# script locally and share aggregate CSVs only.

suppressPackageStartupMessages({
  library(dplyr)
  library(janitor)
  library(jsonlite)
  library(lubridate)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

args <- commandArgs(trailingOnly = TRUE)
arg_or <- function(i, default = NA_character_) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}
repo <- normalizePath(arg_or(1, getwd()), mustWork = FALSE)
config_path <- file.path(repo, "config", "config.json")
if (!file.exists(config_path)) {
  stop("Could not find config/config.json. Usage: Rscript code/resubmission/05_clif_site_mortality_distribution_export.R [repo] [output_dir]")
}

config <- jsonlite::fromJSON(config_path)
site_name <- config$site_name %||% "site"
configured_tables_path <- arg_or(3, config$tables_path)
out_dir <- arg_or(2, file.path(repo, "output", "resubmission", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

has_clif_tables <- function(path) {
  any(file.exists(file.path(path, c("clif_patient.parquet", "patient.parquet", "clif_patient.csv", "patient.csv")))) ||
    any(file.exists(file.path(path, "2.1.0", c("clif_patient.parquet", "patient.parquet"))))
}

resolve_tables_path <- function(path) {
  if (has_clif_tables(path)) return(path)
  nearby_extract <- file.path(dirname(path), "CLIF v2.1")
  if (dir.exists(nearby_extract) && has_clif_tables(nearby_extract)) {
    message("Configured tables_path does not contain CLIF extracts; using nearby CLIF v2.1 extract: ", nearby_extract)
    return(normalizePath(nearby_extract, mustWork = TRUE))
  }
  path
}

tables_path <- resolve_tables_path(configured_tables_path)

message("Site: ", site_name)
message("CLIF tables: ", tables_path)
message("Output: ", out_dir)

safe_ts <- function(x, tz = config$time_zone %||% "UTC") {
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

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(
    ext,
    "csv" = readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    "parquet" = {
      if (!requireNamespace("arrow", quietly = TRUE)) stop("Package `arrow` is required to read parquet files.")
      arrow::read_parquet(path)
    },
    "fst" = {
      if (!requireNamespace("fst", quietly = TRUE)) stop("Package `fst` is required to read fst files.")
      fst::read_fst(path, as.data.table = FALSE)
    },
    stop("Unsupported file extension: ", ext)
  ) %>%
    janitor::clean_names()
}

find_clif_table <- function(table_name) {
  patterns <- c(
    paste0("^clif_", table_name, "\\.(csv|parquet|fst)$"),
    paste0("^", table_name, "\\.(csv|parquet|fst)$")
  )
  files <- list.files(tables_path, recursive = TRUE, full.names = TRUE)
  hits <- files[str_detect(basename(files), regex(paste(patterns, collapse = "|"), ignore_case = TRUE))]
  if (!length(hits)) stop("Could not find CLIF table: ", table_name)
  hits[[1]]
}

get_clif_table <- function(table_name) {
  clif_tables_env <- get0("clif_tables", inherits = TRUE)
  if (!is.null(clif_tables_env)) {
    candidates <- c(paste0("clif_", table_name), table_name)
    key <- intersect(candidates, names(clif_tables_env))
    if (length(key)) return(janitor::clean_names(clif_tables_env[[key[[1]]]]))
  }

  path <- find_clif_table(table_name)
  read_any(path)
}

first_present <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

summarise_mortality <- function(df, stratum_name = "overall", stratum_value = "overall") {
  tibble(
    site = site_name,
    cohort = "All CLIF hospitalizations",
    stratum = stratum_name,
    stratum_value = as.character(stratum_value),
    n_hospitalizations = nrow(df),
    n_patients = n_distinct(df$patient_id),
    n_with_admission_time = sum(!is.na(df$admission_time)),
    n_with_discharge_time = sum(!is.na(df$discharge_time)),
    in_hospital_deaths = sum(df$in_hospital_death == 1L, na.rm = TRUE),
    in_hospital_mortality_rate = in_hospital_deaths / n_hospitalizations,
    death_by_day_28 = sum(df$death_by_day_28 == 1L, na.rm = TRUE),
    death_by_day_28_rate = death_by_day_28 / n_hospitalizations,
    death_by_day_30 = sum(df$death_by_day_30 == 1L, na.rm = TRUE),
    death_by_day_30_rate = death_by_day_30 / n_hospitalizations,
    hospice_discharges = sum(df$hospice_discharge == 1L, na.rm = TRUE),
    hospice_discharge_rate = hospice_discharges / n_hospitalizations,
    discharge_death_or_hospice = sum(df$death_or_hospice_discharge == 1L, na.rm = TRUE),
    discharge_death_or_hospice_rate = discharge_death_or_hospice / n_hospitalizations
  )
}

patient_path <- NA_character_
hospitalization_path <- NA_character_
patient <- get_clif_table("patient")
hospitalization <- get_clif_table("hospitalization")
if (is.null(get0("clif_tables", inherits = TRUE))) {
  patient_path <- find_clif_table("patient")
  hospitalization_path <- find_clif_table("hospitalization")
}

death_col <- first_present(patient, c("death_dttm", "death_datetime", "death_time", "deceased_dttm"))
admit_col <- first_present(hospitalization, c("admission_dttm", "hospital_admission_dttm", "admit_dttm"))
discharge_col <- first_present(hospitalization, c("discharge_dttm", "hospital_discharge_dttm"))
discharge_category_col <- first_present(hospitalization, c("discharge_category", "discharge_disposition", "disposition"))

if (is.na(admit_col) || is.na(discharge_col)) {
  stop("Hospitalization table must contain admission and discharge datetime columns.")
}

patient_min <- patient %>%
  transmute(
    patient_id,
    death_time = if (!is.na(death_col)) safe_ts(.data[[death_col]]) else as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC")
  )

hosp_min <- hospitalization %>%
  transmute(
    patient_id,
    hospitalization_id = if ("hospitalization_id" %in% names(.)) hospitalization_id else NA_character_,
    admission_time = safe_ts(.data[[admit_col]]),
    discharge_time = safe_ts(.data[[discharge_col]]),
    discharge_category = if (!is.na(discharge_category_col)) as.character(.data[[discharge_category_col]]) else NA_character_,
    age_at_admission = if ("age_at_admission" %in% names(.)) suppressWarnings(as.numeric(age_at_admission)) else NA_real_
  )

all_clif <- hosp_min %>%
  left_join(patient_min, by = "patient_id") %>%
  mutate(
    discharge_category_lower = str_to_lower(coalesce(discharge_category, "")),
    expired_discharge = as.integer(str_detect(discharge_category_lower, "expired|deceased|death|died")),
    hospice_discharge = as.integer(str_detect(discharge_category_lower, "hospice")),
    death_time_for_discharge = if_else(expired_discharge == 1L & is.na(death_time), discharge_time, death_time),
    in_hospital_death = as.integer(
      expired_discharge == 1L |
        (!is.na(death_time_for_discharge) &
           !is.na(admission_time) &
           !is.na(discharge_time) &
           death_time_for_discharge >= admission_time &
           death_time_for_discharge <= discharge_time)
    ),
    death_by_day_28 = as.integer(
      !is.na(death_time_for_discharge) &
        !is.na(admission_time) &
        death_time_for_discharge >= admission_time &
        death_time_for_discharge <= admission_time + lubridate::days(28)
    ),
    death_by_day_30 = as.integer(
      !is.na(death_time_for_discharge) &
        !is.na(admission_time) &
        death_time_for_discharge >= admission_time &
        death_time_for_discharge <= admission_time + lubridate::days(30)
    ),
    death_or_hospice_discharge = as.integer(in_hospital_death == 1L | hospice_discharge == 1L),
    admission_year = lubridate::year(admission_time),
    age_group = cut(
      age_at_admission,
      breaks = c(-Inf, 17, 39, 49, 59, 69, 79, Inf),
      labels = c("<18", "18-39", "40-49", "50-59", "60-69", "70-79", "80+"),
      right = TRUE
    )
  )

overall <- summarise_mortality(all_clif)

by_year <- all_clif %>%
  filter(!is.na(admission_year)) %>%
  group_by(admission_year) %>%
  group_modify(~ summarise_mortality(.x, "admission_year", .y$admission_year)) %>%
  ungroup()

by_age_group <- all_clif %>%
  filter(!is.na(age_group)) %>%
  group_by(age_group) %>%
  group_modify(~ summarise_mortality(.x, "age_group", .y$age_group)) %>%
  ungroup()

supplementary_table <- bind_rows(overall, by_year, by_age_group)

metadata <- tibble(
  site = site_name,
  patient_table_path = patient_path,
  hospitalization_table_path = hospitalization_path,
  death_column = death_col,
  admission_column = admit_col,
  discharge_column = discharge_col,
  discharge_category_column = discharge_category_col,
  generated_at = as.character(Sys.time()),
  contains_phi = FALSE,
  note = "Aggregate counts and crude rates only; no row-level patient data exported."
)

readr::write_csv(supplementary_table, file.path(out_dir, "supplement_all_clif_raw_mortality_by_site.csv"))
readr::write_csv(overall, file.path(out_dir, "supplement_all_clif_raw_mortality_overall.csv"))
readr::write_csv(by_year, file.path(out_dir, "supplement_all_clif_raw_mortality_by_admission_year.csv"))
readr::write_csv(by_age_group, file.path(out_dir, "supplement_all_clif_raw_mortality_by_age_group.csv"))
readr::write_csv(metadata, file.path(out_dir, "supplement_all_clif_raw_mortality_metadata.csv"))

message("Wrote:")
message("  ", file.path(out_dir, "supplement_all_clif_raw_mortality_by_site.csv"))
message("  ", file.path(out_dir, "supplement_all_clif_raw_mortality_overall.csv"))
message("  ", file.path(out_dir, "supplement_all_clif_raw_mortality_by_admission_year.csv"))
message("  ", file.path(out_dir, "supplement_all_clif_raw_mortality_by_age_group.csv"))
message("  ", file.path(out_dir, "supplement_all_clif_raw_mortality_metadata.csv"))

print(overall)
