# ================================================================================================
# REFER Resubmission Analysis
#
# Reviewer-responsive estimand:
#   - ZCTA-level fixed pre-ARF exposure interval
#   - ARF onset as time zero
#   - cause-specific Cox proportional hazards models
#   - baseline comorbidity adjustment with Charlson score
#   - unadjusted Aalen-Johansen cumulative incidence by PM2.5/NO2 quartile
# ================================================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(broom)
  library(data.table)
  library(dplyr)
  library(forcats)
  library(fst)
  library(ggplot2)
  library(glue)
  library(janitor)
  library(jsonlite)
  library(lubridate)
  library(purrr)
  library(readr)
  library(rlang)
  library(scales)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(cmd_file)) normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = TRUE) else NA_character_
}
repo_guess <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE) else getwd()

config_path <- Sys.getenv(
  "REFER_RESUBMISSION_CONFIG",
  file.path(repo_guess, "code", "resubmission", "resubmission_config.json")
)
if (!file.exists(config_path)) {
  fallback <- file.path(repo_guess, "config", "config.json")
  if (!file.exists(fallback)) {
    stop(
      "No resubmission config found. Copy code/resubmission/resubmission_config_template.json ",
      "to code/resubmission/resubmission_config.json, or set REFER_RESUBMISSION_CONFIG."
    )
  }
  message("Using main config file because resubmission config was not found: ", fallback)
  config_path <- fallback
}

config <- jsonlite::fromJSON(config_path)
repo <- normalizePath(config$repo %||% repo_guess, mustWork = TRUE)
configured_tables_path <- normalizePath(config$tables_path, mustWork = TRUE)
site_name <- config$site_name %||% "site"
time_zone <- config$time_zone %||% "UTC"
study_start <- as.Date(config$study_start %||% "2018-01-01")
study_end <- as.Date(config$study_end %||% "2024-12-31")
followup_days <- as.integer(config$followup_days %||% 30)
pm25_prior_months <- as.integer(config$pm25_prior_months %||% 12)
no2_prior_months <- as.integer(config$no2_prior_months %||% 12)
no2_lag_years <- as.integer(config$no2_lag_years %||% 1)
charlson_lookback_days <- as.integer(config$charlson_lookback_days %||% 365)
charlson_include_index <- isTRUE(config$charlson_include_index_diagnoses %||% TRUE)
imv_gap_hours <- as.numeric(config$imv_gap_hours %||% 6)
successful_extubation_hours <- as.numeric(config$successful_extubation_hours %||% 48)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_root <- file.path(repo, config$output_dir %||% "output/resubmission", stamp)
fig_dir <- file.path(out_root, "figures")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

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
message("Output: ", out_root)

safe_ts <- function(x, tz = time_zone) {
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

normalize_zip <- function(x) {
  x <- str_replace_all(as.character(x), "[^0-9]", "")
  x <- ifelse(nchar(x) >= 5, substr(x, 1, 5), x)
  ifelse(nchar(x) == 5, x, NA_character_)
}

normalize_fio2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    x > 1 & x <= 100 ~ x / 100,
    x >= 0.21 & x <= 1 ~ x,
    TRUE ~ NA_real_
  )
}

harmonize_sex <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  case_when(
    x %in% c("female", "f", "woman") ~ "Female",
    x %in% c("male", "m", "man") ~ "Male",
    TRUE ~ "Other/Unknown"
  )
}

harmonize_race_ethnicity <- function(race, ethnicity) {
  race <- str_to_lower(str_trim(as.character(race)))
  ethnicity <- str_to_lower(str_trim(as.character(ethnicity)))
  is_hispanic <- ethnicity %in% c("hispanic", "latino", "hispanic or latino") |
    (str_detect(ethnicity, "hispanic|latino") & !str_detect(ethnicity, "^non|not|no "))
  case_when(
    is_hispanic & str_detect(race, "white") ~ "Hispanic White",
    is_hispanic & str_detect(race, "black|african") ~ "Hispanic Black",
    is_hispanic ~ "Hispanic Other/Unknown",
    str_detect(race, "white") ~ "Non-Hispanic White",
    str_detect(race, "black|african") ~ "Non-Hispanic Black",
    str_detect(race, "asian") ~ "Asian",
    TRUE ~ "Other/Unknown"
  )
}

find_table_path <- function(tbl_name) {
  exts <- strsplit(config$file_type %||% "csv/parquet/fst", "[/|,; ]+")[[1]]
  exts <- exts[nzchar(exts)]
  if (!length(exts)) exts <- c("csv", "parquet", "fst")
  all_files <- list.files(
    tables_path,
    pattern = paste0("\\.(", paste(unique(exts), collapse = "|"), ")$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(all_files)) stop("No CLIF table files found under: ", tables_path)
  base <- tools::file_path_sans_ext(tolower(basename(all_files)))
  want <- paste0("clif_", tolower(tbl_name))
  candidates <- all_files[base %in% c(want, tolower(tbl_name))]
  if (!length(candidates)) stop("Could not find CLIF table: ", tbl_name)
  candidates[[1]]
}

read_tbl <- function(tbl_name, required = TRUE) {
  path <- tryCatch(find_table_path(tbl_name), error = function(e) {
    if (required) stop(e$message)
    NA_character_
  })
  if (is.na(path)) return(tibble())
  ext <- tolower(tools::file_ext(path))
  out <- switch(
    ext,
    csv = readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    parquet = arrow::read_parquet(path),
    fst = fst::read_fst(path, as.data.table = FALSE),
    stop("Unsupported CLIF table extension: ", ext)
  )
  janitor::clean_names(out)
}

first_existing <- function(df, candidates, default = NA) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) df[[hit[[1]]]] else rep(default, nrow(df))
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

make_quartile <- function(x) {
  q <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)
  q <- unique(q)
  if (length(q) < 3) return(factor(rep(NA_character_, length(x))))
  cut(x, breaks = q, include.lowest = TRUE, labels = paste0("Q", seq_len(length(q) - 1)))
}

quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
quartile_colors <- c(
  "Q1 lowest" = "#2166AC",
  "Q2" = "#67A9CF",
  "Q3" = "#F4A582",
  "Q4 highest" = "#B2182B"
)

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
      plot.title = element_text(face = "bold", size = 21),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16, color = "grey20"),
      strip.text = element_text(face = "bold", size = 18, color = "grey10"),
      panel.spacing.x = unit(12, "pt"),
      panel.spacing.y = unit(12, "pt"),
      plot.margin = margin(8, 32, 8, 8)
    )
}

message("Reading CLIF tables...")
patient <- read_tbl("patient") %>%
  transmute(
    patient_id,
    birth_date = suppressWarnings(as.Date(first_existing(pick(everything()), c("birth_date")))),
    sex = harmonize_sex(first_existing(pick(everything()), c("sex_category", "sex"))),
    race_ethnicity = harmonize_race_ethnicity(
      first_existing(pick(everything()), c("race_category", "race")),
      first_existing(pick(everything()), c("ethnicity_category", "ethnicity"))
    ),
    death_dttm = safe_ts(first_existing(pick(everything()), c("death_dttm"), default = NA))
  )

hospitalization <- read_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm),
    discharge_category = str_to_lower(as.character(first_existing(pick(everything()), c("discharge_category", "discharge_name"), default = NA))),
    age_at_admission = suppressWarnings(as.numeric(first_existing(pick(everything()), c("age_at_admission"), default = NA))),
    zipcode_five_digit = normalize_zip(first_existing(pick(everything()), c("zipcode_five_digit", "zip_code", "zipcode", "postal_code"), default = NA))
  ) %>%
  left_join(patient, by = "patient_id") %>%
  mutate(
    age = coalesce(
      age_at_admission,
      as.numeric(difftime(admission_dttm, as.POSIXct(birth_date, tz = time_zone), units = "days")) / 365.25
    ),
    index_year = year(admission_dttm)
  )

adt <- read_tbl("adt") %>%
  transmute(
    hospitalization_id,
    in_dttm = safe_ts(in_dttm),
    out_dttm = safe_ts(out_dttm),
    location_category = str_to_lower(as.character(location_category))
  )

icu_bounds <- adt %>%
  filter(str_detect(location_category %||% "", "icu"), !is.na(in_dttm)) %>%
  group_by(hospitalization_id) %>%
  summarise(
    first_icu_in = min(in_dttm, na.rm = TRUE),
    last_icu_out = suppressWarnings(max(out_dttm, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(last_icu_out = if_else(is.finite(last_icu_out), last_icu_out, as.POSIXct(NA, tz = time_zone)))

base <- hospitalization %>%
  inner_join(icu_bounds, by = "hospitalization_id") %>%
  filter(
    !is.na(first_icu_in),
    as.Date(first_icu_in) >= study_start,
    as.Date(first_icu_in) <= study_end,
    !is.na(age),
    age >= 18
  ) %>%
  mutate(
    arf_window_start = first_icu_in - hours(24),
    arf_window_end = first_icu_in + hours(24)
  )

base_ids <- base$hospitalization_id

resp_support <- read_tbl("respiratory_support") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    device_category = str_to_lower(str_trim(as.character(device_category))),
    fio2_set = normalize_fio2(first_existing(pick(everything()), c("fio2_set"), default = NA))
  ) %>%
  filter(!is.na(recorded_dttm))

vitals <- read_tbl("vitals") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    vital_category = str_to_lower(as.character(vital_category)),
    vital_value = suppressWarnings(as.numeric(vital_value))
  ) %>%
  filter(!is.na(recorded_dttm), !is.na(vital_value))

labs <- read_tbl("labs") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    hospitalization_id,
    lab_result_dttm = safe_ts(lab_result_dttm),
    lab_category = str_to_lower(as.character(lab_category)),
    lab_value_numeric = suppressWarnings(as.numeric(first_existing(pick(everything()), c("lab_value_numeric", "lab_value"), default = NA)))
  ) %>%
  filter(!is.na(lab_result_dttm), !is.na(lab_value_numeric))

message("Identifying ARF onset...")
window_tbl <- base %>% select(hospitalization_id, arf_window_start, arf_window_end)

fio2_win <- resp_support %>%
  inner_join(window_tbl, by = "hospitalization_id") %>%
  filter(recorded_dttm >= arf_window_start, recorded_dttm <= arf_window_end, !is.na(fio2_set)) %>%
  select(hospitalization_id, fio2_time = recorded_dttm, fio2_set)

spo2_win <- vitals %>%
  filter(vital_category == "spo2") %>%
  inner_join(window_tbl, by = "hospitalization_id") %>%
  filter(recorded_dttm >= arf_window_start, recorded_dttm <= arf_window_end) %>%
  select(hospitalization_id, spo2_time = recorded_dttm, spo2 = vital_value)

po2_win <- labs %>%
  filter(lab_category == "po2_arterial") %>%
  inner_join(window_tbl, by = "hospitalization_id") %>%
  filter(lab_result_dttm >= arf_window_start, lab_result_dttm <= arf_window_end) %>%
  select(hospitalization_id, po2_time = lab_result_dttm, po2 = lab_value_numeric)

pair_nearest <- function(meas, meas_time, fio2, max_gap_h = 1) {
  if (!nrow(meas) || !nrow(fio2)) return(tibble())
  mdt <- as.data.table(meas)
  fdt <- as.data.table(fio2)
  setnames(mdt, meas_time, "meas_time")
  mdt[, meas_time_keep := meas_time]
  setkey(mdt, hospitalization_id, meas_time)
  setkey(fdt, hospitalization_id, fio2_time)
  out <- fdt[
    mdt,
    roll = "nearest",
    on = .(hospitalization_id, fio2_time = meas_time),
    nomatch = 0L
  ][
    , time_diff_h := abs(as.numeric(difftime(meas_time_keep, fio2_time, units = "hours")))
  ][
    time_diff_h <= max_gap_h
  ]
  as_tibble(out)
}

spo2_fio2 <- pair_nearest(spo2_win, "spo2_time", fio2_win, max_gap_h = 1) %>%
  transmute(
    hospitalization_id,
    arf_onset_time = meas_time_keep,
    criterion = "hypox_spo2_room_air",
    hit = spo2 < 90 & fio2_set <= 0.21 + 1e-6
  ) %>%
  filter(hit)

po2_fio2 <- pair_nearest(po2_win, "po2_time", fio2_win, max_gap_h = 1) %>%
  transmute(
    hospitalization_id,
    arf_onset_time = meas_time_keep,
    criterion = case_when(
      po2 <= 60 & fio2_set <= 0.21 + 1e-6 ~ "hypox_po2_room_air",
      !is.na(po2) & !is.na(fio2_set) & fio2_set > 0 & po2 / fio2_set <= 300 ~ "hypox_pf_ratio",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(criterion))

pco2_win <- labs %>%
  filter(lab_category == "pco2_arterial") %>%
  inner_join(window_tbl, by = "hospitalization_id") %>%
  filter(lab_result_dttm >= arf_window_start, lab_result_dttm <= arf_window_end) %>%
  select(hospitalization_id, pco2_time = lab_result_dttm, pco2 = lab_value_numeric)

ph_win <- labs %>%
  filter(lab_category == "ph_arterial") %>%
  inner_join(window_tbl, by = "hospitalization_id") %>%
  filter(lab_result_dttm >= arf_window_start, lab_result_dttm <= arf_window_end) %>%
  select(hospitalization_id, ph_time = lab_result_dttm, ph = lab_value_numeric)

hyper_hits <- tibble()
if (nrow(pco2_win) && nrow(ph_win)) {
  pdt <- as.data.table(pco2_win)
  hdt <- as.data.table(ph_win)
  pdt[, pco2_time_keep := pco2_time]
  hdt[, ph_time_keep := ph_time]
  setkey(pdt, hospitalization_id, pco2_time)
  setkey(hdt, hospitalization_id, ph_time)
  hyper_hits <- hdt[
    pdt,
    roll = "nearest",
    on = .(hospitalization_id, ph_time = pco2_time),
    nomatch = 0L
  ][
    , time_diff_h := abs(as.numeric(difftime(pco2_time_keep, ph_time_keep, units = "hours")))
  ][
    time_diff_h <= 2 & pco2 >= 45 & ph < 7.35
  ] %>%
    as_tibble() %>%
    transmute(hospitalization_id, arf_onset_time = pco2_time_keep, criterion = "hypercapnic")
}

arf_hits <- bind_rows(spo2_fio2, po2_fio2, hyper_hits)
if (!nrow(arf_hits)) stop("No ARF cases identified with the resubmission ARF onset algorithm.")

arf_onset <- arf_hits %>%
  group_by(hospitalization_id) %>%
  summarise(
    arf_onset = min(arf_onset_time, na.rm = TRUE),
    hypoxemic_arf = any(str_detect(criterion, "^hypox"), na.rm = TRUE),
    hypercapnic_arf = any(criterion == "hypercapnic", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    arf_subtype = case_when(
      hypoxemic_arf & hypercapnic_arf ~ "Mixed",
      hypoxemic_arf ~ "Hypoxemic",
      hypercapnic_arf ~ "Hypercapnic",
      TRUE ~ "Other"
    )
  )

cohort <- base %>%
  inner_join(arf_onset, by = "hospitalization_id") %>%
  filter(!is.na(zipcode_five_digit))

message("Computing Charlson score...")
diagnosis <- read_tbl("hospital_diagnosis") %>%
  transmute(
    hospitalization_id,
    diagnosis_code = as.character(diagnosis_code)
  ) %>%
  left_join(hospitalization %>% select(patient_id, hospitalization_id, admission_dttm), by = "hospitalization_id") %>%
  filter(!is.na(diagnosis_code), !is.na(patient_id))

dx_for_charlson <- diagnosis %>%
  inner_join(
    cohort %>% select(patient_id, index_hospitalization_id = hospitalization_id, arf_onset),
    by = "patient_id",
    relationship = "many-to-many"
  ) %>%
  filter(
    (hospitalization_id == index_hospitalization_id & charlson_include_index) |
      (hospitalization_id != index_hospitalization_id &
         !is.na(admission_dttm) &
         admission_dttm < arf_onset &
         admission_dttm >= arf_onset - days(charlson_lookback_days))
  )

charlson_hosp <- charlson_from_codes(dx_for_charlson) %>%
  group_by(patient_id) %>%
  summarise(
    across(myocardial_infarction:aids_hiv, ~ any(.x, na.rm = TRUE)),
    charlson_score = max(charlson_score, na.rm = TRUE),
    .groups = "drop"
  )

cohort <- cohort %>%
  left_join(charlson_hosp, by = "patient_id") %>%
  mutate(
    across(myocardial_infarction:aids_hiv, ~ coalesce(.x, FALSE)),
    charlson_score = coalesce(charlson_score, 0)
  )

message("Linking ZCTA pollution exposures...")
default_zcta_dir <- file.path(dirname(repo), "CLIF-pollution-microbiome", "data", "exposome_zcta")
default_monthly_no2_dir <- file.path(repo, "data", "resubmission", "no2_zcta_monthly")
default_zcta_acs_path <- file.path(
  dirname(repo),
  "environment_transplant_survival",
  "data",
  "processed",
  "community",
  "zcta_acs_community_covariates_2005_2023.csv.gz"
)
pm25_path <- normalizePath(
  config$zcta_pm25_monthly_path %||% file.path(default_zcta_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet"),
  mustWork = TRUE
)
no2_monthly_dir <- config$zcta_no2_monthly_dir %||% default_monthly_no2_dir
no2_monthly_files <- character()
if (dir.exists(no2_monthly_dir)) {
  no2_monthly_files <- list.files(
    no2_monthly_dir,
    pattern = "^no2_zcta_monthly_[0-9]{4}\\.parquet$",
    full.names = TRUE
  )
}
no2_path <- config$zcta_no2_annual_path %||% file.path(default_zcta_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")
zcta_acs_path <- normalizePath(config$zcta_acs_community_covariates_path %||% default_zcta_acs_path, mustWork = TRUE)

pm25_monthly <- arrow::read_parquet(pm25_path) %>%
  transmute(
    zipcode_five_digit = normalize_zip(zip),
    exposure_month = as.Date(sprintf("%04d-%02d-01", as.integer(year), as.integer(month))),
    pm25_ug_m3 = as.numeric(pm25_ug_m3)
  ) %>%
  filter(!is.na(zipcode_five_digit), !is.na(exposure_month), !is.na(pm25_ug_m3))

pm25_exposure <- cohort %>%
  transmute(
    hospitalization_id,
    zipcode_five_digit,
    exposure_start = floor_date(as.Date(arf_onset), "month") %m-% months(pm25_prior_months),
    exposure_end = floor_date(as.Date(arf_onset), "month") - days(1)
  ) %>%
  left_join(pm25_monthly, by = "zipcode_five_digit", relationship = "many-to-many") %>%
  filter(exposure_month >= exposure_start, exposure_month <= exposure_end) %>%
  group_by(hospitalization_id) %>%
  summarise(
    pm25_12m_zcta = mean(pm25_ug_m3, na.rm = TRUE),
    pm25_months_observed = n_distinct(exposure_month),
    .groups = "drop"
  )

if (length(no2_monthly_files)) {
  message("Using monthly ZCTA NO2 files from: ", normalizePath(no2_monthly_dir, mustWork = TRUE))
  no2_monthly <- purrr::map_dfr(no2_monthly_files, arrow::read_parquet) %>%
    transmute(
      zipcode_five_digit = normalize_zip(zip),
      exposure_month = as.Date(sprintf("%04d-%02d-01", as.integer(year), as.integer(month))),
      no2_ppb = as.numeric(first_existing(pick(everything()), c("no2_ppbv", "no2_ppb", "no2"), default = NA))
    ) %>%
    filter(!is.na(zipcode_five_digit), !is.na(exposure_month), !is.na(no2_ppb)) %>%
    distinct(zipcode_five_digit, exposure_month, .keep_all = TRUE)

  no2_monthly_exposure <- cohort %>%
    transmute(
      hospitalization_id,
      zipcode_five_digit,
      no2_exposure_start = floor_date(as.Date(arf_onset), "month") %m-% months(no2_prior_months),
      no2_exposure_end = floor_date(as.Date(arf_onset), "month") - days(1)
    ) %>%
    left_join(no2_monthly, by = "zipcode_five_digit", relationship = "many-to-many") %>%
    filter(exposure_month >= no2_exposure_start, exposure_month <= no2_exposure_end) %>%
    group_by(hospitalization_id) %>%
    summarise(
      no2_12m_zcta_raw = mean(no2_ppb, na.rm = TRUE),
      no2_months_observed = n_distinct(exposure_month),
      .groups = "drop"
    ) %>%
    mutate(
      no2_monthly_complete_12m = no2_months_observed == no2_prior_months,
      no2_12m_monthly_zcta = if_else(no2_monthly_complete_12m, no2_12m_zcta_raw, NA_real_)
    ) %>%
    select(hospitalization_id, no2_12m_monthly_zcta, no2_months_observed, no2_monthly_complete_12m)

  no2_annual <- arrow::read_parquet(normalizePath(no2_path, mustWork = TRUE)) %>%
    transmute(
      zipcode_five_digit = normalize_zip(zip),
      no2_exposure_year = as.integer(year),
      no2_annual_prior_year_zcta = as.numeric(first_existing(pick(everything()), c("no2", "no2_ppb", "no2_ppbv"), default = NA))
    ) %>%
    filter(!is.na(zipcode_five_digit), !is.na(no2_exposure_year), !is.na(no2_annual_prior_year_zcta)) %>%
    distinct(zipcode_five_digit, no2_exposure_year, .keep_all = TRUE)

  no2_annual_exposure <- cohort %>%
    transmute(
      hospitalization_id,
      zipcode_five_digit,
      no2_exposure_year = year(arf_onset) - no2_lag_years
    ) %>%
    left_join(no2_annual, by = c("zipcode_five_digit", "no2_exposure_year")) %>%
    select(hospitalization_id, no2_annual_prior_year_zcta, no2_exposure_year)

  no2_exposure <- cohort %>%
    select(hospitalization_id) %>%
    left_join(no2_monthly_exposure, by = "hospitalization_id") %>%
    left_join(no2_annual_exposure, by = "hospitalization_id") %>%
    mutate(
      no2_12m_zcta = coalesce(no2_12m_monthly_zcta, no2_annual_prior_year_zcta),
      no2_exposure_source = case_when(
        !is.na(no2_12m_monthly_zcta) ~ "monthly_12m",
        !is.na(no2_annual_prior_year_zcta) ~ "annual_prior_year_fallback",
        TRUE ~ NA_character_
      )
    ) %>%
    select(
      hospitalization_id,
      no2_12m_zcta,
      no2_exposure_source,
      no2_months_observed,
      no2_monthly_complete_12m,
      no2_12m_monthly_zcta,
      no2_annual_prior_year_zcta,
      no2_exposure_year
    )
} else {
  message("Monthly ZCTA NO2 files not found; falling back to annual prior-year NO2: ", no2_path)
  no2_annual <- arrow::read_parquet(normalizePath(no2_path, mustWork = TRUE)) %>%
    transmute(
      zipcode_five_digit = normalize_zip(zip),
      no2_exposure_year = as.integer(year),
      no2_ppb = as.numeric(first_existing(pick(everything()), c("no2", "no2_ppb", "no2_ppbv"), default = NA))
    ) %>%
    filter(!is.na(zipcode_five_digit), !is.na(no2_exposure_year), !is.na(no2_ppb)) %>%
    distinct(zipcode_five_digit, no2_exposure_year, .keep_all = TRUE)

  no2_exposure <- cohort %>%
    transmute(
      hospitalization_id,
      zipcode_five_digit,
      no2_exposure_year = year(arf_onset) - no2_lag_years
    ) %>%
    left_join(no2_annual, by = c("zipcode_five_digit", "no2_exposure_year")) %>%
    transmute(
      hospitalization_id,
      no2_12m_zcta = no2_ppb,
      no2_exposure_source = if_else(!is.na(no2_ppb), "annual_prior_year", NA_character_),
      no2_months_observed = NA_integer_,
      no2_monthly_complete_12m = NA,
      no2_12m_monthly_zcta = NA_real_,
      no2_annual_prior_year_zcta = no2_ppb,
      no2_exposure_year
    )
}

cohort <- cohort %>%
  left_join(pm25_exposure, by = "hospitalization_id") %>%
  left_join(no2_exposure, by = "hospitalization_id")

message("Linking ZCTA ACS social vulnerability covariates...")
zcta_acs <- readr::read_csv(zcta_acs_path, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    zipcode_five_digit = normalize_zip(zip),
    index_year = as.integer(analysis_year),
    acs_community_year = as.integer(community_year),
    acs_total_population = as.numeric(total_population),
    acs_median_household_income = as.numeric(median_household_income),
    acs_pct_poverty = as.numeric(pct_poverty),
    acs_pct_bachelor_plus = as.numeric(pct_bachelor_plus),
    acs_pct_unemployed = as.numeric(pct_unemployed),
    acs_pct_no_vehicle = as.numeric(pct_no_vehicle),
    acs_pct_nonwhite = as.numeric(pct_nonwhite)
  ) %>%
  mutate(
    acs_median_household_income = if_else(acs_median_household_income < 0, NA_real_, acs_median_household_income),
    across(
      c(acs_pct_poverty, acs_pct_bachelor_plus, acs_pct_unemployed, acs_pct_no_vehicle, acs_pct_nonwhite),
      ~ if_else(.x < 0 | .x > 1, NA_real_, .x)
    )
  ) %>%
  filter(!is.na(zipcode_five_digit), !is.na(index_year)) %>%
  distinct(zipcode_five_digit, index_year, .keep_all = TRUE) %>%
  group_by(index_year) %>%
  mutate(
    across(
      c(
        acs_median_household_income,
        acs_pct_poverty,
        acs_pct_bachelor_plus,
        acs_pct_unemployed,
        acs_pct_no_vehicle,
        acs_pct_nonwhite
      ),
      ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x)
    ),
    acs_svi_poverty_rank = percent_rank(acs_pct_poverty),
    acs_svi_unemployed_rank = percent_rank(acs_pct_unemployed),
    acs_svi_no_vehicle_rank = percent_rank(acs_pct_no_vehicle),
    acs_svi_nonwhite_rank = percent_rank(acs_pct_nonwhite),
    acs_svi_low_income_rank = percent_rank(-acs_median_household_income),
    acs_svi_low_education_rank = percent_rank(-acs_pct_bachelor_plus),
    acs_zcta_svi_proxy = rowMeans(
      cbind(
        acs_svi_poverty_rank,
        acs_svi_unemployed_rank,
        acs_svi_no_vehicle_rank,
        acs_svi_nonwhite_rank,
        acs_svi_low_income_rank,
        acs_svi_low_education_rank
      ),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  select(
    zipcode_five_digit,
    index_year,
    acs_community_year,
    acs_total_population,
    acs_median_household_income,
    acs_pct_poverty,
    acs_pct_bachelor_plus,
    acs_pct_unemployed,
    acs_pct_no_vehicle,
    acs_pct_nonwhite,
    acs_zcta_svi_proxy
  )

cohort <- cohort %>%
  left_join(zcta_acs, by = c("zipcode_five_digit", "index_year"))

message("Constructing post-ARF competing-risk outcomes...")
imv_records <- resp_support %>%
  inner_join(cohort %>% select(hospitalization_id, arf_onset, last_icu_out, discharge_dttm), by = "hospitalization_id") %>%
  filter(recorded_dttm >= arf_onset) %>%
  mutate(is_imv = str_detect(device_category %||% "", "imv|invasive")) %>%
  arrange(hospitalization_id, recorded_dttm)

imv_runs <- imv_records %>%
  filter(is_imv) %>%
  group_by(hospitalization_id) %>%
  arrange(recorded_dttm, .by_group = TRUE) %>%
  mutate(
    gap_h = as.numeric(difftime(recorded_dttm, lag(recorded_dttm), units = "hours")),
    run_id = cumsum(if_else(is.na(gap_h) | gap_h > imv_gap_hours, 1L, 0L))
  ) %>%
  group_by(hospitalization_id, run_id) %>%
  summarise(
    imv_start = min(recorded_dttm),
    imv_end = max(recorded_dttm),
    .groups = "drop"
  )

death_times <- cohort %>%
  mutate(
    death_time = case_when(
      !is.na(death_dttm) & death_dttm >= arf_onset ~ death_dttm,
      str_detect(discharge_category %||% "", "expired|death|hospice") ~ discharge_dttm,
      TRUE ~ as.POSIXct(NA, tz = time_zone)
    )
  ) %>%
  select(hospitalization_id, death_time)

event_data <- cohort %>%
  select(hospitalization_id, arf_onset, last_icu_out, discharge_dttm) %>%
  left_join(death_times, by = "hospitalization_id") %>%
  mutate(
    administrative_censor = arf_onset + days(followup_days),
    censor_time = pmin(last_icu_out, discharge_dttm, administrative_censor, na.rm = TRUE)
  )

extubation_events <- imv_runs %>%
  left_join(event_data, by = "hospitalization_id") %>%
  filter(imv_start >= arf_onset) %>%
  group_by(hospitalization_id) %>%
  arrange(imv_start, .by_group = TRUE) %>%
  mutate(next_imv_start = lead(imv_start)) %>%
  filter(
    is.na(next_imv_start) | as.numeric(difftime(next_imv_start, imv_end, units = "hours")) >= successful_extubation_hours
  ) %>%
  filter(imv_end + hours(successful_extubation_hours) <= censor_time) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(hospitalization_id, extubation_time = imv_end)

persistent_rf_events <- imv_runs %>%
  left_join(event_data, by = "hospitalization_id") %>%
  group_by(hospitalization_id) %>%
  slice_max(imv_end, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(last_icu_out), imv_end >= last_icu_out - hours(imv_gap_hours), last_icu_out <= censor_time) %>%
  transmute(hospitalization_id, persistent_rf_time = last_icu_out)

events <- event_data %>%
  left_join(extubation_events, by = "hospitalization_id") %>%
  left_join(persistent_rf_events, by = "hospitalization_id") %>%
  rowwise() %>%
  mutate(
    event_times = list(c(death = death_time, extubation = extubation_time, persistent_rf = persistent_rf_time)),
    event_times = list(event_times[!is.na(event_times) & event_times <= censor_time]),
    event_name = if (length(event_times) == 0) "censor" else names(event_times)[which.min(event_times)],
    event_time = if (length(event_times) == 0) censor_time else min(event_times)
  ) %>%
  ungroup() %>%
  transmute(
    hospitalization_id,
    event_name,
    event_time = as.POSIXct(event_time, origin = "1970-01-01", tz = time_zone),
    ftime_days = as.numeric(difftime(event_time, arf_onset, units = "days")),
    event_code = case_when(
      event_name == "extubation" ~ 1L,
      event_name == "death" ~ 2L,
      event_name == "persistent_rf" ~ 3L,
      TRUE ~ 0L
    )
  ) %>%
  mutate(ftime_days = pmax(ftime_days, 1 / 24))

analysis_df <- cohort %>%
  left_join(events, by = "hospitalization_id") %>%
  mutate(
    site = site_name,
    pm25_q = make_quartile(pm25_12m_zcta),
    no2_q = make_quartile(no2_12m_zcta),
    pm25_per_5 = pm25_12m_zcta / 5,
    no2_per_10 = no2_12m_zcta / 10,
    age_10 = age / 10,
    acs_median_household_income_10k = acs_median_household_income / 10000,
    index_year_f = factor(index_year),
    sex = factor(sex),
    race_ethnicity = forcats::fct_lump_min(factor(race_ethnicity), min = 100, other_level = "Other/Unknown"),
    arf_subtype = factor(arf_subtype)
  ) %>%
  filter(!is.na(ftime_days), is.finite(ftime_days), ftime_days > 0)

readr::write_csv(analysis_df, file.path(out_root, "resubmission_analysis_dataset.csv"))

cohort_summary <- tibble(
  site = site_name,
  n_arf = nrow(analysis_df),
  n_patients = n_distinct(analysis_df$patient_id),
  n_with_pm25 = sum(!is.na(analysis_df$pm25_12m_zcta)),
  n_with_no2 = sum(!is.na(analysis_df$no2_12m_zcta)),
  n_with_zcta_acs = sum(complete.cases(analysis_df[, c(
    "acs_pct_poverty",
    "acs_pct_unemployed",
    "acs_pct_no_vehicle",
    "acs_pct_nonwhite",
    "acs_median_household_income_10k",
    "acs_pct_bachelor_plus"
  )])),
  n_death = sum(analysis_df$event_code == 2, na.rm = TRUE),
  n_extubation = sum(analysis_df$event_code == 1, na.rm = TRUE),
  n_persistent_rf = sum(analysis_df$event_code == 3, na.rm = TRUE),
  n_censored = sum(analysis_df$event_code == 0, na.rm = TRUE),
  mean_charlson = mean(analysis_df$charlson_score, na.rm = TRUE),
  median_charlson = median(analysis_df$charlson_score, na.rm = TRUE)
)
readr::write_csv(cohort_summary, file.path(out_root, "resubmission_cohort_summary.csv"))

message("Fitting cause-specific Cox models...")
fit_cause_cox <- function(data, cause_code, cause_label, exposure_terms, model_label) {
  social_covars <- c(
    "acs_pct_poverty",
    "acs_pct_unemployed",
    "acs_pct_no_vehicle",
    "acs_pct_nonwhite",
    "acs_median_household_income_10k",
    "acs_pct_bachelor_plus"
  )
  covars <- c(
    exposure_terms,
    "age_10",
    "sex",
    "race_ethnicity",
    "charlson_score",
    "arf_subtype",
    "index_year_f",
    social_covars
  )
  model_df <- data %>%
    select(ftime_days, event_code, all_of(covars)) %>%
    drop_na()
  if (!nrow(model_df) || sum(model_df$event_code == cause_code) < 10) return(tibble())
  fml <- as.formula(paste0("Surv(ftime_days, event_code == ", cause_code, ") ~ ", paste(covars, collapse = " + ")))
  fit <- survival::coxph(fml, data = model_df, ties = "efron")
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      site = site_name,
      model = model_label,
      cause = cause_label,
      term,
      n = nrow(model_df),
      events = sum(model_df$event_code == cause_code),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      adjustment_set = paste(
        "age_10",
        "sex",
        "race_ethnicity",
        "charlson_score",
        "arf_subtype",
        "index_year_f",
        paste(social_covars, collapse = " + "),
        sep = " + "
      )
    )
}

cox_results <- bind_rows(
  fit_cause_cox(analysis_df, 1, "Successful extubation", "pm25_per_5", "PM25 single-pollutant"),
  fit_cause_cox(analysis_df, 2, "Death", "pm25_per_5", "PM25 single-pollutant"),
  fit_cause_cox(analysis_df, 3, "Persistent respiratory failure", "pm25_per_5", "PM25 single-pollutant"),
  fit_cause_cox(analysis_df, 1, "Successful extubation", "no2_per_10", "NO2 single-pollutant"),
  fit_cause_cox(analysis_df, 2, "Death", "no2_per_10", "NO2 single-pollutant"),
  fit_cause_cox(analysis_df, 3, "Persistent respiratory failure", "no2_per_10", "NO2 single-pollutant"),
  fit_cause_cox(analysis_df, 1, "Successful extubation", c("pm25_per_5", "no2_per_10"), "PM25 + NO2"),
  fit_cause_cox(analysis_df, 2, "Death", c("pm25_per_5", "no2_per_10"), "PM25 + NO2"),
  fit_cause_cox(analysis_df, 3, "Persistent respiratory failure", c("pm25_per_5", "no2_per_10"), "PM25 + NO2")
)
readr::write_csv(cox_results, file.path(out_root, "resubmission_cause_specific_cox_results.csv"))

message("Building unadjusted Aalen-Johansen CIFs...")
make_cif_data <- function(data, exposure_group, exposure_label) {
  dat <- data %>%
    filter(!is.na(.data[[exposure_group]])) %>%
    mutate(
      status = factor(
        event_code,
        levels = c(0, 1, 2, 3),
        labels = c("censor", "Successful extubation", "Death", "Persistent respiratory failure")
      ),
      exposure_group = .data[[exposure_group]]
    )
  if (!nrow(dat)) return(tibble())
  fit <- survival::survfit(Surv(ftime_days, status) ~ exposure_group, data = dat)
  ss <- summary(fit)
  pstate <- as_tibble(ss$pstate, .name_repair = "minimal")
  se <- as_tibble(ss$std.err, .name_repair = "minimal")
  names(pstate) <- ss$states
  names(se) <- paste0(ss$states, "__se")
  bind_cols(
    tibble(exposure = exposure_label, time = ss$time, strata = ss$strata),
    pstate,
    se
  ) %>%
    pivot_longer(cols = all_of(ss$states), names_to = "state", values_to = "cif") %>%
    pivot_longer(cols = ends_with("__se"), names_to = "state_se", values_to = "std_error") %>%
    mutate(state_se = str_remove(state_se, "__se$")) %>%
    filter(state == state_se) %>%
    select(-state_se) %>%
    mutate(
      quartile = str_replace(strata, "^exposure_group=", ""),
      conf_low = pmax(0, cif - 1.96 * std_error),
      conf_high = pmin(1, cif + 1.96 * std_error)
    ) %>%
    filter(!state %in% c("censor", "(s0)"))
}

cif_plot_data <- bind_rows(
  make_cif_data(analysis_df, "pm25_q", "PM2.5"),
  make_cif_data(analysis_df, "no2_q", "NO2")
)
readr::write_csv(cif_plot_data, file.path(out_root, "resubmission_aalen_johansen_cif_plot_data.csv"))

if (nrow(cif_plot_data)) {
  cif_plot_data <- cif_plot_data %>%
    mutate(
      quartile = recode(
        as.character(quartile),
        "Q1" = "Q1 lowest",
        "Q2" = "Q2",
        "Q3" = "Q3",
        "Q4" = "Q4 highest"
      ),
      quartile = factor(quartile, levels = quartile_labels),
      state = recode(
        as.character(state),
        "Successful extubation" = "Extubation",
        "Persistent respiratory failure" = "Persistent RF"
      ),
      state = factor(state, levels = c("Extubation", "Death", "Persistent RF"))
    )

  p_cif <- ggplot(cif_plot_data, aes(x = time, y = cif, color = quartile)) +
    geom_step(linewidth = 0.86) +
    facet_grid(state ~ exposure, scales = "free_y") +
    scale_color_manual(values = quartile_colors, drop = FALSE) +
    scale_x_continuous(
      breaks = seq(0, followup_days, length.out = 5),
      labels = scales::label_number(accuracy = 0.1, trim = TRUE),
      limits = c(0, followup_days),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(
      title = "Aalen-Johansen Cumulative Incidence by Exposure Quartile",
      x = "Days since ARF onset",
      y = "Cumulative incidence",
      color = "Quartile"
    ) +
    theme_cif_transplant() +
    theme(plot.title = element_text(hjust = 0, margin = margin(b = 8)))
  ggsave(file.path(fig_dir, "resubmission_aalen_johansen_cif_quartiles.png"), p_cif, width = 11, height = 8, dpi = 300)
  ggsave(file.path(fig_dir, "resubmission_aalen_johansen_cif_quartiles.pdf"), p_cif, width = 11, height = 8)
}

message("Done.")
message("Wrote outputs to: ", out_root)
print(cohort_summary)
print(cox_results, n = nrow(cox_results))
