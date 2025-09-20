# ===============================================================================================
# ICU Respiratory Failure Environmental Risk (REFER) Index
# PI: Peter Graffy (graffy@uchicago.edu)
# Run this script after running: 01_REFER_cohort_identification.R
# Purpose: Build cohorts, history, outcomes; link exposome; fit models; save outputs consistently
# ===============================================================================================

# ------------------------------------ 0) Config & Setup ----------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(MASS)       
  library(broom)
  library(patchwork)
  library(ggplot2)
  library(ggeffects)
  library(gt)
  library(purrr)
  library(stringr)
  library(rlang)
  library(tidycensus)
  library(fixest)
  library(marginaleffects)
  library(pscl)
  library(glmmTMB)
  library(digest)
  library(pROC)
  library(glue)
  library(tibble)
  library(forcats)
  library(grid)
  library(jsonlite)
})

# ---- 0.1 Load config (YAML or fallback defaults) ----------------------------------------------
keep_vars <- c("clif_tables", "cohort_min", "cohort_min_periop", "repo")
rm(list = setdiff(ls(envir = .GlobalEnv), keep_vars), envir = .GlobalEnv)

#setwd("~/CLIF-ARFVI") #<------ set your working directory to the main folder path here if needed

# ------------------------------- 0) Config & Setup --------------------------------
# --- Repo-anchored config + path helpers ---
stopifnot(exists("repo"))  # e.g., repo <- "/path/to/your/repo"
rpath <- function(...) file.path(repo, ...)         # join to repo root

source(rpath("utils", "config.R"))

`%||%` <- function(x, y) if (!is.null(x)) x else y  # null-coalesce helper

# print a few key fields the user asked for
site_name   <- config$site_name   %||% "unknown_site"
tables_path <- config$tables_path %||% "data/"
file_type   <- config$file_type   %||% "parquet"

print(paste("Site Name:", site_name))
print(paste("Tables Path:", tables_path))
print(paste("File Type:", file_type))

# build runtime cfg from config + defaults
cfg <- list(
  project_code = config$project_code %||% "refer",
  version      = config$version      %||% "v0_1",
  site_name    = config$site_name    %||% "site",
  date_stamp   = format(Sys.Date(), "%Y%m%d"),
  run_id       = format(Sys.time(), "%Y%m%d_%H%M%S")
)

# directory *names* (keep overridable), but full *paths* anchored to repo
output_basename  <- config$output_dir  %||% "output"
figures_basename <- config$figures_dir %||% "figures"

cfg$output_dir   <- rpath(output_basename)                   # e.g., <repo>/output
cfg$figures_dir  <- figures_basename                         # name only
cfg$figures_path <- rpath(output_basename, figures_basename) # e.g., <repo>/output/figures

# prefix like: refer_sitename_20250908
cfg$prefix <- paste(cfg$project_code, cfg$site_name, cfg$date_stamp, sep = "_")

# ensure output folders exist
dir.create(cfg$output_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$figures_path, recursive = TRUE, showWarnings = FALSE)

# ---- Save helpers (always under repo/output) ----
.path_csv <- function(name)
  file.path(cfg$output_dir, paste0(cfg$prefix, "_", basename(name), "_", cfg$run_id, ".csv"))

.path_rds <- function(name)
  file.path(cfg$output_dir, paste0(cfg$prefix, "_", basename(name), "_", cfg$run_id, ".rds"))

.path_png <- function(name)
  file.path(cfg$figures_path, paste0(cfg$prefix, "_", basename(name), "_", cfg$run_id, ".png"))

.path_pdf <- function(name)
  file.path(cfg$figures_path, paste0(cfg$prefix, "_", basename(name), "_", cfg$run_id, ".pdf"))

# If 'name' includes subdirs (e.g., "models/inhosp"), create them under output/ and output/figures/
.ensure_subdir <- function(name, under) {
  subdir <- dirname(name)
  if (!identical(subdir, "."))
    dir.create(file.path(under, subdir), recursive = TRUE, showWarnings = FALSE)
}

save_tbl <- function(x, name) {
  .ensure_subdir(name, cfg$output_dir)
  readr::write_csv(x, file.path(cfg$output_dir, dirname(name),
                                paste0(cfg$prefix, "_", basename(name), "_", cfg$run_id, ".csv")))
}

save_model <- function(fit, name) {
  .ensure_subdir(name, cfg$output_dir)
  saveRDS(fit, file.path(cfg$output_dir, dirname(name),
                         paste0(cfg$prefix, "_", basename(name), "_", cfg$run_id, ".rds")))
}

save_plot <- function(p, name, w = 11, h = 9, dpi = 300) {
  .ensure_subdir(name, cfg$figures_path)
  ggsave(file.path(cfg$figures_path, dirname(name),
                   paste0(cfg$prefix, "_", basename(name), "_", cfg$run_id, ".png")),
         p, width = w, height = h, dpi = dpi)
  ggsave(file.path(cfg$figures_path, dirname(name),
                   paste0(cfg$prefix, "_", basename(name), "_", cfg$run_id, ".pdf")),
         p, width = w, height = h)
}

# ------------------------------------ 1) Table Access & Helpers ---------------------------------
get_tbl <- function(nm) {
  ct <- get0("clif_tables", inherits = TRUE)
  if (is.null(ct)) stop("Couldn't find 'clif_tables' in your environment.")
  key <- if (nm %in% names(ct)) nm else {
    ci <- names(ct)[tolower(names(ct)) == tolower(nm)]
    if (length(ci) == 1) ci else nm
  }
  if (!key %in% names(ct)) stop(sprintf("Table '%s' not in clif_tables. Available: %s", nm, paste(names(ct), collapse = ", ")))
  janitor::clean_names(ct[[key]])
}

pick_col <- function(df, candidates, required = TRUE) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  if (required) stop(sprintf("None of these columns found: %s", paste(candidates, collapse = ", ")))
  rep(NA, nrow(df))
}

coalesce_any <- function(data, candidates) {
  cols <- dplyr::select(data, dplyr::any_of(candidates))
  if (ncol(cols) == 0) return(rep(NA_character_, nrow(data)))
  dplyr::coalesce(!!!cols)
}

safe_ts <- function(x, tz = "America/Chicago") {
  if (inherits(x, "POSIXt")) return(x)
  if (is.numeric(x)) {
    x2 <- ifelse(x > 1e12, x/1000, x)
    return(lubridate::as_datetime(x2, tz = tz))
  }
  suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c("ymd_HMS","ymd_HM","ymd","ymdTz","ymdT","mdy_HMS","mdy_HM","mdy","dmy_HMS","dmy_HM","dmy","HMS"),
    tz = tz, quiet = TRUE
  ))
}

add_index_fields <- function(df) {
  if (!all(c("admission_dttm","discharge_dttm") %in% names(df))) {
    df <- df |>
      dplyr::left_join(
        hospitalization |>
          dplyr::select(patient_id, hospitalization_id,
                        admission_dttm, discharge_dttm,
                        admitting_service, discharge_service, zip_code),
        by = c("patient_id","hospitalization_id")
      )
  }
  admit_raw_vec <- coalesce_any(df, c("admission_dttm","admit_dttm","admit_time","admission_time"))
  disch_raw_vec <- coalesce_any(df, c("discharge_dttm","discharge_time","dc_dttm","disposition_time"))
  df |>
    dplyr::mutate(
      index_admit     = safe_ts(admit_raw_vec),
      index_discharge = safe_ts(disch_raw_vec),
      index_year      = lubridate::year(index_admit),
      index_date      = as.Date(index_admit)
    )
}

# Core CLIF tables
patient         <- get_tbl("clif_patient")
hospitalization <- get_tbl("clif_hospitalization")
diagnosis       <- get_tbl("clif_hospital_diagnosis")
support         <- get_tbl("clif_respiratory_support")
med_admin       <- get_tbl("clif_medication_admin_continuous")
icu_stay        <- get_tbl("clif_adt")
vitals          <- get_tbl("clif_vitals")
labs_df         <- get_tbl("clif_labs")   # avoid name clash with ggplot2::labs()

# ------------------------------------ 2) Build Cohorts ------------------------------------------
arf_idx    <- cohort_min        |> clean_names() |> add_index_fields() |> mutate(cohort = "ARF")
periop_idx <- cohort_min_periop |> clean_names() |> add_index_fields() |> mutate(cohort = "PERIOP")
cohort_all <- bind_rows(arf_idx, periop_idx) |> filter(!is.na(index_admit))
#save_tbl(cohort_all, "cohort_all")

# ------------------------------------ 3) Lookback / History -------------------------------------
lookback_days <- 365
cohort_lb <- cohort_all |>
  transmute(patient_id, hospitalization_id, index_admit, index_year, cohort,
            lb_start = index_admit - days(lookback_days))

adt_tmp <- icu_stay |>
  dplyr::left_join(
    hospitalization |> dplyr::select(hospitalization_id, patient_id),
    by = "hospitalization_id"
  )

icu_segs <- adt_tmp |>
  mutate(
    in_raw  = pick_col(adt_tmp,  c("in_dttm","adt_in_dttm","icu_in","arrival_dttm","in_time")),
    out_raw = pick_col(adt_tmp,  c("out_dttm","adt_out_dttm","icu_out","departure_dttm","out_time")),
    in_ts   = safe_ts(in_raw),
    out_ts  = safe_ts(out_raw),
    loccat  = tolower(pick_col(adt_tmp, c("location_category","loc_category","unit_category")))
  ) |>
  filter(str_detect(loccat, "\\bicu\\b|\\bccu\\b|\\bmicu\\b|\\bsicu\\b|\\bcicu\\b|\\bnicu\\b|\\bpicu\\b")) |>
  filter(!is.na(patient_id), !is.na(in_ts), !is.na(out_ts), out_ts > in_ts)

# Prior ICU stays
icu_hist <- icu_segs |>
  inner_join(cohort_lb, by = join_by(patient_id), relationship = "many-to-many") |>
  filter(out_ts < index_admit, out_ts >= lb_start, hospitalization_id.x != hospitalization_id.y) |>
  group_by(patient_id, hospitalization_id.y) |>
  summarise(prior_icu_stays = n_distinct(hospitalization_id.x), .groups = "drop") |>
  rename(hospitalization_id = hospitalization_id.y)

# Baseline meds in lookback
med_tmp <- med_admin |>
  dplyr::left_join(
    hospitalization |> dplyr::select(hospitalization_id, patient_id),
    by = "hospitalization_id"
  )

baseline_meds <- med_tmp |>
  mutate(
    admin_ts = safe_ts(pick_col(med_tmp, c("admin_dttm","start_dttm","start_time","start_ts"), TRUE)),
    med_low  = tolower(pick_col(med_tmp, c("medication_name","medication","drug_name","med_name"), TRUE))
  ) |>
  inner_join(cohort_lb, by = "patient_id", relationship = "many-to-many") |>
  filter(!is.na(admin_ts), admin_ts >= lb_start, admin_ts < index_admit) |>
  mutate(
    is_acei_arb = str_detect(med_low, "lisinopril|losartan|valsartan|enalapril|olmesartan|candesartan"),
    is_diuretic = str_detect(med_low, "furosemide|bumetanide|torsemide|hydrochlorothiazide"),
    is_bb       = str_detect(med_low, "metoprolol|carvedilol|atenolol|propranolol")
  ) |>
  group_by(patient_id, hospitalization_id.y) |>
  summarise(
    any_acei_arb = as.integer(any(is_acei_arb, na.rm = TRUE)),
    any_diuretic = as.integer(any(is_diuretic, na.rm = TRUE)),
    any_bb       = as.integer(any(is_bb, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  rename(hospitalization_id = hospitalization_id.y)

history_features <- cohort_all |>
  dplyr::select(patient_id, hospitalization_id) |>
  dplyr::distinct() |>
  dplyr::left_join(icu_hist,      by = c("patient_id","hospitalization_id")) |>
  dplyr::left_join(baseline_meds, by = c("patient_id","hospitalization_id")) |>
  dplyr::mutate(dplyr::across(c(prior_icu_stays, any_acei_arb, any_diuretic, any_bb), ~tidyr::replace_na(., 0)))

#save_tbl(history_features, "history_features")

# ------------------------------------ 4) Outcomes ------------------------------------------------
# ICU LOS
icu_los <- icu_segs |>
  semi_join(cohort_all, by = "hospitalization_id") |>
  mutate(seg_days = as.numeric(difftime(out_ts, in_ts, units = "days"))) |>
  group_by(hospitalization_id) |>
  summarise(icu_los_days = sum(seg_days, na.rm = TRUE), .groups = "drop")

# Hospital LOS
hosp_los <- cohort_all |>
  transmute(hospitalization_id,
            hosp_los_days = as.numeric(difftime(index_discharge, index_admit, units = "days"))) |>
  left_join(hospitalization |> dplyr::select(hospitalization_id, discharge_category, county_code),
            by = "hospitalization_id")

# Mortality
mortality_instay <- cohort_all |>
  left_join(patient |> dplyr::select(patient_id, death_dttm), by = "patient_id") |>
  mutate(
    death_ts      = safe_ts(death_dttm),
    in_hosp_death = as.integer(!is.na(death_ts) & death_ts >= index_admit & death_ts <= index_discharge),
    death_30d     = as.integer(!is.na(death_ts) & death_ts <= (index_admit + days(30)))
  ) |>
  dplyr::select(hospitalization_id, in_hosp_death, death_30d)

# Vent flag + durations
vent_flag <- support |>
  mutate(dev_low = tolower(coalesce(device_name, ""))) |>
  filter(str_detect(dev_low, "\\bvent\\b|vent;|vent ")) |>
  filter(!str_detect(dev_low, "bipap|cpap|high flow|nasal cannula|\\bnc\\b|trach collar|oxytrach|room air")) |>
  semi_join(cohort_all, by = "hospitalization_id") |>
  distinct(hospitalization_id) |>
  mutate(vent_proc_flag = 1L)

support_tmp <- support |>
  left_join(hospitalization |> dplyr::select(hospitalization_id, patient_id), by = "hospitalization_id") |>
  mutate(
    rec_time = safe_ts(coalesce(.data$recorded_time, .data$recorded_dttm)),
    dev_low  = tolower(coalesce(.data$device_name, ""))
  ) |>
  filter(!is.na(rec_time)) |>
  semi_join(cohort_all, by = "hospitalization_id")

support_class <- support_tmp |>
  mutate(
    is_niv = str_detect(dev_low, "bipap|cpap|high flow|hf vent|nasal cannula|\\bnc\\b|venturi|face mask|face tent|trach collar|oxytrach|room air|t-piece|ram cannula|aerosol mask|o2 hood"),
    has_vent_token = str_detect(dev_low, "(^|[ ;])vent([ ;]|$)"),
    is_invasive_vent = has_vent_token & !is_niv
  )

gap_hours <- 6
vent_durations <- support_class |>
  arrange(hospitalization_id, rec_time) |>
  group_by(hospitalization_id) |>
  mutate(
    next_time   = lead(rec_time),
    next_invas  = lead(is_invasive_vent),
    gap_hr      = as.numeric(difftime(next_time, rec_time, units = "hours")),
    add_hours   = if_else(is_invasive_vent & next_invas & !is.na(gap_hr) & gap_hr > 0 & gap_hr <= gap_hours, gap_hr, 0),
    next_niv    = lead(is_niv),
    add_niv_hrs = if_else(is_niv & next_niv & !is.na(gap_hr) & gap_hr > 0 & gap_hr <= gap_hours, gap_hr, 0)
  ) |>
  summarise(vent_hours = sum(add_hours, na.rm = TRUE),
            niv_hours  = sum(add_niv_hrs, na.rm = TRUE), .groups = "drop") |>
  mutate(vent_proc_flag = as.integer(vent_hours > 0))

# AKI (creatinine swing)
aki_flag <- labs_df |>
  mutate(name_low = tolower(pick_col(labs_df, c("lab_name","test_name","component","loinc_name")))) |>
  filter(str_detect(name_low, "creatinine")) |>
  semi_join(cohort_all, by = "hospitalization_id") |>
  group_by(hospitalization_id) |>
  summarise(aki_flag = as.integer((max(lab_value_numeric, na.rm = TRUE) - min(lab_value_numeric, na.rm = TRUE)) >= 0.3),
            .groups = "drop")

vaso_flag <- med_admin |>
  dplyr::mutate(med_low = tolower(pick_col(med_admin, c("medication_name","medication","drug_name","med_name")))) |>
  dplyr::filter(stringr::str_detect(med_low, "norepinephrine|epinephrine|phenylephrine|vasopressin|dopamine")) |>
  dplyr::semi_join(cohort_all, by="hospitalization_id") |>
  dplyr::distinct(hospitalization_id) |>
  dplyr::mutate(vaso_flag = 1L)

# Assemble outcomes
outcomes <- cohort_all |>
  dplyr::select(patient_id, hospitalization_id, cohort, index_admit, index_discharge, index_year, census_tract, county_code) |>
  left_join(icu_los,  by = "hospitalization_id") |>
  left_join(hosp_los, by = "hospitalization_id") |>
  left_join(mortality_instay, by = "hospitalization_id") |>
  left_join(vaso_flag, by = "hospitalization_id") |>
  left_join(vent_flag, by = "hospitalization_id") |>
  left_join(vent_durations, by = "hospitalization_id") |>
  left_join(aki_flag,  by = "hospitalization_id") |>
  mutate(across(c(aki_flag, in_hosp_death, death_30d), ~ replace_na(., 0L)),
         vent_hours = coalesce(vent_hours, 0), niv_hours = coalesce(niv_hours, 0))

#save_tbl(outcomes, "outcomes")

# ------------------------------------ 5) Link Exposome ------------------------------------------
outcomes <- outcomes |>
  mutate(fips_county = str_pad(as.character(county_code.x), width = 5, pad = "0"),
         GEOID = fips_county)

svi    <- readr::read_csv(rpath("exposome", "SVI_county_year.csv"))
pm25   <- readr::read_csv(rpath("exposome", "pm25_county_year.csv"))
no2    <- readr::read_csv(rpath("exposome", "no2_county_year.csv"))
daymet <- readr::read_csv(rpath("exposome", "daymet_county_year_allvars.csv"))

exposome <- svi |>
  left_join(pm25,   by = c("GEOID","year")) |>
  left_join(no2,    by = c("GEOID","year")) |>
  left_join(daymet, by = c("GEOID","year"))

# --- NEW: add tract-level ACS ----------------------------------------------

# Path to the zip in exposome/
zip_path <- rpath("exposome", "acs_estimates.csv.zip")  # change name if needed

# Inspect ZIP contents and pick the CSV you want
zip_files <- unzip(zip_path, list = TRUE)$Name
csv_inside <- zip_files[grepl("\\.csv$", zip_files, ignore.case = TRUE)]

if (length(csv_inside) == 0) stop("No CSV found inside the ZIP.")
if (length(csv_inside) > 1) {
  # prefer a file that looks like ACS, otherwise take the first CSV
  pick <- csv_inside[grepl("acs|estimate", csv_inside, ignore.case = TRUE)]
  csv_inside <- if (length(pick)) pick[1] else csv_inside[1]
}

# Read directly via a connection without extracting
acs <- readr::read_csv(unz(zip_path, csv_inside)) |>
  transmute(
    year  = as.integer(year),
    geoid = str_pad(as.character(geoid), width = 11, pad = "0"),
    median_income,
    pov_rate_pct,
    pct_lt_hs,
    unemp_rate_pct,
    pct_insured
  ) |>
  distinct(geoid, year, .keep_all = TRUE)

# Sanity: ensure year for linkage is integer
outcomes <- outcomes |>
  mutate(index_year = as.integer(index_year))

outcomes_exp <- outcomes |>
  left_join(exposome, by = c("county_code.x" = "GEOID","index_year" = "year"))

# Join tract ACS (adds 4 vars with acs_ prefix)
outcomes_exp <- outcomes_exp |>
  left_join(
    acs |> rename_with(~ paste0("acs_", .x), -c(geoid, year)),
    by = c("census_tract" = "geoid", "index_year" = "year")
  )

#save_tbl(outcomes_exp, "outcomes_exposome")

# ------------------------------------ 6) ARF Analytic Frame + Demographics ----------------------
arf_exp <- outcomes_exp |> filter(cohort == "ARF") |>
  left_join(patient |> dplyr::select(patient_id, race_name, ethnicity_category, sex_category, birth_date),
            by = "patient_id") |>
  mutate(
    race_ethnicity = case_when(
      str_to_lower(ethnicity_category) %in% c("hispanic", "latino", "latinx") ~ paste0("Hispanic ", race_name),
      TRUE ~ paste0("Non-Hispanic ", race_name)
    ),
    age = as.numeric(difftime(index_admit, birth_date, units = "days")) / 365.25
  ) |>
  mutate(
    re_low = str_to_lower(race_ethnicity),
    is_nonhisp = str_detect(re_low, "\\bnon[- ]?hispanic\\b"),
    is_hisp    = str_detect(re_low, "\\bhispanic\\b") & !is_nonhisp,
    is_white   = str_detect(re_low, "\\bwhite\\b"),
    is_black   = str_detect(re_low, "black"),
    is_asian_any = str_detect(re_low, "asian|mideast|filipino|chinese|korean|vietnamese|pacific islander|samoan"),
    race_ethnicity_simple = case_when(
      is_white & is_hisp    ~ "Hispanic White",
      is_white & is_nonhisp ~ "Non-Hispanic White",
      is_black & is_hisp    ~ "Hispanic Black",
      is_black & is_nonhisp ~ "Non-Hispanic Black",
      is_asian_any          ~ "Asian",
      TRUE                  ~ "Other"
    )
  ) |>
  dplyr::select(-re_low, -is_nonhisp, -is_hisp, -is_white, -is_black, -is_asian_any) |>
  mutate(
    sex_category = factor(sex_category),
    race_ethnicity_simple = factor(race_ethnicity_simple,
                                   levels = c("Non-Hispanic White","Hispanic White","Non-Hispanic Black","Hispanic Black","Asian","Other"))
  )

#save_tbl(arf_exp, "arf_exp")

# ------------------------------------ 7) Models (Adjusted) --------------------------------------
# Helper to tidy-save any model
tidy_and_save <- function(fit, name, exponentiate = FALSE, folder = "models") {
  tt <- broom::tidy(fit, exponentiate = exponentiate, conf.int = TRUE)
  
  # subfolder under <repo>/output
  out_dir <- file.path(cfg$output_dir, folder)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  sanitize_filename <- function(x) {
    x <- gsub("[^[:alnum:]_-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    tolower(x)
  }
  clean_name <- sanitize_filename(name)
  
  out_file <- file.path(out_dir, paste0(cfg$prefix, "_model_tidy_", clean_name, "_", cfg$run_id, ".csv"))
  readr::write_csv(tt, out_file)
  tt
}


# Logistic: in-hospital death
fit_mort_adj <- glm(in_hosp_death ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                    data = arf_exp, family = binomial())
tt_mort <- tidy_and_save(fit_mort_adj, "inhosp_death_adj", exponentiate = TRUE)

# Logistic: 30-day death
fit_mort30_adj <- glm(death_30d ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                      data = arf_exp, family = binomial())
tt_mort30 <- tidy_and_save(fit_mort30_adj, "death30d_adj", exponentiate = TRUE)

# NegBin: ICU LOS
fit_icu_nb <- glm.nb(icu_los_days ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                     data = arf_exp)
tt_icu <- tidy_and_save(fit_icu_nb, "icu_los_adj", exponentiate = TRUE)

# NegBin: Vent hours
fit_vent_nb <- glm.nb(vent_hours ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                      data = arf_exp)
tt_vent <- tidy_and_save(fit_vent_nb, "vent_hours_adj", exponentiate = TRUE)

# Optional: AKI & vaso (logistic)
fit_aki  <- glm(aki_flag  ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                data = arf_exp, family = binomial())

tidy_and_save(fit_aki,  "aki_adj",  exponentiate = TRUE)


get_model_auc <- function(fit, data, outcome, label = NULL) {
  # Build the exact complete-case frame used by the model
  df <- model.frame(formula(fit), data = data, na.action = na.omit)
  
  # Preds on those rows
  p  <- predict(fit, newdata = df, type = "response")
  y  <- df[[outcome]]
  
  # Ensure binary 0/1 with 0 = control, 1 = case
  if (is.logical(y)) y <- as.integer(y)
  if (is.numeric(y)) y <- factor(y, levels = c(0,1))
  if (is.factor(y) && !identical(levels(y), c("0","1"))) {
    y <- forcats::fct_relabel(y, as.character)
    y <- forcats::fct_drop(y)
  }
  
  r   <- pROC::roc(response = y, predictor = p, levels = c("0","1"), quiet = TRUE, direction = "<")
  ci  <- pROC::ci.auc(r)
  
  tibble(
    model   = label %||% deparse1(formula(fit)),
    outcome = outcome,
    n       = nrow(df),
    auc     = as.numeric(pROC::auc(r)),
    auc_lo  = as.numeric(ci[1]),
    auc_hi  = as.numeric(ci[3])
  )
}

auc_tbl <- bind_rows(
  get_model_auc(fit_mort_adj,    arf_exp, "in_hosp_death", "In-hospital death (adjusted)"),
  get_model_auc(fit_mort30_adj,  arf_exp, "death_30d",     "30-day death (adjusted)")
)

print(auc_tbl)

# Save alongside other outputs (repo-anchored helper)
save_tbl(auc_tbl, "metrics/auc_adjusted_models")

# ------------------------------------ 8) Plots: NO2 main effects --------------------------------
ref_sex <- arf_exp %>% count(sex_category, sort = TRUE) %>% slice(1) %>% pull(sex_category)
ref_re  <- arf_exp %>% count(race_ethnicity_simple, sort = TRUE) %>% slice(1) %>% pull(race_ethnicity_simple)

make_grid <- function(df, n = 100) {
  rng <- df %>% summarize(no2_min = min(no2_mean, na.rm = TRUE),
                          no2_max = max(no2_mean, na.rm = TRUE))
  tibble(
    no2_mean = seq(rng$no2_min, rng$no2_max, length.out = n),
    pm25_mean = mean(df$pm25_mean, na.rm = TRUE),
    age = mean(df$age, na.rm = TRUE),
    sex_category = ref_sex,
    race_ethnicity_simple = ref_re,
    svi_overall = mean(df$svi_overall, na.rm = TRUE),
    # NEW tract-level ACS covariates
    acs_median_income   = mean(df$acs_median_income, na.rm = TRUE),
    acs_pct_lt_hs       = mean(df$acs_pct_lt_hs, na.rm = TRUE),
    acs_unemp_rate_pct  = mean(df$acs_unemp_rate_pct, na.rm = TRUE),
    acs_pct_insured     = mean(df$acs_pct_insured, na.rm = TRUE),
  ) %>% mutate(no2_10 = no2_mean / 10)
}
grid <- make_grid(arf_exp)

pred_ci_logistic <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>% mutate(link = pr$fit, se = pr$se.fit,
                     pred = plogis(link),
                     lo = plogis(link - 1.96*se),
                     hi = plogis(link + 1.96*se))
}
pred_ci_negbin <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>% mutate(link = pr$fit, se = pr$se.fit,
                     pred = exp(link),
                     lo = exp(link - 1.96*se),
                     hi = exp(link + 1.96*se))
}

df_mort_inhosp <- pred_ci_logistic(fit_mort_adj, grid)
df_mort_30d    <- pred_ci_logistic(fit_mort30_adj, grid)
df_vent        <- pred_ci_negbin (fit_vent_nb,     grid)

p1 <- ggplot(df_mort_inhosp, aes(x = no2_10, y = pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(size = 1.2) +
  labs(x = "NO2 (per 10 ppb)", y = "Predicted probability",
       title = "In-hospital mortality vs NO2",
       subtitle = paste("Adj: PM2.5, age, sex, race/ethnicity, SVI; \nrefs:", ref_sex, "/", ref_re)) +
  theme_classic(base_size = 14)

p2 <- ggplot(df_mort_30d, aes(x = no2_10, y = pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(size = 1.2) +
  labs(x = "NO2 (per 10 ppb)", y = "Predicted probability",
       title = "30-day mortality vs NO2",
       subtitle = paste("Adj: PM2.5, age, sex, race/ethnicity, SVI; \nrefs:", ref_sex, "/", ref_re)) +
  theme_classic(base_size = 14)

p3 <- ggplot(df_vent, aes(x = no2_10, y = pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(size = 1.2) +
  labs(x = "NO2 (per 10 ppb)", y = "Predicted mean ventilation hours",
       title = "Ventilation hours vs NO2",
       subtitle = paste("NB model; Adj: PM2.5, age, sex, race/ethnicity, SVI; \nrefs:", ref_sex, "/", ref_re)) +
  theme_classic(base_size = 14)

combo_main <- (p1 | p2) / p3
save_plot(combo_main, "no2_outcomes_combined")

# ------------------------------------ 9) Effect Modification: ARF Subtype -----------------------
arf_exp <- arf_exp |>
  left_join(cohort_min |> dplyr::select(patient_id, hospitalization_id, hypoxemic_arf, hypercapnic_arf, mixed_arf),
            by = c("patient_id","hospitalization_id")) |>
  mutate(
    arf_subtype = case_when(
      mixed_arf == 1 ~ "Mixed",
      hypoxemic_arf == 1 ~ "Hypoxemic",
      hypercapnic_arf == 1 ~ "Hypercapnic",
      TRUE ~ "Other"
    ),
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic","Hypercapnic","Mixed","Other")),
    no2_10 = no2_mean/10
  )

fit_mort_inhosp_sub <- glm(in_hosp_death ~ pm25_mean + no2_10 * arf_subtype + age + sex_category + race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                           data = arf_exp, family = binomial())
fit_mort_30d_sub    <- glm(death_30d     ~ pm25_mean + no2_10 * arf_subtype + age + sex_category + race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                           data = arf_exp, family = binomial())
fit_vent_nb_sub     <- glm.nb(vent_hours  ~ pm25_mean + no2_10 * arf_subtype + age + sex_category + race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                              data = arf_exp)

tidy_and_save(fit_mort_inhosp_sub, "inhosp_death_x_subtype", exponentiate = TRUE)
tidy_and_save(fit_mort_30d_sub,    "death30d_x_subtype",    exponentiate = TRUE)
tidy_and_save(fit_vent_nb_sub,     "vent_hours_x_subtype",  exponentiate = TRUE)

# Pred grids
make_grid_sub <- function(df, n = 150) {
  rng <- df %>% summarize(lo = quantile(no2_10, 0.01, na.rm = TRUE),
                          hi = quantile(no2_10, 0.99, na.rm = TRUE))
  expand.grid(no2_10 = seq(rng$lo, rng$hi, length.out = n),
              arf_subtype = levels(df$arf_subtype)) %>%
    as_tibble() %>%
    mutate(
      pm25_mean = mean(df$pm25_mean, na.rm = TRUE),
      age = mean(df$age, na.rm = TRUE),
      sex_category = ref_sex,
      race_ethnicity_simple = ref_re,
      svi_overall = mean(df$svi_overall, na.rm = TRUE),
      # NEW tract-level ACS covariates
      acs_median_income   = mean(df$acs_median_income, na.rm = TRUE),
      acs_pct_lt_hs       = mean(df$acs_pct_lt_hs, na.rm = TRUE),
      acs_unemp_rate_pct  = mean(df$acs_unemp_rate_pct, na.rm = TRUE),
      acs_pct_insured     = mean(df$acs_pct_insured, na.rm = TRUE),
    )
}
grid_sub <- make_grid_sub(arf_exp)

pred_ci_log <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>% mutate(pred = plogis(pr$fit), lo = plogis(pr$fit - 1.96*pr$se.fit), hi = plogis(pr$fit + 1.96*pr$se.fit))
}
pred_ci_nb  <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>% mutate(pred = exp(pr$fit), lo = exp(pr$fit - 1.96*pr$se.fit), hi = exp(pr$fit + 1.96*pr$se.fit))
}

# --- Drop unused levels (subtype + covariates just in case) ---
arf_exp <- arf_exp %>%
  dplyr::mutate(
    arf_subtype = droplevels(arf_subtype),
    sex_category = droplevels(sex_category),
    race_ethnicity_simple = droplevels(race_ethnicity_simple)
  )

# --- Fit interaction models (NO2 per 10 ppb already in no2_10) ---
fit_mort_inhosp_sub <- glm(
  in_hosp_death ~ pm25_mean + no2_10 * arf_subtype + age + sex_category +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)

fit_mort_30d_sub <- glm(
  death_30d ~ pm25_mean + no2_10 * arf_subtype + age + sex_category +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)

fit_vent_nb_sub <- MASS::glm.nb(
  vent_hours ~ pm25_mean + no2_10 * arf_subtype + age + sex_category +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp
)

# --- Build grid from the model's levels and coerce factors before predict ---
arf_lvls <- fit_mort_inhosp_sub$xlevels$arf_subtype  # training levels used by model

make_grid_sub_safe <- function(df, arf_lvls, n = 150) {
  # helpers
  safe_mean <- function(x) {
    m <- mean(x, na.rm = TRUE)
    if (is.nan(m)) NA_real_ else m
  }
  top_level <- function(x) {
    if (is.factor(x)) x <- droplevels(x)
    tab <- sort(table(x), decreasing = TRUE)
    if (length(tab) == 0) NA_character_ else names(tab)[1]
  }

  rng <- df %>%
    dplyr::summarize(
      lo = stats::quantile(no2_10, 0.01, na.rm = TRUE),
      hi = stats::quantile(no2_10, 0.99, na.rm = TRUE)
    )

  expand.grid(
    no2_10 = seq(rng$lo, rng$hi, length.out = n),
    arf_subtype = arf_lvls,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      pm25_mean             = safe_mean(df$pm25_mean),
      age                   = safe_mean(df$age),
      sex_category          = top_level(df$sex_category),
      race_ethnicity_simple = top_level(df$race_ethnicity_simple),
      svi_overall           = safe_mean(df$svi_overall),

      # NEW tract-level ACS covariates
      acs_median_income   = safe_mean(df$acs_median_income),
      acs_pct_lt_hs       = safe_mean(df$acs_pct_lt_hs),
      acs_unemp_rate_pct  = safe_mean(df$acs_unemp_rate_pct),
      acs_pct_insured     = safe_mean(df$acs_pct_insured),

      # coerce to model's factor levels
      arf_subtype = factor(arf_subtype, levels = arf_lvls),
      sex_category = factor(sex_category, levels = levels(df$sex_category)),
      race_ethnicity_simple = factor(race_ethnicity_simple, levels = levels(df$race_ethnicity_simple))
    )
}

grid_sub <- make_grid_sub_safe(arf_exp, arf_lvls)

# --- Predict with your existing helpers pred_ci_log() / pred_ci_nb() ---
df_mort_inhosp_sub <- pred_ci_log(fit_mort_inhosp_sub, grid_sub)
df_mort_30d_sub    <- pred_ci_log(fit_mort_30d_sub,    grid_sub)
df_vent_sub        <- pred_ci_nb (fit_vent_nb_sub,     grid_sub)

theme_pub <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

pal <- scale_color_brewer(palette = "Dark2")
fill_pal <- scale_fill_brewer(palette = "Dark2")

p1s <- ggplot(df_mort_inhosp_sub, aes(no2_10, pred, color = arf_subtype, fill = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) + pal + fill_pal +
  labs(title = "In-hospital mortality vs NO2 by ARF subtype",
       x = "NO2 (per 10 ppb)", y = "Predicted probability",
       color = "ARF subtype", fill = "ARF subtype",
       subtitle = paste("Adjusted for PM2.5, age, sex, race/ethnicity, SVI; \nrefs:",
                        ref_sex, "/", ref_re)) + theme_pub

p2s <- ggplot(df_mort_30d_sub, aes(no2_10, pred, color = arf_subtype, fill = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) + pal + fill_pal +
  labs(title = "30-day mortality vs NO2 by ARF subtype",
       x = "NO2 (per 10 ppb)", y = "Predicted probability",
       color = "ARF subtype", fill = "ARF subtype",
       subtitle = paste("Adjusted for PM2.5, age, sex, race/ethnicity, SVI; \nrefs:",
                        ref_sex, "/", ref_re)) + theme_pub

p3s <- ggplot(df_vent_sub, aes(no2_10, pred, color = arf_subtype, fill = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) + pal + fill_pal +
  labs(title = "Ventilation hours vs NO2 by ARF subtype",
       x = "NO2 (per 10 ppb)", y = "Predicted mean hours",
       color = "ARF subtype", fill = "ARF subtype",
       subtitle = paste("Negative binomial; adjusted for PM2.5, age, sex, race/ethnicity, SVI; \nrefs:",
                        ref_sex, "/", ref_re)) + theme_pub

combo_sub <- (p1s | p2s) / p3s + plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_plot(combo_sub, "no2_by_subtype_combined")

# --- drop unused levels just once ---
arf_exp <- arf_exp %>%
  mutate(
    arf_subtype = droplevels(arf_subtype),
    sex_category = droplevels(sex_category),
    race_ethnicity_simple = droplevels(race_ethnicity_simple)
  )

# Logistic: in-hospital and 30-day mortality
fit_mort_inhosp_sex <- glm(
  in_hosp_death ~ pm25_mean + no2_10*sex_category + arf_subtype + age +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)

fit_mort_30d_sex <- glm(
  death_30d ~ pm25_mean + no2_10*sex_category + arf_subtype + age +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)

# Negative binomial: vent hours
fit_vent_nb_sex <- MASS::glm.nb(
  vent_hours ~ pm25_mean + no2_10*sex_category + arf_subtype + age +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp
)

# levels from fitted model to avoid new levels at predict-time
sx_lvls  <- fit_mort_inhosp_sex$xlevels$sex_category
sub_lvls <- fit_mort_inhosp_sex$xlevels$arf_subtype

make_grid_sex <- function(df, sx_lvls, sub_lvls, n = 150) {
  safe_mean <- function(x) { m <- mean(x, na.rm = TRUE); if (is.nan(m)) NA_real_ else m }
  top_level <- function(x) { if (is.factor(x)) x <- droplevels(x); names(sort(table(x), TRUE))[1] }
  
  rng <- df %>% summarize(lo = quantile(no2_10, 0.01, na.rm = TRUE),
                          hi = quantile(no2_10, 0.99, na.rm = TRUE))
  
  expand.grid(
    no2_10 = seq(rng$lo, rng$hi, length.out = n),
    sex_category = sx_lvls,
    arf_subtype = sub_lvls,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) |> 
    as_tibble() |>
    mutate(
      pm25_mean             = safe_mean(df$pm25_mean),
      age                   = safe_mean(df$age),
      race_ethnicity_simple = top_level(df$race_ethnicity_simple),
      svi_overall           = safe_mean(df$svi_overall),
      acs_median_income     = safe_mean(df$acs_median_income),
      acs_pct_lt_hs         = safe_mean(df$acs_pct_lt_hs),
      acs_unemp_rate_pct    = safe_mean(df$acs_unemp_rate_pct),
      acs_pct_insured       = safe_mean(df$acs_pct_insured),
      sex_category          = factor(sex_category, levels = levels(df$sex_category)),
      arf_subtype           = factor(arf_subtype,  levels = levels(df$arf_subtype)),
      race_ethnicity_simple = factor(race_ethnicity_simple, levels = levels(df$race_ethnicity_simple))
    )
}

grid_sex <- make_grid_sex(arf_exp, sx_lvls, sub_lvls)

df_mort_inhosp_sex <- pred_ci_log(fit_mort_inhosp_sex, grid_sex)
df_mort_30d_sex    <- pred_ci_log(fit_mort_30d_sex,    grid_sex)
df_vent_sex        <- pred_ci_nb (fit_vent_nb_sex,     grid_sex)

p1_sex <- ggplot(df_mort_inhosp_sex, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ sex_category) +
  pal + theme_pub +
  labs(title = "In-hospital mortality vs NO₂, stratified by sex",
       x = "NO₂ (per 10 ppb)", y = "Predicted probability",
       color = "ARF subtype")

p2_sex <- ggplot(df_mort_30d_sex, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ sex_category) +
  pal + theme_pub +
  labs(title = "30-day mortality vs NO₂, stratified by sex",
       x = "NO₂ (per 10 ppb)", y = "Predicted probability",
       color = "ARF subtype")

p3_sex <- ggplot(df_vent_sex, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ sex_category) +
  pal + theme_pub +
  labs(title = "Ventilation hours vs NO₂, stratified by sex",
       x = "NO₂ (per 10 ppb)", y = "Predicted mean hours",
       color = "ARF subtype")


combo_sex <- (p1_sex | p2_sex) / p3_sex + plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_plot(combo_sex, "no2_by_sex_subtype_combined")


fit_mort_inhosp_re <- glm(
  in_hosp_death ~ pm25_mean + no2_10*race_ethnicity_simple + arf_subtype + age +
    sex_category + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)

fit_mort_30d_re <- glm(
  death_30d ~ pm25_mean + no2_10*race_ethnicity_simple + arf_subtype + age +
    sex_category + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)

fit_vent_nb_re <- MASS::glm.nb(
  vent_hours ~ pm25_mean + no2_10*race_ethnicity_simple + arf_subtype + age +
    sex_category + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp
)

re_lvls  <- fit_mort_inhosp_re$xlevels$race_ethnicity_simple
sub_lvls <- fit_mort_inhosp_re$xlevels$arf_subtype

make_grid_re <- function(df, re_lvls, sub_lvls, n = 150) {
  safe_mean <- function(x) { m <- mean(x, na.rm = TRUE); if (is.nan(m)) NA_real_ else m }
  rng <- df %>% summarize(lo = quantile(no2_10, 0.01, na.rm = TRUE),
                          hi = quantile(no2_10, 0.99, na.rm = TRUE))
  expand.grid(
    no2_10 = seq(rng$lo, rng$hi, length.out = n),
    race_ethnicity_simple = re_lvls,
    arf_subtype = sub_lvls,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) |>
    as_tibble() |>
    mutate(
      pm25_mean          = safe_mean(df$pm25_mean),
      age                = safe_mean(df$age),
      sex_category       = names(sort(table(df$sex_category), TRUE))[1],
      svi_overall        = safe_mean(df$svi_overall),
      acs_median_income  = safe_mean(df$acs_median_income),
      acs_pct_lt_hs      = safe_mean(df$acs_pct_lt_hs),
      acs_unemp_rate_pct = safe_mean(df$acs_unemp_rate_pct),
      acs_pct_insured    = safe_mean(df$acs_pct_insured),
      race_ethnicity_simple = factor(race_ethnicity_simple, levels = levels(df$race_ethnicity_simple)),
      sex_category          = factor(sex_category, levels = levels(df$sex_category)),
      arf_subtype           = factor(arf_subtype, levels = levels(df$arf_subtype))
    )
}

grid_re <- make_grid_re(arf_exp, re_lvls, sub_lvls)

df_mort_inhosp_re <- pred_ci_log(fit_mort_inhosp_re, grid_re)
df_mort_30d_re    <- pred_ci_log(fit_mort_30d_re,    grid_re)
df_vent_re        <- pred_ci_nb (fit_vent_nb_re,     grid_re)

# set y-axis limits for all races except Hispanic Black
mort_y_lim <- c(0, 0.3)
vent_y_lim <- c(0, 20)

# mortality: in-hospital
p1_re <- ggplot(df_mort_inhosp_re, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ race_ethnicity_simple, nrow = 2) +
  pal + theme_pub +
  labs(title = "In-hospital mortality vs NO₂, \nstratified by race/ethnicity",
       x = "NO₂ (per 10 ppb)", y = "Predicted probability",
       color = "ARF subtype") +
  scale_y_continuous(limits = mort_y_lim) +
  # allow Hispanic Black to exceed limits
  facet_wrap(~ race_ethnicity_simple, nrow = 2, scales = "free_y")

# mortality: 30-day
p2_re <- ggplot(df_mort_30d_re, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ race_ethnicity_simple, nrow = 2) +
  pal + theme_pub +
  labs(title = "30-day mortality vs NO₂, \nstratified by race/ethnicity",
       x = "NO₂ (per 10 ppb)", y = "Predicted probability",
       color = "ARF subtype") +
  scale_y_continuous(limits = mort_y_lim) +
  facet_wrap(~ race_ethnicity_simple, nrow = 2, scales = "free_y")

# ventilation hours
p3_re <- ggplot(df_vent_re, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ race_ethnicity_simple, nrow = 2) +
  pal + theme_pub +
  labs(title = "Ventilation hours vs NO₂, stratified by race/ethnicity",
       x = "NO₂ (per 10 ppb)", y = "Predicted mean hours",
       color = "ARF subtype") +
  scale_y_continuous(limits = vent_y_lim) +
  facet_wrap(~ race_ethnicity_simple, nrow = 2, scales = "free_y")

combo_re <- (p1_re | p2_re) / p3_re + plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_plot(combo_re, "no2_by_race_subtype_combined")


# ------------------------------------ 10) Descriptive Tables ------------------------------------
arf_summary <- outcomes_exp %>%
  filter(cohort == "ARF") %>%
  summarise(
    n            = n(),
    mort_inhosp  = mean(in_hosp_death, na.rm = TRUE),
    mort_30d     = mean(death_30d, na.rm = TRUE),
    mean_ICU_los = mean(icu_los_days, na.rm = TRUE),
    mean_PM25    = mean(pm25_mean, na.rm = TRUE),
    mean_NO2     = mean(no2_mean, na.rm = TRUE),
    mean_SVI     = mean(svi_overall, na.rm = TRUE)
  )
save_tbl(arf_summary, "arf_summary")

arf_strat <- outcomes_exp %>%
  filter(cohort == "ARF") %>%
  mutate(svi_tertile = ntile(svi_overall, 3),
         pm25_quint  = ntile(pm25_mean, 5)) %>%
  group_by(svi_tertile) %>%
  summarise(n = n(),
            mort_inhosp = mean(in_hosp_death, na.rm = TRUE),
            mort_30d    = mean(death_30d, na.rm = TRUE),
            mean_ICU_los= mean(icu_los_days, na.rm = TRUE), .groups = "drop")
save_tbl(arf_strat, "arf_strat")


# ------------- 1) Variable catalog & labels (edit as needed) --------------------
tbl1_vars <- list(
  cont = c(
    age              = "Age (years)",
    pm25_mean        = "PM2.5 (annual mean, μg/m³)",
    no2_mean         = "NO₂ (annual mean, ppb)",
    svi_overall      = "SVI (overall)",
    icu_los_days     = "ICU length of stay (days)",
    hosp_los_days    = "Hospital length of stay (days)",
    vent_hours       = "Invasive ventilation (hours)"
  ),
  cat = c(
    sex_category            = "Sex",
    race_ethnicity_simple   = "Race/Ethnicity",
    in_hosp_death           = "In-hospital death",
    death_30d               = "30-day death",
    vaso_flag               = "Vasopressor use",
    aki_flag                = "AKI"
  )
)

# Include optional history features if present
opt_hist <- c("prior_icu_stays","any_acei_arb","any_diuretic","any_bb")
opt_hist_present <- intersect(opt_hist, names(arf_exp))
if ("prior_icu_stays" %in% opt_hist_present) {
  tbl1_vars$cont <- c(tbl1_vars$cont, prior_icu_stays = "Prior ICU stays (count)")
}
bin_map <- c(any_acei_arb = "ACEi/ARB during lookback",
             any_diuretic = "Diuretic during lookback",
             any_bb       = "Beta-blocker during lookback")
add_bins <- intersect(names(bin_map), opt_hist_present)
if (length(add_bins)) {
  tbl1_vars$cat <- c(tbl1_vars$cat, stats::setNames(as.list(bin_map[add_bins]), add_bins))
}

# ------------- 2) Clean factor levels for reproducible categories ---------------
arf_exp_tbl <- arf_exp %>%
  mutate(
    sex_category = forcats::fct_explicit_na(sex_category, na_level = "(Missing)"),
    race_ethnicity_simple = forcats::fct_explicit_na(race_ethnicity_simple, na_level = "(Missing)")
  )

# ------------- 3) Pretty Table 1 (gtsummary) -----------------------------------
# Formatting for continuous & categorical
# 0) Combine your label dictionaries into ONE named character vector
all_labels <- c(tbl1_vars$cont, tbl1_vars$cat)  # named chr: var -> "Pretty Label"

# 1) Build a pure list of formulas: var ~ "Pretty Label"
labels_list <- lapply(names(all_labels), function(v) {
  rlang::new_formula(rlang::sym(v), as.character(all_labels[[v]]))
})
# ensure there are no names on the outer list (gtsummary prefers an unnamed list of formulas)
names(labels_list) <- NULL

# 2) (Optional) sanity check: make sure all labeled vars exist in data you select
vars_keep <- intersect(names(all_labels), names(arf_exp_tbl))

# 3) Build Table 1
tbl1 <- arf_exp_tbl %>%
  dplyr::select(dplyr::any_of(vars_keep)) %>%
  gtsummary::tbl_summary(
    type = list(
      gtsummary::all_continuous()  ~ "continuous",
      gtsummary::all_categorical() ~ "categorical"
    ),
    statistic = c(
      list(gtsummary::all_continuous()  ~ "{mean} ± {sd} ({median}; {p25}, {p75})"),
      list(gtsummary::all_categorical() ~ "{n} ({p}%)")
    ),
    label   = labels_list,
    missing = "ifany",
    digits  = gtsummary::all_continuous() ~ 1
  ) %>%
  # gtsummary::add_overall(last = TRUE) %>%  # <- remove; no 'by' stratifier
  gtsummary::bold_labels() %>%
  gtsummary::modify_caption(
    paste0("**Table 1. Baseline characteristics and outcomes — ARF cohort** (Site: ", site_name, ")")
  )

# Save (HTML + RTF) — convert with gtsummary::as_gt()
tbl1_gt <- gtsummary::as_gt(tbl1)

gt::gtsave(tbl1_gt, filename = file.path(cfg$figures_path,
                                         paste0(cfg$prefix, "_table1_", cfg$run_id, ".html")))
# If your gt version supports RTF (≥ 0.10.0 typically); otherwise skip:
try(
  gt::gtsave(tbl1_gt, filename = file.path(cfg$figures_path,
                                           paste0(cfg$prefix, "_table1_", cfg$run_id, ".rtf"))),
  silent = TRUE
)

# ------------- 4) Federated-friendly machine-readable export -------------------
# We’ll export two tidy frames:
#   A) continuous_stats: variable, n, mean, sd, q25, median, q75
#   B) categorical_stats: variable, level, n, pct
# This allows exact recomputation across nodes (pooled mean/SD, and summed counts).

# A) Continuous
continuous_vars <- intersect(names(tbl1_vars$cont), names(arf_exp_tbl))
continuous_stats <- map_dfr(continuous_vars, function(v) {
  x <- arf_exp_tbl[[v]]
  x <- x[is.finite(as.numeric(x))]  # drop NA / non-numeric safely
  tibble(
    site      = site_name,
    variable  = v,
    label     = unname(tbl1_vars$cont[v]),
    n         = length(x),
    mean      = if (length(x)) mean(x, na.rm = TRUE) else NA_real_,
    sd        = if (length(x)) stats::sd(x, na.rm = TRUE) else NA_real_,
    q25       = if (length(x)) as.numeric(quantile(x, 0.25, na.rm = TRUE)) else NA_real_,
    median    = if (length(x)) as.numeric(quantile(x, 0.50, na.rm = TRUE)) else NA_real_,
    q75       = if (length(x)) as.numeric(quantile(x, 0.75, na.rm = TRUE)) else NA_real_
  )
})

# B) Categorical (ensure factors)
categorical_vars <- intersect(names(tbl1_vars$cat), names(arf_exp_tbl))
categorical_stats <- map_dfr(categorical_vars, function(v) {
  x <- arf_exp_tbl[[v]]
  x_fac <- if (is.factor(x)) x else factor(x)
  cnt <- as.data.frame(table(x_fac, useNA = "ifany"), stringsAsFactors = FALSE)
  names(cnt) <- c("level", "n")
  cnt %>%
    mutate(
      site     = site_name,
      variable = v,
      label    = unname(tbl1_vars$cat[v]),
      total_n  = sum(n),
      pct      = ifelse(total_n > 0, 100 * n / total_n, NA_real_)
    ) %>%
    dplyr::select(site, variable, label, level, n, pct, total_n)
})

# Save site-level tidy stats
save_tbl(continuous_stats, "table1_continuous_site")
save_tbl(categorical_stats, "table1_categorical_site")

message("All done. Outputs written to: ", normalizePath(cfg$output_dir))


# =========================
# A) ARF Incidence
# =========================


# 1. Estimate county-specific slope (no2_mean ~ year)
no2_slopes <- no2 %>%
  group_by(GEOID) %>%
  do({
    fit <- lm(no2_mean ~ year, data = .)
    tibble(
      slope = coef(fit)[["year"]],
      intercept = coef(fit)[["(Intercept)"]]
    )
  })

# 2. Predict 2018 from slope + intercept
no2_2018 <- no2_slopes %>%
  mutate(
    year = 2018,
    no2_mean = intercept + slope * 2018
  ) %>%
  dplyr::select(GEOID, year, no2_mean)

# 3. Append backfilled 2018 to your original no2 dataset
no2_complete <- no2 %>%
  bind_rows(no2_2018) %>%
  arrange(GEOID, year)

years <- 2018:2024

pop_data <- map_dfr(years, function(y) {
  get_estimates(
    geography = "county",
    product = "population", 
    year = y,
    geometry = FALSE
  ) %>%
    mutate(year = y)
})

pops <- c("POP", "POPESTIMATE")

pop_data_clean <- pop_data %>%
  filter(variable %in% pops) %>%
  dplyr::select(GEOID, NAME, year, population = value)


arf_counts <- arf_exp %>%
  mutate(year = year(index_admit)) %>% 
  group_by(fips_county, year) %>%
  summarise(arf_cases = n(), .groups = "drop")

# --- 0) Ensure keys are consistent types/widths ---
pop_data_clean <- pop_data_clean %>%
  mutate(
    GEOID = str_pad(as.character(GEOID), width = 5, side = "left", pad = "0"),
    year  = as.integer(year)
  )

no2_complete <- no2_complete %>%
  mutate(
    GEOID = str_pad(as.character(GEOID), width = 5, side = "left", pad = "0"),
    year  = as.integer(year)
  )

# If arf_exp has county as numeric or short FIPS, standardize to 5 digits
# Replace `county_fips` and `admit_date` with the actual column names in arf_exp if different
arf_counts <- arf_exp %>%
  mutate(
    GEOID = str_pad(as.character(fips_county), width = 5, side = "left", pad = "0"),
    year  = year(index_admit)
  ) %>%
  group_by(GEOID, year) %>%
  summarise(arf_cases = n(), .groups = "drop")

# --- 1) (Optional) if there are duplicate NO2 rows per GEOID-year, collapse to mean ---
no2_collapse <- no2_complete %>%
  group_by(GEOID, year) %>%
  summarise(no2_mean = mean(no2_mean, na.rm = TRUE), .groups = "drop")

# --- 2) Build base panel with pop + NO2 (keep only county-years where both exist) ---
base_panel <- pop_data_clean %>%
  inner_join(no2_collapse, by = c("GEOID", "year"))

# --- 3) Add ARF counts; fill missing with zeros (no observed cases that year) ---
analysis_df <- base_panel %>%
  left_join(arf_counts, by = c("GEOID", "year")) %>%
  mutate(
    arf_cases = coalesce(arf_cases, 0L)
  )

# --- 4) Restrict to target years (if desired) ---
analysis_df <- analysis_df %>%
  filter(year >= 2018, year <= 2024)

# --- 5) Quick sanity checks ---
# Any county-years missing pop or NO2? (should be none after inner_join)
stopifnot(!any(is.na(analysis_df$population)))
stopifnot(!any(is.na(analysis_df$no2_mean)))

# Any ARF cases that didn't join to pop/NO2? (should be 0 after inner_join base)
arf_unmatched <- arf_counts %>%
  anti_join(base_panel, by = c("GEOID", "year"))
if (nrow(arf_unmatched) > 0) {
  message("Heads up: Some ARF counts had no matching pop+NO2 rows. Examples:")
  print(head(arf_unmatched, 10))
}

# Peek
dplyr::glimpse(analysis_df)
dplyr::count(analysis_df, year)
summary(analysis_df$arf_cases)

# Output dir for figures
fig_dir <- cfg$figures_path
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# -------------------------
# Helper: save plots
# -------------------------
save_fig <- function(p, name, width = 7, height = 5, dpi = 320) {
  ggsave(filename = file.path(fig_dir, name), plot = p,
         width = width, height = height, dpi = dpi, device = ragg::agg_png)
}

# =========================
# Data prep
# Assumes: analysis_df has GEOID, NAME, year, population, no2_mean, arf_cases
# =========================
analysis_df <- analysis_df %>%
  mutate(population = as.numeric(population),
         log_pop    = log(population))

stopifnot(all(is.finite(analysis_df$log_pop)))

# =========================
# A) Main model: FE Negative Binomial (NO2)
# =========================
fit_nb <- fenegbin(
  arf_cases ~ no2_mean | GEOID + year,
  data   = analysis_df,
  offset = ~ log_pop
)

# Clustered SE table
etable(fit_nb, se = "cluster", cluster = ~ GEOID)

# IRR point estimate (per 1 ppb NO2)
irr_no2 <- exp(coef(fit_nb)["no2_mean"])

# Tidy + plot
no2_res <- tidy(fit_nb, conf.int = TRUE, se.type = "cluster", cluster = ~ GEOID) %>%
  filter(term == "no2_mean") %>%
  transmute(term, IRR = exp(estimate),
            IRR_low = exp(conf.low), IRR_high = exp(conf.high))

p_nb_point <- ggplot(no2_res,
                     aes(x = term, y = IRR, ymin = IRR_low, ymax = IRR_high)) +
  geom_pointrange(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = expression("Incidence Rate Ratio for NO"[2]*" and ARF"),
       x = NULL, y = "Incidence Rate Ratio (95% CI)") +
  theme_minimal(base_size = 14)
save_fig(p_nb_point, "irr_no2_nb.png")

# =========================
# B) Historic burden tertiles (descriptive + FE Poisson by stratum)
# =========================
# Long-term mean NO2 (2019+; keep if you included 2018 extrapolation)
county_burden <- analysis_df %>%
  filter(year >= 2019) %>%
  group_by(GEOID) %>%
  summarise(no2_historic_mean = mean(no2_mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(no2_tertile = ntile(no2_historic_mean, 3),
         no2_tertile = factor(no2_tertile, levels = 1:3,
                              labels = c("Low NO\u2082","Medium NO\u2082","High NO\u2082")))

analysis_df2 <- analysis_df %>% left_join(county_burden, by = "GEOID")

# Stratified FE Poisson (more stable than NB for interactions)
fits_tertile <- analysis_df2 %>%
  group_by(no2_tertile) %>%
  group_map(~ fepois(arf_cases ~ no2_mean | GEOID + year, data = .x, offset = ~ log_pop))

names(fits_tertile) <- levels(analysis_df2$no2_tertile)

tidy_tertiles <- imap_dfr(fits_tertile, ~ tidy(.x, conf.int = TRUE,
                                               se.type = "cluster", cluster = ~ GEOID) %>%
                            filter(term == "no2_mean") %>%
                            transmute(tertile = .y,
                                      IRR = exp(estimate),
                                      IRR_low = exp(conf.low),
                                      IRR_high = exp(conf.high)))

p_tertile_pts <- ggplot(tidy_tertiles,
                        aes(x = tertile, y = IRR, ymin = IRR_low, ymax = IRR_high)) +
  geom_pointrange(color = "steelblue", linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Effect of NO\u2082 on ARF Incidence by County Burden Tertile",
       x = "County NO\u2082 Burden Tertile", y = "Incidence Rate Ratio (95% CI)") +
  theme_minimal(base_size = 14)
save_fig(p_tertile_pts, "irr_no2_by_tertile.png")

# =========================
# C) Continuous interaction: NO2 × historic burden (FE Poisson)
#   + Marginal effects curve (IRR per +10 ppb)
# =========================
fit_nb_interact <- fepois(
  arf_cases ~ no2_mean * scale(no2_historic_mean) | GEOID + year,
  data = analysis_df2,
  offset = ~ log_pop
)

per <- 1

# 1) Grid over the FULL observed range (no trimming)
burden_grid_full <- seq(
  min(analysis_df2$no2_historic_mean, na.rm = TRUE),
  max(analysis_df2$no2_historic_mean, na.rm = TRUE),
  length.out = 200
)

# 2) Coefs & clustered vcov (from fit_nb_interact)
b <- coef(fit_nb_interact)
V <- vcov(fit_nb_interact, cluster = ~ GEOID)

b_no2 <- b[["no2_mean"]]
b_int <- b[["no2_mean:scale(no2_historic_mean)"]]

mu    <- mean(analysis_df2$no2_historic_mean, na.rm = TRUE)
sdv   <- sd(analysis_df2$no2_historic_mean,   na.rm = TRUE)
z     <- (burden_grid_full - mu) / sdv

slope     <- as.numeric(b_no2 + b_int * z)
var_slope <- V["no2_mean","no2_mean"] +
  (z^2) * V["no2_mean:scale(no2_historic_mean)","no2_mean:scale(no2_historic_mean)"] +
  2 * z * V["no2_mean","no2_mean:scale(no2_historic_mean)"]
se_slope  <- sqrt(pmax(0, var_slope))

curve_full <- tibble(
  burden  = burden_grid_full,
  IRR     = exp(slope * per),
  IRR_low = exp(slope * per - 1.96 * se_slope * per),
  IRR_high= exp(slope * per + 1.96 * se_slope * per)
)

# Tertile cut lines
cuts <- quantile(analysis_df2$no2_historic_mean, c(1/3, 2/3), na.rm = TRUE)

# 3) Plot: full range + log y-scale
p_me_full <- ggplot(curve_full, aes(x = burden, y = IRR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = IRR_low, ymax = IRR_high), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = cuts, linetype = "dotted", linewidth = 0.6) +
  geom_rug(data = analysis_df2, aes(x = no2_historic_mean, y = NULL),
           inherit.aes = FALSE, sides = "b", alpha = 0.25) +
  scale_y_log10() +  # show entire spectrum cleanly
  labs(
    title    = expression("Marginal Effect of NO"[2]*" on ARF Incidence vs Historic NO"[2]*" Burden"),
    subtitle = bquote("IRR per +"*.(per)*" ppb NO"[2]*"; county & year fixed effects, clustered SEs"),
    x        = expression("Historic NO"[2]*" burden (county mean over study period)"),
    y        = "Incidence Rate Ratio (log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

save_fig(p_me_full, "marginal_effect_no2_by_burden_FULL.png", width = 9, height = 6)

# =========================
# D) Year-over-year (YOY) change models: lag 0/1/2
# =========================
analysis_df_lag <- analysis_df %>%
  arrange(GEOID, year) %>%
  group_by(GEOID) %>%
  mutate(
    no2_change      = no2_mean - lag(no2_mean),  # ΔNO2 (t-1 -> t)
    no2_change_lag0 = no2_change,                # concurrent
    no2_change_lag1 = lag(no2_change, 1),        # 1-year lag
    no2_change_lag2 = lag(no2_change, 2)         # 2-year lag
  ) %>%
  ungroup()

fit_lag0 <- fepois(arf_cases ~ no2_change_lag0 | GEOID + year, data = analysis_df_lag, offset = ~ log_pop)
fit_lag1 <- fepois(arf_cases ~ no2_change_lag1 | GEOID + year, data = analysis_df_lag, offset = ~ log_pop)
fit_lag2 <- fepois(arf_cases ~ no2_change_lag2 | GEOID + year, data = analysis_df_lag, offset = ~ log_pop)

tidy_lags <- list(Lag0 = fit_lag0, Lag1 = fit_lag1, Lag2 = fit_lag2) %>%
  imap_dfr(~ tidy(.x, conf.int = TRUE, se.type = "cluster", cluster = ~ GEOID) %>%
             filter(term == paste0("no2_change_", tolower(.y))) %>%
             transmute(Lag = .y,
                       IRR = exp(estimate),
                       IRR_low = exp(conf.low),
                       IRR_high = exp(conf.high)))

p_lags <- ggplot(tidy_lags, aes(x = Lag, y = IRR, ymin = IRR_low, ymax = IRR_high)) +
  geom_pointrange(size = 1.1, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = expression("Effect of Year-to-Year Change in NO"[2]*" on ARF Incidence"),
       subtitle = "Comparing concurrent (Lag 0), 1-year lag, and 2-year lag effects",
       x = "Lag (years)", y = "Incidence Rate Ratio (95% CI)") +
  theme_minimal(base_size = 14)
save_fig(p_lags, "lagged_delta_no2_effects.png")

# =========================
# E) PM2.5 models (level effects + side-by-side comparison)
# =========================
# Assumes data frame 'pm25' with GEOID, year, pm25_mean
analysis_df_pm <- analysis_df %>% left_join(pm25, by = c("GEOID","year"))

fit_no2_pois  <- fepois(arf_cases ~ no2_mean  | GEOID + year, data = analysis_df_pm, offset = ~ log_pop)
fit_pm25_pois <- fepois(arf_cases ~ pm25_mean | GEOID + year, data = analysis_df_pm, offset = ~ log_pop)
fit_both_pois <- fepois(arf_cases ~ no2_mean + pm25_mean | GEOID + year, data = analysis_df_pm, offset = ~ log_pop)

etable(fit_no2_pois, fit_pm25_pois, fit_both_pois,
       se = "cluster", cluster = ~ GEOID,
       dict = c(no2_mean = "NO\u2082 (per 1 ppb)",
                pm25_mean = "PM\u2082\u2022\u2085 (per 1 \u00B5g/m\u00B3)"))

tidy_models <- list(NO2 = fit_no2_pois, PM25 = fit_pm25_pois, Both = fit_both_pois) %>%
  imap_dfr(~ tidy(.x, conf.int = TRUE, se.type = "cluster", cluster = ~ GEOID) %>%
             filter(term %in% c("no2_mean", "pm25_mean")) %>%
             transmute(model = .y, term,
                       IRR = exp(estimate),
                       IRR_low = exp(conf.low),
                       IRR_high = exp(conf.high)))

p_side <- ggplot(tidy_models,
                 aes(x = term, y = IRR, ymin = IRR_low, ymax = IRR_high, color = model)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("no2_mean" = expression("NO"[2]),
                              "pm25_mean" = expression("PM"[2.5]))) +
  labs(title = "Pollutant Effects on ARF Incidence",
       subtitle = "Fixed-effects Poisson models (county + year), clustered SEs",
       x = "Pollutant", y = "Incidence Rate Ratio (95% CI)", color = "Model") +
  theme_minimal(base_size = 14)
save_fig(p_side, "pollutant_effects_side_by_side.png")

# output folder
out_dir <- file.path(cfg$output_dir, "data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

arf_counts_export <- analysis_df %>%
  mutate(
    GEOID = str_pad(as.character(GEOID), 5, side = "left", pad = "0"),
    year = as.integer(year),
    arf_cases = as.integer(arf_cases)
  ) %>%
  filter(year >= 2018, year <= 2024) %>%
  dplyr::select(GEOID, dplyr::any_of("NAME"), year, arf_cases) %>%
  arrange(GEOID, year)

write_csv(arf_counts_export, file.path(out_dir, "arf_counts_by_county_year.csv"))



# ---------------- Federated Output Script ----------------
# Each site runs this locally. Only aggregated predictions are exported.
# No line-level or patient data leaves the site.
# ---------------------------------------------------------

# ---- 0) Site config ----
site_plain_name <- config$site_name
out_dir         <- file.path(cfg$output_dir, "federated_preds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

site_id <- digest(paste(site_plain_name, sep="::"), algo = "sha1")

# ---- 1) Data prep ----
stopifnot(exists("arf_exp"))

arf_exp <- arf_exp %>%
  mutate(
    arf_subtype = case_when(
      mixed_arf == 1 ~ "Mixed",
      hypoxemic_arf == 1 ~ "Hypoxemic",
      hypercapnic_arf == 1 ~ "Hypercapnic",
      TRUE ~ "Other"
    ),
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic","Hypercapnic","Mixed","Other")),
    no2_10 = no2_mean/10,
    sex_category = droplevels(factor(sex_category)),
    race_ethnicity_simple = droplevels(factor(race_ethnicity_simple))
  )

# ---- 2) Shared NO2 grid ----
no2_grid <- tibble(no2_10 = seq(0, 6, by = 0.02))   # 0–60 ppb

# --- Safe helpers (unchanged) ---
safe_mean <- function(x) { m <- mean(x, na.rm = TRUE); if (is.nan(m)) NA_real_ else m }
top_level <- function(x) {
  if (is.factor(x)) x <- droplevels(x)
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0) NA_character_ else names(tab)[1]
}

# --- Shared NO2 grid (unchanged) ---
no2_grid <- tibble::tibble(no2_10 = seq(0, 6, by = 0.02))

# --- FIXED: make_grid that avoids duplicate names ---
make_grid <- function(df, by_vars = character()) {
  # Start with NO2 grid
  base <- tibble::tibble(no2_10 = no2_grid$no2_10)
  
  # Add anchors ONLY if that variable is NOT being stratified
  if (!"pm25_mean" %in% by_vars)
    base$pm25_mean <- safe_mean(df$pm25_mean)
  if (!"age" %in% by_vars)
    base$age <- safe_mean(df$age)
  if (!"sex_category" %in% by_vars)
    base$sex_category <- top_level(df$sex_category)
  if (!"race_ethnicity_simple" %in% by_vars)
    base$race_ethnicity_simple <- top_level(df$race_ethnicity_simple)
  if (!"svi_overall" %in% by_vars)
    base$svi_overall <- safe_mean(df$svi_overall)
  if (!"acs_median_income" %in% by_vars)
    base$acs_median_income <- safe_mean(df$acs_median_income)
  if (!"acs_pct_lt_hs" %in% by_vars)
    base$acs_pct_lt_hs <- safe_mean(df$acs_pct_lt_hs)
  if (!"acs_unemp_rate_pct" %in% by_vars)
    base$acs_unemp_rate_pct <- safe_mean(df$acs_unemp_rate_pct)
  if (!"acs_pct_insured" %in% by_vars)
    base$acs_pct_insured <- safe_mean(df$acs_pct_insured)
  
  # Expand over requested stratifiers
  if ("arf_subtype" %in% by_vars) {
    base <- tidyr::expand_grid(base, arf_subtype = levels(df$arf_subtype))
  }
  if ("sex_category" %in% by_vars) {
    base <- tidyr::expand_grid(base, sex_category = levels(df$sex_category))
  }
  if ("race_ethnicity_simple" %in% by_vars) {
    base <- tidyr::expand_grid(base, race_ethnicity_simple = levels(df$race_ethnicity_simple))
  }
  
  # Coerce to factors when present
  if ("arf_subtype" %in% names(base)) {
    base$arf_subtype <- factor(base$arf_subtype, levels = levels(df$arf_subtype))
  }
  if ("sex_category" %in% names(base)) {
    base$sex_category <- factor(base$sex_category, levels = levels(df$sex_category))
  }
  if ("race_ethnicity_simple" %in% names(base)) {
    base$race_ethnicity_simple <- factor(base$race_ethnicity_simple,
                                         levels = levels(df$race_ethnicity_simple))
  }
  
  tibble::as_tibble(base)
}



# ---- 3) Prediction helpers ----
pred_ci_log <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  tibble(pred = plogis(pr$fit),
         lo   = plogis(pr$fit - 1.96*pr$se.fit),
         hi   = plogis(pr$fit + 1.96*pr$se.fit))
}
pred_ci_nb <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  tibble(pred = exp(pr$fit),
         lo   = exp(pr$fit - 1.96*pr$se.fit),
         hi   = exp(pr$fit + 1.96*pr$se.fit))
}

# ---- 4) Fit models ----
fit_inhosp <- glm(in_hosp_death ~ pm25_mean + no2_10*arf_subtype + age + sex_category +
                    race_ethnicity_simple + svi_overall + acs_median_income +
                    acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                  data = arf_exp, family = binomial())

fit_30d <- glm(death_30d ~ pm25_mean + no2_10*arf_subtype + age + sex_category +
                 race_ethnicity_simple + svi_overall + acs_median_income +
                 acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
               data = arf_exp, family = binomial())

fit_vent <- MASS::glm.nb(vent_hours ~ pm25_mean + no2_10*arf_subtype + age + sex_category +
                           race_ethnicity_simple + svi_overall + acs_median_income +
                           acs_pct_lt_hs + acs_unemp_rate_pct + acs_pct_insured,
                         data = arf_exp)

# ---- 5) Predictions by sex ----
make_grid_from_fit <- function(df, fit, by_vars = character(), no2_vals = seq(0, 6, by = 0.02)) {
  xlv <- fit$xlevels
  # Helper to get levels from model if available, else from df
  levs <- function(var) if (!is.null(xlv[[var]])) xlv[[var]] else levels(droplevels(factor(df[[var]])))
  
  base <- tibble::tibble(
    no2_10               = no2_vals,
    pm25_mean            = mean(df$pm25_mean, na.rm = TRUE),
    age                  = mean(df$age, na.rm = TRUE),
    svi_overall          = mean(df$svi_overall, na.rm = TRUE),
    acs_median_income    = mean(df$acs_median_income, na.rm = TRUE),
    acs_pct_lt_hs        = mean(df$acs_pct_lt_hs, na.rm = TRUE),
    acs_unemp_rate_pct   = mean(df$acs_unemp_rate_pct, na.rm = TRUE),
    acs_pct_insured      = mean(df$acs_pct_insured, na.rm = TRUE)
  )
  
  if (!"sex_category" %in% by_vars)           base$sex_category <- levs("sex_category")[1]
  if (!"race_ethnicity_simple" %in% by_vars)  base$race_ethnicity_simple <- levs("race_ethnicity_simple")[1]
  
  if ("arf_subtype" %in% by_vars)            base <- tidyr::expand_grid(base, arf_subtype = levs("arf_subtype"))
  if ("sex_category" %in% by_vars)           base <- tidyr::expand_grid(base, sex_category = levs("sex_category"))
  if ("race_ethnicity_simple" %in% by_vars)  base <- tidyr::expand_grid(base, race_ethnicity_simple = levs("race_ethnicity_simple"))
  
  if ("arf_subtype" %in% names(base))            base$arf_subtype <- factor(base$arf_subtype, levels = levs("arf_subtype"))
  if ("sex_category" %in% names(base))           base$sex_category <- factor(base$sex_category, levels = levs("sex_category"))
  if ("race_ethnicity_simple" %in% names(base))  base$race_ethnicity_simple <- factor(base$race_ethnicity_simple, levels = levs("race_ethnicity_simple"))
  base
}

# Build grids directly from the fitted model’s levels
grid_sex  <- make_grid_from_fit(arf_exp, fit_inhosp, by_vars = c("arf_subtype","sex_category"))
grid_race <- make_grid_from_fit(arf_exp, fit_inhosp, by_vars = c("arf_subtype","race_ethnicity_simple"))

# Apply to your grids
pred_inhosp_sex <- dplyr::bind_cols(grid_sex,  pred_ci_log(fit_inhosp, grid_sex)) |>
  dplyr::mutate(site_id = site_id, outcome = "in_hosp_death")

pred_30d_sex <- dplyr::bind_cols(grid_sex, pred_ci_log(fit_30d, grid_sex)) |>
  dplyr::mutate(site_id = site_id, outcome = "death_30d")

pred_vent_sex <- dplyr::bind_cols(grid_sex, pred_ci_nb(fit_vent, grid_sex)) |>
  dplyr::mutate(site_id = site_id, outcome = "vent_hours")


# ---- 6) Predictions by race ----
pred_inhosp_race <- bind_cols(grid_race, pred_ci_log(fit_inhosp, grid_race)) %>%
  mutate(site_id = site_id, outcome = "in_hosp_death")

pred_30d_race <- bind_cols(grid_race, pred_ci_log(fit_30d, grid_race)) %>%
  mutate(site_id = site_id, outcome = "death_30d")

pred_vent_race <- bind_cols(grid_race, pred_ci_nb(fit_vent, grid_race)) %>%
  mutate(site_id = site_id, outcome = "vent_hours")

# ---- 7) Write outputs ----
safe_cols <- c("site_id","outcome","no2_10","arf_subtype","sex_category",
               "race_ethnicity_simple","pred","lo","hi")

write_csv(dplyr::select(pred_inhosp_sex, all_of(safe_cols)), file.path(out_dir, "site_preds_inhosp_by_sex.csv"))
write_csv(dplyr::select(pred_30d_sex,    all_of(safe_cols)), file.path(out_dir, "site_preds_30d_by_sex.csv"))
write_csv(dplyr::select(pred_vent_sex,   all_of(safe_cols)), file.path(out_dir, "site_preds_vent_by_sex.csv"))

write_csv(dplyr::select(pred_inhosp_race, all_of(safe_cols)), file.path(out_dir, "site_preds_inhosp_by_race.csv"))
write_csv(dplyr::select(pred_30d_race,    all_of(safe_cols)), file.path(out_dir, "site_preds_30d_by_race.csv"))
write_csv(dplyr::select(pred_vent_race,   all_of(safe_cols)), file.path(out_dir, "site_preds_vent_by_race.csv"))

message("Done. Outputs written to: ", normalizePath(out_dir))


#===================================================
# ROC Model Comparisons for SVI Improvement
#===================================================

scale_exposures <- FALSE

# ==== HELPERS (same as before, with minor robustness tweaks) ====
build_seq_formulas <- function(outcome, base_term, add_terms) {
  terms <- c()
  set_names(
    map(add_terms, function(x) {
      terms <<- c(terms, x)
      as.formula(glue("{outcome} ~ {base_term} + {paste(terms, collapse = ' + ')}"))
    }) |> (\(lst) c(list(as.formula(glue("{outcome} ~ {base_term}"))), lst))(),
    nm = c("SVI", paste0("SVI + ", accumulate(add_terms, ~paste(.x, .y, sep = " + "))))
  )
}

get_auc <- function(formula, data) {
  df <- model.frame(formula, data = data, na.action = na.omit)
  fit <- glm(formula, data = df, family = binomial())
  preds <- predict(fit, type = "response")
  y <- df[[all.vars(formula)[1]]]
  
  # coerce outcome to {0,1} factor with "0" as control, "1" as case
  if (is.logical(y)) y <- as.integer(y)
  if (is.numeric(y)) y <- factor(y, levels = c(0, 1))
  if (is.factor(y) && !identical(levels(y), c("0","1"))) {
    y <- fct_relabel(y, as.character)
    if (!all(levels(y) %in% c("0","1"))) stop("Outcome must be binary (0/1).")
    y <- fct_drop(y)
  }
  
  roc_obj <- pROC::roc(response = y, predictor = preds, levels = c("0","1"), quiet = TRUE)
  list(model = fit, roc = roc_obj, auc = as.numeric(pROC::auc(roc_obj)), n = nrow(df))
}

compare_sequential <- function(roc_list) {
  tibble(
    model_prev = names(roc_list)[-length(roc_list)],
    model_curr = names(roc_list)[-1],
    p_value    = map2_dbl(roc_list[-length(roc_list)], roc_list[-1],
                          ~ pROC::roc.test(.x$roc, .y$roc, method = "delong")$p.value),
    auc_prev   = map_dbl(roc_list[-length(roc_list)], "auc"),
    auc_curr   = map_dbl(roc_list[-1], "auc"),
    delta_auc  = auc_curr - auc_prev
  )
}

# ==== DATA PREP (optional scaling) ====
env_terms_all <- c("no2_mean", "pm25_mean", "tmax_mean", "tmin_mean", "prcp_mean", "vp_mean")

present_terms <- env_terms_all[env_terms_all %in% names(arf_exp)]
missing_terms <- setdiff(env_terms_all, present_terms)
if (length(missing_terms)) message("Skipping missing variables: ", paste(missing_terms, collapse = ", "))

arf_exp_use <- arf_exp
if (scale_exposures && length(present_terms)) {
  arf_exp_use <- arf_exp_use %>%
    mutate(across(all_of(present_terms), ~ as.numeric(scale(.x)), .names = "{.col}"))
}

# ==== IN-HOSPITAL MORTALITY ====
formulas_inhosp <- build_seq_formulas(
  outcome   = "in_hosp_death",
  base_term = "svi_overall",
  add_terms = present_terms
)

aucs_inhosp <- imap(formulas_inhosp, ~ get_auc(.x, arf_exp_use))
aucs_inhosp_df <- tibble(
  model = names(aucs_inhosp),
  auc   = map_dbl(aucs_inhosp, "auc"),
  n     = map_dbl(aucs_inhosp, "n")
)

inhosp_comp <- compare_sequential(aucs_inhosp)

# ==== 30-DAY MORTALITY ====
formulas_30d <- build_seq_formulas(
  outcome   = "death_30d",
  base_term = "svi_overall",
  add_terms = present_terms
)

aucs_30d <- imap(formulas_30d, ~ get_auc(.x, arf_exp_use))
aucs_30d_df <- tibble(
  model = names(aucs_30d),
  auc   = map_dbl(aucs_30d, "auc"),
  n     = map_dbl(aucs_30d, "n")
)

mort30_comp <- compare_sequential(aucs_30d)

# ==== COMBINE & (OPTIONAL) SAVE ====
seq_label <- if (scale_exposures) "zscored" else "raw"
summary_auc <- bind_rows(
  aucs_inhosp_df  %>% mutate(outcome = "in_hosp_death"),
  aucs_30d_df     %>% mutate(outcome = "death_30d")
) %>%
  relocate(outcome)

summary_comp <- bind_rows(
  inhosp_comp %>% mutate(outcome = "in_hosp_death"),
  mort30_comp %>% mutate(outcome = "death_30d")
) %>% relocate(outcome)

print(summary_auc)
print(summary_comp)

# Uncomment to write out
# Correct paths + proper interpolation
readr::write_csv(
  summary_auc,
  file.path(cfg$output_dir, sprintf("auc_sequential_%s.csv", seq_label))
)

readr::write_csv(
  summary_comp,
  file.path(cfg$output_dir, sprintf("auc_diffs_delong_%s.csv", seq_label))
)

# ---- IN-HOSPITAL DEATH ----
# Base plot with first model
plot(aucs_inhosp[[1]]$roc, legacy.axes = TRUE, lwd = 2,
     main = "ROC: In-hospital Death", col = "#1b9e77")
# Add the rest
cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")
for (i in 2:length(aucs_inhosp)) {
  plot(aucs_inhosp[[i]]$roc, add = TRUE, lwd = 2, col = cols[(i - 1) %% length(cols) + 1])
}
legend("topleft",
       legend = paste0(names(aucs_inhosp), " (AUC=", sprintf("%.3f", sapply(aucs_inhosp, `[[`, "auc")), ")"),
       col = cols[seq_along(aucs_inhosp)], lwd = 2, cex = 0.85, bty = "n")

# ---- 30-DAY DEATH ----
plot(aucs_30d[[1]]$roc, legacy.axes = TRUE, lwd = 2,
     main = "ROC: 30-day Death", col = "#1b9e77")
for (i in 2:length(aucs_30d)) {
  plot(aucs_30d[[i]]$roc, add = TRUE, lwd = 2, col = cols[(i - 1) %% length(cols) + 1])
}
legend("topleft",
       legend = paste0(names(aucs_30d), " (AUC=", sprintf("%.3f", sapply(aucs_30d, `[[`, "auc")), ")"),
       col = cols[seq_along(aucs_30d)], lwd = 2, cex = 0.85, bty = "n")


# Helper: convert pROC roc to a tidy df
roc_to_df <- function(roc_obj, label) {
  tibble(
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities)
  ) %>%
    mutate(model = label)
}

# Build tidy data with model names + AUC in label
fmt_lab <- function(name, auc) paste0(name, " (AUC=", sprintf("%.3f", auc), ")")

inhosp_df <- map2(aucs_inhosp, names(aucs_inhosp), ~ roc_to_df(.x$roc, fmt_lab(.y, .x$auc))) %>%
  list_rbind() %>%
  mutate(outcome = "In-hospital Death")

m30_df <- map2(aucs_30d, names(aucs_30d), ~ roc_to_df(.x$roc, fmt_lab(.y, .x$auc))) %>%
  list_rbind() %>%
  mutate(outcome = "30-day Death")

roc_long <- bind_rows(inhosp_df, m30_df)

# Order legend by model complexity (SVI → … → +vp_mean)
model_order <- unique(roc_long$model)
roc_long$model <- factor(roc_long$model, levels = model_order)

# --- Build short labels with proper plotmath tokens
shorten_name <- function(x) {
  x <- str_replace(x, "^SVI \\+ ", "+")
  x <- str_replace_all(x, "no2_mean",  "NO[2]")
  x <- str_replace_all(x, "pm25_mean", "PM[2.5]")
  x <- str_replace_all(x, "tmax_mean", "T[max]")
  x <- str_replace_all(x, "tmin_mean", "T[min]")
  x <- str_replace_all(x, "prcp_mean", "Prcp")
  x <- str_replace_all(x, "vp_mean",   "VP")
  x
}

# Build code → label map, plus parsed plotmath expressions
mk_labels <- function(roc_list) {
  nms   <- names(roc_list)
  codes <- paste0("M", seq_along(nms) - 1)
  aucs  <- map_dbl(roc_list, "auc")
  short <- ifelse(nms == "SVI", "SVI", shorten_name(nms))
  
  # Plotmath expression string: e.g., 'M1:~+NO[2]~plain("(0.536)")'
  expr_str <- sprintf('%s:~%s~plain("(%.3f)")', codes, short, aucs)
  tibble(
    name       = nms,
    code       = codes,
    auc        = aucs,
    label_expr = map(expr_str, ~ parse(text = .x)[[1]])
  )
}

labs_inhosp <- mk_labels(aucs_inhosp)
labs_30d    <- mk_labels(aucs_30d)

# --- Tidy ROC data; use model CODE (M0..M6) as the legend key
roc_to_df <- function(roc_obj, code) {
  tibble(
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities),
    model       = code
  )
}

inhosp_df <- map2(aucs_inhosp, labs_inhosp$code, ~ roc_to_df(.x$roc, .y)) |> list_rbind()
m30_df    <- map2(aucs_30d,    labs_30d$code,    ~ roc_to_df(.x$roc, .y)) |> list_rbind()

# --- Styling (classic theme)
style_plot <- function(p) {
  p +
    theme_classic(base_size = 12) +
    theme(
      legend.position   = "bottom",
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 9),
      legend.key.width  = unit(1.1, "lines"),
      legend.margin     = margin(t = 0, r = 0, b = 0, l = 0),
      plot.margin       = margin(t = 5, r = 5, b = 35, l = 5)
    ) +
    guides(color = guide_legend(ncol = 2, byrow = TRUE))
}

# --- Plots (with parsed legend labels so NO2, PM2.5, T[max], T[min] render correctly)
p_inhosp <- style_plot(
  ggplot(inhosp_df,
         aes(x = 1 - specificity, y = sensitivity, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5) +
    geom_line(linewidth = 1) +
    coord_equal() +
    labs(
      title = "ROC Curves: In-hospital Death",
      x = "1 - Specificity (FPR)",
      y = "Sensitivity (TPR)",
      color = "Models"
    ) +
    scale_color_discrete(
      breaks = labs_inhosp$code,
      labels = labs_inhosp$label_expr
    )
)

p_30d <- style_plot(
  ggplot(m30_df,
         aes(x = 1 - specificity, y = sensitivity, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5) +
    geom_line(linewidth = 1) +
    coord_equal() +
    labs(
      title = "ROC Curves: 30-day Death",
      x = "1 - Specificity (FPR)",
      y = "Sensitivity (TPR)",
      color = "Models"
    ) +
    scale_color_discrete(
      breaks = labs_30d$code,
      labels = labs_30d$label_expr
    )
)

p_inhosp
p_30d

ggsave(
  file.path(cfg$figures_path, "roc_inhosp_classic.png"),
  p_inhosp,
  width = 7, height = 6, dpi = 300
)

ggsave(
  file.path(cfg$figures_path, "roc_30d_classic.png"),
  p_30d,
  width = 7, height = 6, dpi = 300
)


# ============================================================
# Best-performing subtype models with FULL adjustment set
# Base terms are always included; we add Tmax, Tmin, Prcp, VP sequentially
# ============================================================


# ----- Fixed adjustment set (use only those present) -----
base_terms_all <- c(
  "pm25_mean","no2_mean","age","sex_category","race_ethnicity_simple",
  "svi_overall","acs_median_income","acs_pct_lt_hs","acs_unemp_rate_pct","acs_pct_insured"
)
base_terms <- intersect(base_terms_all, names(arf_exp))

# Variables to add sequentially (use only those present)
add_terms_all <- c("tmax_mean","tmin_mean","prcp_mean","vp_mean")
add_terms <- intersect(add_terms_all, names(arf_exp))
if (length(setdiff(add_terms_all, add_terms)))
  message("Skipping missing env vars: ", paste(setdiff(add_terms_all, add_terms), collapse = ", "))

# ----- Helpers -----
get_auc_cc <- function(formula, data, outcome) {
  df <- model.frame(formula, data = data, na.action = na.omit)
  if (!nrow(df)) return(list(roc=NULL, auc=NA_real_, n=0))
  fit <- glm(formula, data = df, family = binomial())
  p   <- predict(fit, type = "response")
  y   <- df[[outcome]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.numeric(y)) y <- factor(y, levels = c(0,1))
  if (is.factor(y) && !identical(levels(y), c("0","1"))) {
    y <- forcats::fct_relabel(y, as.character) |> forcats::fct_drop()
  }
  if (length(levels(y)) < 2) return(list(roc=NULL, auc=NA_real_, n=nrow(df)))
  r <- pROC::roc(response = y, predictor = p, levels = c("0","1"), quiet = TRUE, direction = "<")
  list(roc = r, auc = as.numeric(pROC::auc(r)), n = nrow(df))
}

# Build cumulative forms: outcome ~ (BASE) + sequential add-ons
build_seq_forms_adj <- function(outcome, base_terms, add_terms) {
  rhs_base <- paste(base_terms, collapse = " + ")
  fmls  <- list(as.formula(glue("{outcome} ~ {rhs_base}")))
  nms   <- c("Adj")
  terms <- character()
  for (t in add_terms) {
    terms <- c(terms, t)
    rhs <- paste(c(rhs_base, terms), collapse = " + ")
    fmls[[length(fmls)+1]] <- as.formula(glue("{outcome} ~ {rhs}"))
    nms <- c(nms, paste0("Adj + ", paste(terms, collapse = " + ")))
  }
  names(fmls) <- nms
  fmls
}

# Legend labels with plotmath subscripts
shorten_name <- function(x) {
  x <- str_replace(x, "^Adj \\+ ", "+")
  x <- str_replace_all(x, "no2_mean",  "NO[2]")
  x <- str_replace_all(x, "pm25_mean", "PM[2.5]")
  x <- str_replace_all(x, "tmax_mean", "T[max]")
  x <- str_replace_all(x, "tmin_mean", "T[min]")
  x <- str_replace_all(x, "prcp_mean", "Prcp")
  x <- str_replace_all(x, "vp_mean",   "VP")
  x
}
mk_labels <- function(roc_list) {
  keep <- !vapply(roc_list, function(x) is.null(x$roc) || is.na(x$auc), logical(1))
  rl   <- roc_list[keep]
  nms  <- names(rl)
  codes <- paste0("M", seq_along(nms) - 1)
  aucs  <- map_dbl(rl, "auc")
  short <- ifelse(nms == "Adj", "Adj", shorten_name(nms))
  expr_str <- sprintf('%s:~%s~plain("(%.3f)")', codes, short, aucs)
  tibble(name = nms, code = codes, auc = aucs,
         label_expr = map(expr_str, ~ parse(text = .x)[[1]]))
}
roc_to_df <- function(roc_obj, code) {
  tibble(
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities),
    model = code
  )
}
style_roc <- function(p, ncol_legend = 3) {
  p +
    theme_classic(base_size = 12) +
    theme(
      legend.position   = "bottom",
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 9),
      legend.key.width  = unit(0.9, "lines"),
      legend.margin     = margin(0, 0, 0, 0),
      plot.margin       = margin(6, 8, 46, 8)
    ) +
    guides(color = guide_legend(ncol = ncol_legend, byrow = TRUE))
}

best_by_subtype <- function(outcome_var, pretty_outcome, out_stub) {
  stopifnot("arf_subtype" %in% names(arf_exp))
  subtypes <- levels(droplevels(factor(arf_exp$arf_subtype)))
  forms <- build_seq_forms_adj(outcome_var, base_terms, add_terms)
  
  all_auc_rows <- list(); best_rows <- list()
  roc_dir <- file.path("roc_by_subtype", out_stub)
  
  for (st in subtypes) {
    dat <- arf_exp |> filter(arf_subtype == st)
    aucs <- imap(forms, ~ get_auc_cc(.x, dat, outcome_var))
    auc_df <- tibble(
      subtype = st,
      model   = names(forms),
      auc     = map_dbl(aucs, "auc"),
      n       = map_dbl(aucs, "n")
    )
    all_auc_rows[[st]] <- auc_df
    
    valid <- which(!is.na(auc_df$auc))
    if (!length(valid)) {
      best_rows[[st]] <- tibble(subtype = st, best_model = NA_character_, auc = NA_real_, n = 0L)
      next
    }
    max_auc <- max(auc_df$auc[valid])
    cand    <- valid[which(abs(auc_df$auc[valid] - max_auc) < 1e-12)]
    best_ix <- min(cand) # tie -> simpler model
    best_row <- auc_df[best_ix, , drop = FALSE]
    best_rows[[st]] <- best_row |> transmute(subtype, best_model = model, auc, n)
    
    # --- Build ROC overlay for this subtype (highlight best) ---
    # Keep only models with defined ROC
    keep_names <- names(forms)[valid]
    roc_list   <- aucs[valid]; names(roc_list) <- keep_names
    
    # Legend labels w/ AUC embedded
    labs_sub <- mk_labels(roc_list)
    
    # Tidy ROC data
    df_roc <- map2(roc_list, labs_sub$code, ~ roc_to_df(.x$roc, .y)) |> list_rbind()
    
    # Correct: best model name is in `best_row$model` (not best_row$best_model)
    best_model_name <- best_row$model
    idx <- match(best_model_name, labs_sub$name)
    best_code <- if (!is.na(idx)) labs_sub$code[idx] else NA_character_
    
    if (nrow(df_roc) == 0L) {
      next  # nothing to plot for this subtype
    }
    
    df_roc <- df_roc |>
      dplyr::mutate(
        is_best = if (!is.na(best_code)) model == best_code else FALSE,
        lw      = ifelse(is_best, 1.6, 0.9),
        alph    = ifelse(is_best, 1.0, 0.85)
      )
    
    p <- style_roc(
      ggplot(df_roc, aes(x = 1 - specificity, y = sensitivity, color = model)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5) +
        geom_line(aes(linewidth = lw, alpha = alph)) +
        scale_linewidth_identity() +
        scale_alpha_identity() +
        coord_equal() +
        labs(
          title = glue("{pretty_outcome}: ROC — {st}"),
          x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = "Models"
        ) +
        scale_color_discrete(
          breaks = labs_sub$code,
          labels = labs_sub$label_expr
        ),
      ncol_legend = 3
    )
    
    # save under <repo>/output/figures/roc_by_subtype/<out_stub>/...
    save_plot(p, file.path(roc_dir, paste0("roc_", out_stub, "_", gsub("[^[:alnum:]_]+", "", st))))
    
  }
  
  auc_all  <- bind_rows(all_auc_rows)  |> arrange(subtype, match(model, names(forms)))
  best_all <- bind_rows(best_rows)     |> arrange(subtype)
  
  # plain ASCII model names for CSV
  plainify <- function(x) x |>
    str_replace("^Adj \\+ ", "+") |>
    str_replace_all("no2_mean","NO2") |>
    str_replace_all("pm25_mean","PM2.5") |>
    str_replace_all("tmax_mean","Tmax") |>
    str_replace_all("tmin_mean","Tmin") |>
    str_replace_all("prcp_mean","Prcp") |>
    str_replace_all("vp_mean","VP")
  
  best_all <- best_all |> mutate(best_model_pretty = ifelse(best_model=="Adj","Adj", plainify(best_model)))
  
  save_tbl(auc_all,  file.path("metrics", paste0("auc_all_models_by_subtype_", out_stub)))
  save_tbl(best_all, file.path("metrics", paste0("best_models_by_subtype_",   out_stub)))
  
  invisible(list(auc_all = auc_all, best = best_all))
}

# =========================
# Run (both outcomes)
# =========================
res_inhosp <- best_by_subtype("in_hosp_death", "In-hospital Death", "inhosp")
res_30d    <- best_by_subtype("death_30d",     "30-day Death",      "30d")

print(res_inhosp$best)
print(res_30d$best)

# ============================================================
# Forward stepwise AUC search per ARF subtype (mortality)
# Finds a compact, high-AUC model for each subtype & outcome
# ============================================================

# ---------- KNOBS ----------
# Start from SVI-only (change to NULL to start from intercept-only)
base_term <- "no2_mean"

# Full slate of potential predictors (we'll intersect with columns present)
candidates_all <- c(
  # demographics / social
  "age", "sex_category", "race_ethnicity_simple",
  "svi_overall", "acs_median_income", "acs_pct_lt_hs",
  "acs_unemp_rate_pct", "acs_pct_insured",
  # pollutants
  "no2_mean", "pm25_mean",
  # daymet
  "tmax_mean", "tmin_mean", "prcp_mean", "vp_mean"
)

# Stepwise limits / stopping rules
max_terms <- 5        # max number of terms to add beyond the base (keep it clinical)
min_delta <- 0.002    # minimum AUC gain to consider meaningful
p_thresh  <- 0.05     # DeLong p-value threshold vs previous step

# ---------- AUC helper on complete cases ----------
get_auc_cc <- function(formula, data, outcome) {
  df <- model.frame(formula, data = data, na.action = na.omit)
  if (!nrow(df)) return(list(roc=NULL, auc=NA_real_, n=0, df=df))
  fit <- glm(formula, data = df, family = binomial())
  p   <- predict(fit, type = "response")
  y   <- df[[outcome]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.numeric(y)) y <- factor(y, levels = c(0,1))
  if (is.factor(y) && !identical(levels(y), c("0","1"))) {
    y <- forcats::fct_relabel(y, as.character) |> forcats::fct_drop()
  }
  if (length(levels(y)) < 2) return(list(roc=NULL, auc=NA_real_, n=nrow(df), df=df))
  r <- pROC::roc(response = y, predictor = p, levels = c("0","1"), quiet = TRUE, direction = "<")
  list(roc = r, auc = as.numeric(pROC::auc(r)), n = nrow(df), df = df)
}

# ---------- forward stepwise search by AUC (greedy) ----------
stepwise_search <- function(data, outcome, base_term, candidates, max_terms, min_delta, p_thresh) {
  # start terms
  terms <- character(0)
  if (!is.null(base_term) && base_term %in% names(data)) terms <- base_term
  
  # compute base
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
  f_curr <- as.formula(glue("{outcome} ~ {rhs}"))
  res_curr <- get_auc_cc(f_curr, data, outcome)
  if (is.na(res_curr$auc)) return(list(path = tibble(), best = NULL, rocs = list()))
  
  steps <- tibble(step = 0L, added = if (length(terms)) "(base)" else "(intercept)",
                  formula = deparse1(f_curr), auc = res_curr$auc, n = res_curr$n,
                  delta_auc = NA_real_, p_value = NA_real_)
  
  remaining <- setdiff(candidates, terms)
  rocs <- list(`M0` = res_curr$roc)
  names_map <- list(`M0` = steps$formula[1])
  
  # iterate
  k <- 0L
  while (length(remaining) > 0 && k < max_terms) {
    # try adding each remaining term
    trial <- map(remaining, function(add) {
      rhs2 <- paste(c(terms, add), collapse = " + ")
      f <- as.formula(glue("{outcome} ~ {rhs2}"))
      list(term = add, res = get_auc_cc(f, data, outcome), fml = f)
    })
    
    # keep only valid AUCs
    trial <- keep(trial, ~ !is.na(.x$res$auc))
    if (!length(trial)) break
    
    # pick best by AUC
    aucs <- map_dbl(trial, ~ .x$res$auc)
    best_idx <- which.max(aucs)
    cand     <- trial[[best_idx]]
    
    # compare with DeLong vs current (handle differing N)
    p_val <- tryCatch(pROC::roc.test(res_curr$roc, cand$res$roc, method = "delong")$p.value,
                      error = function(e) NA_real_)
    d_auc <- cand$res$auc - res_curr$auc
    
    # stopping rule
    if (!is.na(p_val) && p_val > p_thresh && (is.na(d_auc) || d_auc < min_delta)) break
    
    # accept
    terms <- c(terms, cand$term)
    res_curr <- cand$res
    k <- k + 1L
    
    steps <- add_row(
      steps,
      step = k,
      added = cand$term,
      formula = deparse1(cand$fml),
      auc = cand$res$auc,
      n = cand$res$n,
      delta_auc = d_auc,
      p_value = p_val
    )
    
    rocs[[paste0("M", k)]] <- cand$res$roc
    names_map[[paste0("M", k)]] <- deparse1(cand$fml)
    
    remaining <- setdiff(remaining, cand$term)
  }
  
  list(path = steps, best = tail(steps, 1), rocs = rocs, names_map = names_map)
}

# ---------- pretty legend labels (plotmath) ----------
shorten_name <- function(x) {
  x %>%
    str_replace("^.*~", "") %>%                 # strip outcome LHS
    str_replace("^\\s*", "") %>%
    str_replace("^1$", "Intercept") %>%
    str_replace_all("no2_mean",  "NO[2]") %>%
    str_replace_all("pm25_mean", "PM[2.5]") %>%
    str_replace_all("tmax_mean", "T[max]") %>%
    str_replace_all("tmin_mean", "T[min]") %>%
    str_replace_all("prcp_mean", "Prcp") %>%
    str_replace_all("vp_mean",   "VP") %>%
    str_replace_all("\\s*\\+\\s*", " + ")
}
mk_labels_from_map <- function(rocs, names_map) {
  codes <- names(rocs)
  nm    <- unlist(names_map[codes], use.names = FALSE)
  # extract RHS and compress for label (“Adj path” style)
  rhs_short <- ifelse(grepl("~ 1$", nm), "Intercept",
                      shorten_name(sub(".*~", "", nm)))
  # compute AUCs for printing
  aucs <- map_dbl(rocs, ~ as.numeric(pROC::auc(.x)))
  expr_str <- sprintf('%s:~%s~plain("(%.3f)")', codes, rhs_short, aucs)
  tibble(code = codes, expr = map(expr_str, ~ parse(text = .x)[[1]]))
}
roc_to_df <- function(roc_obj, code) {
  tibble(
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities),
    model = code
  )
}
style_roc <- function(p, ncol_legend = 3) {
  p +
    theme_classic(base_size = 12) +
    theme(
      legend.position   = "bottom",
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 9),
      legend.key.width  = unit(0.9, "lines"),
      legend.margin     = margin(0, 0, 0, 0),
      plot.margin       = margin(6, 8, 46, 8)
    ) +
    guides(color = guide_legend(ncol = ncol_legend, byrow = TRUE))
}

# ---------- Runner per outcome across subtypes ----------
stepwise_by_subtype <- function(outcome_var, pretty_outcome, out_stub) {
  stopifnot("arf_subtype" %in% names(arf_exp))
  subtypes <- levels(droplevels(factor(arf_exp$arf_subtype)))
  
  # available candidates in data
  present <- intersect(candidates_all, names(arf_exp))
  # remove base from candidates so we don't add it twice
  candidates <- if (!is.null(base_term)) setdiff(present, base_term) else present
  
  all_paths <- list()
  best_rows <- list()
  
  roc_dir <- file.path("roc_by_subtype_stepwise", out_stub)
  
  for (st in subtypes) {
    dat <- arf_exp |> filter(arf_subtype == st)
    sw  <- stepwise_search(dat, outcome_var, base_term, candidates, max_terms, min_delta, p_thresh)
    
    if (nrow(sw$path) == 0) {
      best_rows[[st]] <- tibble(subtype = st, best_step = NA_integer_,
                                best_formula = NA_character_, auc = NA_real_, n = 0L)
      next
    }
    
    # record path & best
    path <- sw$path |> mutate(subtype = st, outcome = outcome_var, .before = 1)
    all_paths[[st]] <- path
    best <- sw$best
    best_rows[[st]] <- tibble(
      subtype = st, outcome = outcome_var,
      best_step = best$step, best_formula = best$formula,
      auc = best$auc, n = best$n
    )
    
    # ----- ROC overlay for the stepwise path -----
    # convert ROCs to tidy for plotting
    labs_df <- mk_labels_from_map(sw$rocs, sw$names_map)
    df_roc <- imap(sw$rocs, ~ roc_to_df(.x, .y)) |> list_rbind()
    best_code <- paste0("M", best$step)
    
    df_roc <- df_roc |>
      mutate(lw = ifelse(model == best_code, 1.6, 0.9),
             alph = ifelse(model == best_code, 1.0, 0.85))
    
    p <- style_roc(
      ggplot(df_roc, aes(x = 1 - specificity, y = sensitivity, color = model)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5) +
        geom_line(aes(linewidth = lw, alpha = alph)) +
        scale_linewidth_identity() + scale_alpha_identity() +
        coord_equal() +
        labs(
          title = glue("{pretty_outcome}: ROC - {st}"),
          x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = "Models"
        ) +
        scale_color_discrete(
          breaks = labs_df$code,
          labels = labs_df$expr
        ),
      ncol_legend = 1
    )
    
    # save figure: <repo>/output/figures/roc_by_subtype_stepwise/<out_stub>/...
    save_plot(p, file.path(roc_dir, paste0("stepwise_roc_", out_stub, "_", gsub("[^[:alnum:]_]+", "", st))))
  }
  
  # bind & save tables
  paths_tbl <- bind_rows(all_paths)
  best_tbl  <- bind_rows(best_rows)
  
  save_tbl(paths_tbl, file.path("metrics", paste0("stepwise_paths_", out_stub)))
  save_tbl(best_tbl,  file.path("metrics", paste0("stepwise_best_",  out_stub)))
  
  invisible(list(paths = paths_tbl, best = best_tbl))
}

# =========================
# Run for both outcomes
# =========================
res_inhosp_sw <- stepwise_by_subtype("in_hosp_death", "In-hospital Death", "inhosp")
res_30d_sw    <- stepwise_by_subtype("death_30d",     "30-day Death",      "30d")

# Quick peek
print(res_inhosp_sw$best)
print(res_30d_sw$best)


# ------------------------------------ Stratify by AGE CATEGORY ------------------------------------

pred_ci_log <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  tibble::as_tibble(dplyr::bind_cols(
    newdata,
    pred = plogis(pr$fit),
    lo   = plogis(pr$fit - 1.96 * pr$se.fit),
    hi   = plogis(pr$fit + 1.96 * pr$se.fit)
  ))
}
pred_ci_nb <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  tibble::as_tibble(dplyr::bind_cols(
    newdata,
    pred = exp(pr$fit),
    lo   = exp(pr$fit - 1.96 * pr$se.fit),
    hi   = exp(pr$fit + 1.96 * pr$se.fit)
  ))
}

# 1) Create age_category (adjust cutpoints as you like)
arf_exp <- arf_exp |>
  mutate(
    age_category = cut(
      age,
      breaks = c(-Inf, 49, 64, 79, Inf),          # <50, 50–64, 65–79, ≥80
      labels = c("<50", "50–64", "65–79", "≥80"),
      right = TRUE, include.lowest = TRUE
    ) |> droplevels()
  )

# Ensure factors are clean
arf_exp <- arf_exp |>
  mutate(
    arf_subtype = droplevels(arf_subtype),
    sex_category = droplevels(sex_category),
    race_ethnicity_simple = droplevels(race_ethnicity_simple),
    age_category = droplevels(age_category)
  )

# 2) Fit interaction models: no2 (per 10 ppb) x age_category
fit_mort_inhosp_age <- glm(
  in_hosp_death ~ pm25_mean + no2_10*age_category + arf_subtype + sex_category +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)

fit_mort_30d_age <- glm(
  death_30d ~ pm25_mean + no2_10*age_category + arf_subtype + sex_category +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)

fit_vent_nb_age <- MASS::glm.nb(
  vent_hours ~ pm25_mean + no2_10*age_category + arf_subtype + sex_category +
    race_ethnicity_simple + svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp
)

# (optional) Save tidy model tables
tidy_and_save(fit_mort_inhosp_age, "inhosp_death_x_agecat", exponentiate = TRUE)
tidy_and_save(fit_mort_30d_age,    "death30d_x_agecat",    exponentiate = TRUE)
tidy_and_save(fit_vent_nb_age,     "vent_hours_x_agecat",  exponentiate = TRUE)

# 3) Build prediction grid (vary no2_10 across observed range; facet by age_category)
age_lvls <- fit_mort_inhosp_age$xlevels$age_category
sub_lvls <- fit_mort_inhosp_age$xlevels$arf_subtype

make_grid_age <- function(df, age_lvls, sub_lvls, n = 150) {
  safe_mean <- function(x) { m <- mean(x, na.rm = TRUE); if (is.nan(m)) NA_real_ else m }
  top_level <- function(x) {
    if (is.factor(x)) x <- droplevels(x)
    tab <- sort(table(x), decreasing = TRUE)
    if (length(tab) == 0) NA_character_ else names(tab)[1]
  }
  rng <- df %>% summarize(
    lo = quantile(no2_10, 0.01, na.rm = TRUE),
    hi = quantile(no2_10, 0.99, na.rm = TRUE)
  )
  expand.grid(
    no2_10 = seq(rng$lo, rng$hi, length.out = n),
    age_category = age_lvls,
    arf_subtype  = sub_lvls,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) |>
    as_tibble() |>
    mutate(
      pm25_mean             = safe_mean(df$pm25_mean),
      sex_category          = top_level(df$sex_category),
      race_ethnicity_simple = top_level(df$race_ethnicity_simple),
      svi_overall           = safe_mean(df$svi_overall),
      acs_median_income     = safe_mean(df$acs_median_income),
      acs_pct_lt_hs         = safe_mean(df$acs_pct_lt_hs),
      acs_unemp_rate_pct    = safe_mean(df$acs_unemp_rate_pct),
      acs_pct_insured       = safe_mean(df$acs_pct_insured),
      
      # coerce to model factor levels
      age_category          = factor(age_category, levels = levels(df$age_category)),
      arf_subtype           = factor(arf_subtype,  levels = levels(df$arf_subtype)),
      sex_category          = factor(sex_category, levels = levels(df$sex_category)),
      race_ethnicity_simple = factor(race_ethnicity_simple, levels = levels(df$race_ethnicity_simple))
    )
}

grid_age <- make_grid_age(arf_exp, age_lvls, sub_lvls)

# rebuild the three frames using the updated helpers (or keep the quick patch above)
df_mort_inhosp_age <- pred_ci_log(fit_mort_inhosp_age, grid_age)
df_mort_30d_age    <- pred_ci_log(fit_mort_30d_age,    grid_age)
df_vent_age        <- pred_ci_nb (fit_vent_nb_age,     grid_age)

# sanity
stopifnot(all(c("age_category","arf_subtype","no2_10") %in% names(df_mort_inhosp_age)))

p1_age <- ggplot(df_mort_inhosp_age, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ age_category) + pal + theme_pub +
  labs(title = "In-hospital mortality vs NO₂, stratified by age",
       x = "NO₂ (per 10 ppb)", y = "Predicted probability", color = "ARF subtype")

p2_age <- ggplot(df_mort_30d_age, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ age_category) + pal + theme_pub +
  labs(title = "30-day mortality vs NO₂, stratified by age",
       x = "NO₂ (per 10 ppb)", y = "Predicted probability", color = "ARF subtype")

p3_age <- ggplot(df_vent_age, aes(no2_10, pred, color = arf_subtype)) +
  geom_line(size = 1.2) +
  geom_line(aes(y = lo), linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = hi), linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ age_category) + pal + theme_pub +
  labs(title = "Ventilation hours vs NO₂, stratified by age",
       x = "NO₂ (per 10 ppb)", y = "Predicted mean hours", color = "ARF subtype")

combo_age <- (p1_age | p2_age) / p3_age + plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_plot(combo_age, "no2_by_agecat_subtype_combined")

export_glm_spec <- function(fit, name, out_dir = file.path(cfg$output_dir, "federated_specs")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fam <- family(fit)
  spec <- list(
    model_name = name,
    formula    = deparse1(formula(fit)),
    family     = fam$family,
    link       = fam$link,
    coef       = coef(fit),                 # named vector w/ model.matrix colnames
    xlevels    = fit$xlevels,               # factor levels used during fit
    contrasts  = fit$contrasts,             # usually NULL if defaults
    theta      = if (!is.null(fit$theta)) unname(fit$theta) else NULL
  )
  fp <- file.path(out_dir, paste0(cfg$prefix, "_", name, "_spec_", cfg$run_id, ".json"))
  writeLines(toJSON(spec, pretty = TRUE, auto_unbox = TRUE), fp)
  invisible(fp)
}

# export specs for your three adjusted models (examples)
export_glm_spec(fit_mort_adj,   "inhosp_death_adj")
export_glm_spec(fit_mort30_adj, "death30d_adj")
export_glm_spec(fit_vent_nb,    "vent_hours_adj")

arf_exp <- arf_exp |>
  mutate(
    covid_period = case_when(
      index_admit < as.Date("2020-03-01")              ~ "Pre-COVID",
      index_admit >= as.Date("2020-03-01") & index_admit <= as.Date("2021-12-31") ~ "During COVID",
      index_admit >= as.Date("2022-01-01")             ~ "Post-COVID",
      TRUE ~ NA_character_
    ) |> factor(levels = c("Pre-COVID", "During COVID", "Post-COVID"))
  )


# ========================== COVID ANALYSIS PIPELINE ===========================


# --- 0) Helpers this block relies on (safe fallbacks if not already defined) --
`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!exists("pred_ci_log", mode = "function")) {
  pred_ci_log <- function(fit, newdata) {
    pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
    tibble::as_tibble(dplyr::bind_cols(
      newdata,
      pred = plogis(pr$fit),
      lo   = plogis(pr$fit - 1.96 * pr$se.fit),
      hi   = plogis(pr$fit + 1.96 * pr$se.fit)
    ))
  }
}
if (!exists("pred_ci_nb", mode = "function")) {
  pred_ci_nb <- function(fit, newdata) {
    pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
    tibble::as_tibble(dplyr::bind_cols(
      newdata,
      pred = exp(pr$fit),
      lo   = exp(pr$fit - 1.96 * pr$se.fit),
      hi   = exp(pr$fit + 1.96 * pr$se.fit)
    ))
  }
}
if (!exists("theme_pub")) {
  theme_pub <- theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))
}
if (!exists("pal"))      pal      <- scale_color_brewer(palette = "Dark2")
if (!exists("fill_pal")) fill_pal <- scale_fill_brewer(palette = "Dark2")

# --- 1) Derive covid_period robustly from available date columns --------------
find_index_date <- function(df) {
  cand <- c("admission_date", "admission_dttm", "admit_date",
            "hospitalization_admission_date", "index_admit",
            "encounter_start_date")
  have <- cand[cand %in% names(df)]
  if (length(have) == 0) stop("No admission/index date column found.")
  x <- df[[have[1]]]
  if (inherits(x, "POSIXt")) as.Date(x) else as.Date(x)
}

arf_exp <- arf_exp |>
  mutate(
    .index_date = find_index_date(cur_data_all()),
    covid_period = case_when(
      .index_date <  as.Date("2020-03-01")                     ~ "Pre-COVID",
      .index_date >= as.Date("2020-03-01") &
        .index_date <= as.Date("2021-12-31")                     ~ "During COVID",
      .index_date >= as.Date("2022-01-01")                     ~ "Post-COVID",
      TRUE                                                     ~ NA_character_
    ),
    covid_period = factor(covid_period, levels = c("Pre-COVID","During COVID","Post-COVID"))
  )

# Ensure NO2 per 10 ppb exists
arf_exp <- arf_exp %>% mutate(no2_10 = if (!"no2_10" %in% names(.)) no2_mean/10 else no2_10)

# Reference categories used in plots/grids
ref_sex <- arf_exp %>% count(sex_category, sort = TRUE) %>% slice(1) %>% pull(sex_category)
ref_re  <- arf_exp %>% count(race_ethnicity_simple, sort = TRUE) %>% slice(1) %>% pull(race_ethnicity_simple)

# Drop unused levels to stabilize model.matrix
arf_exp <- arf_exp %>%
  mutate(
    arf_subtype = droplevels(arf_subtype),
    sex_category = droplevels(sex_category),
    race_ethnicity_simple = droplevels(race_ethnicity_simple),
    covid_period = droplevels(covid_period)
  )

# --- 2) Models that ADJUST for covid_period (covariate) -----------------------
fit_mort_inhosp_covid <- glm(
  in_hosp_death ~ pm25_mean + no2_10 + covid_period +
    age + sex_category + race_ethnicity_simple + arf_subtype +
    svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)
fit_mort_30d_covid <- glm(
  death_30d ~ pm25_mean + no2_10 + covid_period +
    age + sex_category + race_ethnicity_simple + arf_subtype +
    svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)
fit_icu_nb_covid <- glm.nb(
  icu_los_days ~ pm25_mean + no2_10 + covid_period +
    age + sex_category + race_ethnicity_simple + arf_subtype +
    svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp
)
fit_vent_nb_covid <- glm.nb(
  vent_hours ~ pm25_mean + no2_10 + covid_period +
    age + sex_category + race_ethnicity_simple + arf_subtype +
    svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp
)

tidy_and_save(fit_mort_inhosp_covid, "inhosp_death_adj_plus_covid", exponentiate = TRUE)
tidy_and_save(fit_mort_30d_covid,    "death30d_adj_plus_covid",    exponentiate = TRUE)
tidy_and_save(fit_icu_nb_covid,      "icu_los_adj_plus_covid",     exponentiate = TRUE)
tidy_and_save(fit_vent_nb_covid,     "vent_hours_adj_plus_covid",  exponentiate = TRUE)

# --- 3) Models with NO2 × covid_period INTERACTION ---------------------------
fit_mort_inhosp_covid_int <- glm(
  in_hosp_death ~ pm25_mean + no2_10*covid_period +
    age + sex_category + race_ethnicity_simple + arf_subtype +
    svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)
fit_mort_30d_covid_int <- glm(
  death_30d ~ pm25_mean + no2_10*covid_period +
    age + sex_category + race_ethnicity_simple + arf_subtype +
    svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp, family = binomial()
)
fit_icu_nb_covid_int <- glm.nb(
  icu_los_days ~ pm25_mean + no2_10*covid_period +
    age + sex_category + race_ethnicity_simple + arf_subtype +
    svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp
)
fit_vent_nb_covid_int <- glm.nb(
  vent_hours ~ pm25_mean + no2_10*covid_period +
    age + sex_category + race_ethnicity_simple + arf_subtype +
    svi_overall + acs_median_income + acs_pct_lt_hs +
    acs_unemp_rate_pct + acs_pct_insured,
  data = arf_exp
)

tidy_and_save(fit_mort_inhosp_covid_int, "inhosp_death_no2_x_covid", exponentiate = TRUE)
tidy_and_save(fit_mort_30d_covid_int,    "death30d_no2_x_covid",    exponentiate = TRUE)
tidy_and_save(fit_icu_nb_covid_int,      "icu_los_no2_x_covid",     exponentiate = TRUE)
tidy_and_save(fit_vent_nb_covid_int,     "vent_hours_no2_x_covid",  exponentiate = TRUE)

# --- 4) Prediction grid & plots by COVID period (using the *interaction* fits) -
safe_mean <- function(x) { m <- mean(x, na.rm = TRUE); if (is.nan(m)) NA_real_ else m }
top_level <- function(x) { if (is.factor(x)) x <- droplevels(x); nm <- names(sort(table(x), TRUE)); if (length(nm)) nm[1] else NA_character_ }

cp_lvls <- fit_mort_inhosp_covid_int$xlevels$covid_period
sub_lvls <- fit_mort_inhosp_covid_int$xlevels$arf_subtype

rng_no2 <- arf_exp %>% summarize(lo = quantile(no2_10, 0.01, na.rm = TRUE),
                                 hi = quantile(no2_10, 0.99, na.rm = TRUE))

make_grid_covid <- function(df, cp_lvls, sub_lvls, n = 150) {
  expand.grid(
    no2_10 = seq(rng_no2$lo, rng_no2$hi, length.out = n),
    covid_period = cp_lvls,
    arf_subtype  = sub_lvls,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) %>%
    tibble::as_tibble() %>%
    mutate(
      pm25_mean             = safe_mean(df$pm25_mean),
      age                   = safe_mean(df$age),
      sex_category          = top_level(df$sex_category),
      race_ethnicity_simple = top_level(df$race_ethnicity_simple),
      svi_overall           = safe_mean(df$svi_overall),
      acs_median_income     = safe_mean(df$acs_median_income),
      acs_pct_lt_hs         = safe_mean(df$acs_pct_lt_hs),
      acs_unemp_rate_pct    = safe_mean(df$acs_unemp_rate_pct),
      acs_pct_insured       = safe_mean(df$acs_pct_insured),
      covid_period          = factor(covid_period, levels = levels(df$covid_period)),
      arf_subtype           = factor(arf_subtype,  levels = levels(df$arf_subtype)),
      sex_category          = factor(sex_category, levels = levels(df$sex_category)),
      race_ethnicity_simple = factor(race_ethnicity_simple, levels = levels(df$race_ethnicity_simple))
    )
}

grid_covid <- make_grid_covid(arf_exp, cp_lvls, sub_lvls)

# Predictions
df_mort_inhosp_covid <- pred_ci_log(fit_mort_inhosp_covid_int, grid_covid)
df_mort_30d_covid    <- pred_ci_log(fit_mort_30d_covid_int,    grid_covid)
df_icu_covid         <- pred_ci_nb (fit_icu_nb_covid_int,      grid_covid)
df_vent_covid        <- pred_ci_nb (fit_vent_nb_covid_int,     grid_covid)

# Plots: facet by covid_period, color by ARF subtype
p1_cov <- ggplot(df_mort_inhosp_covid, aes(no2_10, pred, color = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = arf_subtype), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ covid_period) + pal + fill_pal + theme_pub +
  labs(title = "In-hospital mortality vs NO₂ by COVID period",
       x = "NO₂ (per 10 ppb)", y = "Predicted probability",
       color = "ARF subtype", fill = "ARF subtype")

p2_cov <- ggplot(df_mort_30d_covid, aes(no2_10, pred, color = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = arf_subtype), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ covid_period) + pal + fill_pal + theme_pub +
  labs(title = "30-day mortality vs NO₂ by COVID period",
       x = "NO₂ (per 10 ppb)", y = "Predicted probability",
       color = "ARF subtype", fill = "ARF subtype")

p3_cov <- ggplot(df_icu_covid, aes(no2_10, pred, color = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = arf_subtype), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ covid_period) + pal + fill_pal + theme_pub +
  labs(title = "ICU LOS vs NO₂ by COVID period",
       x = "NO₂ (per 10 ppb)", y = "Predicted mean days",
       color = "ARF subtype", fill = "ARF subtype")

p4_cov <- ggplot(df_vent_covid, aes(no2_10, pred, color = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = arf_subtype), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ covid_period) + pal + fill_pal + theme_pub +
  labs(title = "Ventilation hours vs NO₂ by COVID period",
       x = "NO₂ (per 10 ppb)", y = "Predicted mean hours",
       color = "ARF subtype", fill = "ARF subtype")

combo_covid <- (p1_cov | p2_cov) / (p3_cov | p4_cov) + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot(combo_covid, "no2_by_covid_period_combined")

# --- 5) AUCs overall and period-stratified (for adjusted + COVID covariate) ---
get_model_auc <- get_model_auc  # use your existing version

# Overall AUC (adjusted + covid covariate)
auc_covid_overall <- dplyr::bind_rows(
  get_model_auc(fit_mort_inhosp_covid, arf_exp, "in_hosp_death", "In-hospital death | + COVID period"),
  get_model_auc(fit_mort_30d_covid,    arf_exp, "death_30d",     "30-day death | + COVID period")
)

# Period-specific AUCs
auc_covid_by_period <- purrr::map_dfr(levels(arf_exp$covid_period), function(cp) {
  df_cp <- arf_exp %>% filter(covid_period == cp)
  dplyr::bind_rows(
    get_model_auc(fit_mort_inhosp_covid, df_cp, "in_hosp_death", paste0("In-hospital death | ", cp)),
    get_model_auc(fit_mort_30d_covid,    df_cp, "death_30d",     paste0("30-day death | ", cp))
  ) %>%
    mutate(covid_period = cp)
})

save_tbl(auc_covid_overall,   "metrics/auc_adjusted_plus_covid_overall")
save_tbl(auc_covid_by_period, "metrics/auc_adjusted_plus_covid_by_period")



# =========================
# Wrap up
# =========================

message("\n🎯 Cohort identification, linkage, and analyses complete!")
message("📂 All outputs saved under: ", normalizePath(cfg$output_dir))
message("   Figures: ", normalizePath(cfg$figures_path))
message("🔒 Reminder: Remove any local CSVs with patient-level data outside the repo, and upload the 'output' folder to Box.\n")


