# ================================================================================================
# ICU REspiratory Failure Environmental Risk (REFER) Index | PI: Peter Graffy (graffy@uchicago.edu)
# Minimal Cohort Builder (ICU REFER) — Inclusion/Exclusion Only
# Years: 2018–2024; Adults (≥18); ARF criteria within ±24h of ICU admit
# Outputs:
#   cohort_inclusion         : 1 row / hospitalization that meets inclusion (used for 02-analysis)
#   exclusion_breakdown      : hospitalizations with reason for exclusion
#   exclusions_raw           : raw file for spot checking exclusions
#   selection_flow_counts    : flow counts for inclusion criteria
#   2 figures                : selection and exclusions bar charts
# =================================================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(knitr)
  library(fst)
  library(here)
  library(tidyverse)
  library(arrow)
  library(gtsummary)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(fuzzyjoin)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(glue)
  library(scales)
})

# Objective: identify a cohort of hospitalizations from CLIF tables for acute respiratory failure analysis

# Specify inpatient cohort parameters
start_date <- "2018-01-01"
end_date <- "2024-12-31"
include_pediatric <- FALSE
include_er_deaths <- TRUE

# Specify required CLIF tables (updated for CLIF v2.1 and this project)
tables <- c("patient", "hospitalization", "vitals", "labs", 
            "medication_admin_continuous", "adt", 
            "respiratory_support", "hospital_diagnosis", 
            "microbiology_culture")

# Load configuration utility
source("utils/config.R")
site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type

print(paste("Site Name:", site_name))
print(paste("Tables Path:", tables_path))
print(paste("File Type:", file_type))

# --- Config sanity checks ---
stopifnot(exists("config"))
tables_path <- normalizePath(config$tables_path, mustWork = TRUE)

# Allow multiple extensions from config, e.g. "csv/parquet/fst" or "csv"
exts <- strsplit(config$file_type, "[/|,; ]+")[[1]]
exts <- exts[nzchar(exts)]
if (length(exts) == 0) exts <- c("csv","parquet","fst")

# Build a pattern that matches any of the extensions
ext_pat <- paste0("\\.(", paste(unique(exts), collapse = "|"), ")$")

# Look for CLIF-ish filenames in this folder OR subfolders
all_files <- list.files(
  path = tables_path,
  pattern = ext_pat,
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(all_files) == 0) {
  stop("No files with extensions {", paste(exts, collapse = ", "), 
       "} were found under: ", tables_path)
}

# Keep only files whose base name looks like "clif_*.ext" OR matches "*patient.*" for safety
bn <- basename(all_files)
looks_clif <- grepl("^clif_.*", bn, ignore.case = TRUE)

# For matching to required names, remove extension and (optionally) add clif_ prefix
base_no_ext <- tools::file_path_sans_ext(tolower(bn))
# If a file is "patient.csv", treat it as "clif_patient"
base_norm <- ifelse(looks_clif, base_no_ext,
                    paste0("clif_", base_no_ext))

# Map normalized basenames to full paths
found_map <- stats::setNames(all_files, base_norm)

# ---- Required tables (from your table_flags or a default list) ----
if (!exists("table_flags")) {
  # fallback if table_flags isn't defined yet
  required_raw <- c("patient","hospitalization","vitals","labs",
                    "medication_admin_continuous","adt","respiratory_support",
                    "hospital_diagnosis")
} else {
  required_raw <- names(table_flags)[table_flags]
}

required_files <- paste0("clif_", tolower(required_raw))

# What do we have?
cat("Detected CLIF-like files (normalized names):\n")
print(sort(unique(names(found_map))))

# Compute missing
missing <- setdiff(required_files, names(found_map))

if (length(missing) > 0) {
  cat("\nSearch summary:\n")
  cat(" - Searched path: ", tables_path, "\n", sep = "")
  cat(" - Recursive: TRUE\n")
  cat(" - Extensions: ", paste(exts, collapse = ", "), "\n", sep = "")
  cat(" - Total files found: ", length(all_files), "\n", sep = "")
  
  # Help user see near-misses by ignoring the clif_ prefix
  have_core <- sub("^clif_", "", names(found_map))
  need_core <- sub("^clif_", "", required_files)
  maybe_present <- intersect(need_core, have_core)
  
  if (length(maybe_present)) {
    cat("\nFiles that exist but may be missing the 'clif_' prefix or expected case:\n")
    print(maybe_present)
  }
  
  stop("Missing required tables: ", paste(missing, collapse = ", "))
}

# If we made it here, collect the paths in a named list for easy reading
clif_paths <- found_map[required_files]

# Optionally read them (example loaders by extension)
read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "parquet" = arrow::read_parquet(path),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported extension: ", ext))
}

# Example: load into a named list of tibbles/data.frames
clif_tables <- lapply(clif_paths, read_any)

# Access like:
# clif_tables[["clif_patient"]]
# clif_tables[["clif_hospitalization"]]
cat("\nLoaded required CLIF tables: ", paste(names(clif_tables), collapse = ", "), "\n", sep = "")
 
 # ---- Fast, low-copy datetime parser ----
 safe_posix <- function(x) {
   if (inherits(x, "POSIXct")) return(x)
   if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
   fasttime::fastPOSIXct(as.character(x), tz = "UTC")
 }
 
 safe_date <- function(x) {
   if (inherits(x, "Date")) return(x)
   suppressWarnings(as.Date(x))
 }
 
 # ---- Load ONLY the minimal tables/columns we need ----
 get_min <- function(tbl_name, cols) {
   nm <- paste0("clif_", tbl_name)
   if (!exists("clif_tables") || is.null(clif_tables[[nm]])) {
     stop("Missing table in clif_tables: ", nm)
   }
   out <- clif_tables[[nm]]
   # standardize names
   out <- out %>% rename_with(tolower)
   # keep only needed columns that exist
   cols_keep <- intersect(tolower(cols), names(out))
   out %>% select(any_of(cols_keep))
 }
 
 patient <- get_min("patient",
                    c("patient_id","birth_date","sex_category","race_category","ethnicity_category","preferred_language")
 ) %>% mutate(
   birth_date = safe_date(birth_date)
 )
 
 hospitalization <- get_min("hospitalization",
                            c("patient_id","hospitalization_id","admission_dttm","discharge_dttm","age_at_admission",
                              "zipcode_nine_digit", "zipcode_five_digit", "census_tract", "county_code")
 ) %>% mutate(
   admission_dttm = safe_posix(admission_dttm),
   discharge_dttm = safe_posix(discharge_dttm)
 )
 
 adt <- get_min("adt",
                c("hospitalization_id","in_dttm","out_dttm","location_category","location_type")
 ) %>% mutate(
   in_dttm  = safe_posix(in_dttm),
   out_dttm = safe_posix(out_dttm)
 )
 
 # For ARF logic + room air + P/F we need:
 resp_support <- get_min("respiratory_support",
                         c("hospitalization_id","recorded_dttm","device_category","mode_category","fio2_set")
 ) %>% mutate(
   recorded_dttm = safe_posix(recorded_dttm),
   fio2_set = suppressWarnings(as.numeric(fio2_set))
 )
 
 vitals <- get_min("vitals",
                   c("hospitalization_id","recorded_dttm","vital_category","vital_value")
 ) %>% mutate(
   recorded_dttm = safe_posix(recorded_dttm),
   vital_value = suppressWarnings(as.numeric(vital_value))
 )
 
 labs <- get_min("labs",
                 c("hospitalization_id","lab_result_dttm","lab_category","lab_value_numeric")
 ) %>% mutate(
   lab_result_dttm  = safe_posix(lab_result_dttm),
   lab_value_numeric = suppressWarnings(as.numeric(lab_value_numeric))
 )
 
 # ----------------------------
 # Parameters (tune as desired)
 # ----------------------------
 START_DATE <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
 END_DATE   <- as.POSIXct("2024-12-31 23:59:59", tz = "UTC")
 WINDOW_H   <- 24               # ± hours around ICU admit
 N_SPO2_MIN <- 6                # proxy for "continuous" pulse ox
 ROOM_AIR_FIO2 <- 0.21
 JOIN_NEAR_H  <- 1              # max time gap to pair with FiO2
 HYPERPAIR_H  <- 2              # pCO2–pH pairing window
 
 # ----------------------------
 # Identify ICU stays per hospitalization
 # ----------------------------
 icu_segments <- adt %>%
   mutate(
     is_icu = str_detect(tolower(coalesce(location_category, "")), "icu")
   ) %>%
   filter(is_icu)
 
 icu_bounds <- icu_segments %>%
   group_by(hospitalization_id) %>%
   summarize(
     first_icu_in = suppressWarnings(min(in_dttm, na.rm = TRUE)),
     last_icu_out = suppressWarnings(max(out_dttm, na.rm = TRUE)),
     .groups = "drop"
   ) %>%
   mutate(
     first_icu_in = ifelse(is.infinite(first_icu_in), NA, first_icu_in),
     last_icu_out = ifelse(is.infinite(last_icu_out), NA, last_icu_out),
     icu_los_hours = as.numeric(difftime(last_icu_out, first_icu_in, units = "hours"))
   )
 
 # ----------------------------
 # Build base set of candidate hospitalizations
 # ----------------------------
 base <- hospitalization %>%
   inner_join(icu_bounds, by = "hospitalization_id") %>%
   # ICU entry within 2018–2024
   filter(!is.na(first_icu_in),
          first_icu_in >= START_DATE,
          first_icu_in <= END_DATE) %>%
   # Adults: prefer age_at_admission, fallback to birth_date
   left_join(patient %>% select(patient_id, birth_date, sex_category, race_category, ethnicity_category),
             by = "patient_id") %>%
   mutate(
     age_years = coalesce(
       suppressWarnings(as.numeric(age_at_admission)),
       ifelse(!is.na(birth_date),
              as.numeric(floor((as.Date(admission_dttm) - birth_date)/365.25)), NA_real_)
     )
   )
 
 # ----------------------------
 # Filter to the ±24h analysis window for signals
 # ----------------------------
 base <- base %>%
   mutate(
     first_icu_in  = as.POSIXct(first_icu_in,  origin = "1970-01-01", tz = "UTC"),
     last_icu_out  = as.POSIXct(last_icu_out,  origin = "1970-01-01", tz = "UTC"),
     admission_dttm  = as.POSIXct(admission_dttm,  origin = "1970-01-01", tz = "UTC"),
     discharge_dttm  = as.POSIXct(discharge_dttm,  origin = "1970-01-01", tz = "UTC")
   )
 
 stopifnot(length(WINDOW_H) == 1L)  # guardrail
 
 # 2) use durations (dhours) instead of periods (hours)
 win <- base %>%
   transmute(
     hospitalization_id,
     win_start = first_icu_in - dhours(WINDOW_H),
     win_end   = first_icu_in + dhours(WINDOW_H)
   )
 
 # Keep only records in that window and for candidate hospitalizations
 vitals_win <- vitals %>%
   filter(vital_category == "spo2") %>%
   inner_join(win, by = "hospitalization_id") %>%
   filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
   select(hospitalization_id, recorded_dttm, spo2 = vital_value)
 
 labs_win <- labs %>%
   filter(lab_category %in% c("po2_arterial","pco2_arterial","ph_arterial")) %>%
   inner_join(win, by = "hospitalization_id") %>%
   filter(lab_result_dttm >= win_start, lab_result_dttm <= win_end) %>%
   select(hospitalization_id, lab_result_dttm, lab_category, val = lab_value_numeric)
 
 fio2_win <- resp_support %>%
   inner_join(win, by = "hospitalization_id") %>%
   filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
   select(hospitalization_id, recorded_dttm, fio2_set)
 
 # ----------------------------
 # Pair SpO2 & PaO2 to FiO2 to assess room air / P/F
 # ----------------------------
 H_NEAR <- as.numeric(JOIN_NEAR_H)     # hours
 H_HYP  <- as.numeric(HYPERPAIR_H)     # hours
 
 # --- Prep as data.table and key on (hospitalization_id, time) ---
 setDT(vitals_win);  setDT(fio2_win);  setDT(labs_win)
 
 setkey(vitals_win, hospitalization_id, recorded_dttm)
 setkey(fio2_win,   hospitalization_id, recorded_dttm)
 setkey(labs_win,   hospitalization_id, lab_result_dttm)
 
 # 1) SpO2 ↔ FiO2 (nearest within H_NEAR)
 # Ensure POSIXct
 vitals_win$recorded_dttm <- as.POSIXct(vitals_win$recorded_dttm, tz = "UTC")
 fio2_win$recorded_dttm   <- as.POSIXct(fio2_win$recorded_dttm,   tz = "UTC")
 
 # Make data.tables and give distinct time names
 # vdt: SpO2; fdt: FiO2  (ensure POSIXct already)
 vdt <- as.data.table(vitals_win)[
   , .(hospitalization_id, spo2_time = recorded_dttm, spo2 = spo2)
 ]
 fdt <- as.data.table(fio2_win)[
   , .(hospitalization_id, fio2_time = recorded_dttm, fio2_set)
 ]
 setkey(vdt, hospitalization_id, spo2_time)
 setkey(fdt, hospitalization_id, fio2_time)
 
 # KEEP a copy of the SpO2 time so it survives the join
 vdt[, spo2_time_keep := spo2_time]
 
 spo2_fio2 <- fdt[
   vdt, roll = "nearest", on = .(hospitalization_id, fio2_time = spo2_time), nomatch = 0L
 ][
   , timediff_h := abs(as.numeric(difftime(spo2_time_keep, fio2_time, units = "hours")))
 ][
   timediff_h <= as.numeric(JOIN_NEAR_H),
   .(hospitalization_id, spo2_time = spo2_time_keep, fio2_time, spo2, fio2_set,
     on_room_air = !is.na(fio2_set) & fio2_set <= ROOM_AIR_FIO2 + 1e-6,
     timediff_h)
 ]
 
 # 2) PaO2 ↔ FiO2 (nearest within H_NEAR) for P/F ratio
 po2dt <- as.data.table(labs_win)[lab_category == "po2_arterial",
                                  .(hospitalization_id, po2_time = lab_result_dttm, po2 = val)
 ]
 setkey(po2dt, hospitalization_id, po2_time)
 
 po2dt[, po2_time_keep := po2_time]
 
 po2_fio2 <- fdt[
   po2dt, roll = "nearest", on = .(hospitalization_id, fio2_time = po2_time), nomatch = 0L
 ][
   , timediff_h := abs(as.numeric(difftime(po2_time_keep, fio2_time, units = "hours")))
 ][
   timediff_h <= as.numeric(JOIN_NEAR_H),
   .(hospitalization_id, po2_time = po2_time_keep, fio2_time, po2, fio2_set,
     pf_ratio = fifelse(!is.na(fio2_set) & fio2_set > 0, po2 / fio2_set, as.numeric(NA)),
     timediff_h)
 ]
 
 # 3) pCO2 ↔ pH (nearest within H_HYP) for hypercapnia pair
 pco2dt <- as.data.table(labs_win)[lab_category == "pco2_arterial",
                                   .(hospitalization_id, pco2_time = lab_result_dttm, pco2 = val)
 ]
 phdt <- as.data.table(labs_win)[lab_category == "ph_arterial",
                                 .(hospitalization_id, ph_time = lab_result_dttm, ph = val)
 ]
 setkey(pco2dt, hospitalization_id, pco2_time)
 setkey(phdt,   hospitalization_id, ph_time)
 
 # Keep copies so both times are available post-join
 pco2dt[, pco2_time_keep := pco2_time]
 phdt[,   ph_time_keep   := ph_time]
 
 hyper_pairs <- phdt[
   pco2dt, roll = "nearest", on = .(hospitalization_id, ph_time = pco2_time), nomatch = 0L
 ][
   , timediff_h := abs(as.numeric(difftime(pco2_time_keep, ph_time, units = "hours")))
 ][
   timediff_h <= as.numeric(HYPERPAIR_H),
   .(hospitalization_id,
     pco2_time = pco2_time_keep, ph_time = ph_time_keep,
     pco2, ph,
     hyper_pair_hit = (pco2 >= 45 & ph < 7.35),
     timediff_h)
 ]
 
 
 # ----------------------------
 # Compute per-hospitalization ARF criteria flags
 # ----------------------------
 hypox_roomair_spo2 <- spo2_fio2 %>%
   mutate(hit = (spo2 < 90 & on_room_air)) %>%
   group_by(hospitalization_id) %>%
   summarize(any_spo2_roomair_hit = any(hit, na.rm = TRUE),
             spo2_n = n(),
             .groups = "drop")
 
 hypox_roomair_po2 <- po2_fio2 %>%
   mutate(hit = (po2 <= 60 & !is.na(fio2_set) & fio2_set <= ROOM_AIR_FIO2 + 1e-6)) %>%
   group_by(hospitalization_id) %>%
   summarize(any_po2_roomair_hit = any(hit, na.rm = TRUE),
             .groups = "drop")
 
 hypox_pf <- po2_fio2 %>%
   mutate(hit = (!is.na(pf_ratio) & pf_ratio <= 300)) %>%
   group_by(hospitalization_id) %>%
   summarize(any_pf_hit = any(hit, na.rm = TRUE), .groups = "drop")
 
 hyper_flags <- hyper_pairs %>%
   group_by(hospitalization_id) %>%
   summarize(any_hyper_pair = any(hyper_pair_hit, na.rm = TRUE), .groups = "drop")
 
 # Data availability flags (ABG or “continuous” SpO2 within window)
 abg_avail <- labs_win %>%
   filter(lab_category %in% c("po2_arterial","pco2_arterial","ph_arterial")) %>%
   distinct(hospitalization_id) %>%
   mutate(has_abg = TRUE)
 
 spo2_density <- vitals_win %>%
   group_by(hospitalization_id) %>%
   summarize(n_spo2 = n(), .groups = "drop") %>%
   mutate(has_cont_spo2 = n_spo2 >= N_SPO2_MIN)
 
 data_avail <- base %>%
   select(hospitalization_id) %>%
   left_join(abg_avail, by = "hospitalization_id") %>%
   left_join(spo2_density, by = "hospitalization_id") %>%
   mutate(
     has_abg = coalesce(has_abg, FALSE),
     has_cont_spo2 = coalesce(n_spo2 >= N_SPO2_MIN, FALSE),
     meets_data_rule = has_abg | has_cont_spo2
   )
 
 # ----------------------------
 # Join all flags and apply inclusion/exclusion
 # ----------------------------
 
 # ----  Bring geography into `base` safely  ----
 # Helper: keep only columns that actually exist
 keep_any <- function(df, cols) dplyr::select(df, dplyr::any_of(intersect(cols, names(df))))
 
 geo_cols <- c("zipcode_nine_digit", "zipcode_five_digit", "census_tract", "county_code" ,"census_block_group","latitude","longitude")
 
 flags <- base %>%
   # core demographics presence
   mutate(
     has_demo = !(is.na(age_years) | is.na(sex_category) | is.na(race_category)),
     adult = !is.na(age_years) & age_years >= 18,
     has_geo = !is.na(census_tract) | !is.na(zipcode_nine_digit) | !is.na(zipcode_five_digit) | !is.na(county_code),
     icu_24h = !is.na(icu_los_hours) & icu_los_hours >= 24
   ) %>%
   # keep first ICU stay per hospitalization by definition (already collapsed)
   left_join(hypox_roomair_spo2, by = "hospitalization_id") %>%
   left_join(hypox_roomair_po2, by = "hospitalization_id") %>%
   left_join(hypox_pf,          by = "hospitalization_id") %>%
   left_join(hyper_flags,       by = "hospitalization_id") %>%
   left_join(data_avail %>% select(hospitalization_id, meets_data_rule), by = "hospitalization_id") %>%
   mutate(
     any_hypox = coalesce(any_spo2_roomair_hit, FALSE) |
       coalesce(any_po2_roomair_hit, FALSE) |
       coalesce(any_pf_hit, FALSE),
     any_hypercap = coalesce(any_hyper_pair, FALSE),
     arf_criterion_met = any_hypox | any_hypercap,
     mixed_arf = any_hypox & any_hypercap
   )
 
 # Build inclusion mask
 incl <- flags %>%
   mutate(
     include =
       adult &
       icu_24h &
       has_demo &
       has_geo &
       meets_data_rule &
       arf_criterion_met
   )
 
 # ----------------------------
 # Final cohort and exclusion table
 # ----------------------------
 cohort_min <- incl %>%
   filter(include) %>%
   transmute(
     patient_id, hospitalization_id,
     admission_dttm, discharge_dttm,
     first_icu_in, last_icu_out,
     icu_los_hours,
     age_years, sex_category, race_category, ethnicity_category,
     census_tract, county_code, # <- edit as needed to reflect geo subunit at site
     # ARF subtype flags
     hypoxemic_arf = any_hypox,
     hypercapnic_arf = any_hypercap,
     mixed_arf,
     # provenance / data checks
     data_window_start = first_icu_in - hours(WINDOW_H),
     data_window_end   = first_icu_in + hours(WINDOW_H)
   )
 
 # Tabulate reasons for exclusion (first failing reason)
 exclusion_reasons <- incl %>%
   filter(!include) %>%
   mutate(reason = case_when(
     is.na(age_years) | is.na(sex_category) | is.na(race_category) ~ "Missing demographics",
     is.na(age_years) | age_years < 18 ~ "Under 18",
     is.na(first_icu_in) | is.na(last_icu_out) ~ "Missing ICU timing",
     icu_los_hours < 24 ~ "ICU stay < 24h",
     !has_geo ~ "Missing geo code",
     !meets_data_rule ~ "No ABG or continuous SpO2 in ±24h",
     !arf_criterion_met ~ "No ARF criteria in ±24h",
     TRUE ~ "Other"
   )) %>%
   select(patient_id, hospitalization_id, reason)
 
 excluded_tbl <- exclusion_reasons
 
 # ----------------------------
 # Quick counts
 # ----------------------------
 cat("Cohort selection summary:\n",
     "  Candidates (ICU 2018–2024): ", nrow(flags), "\n",
     "  Included:                   ", nrow(cohort_min), "\n",
     "  Excluded:                   ", nrow(excluded_tbl), "\n", sep = "")
 
 cat("\nSubtype counts among included:\n")
 print(cohort_min %>%
         summarize(
           hypoxemic_n = sum(hypoxemic_arf, na.rm = TRUE),
           hypercapnic_n = sum(hypercapnic_arf, na.rm = TRUE),
           mixed_n = sum(mixed_arf, na.rm = TRUE)
         ))
 
 # ----------------------------
 # Save output
 # ----------------------------

 # ---------- Output folder ----------
 ts <- format(Sys.time(), "%Y%m%d")
 out_dir <- file.path("output", glue("arf_cohort_{ts}"))
 if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
 
 # ---------- 1) Save INCLUDED COHORT ----------
 stopifnot(exists("cohort_min"))
 cohort_path <- file.path(out_dir, glue("cohort_inclusion_{ts}.csv"))
 write_csv(cohort_min, cohort_path)
 message("Saved included cohort: ", cohort_path)
 
 # ---------- 2a) Exclusion breakdown (table) ----------
 stopifnot(exists("excluded_tbl"), exists("flags"))
 
 # Counts by first failing reason (from excluded_tbl you created earlier)
 excl_breakdown <- excluded_tbl %>%
   count(reason, sort = TRUE) %>%
   mutate(
     total_candidates = nrow(flags),
     percent = round(100 * n / total_candidates, 1)
   )
 
 excl_tbl_path <- file.path(out_dir, glue("exclusion_breakdown_{ts}.csv"))
 write_csv(excl_breakdown, excl_tbl_path)
 message("Saved exclusion breakdown table: ", excl_tbl_path)
 
 # Also save the raw excluded rows for auditing
 excl_raw_path <- file.path(out_dir, glue("exclusions_raw_{ts}.csv"))
 write_csv(excluded_tbl, excl_raw_path)
 message("Saved raw exclusions: ", excl_raw_path)
 
 # ---------- 2b) Sequential flow (funnel-style counts) ----------
 # Compute stepwise remaining / excluded using your logical flags
 # (coalesce() guards against any stray NAs)
 cand_n <- nrow(flags)
 
 step1 <- flags %>% filter(coalesce(adult, FALSE))
 step2 <- step1 %>% filter(coalesce(has_demo, FALSE))
 step3 <- step2 %>% filter(coalesce(icu_24h, FALSE))
 step4 <- step3 %>% filter(coalesce(has_geo, FALSE))
 step5 <- step4 %>% filter(coalesce(meets_data_rule, FALSE))
 step6 <- step5 %>% filter(coalesce(arf_criterion_met, FALSE))  # final included
 
 flow_df <- tibble::tibble(
   step = c(
     "ICU 2018–2024 (candidates)",
     "≥18 years",
     "Demographics present",
     "ICU stay ≥24h",
     "Geography present (tract/bgrp/ZIP)",
     "ABG or continuous SpO₂ in ±24h",
     "Meets ARF criteria in ±24h"
   ),
   remaining = c(
     cand_n,
     nrow(step1),
     nrow(step2),
     nrow(step3),
     nrow(step4),
     nrow(step5),
     nrow(step6)
   )
 ) %>%
   mutate(excluded_at_step = dplyr::lag(remaining, default = remaining[1]) - remaining)
 
 flow_df <- flow_df %>% mutate(step_f = factor(step, levels = rev(step)))
 
 # Save the flow table
 flow_tbl_path <- file.path(out_dir, glue("selection_flow_counts_{ts}.csv"))
 write_csv(flow_df, flow_tbl_path)
 message("Saved sequential flow counts: ", flow_tbl_path)
 
 # ---------- 2c) Quick visuals ----------
 # Exclusion reasons bar (first failing reason)
 p_excl <- ggplot(excl_breakdown,
                  aes(x = n, y = reorder(reason, n))) +
   geom_col() +
   labs(
     title = "Exclusion breakdown (first failing reason)",
     x = "Excluded (n)", y = NULL,
     caption = glue("Total candidates: {cand_n} | Included: {nrow(step6)}")
   ) +
   theme_classic(base_size = 11)
 
 excl_png <- file.path(out_dir, glue("exclusion_breakdown_{ts}.png"))
 ggsave(excl_png, plot = p_excl, width = 8, height = 5, dpi = 300)
 message("Saved exclusion bar chart: ", excl_png)
 
 # Funnel-style remaining counts per step
 p_flow <- ggplot(flow_df, aes(x = remaining, y = step_f)) +
   geom_col() +
   geom_text(aes(label = comma(remaining)), hjust = -0.1, size = 3) +
   scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
   labs(
     title = "Cohort selection flow",
     x = "Remaining (n)", y = NULL,
     caption = glue("Final included: {nrow(step6)}")
   ) +
   theme_minimal(base_size = 11)
 
 flow_png <- file.path(out_dir, glue("selection_flow_{ts}.png"))
 ggsave(flow_png, plot = p_flow, width = 8, height = 5, dpi = 300)
 message("Saved funnel-style flow chart: ", flow_png)
 














