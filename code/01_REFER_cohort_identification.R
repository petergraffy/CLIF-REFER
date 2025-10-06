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
#   cohort_inclusion_periop  : control cohort of perioperative ARF
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
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
})

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
repo <- config$repo
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
                    "hospital_diagnosis", "patient_assessments")
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
   out %>% dplyr::select(any_of(cols_keep))
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
   left_join(patient %>% dplyr::select(patient_id, birth_date, sex_category, race_category, ethnicity_category),
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
   dplyr::select(hospitalization_id, recorded_dttm, spo2 = vital_value)
 
 labs_win <- labs %>%
   filter(lab_category %in% c("po2_arterial","pco2_arterial","ph_arterial")) %>%
   inner_join(win, by = "hospitalization_id") %>%
   filter(lab_result_dttm >= win_start, lab_result_dttm <= win_end) %>%
   dplyr::select(hospitalization_id, lab_result_dttm, lab_category, val = lab_value_numeric)
 
 fio2_win <- resp_support %>%
   inner_join(win, by = "hospitalization_id") %>%
   filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
   dplyr::select(hospitalization_id, recorded_dttm, fio2_set)
 
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
   dplyr::select(hospitalization_id) %>%
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
   left_join(data_avail %>% dplyr::select(hospitalization_id, meets_data_rule), by = "hospitalization_id") %>%
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
   dplyr::select(patient_id, hospitalization_id, reason)
 
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
 
 # ---- Naming helpers ----
 sanitize_tag <- function(x) {
   x <- if (is.null(x)) "SITE" else as.character(x)
   x <- iconv(x, to = "ASCII//TRANSLIT")        # remove accents
   x <- gsub("[^A-Za-z0-9]+", "_", x)           # non-alnum -> underscore
   x <- gsub("^_+|_+$", "", x)                  # trim underscores
   if (!nzchar(x)) "SITE" else x
 }
 
 SITE_NAME  <- sanitize_tag(if (exists("config")) config$site_name else Sys.getenv("SITE_NAME", "SITE"))
 SYSTEM_TIME <- format(Sys.Date(), "%Y%m%d")  # e.g., 20250905
 
 # Construct "[RESULT]_[SITE]_[TIME][.ext]"
 make_name <- function(result_name, ext = NULL) {
   base <- paste(result_name, SITE_NAME, SYSTEM_TIME, sep = "_")
   if (is.null(ext)) base else paste0(base, ".", ext)
 }

 # ---------- Output folder ----------
 ts <- format(Sys.time(), "%Y%m%d")
 out_dir <- file.path(repo, "output", paste0("run_", SITE_NAME, "_", SYSTEM_TIME))
 if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
 
 #  ---------- 1) INCLUDED COHORT ----------
 stopifnot(exists("cohort_min"))
 #write_csv(cohort_min, file.path(out_dir, make_name("cohort_inclusion", "csv")))
 
 message("Saved ARF cohort CSV: ",
         file.path(out_dir, make_name("cohort_inclusion", "csv")))
 
 # ---------- 2) EXCLUSION BREAKDOWNS ----------
 stopifnot(exists("excluded_tbl"), exists("flags"))
 
 excl_breakdown <- excluded_tbl %>%
   count(reason, sort = TRUE) %>%
   mutate(
     total_candidates = nrow(flags),
     percent = round(100 * n / total_candidates, 1)
   )
 
 write_csv(excl_breakdown, file.path(out_dir, make_name("exclusion_breakdown", "csv")))
 #write_csv(excluded_tbl,   file.path(out_dir, make_name("exclusions_raw", "csv")))
 
 message("Saved exclusion breakdown CSV: ",
         file.path(out_dir, make_name("exclusion_breakdown", "csv")))
 
 message("Saved raw exclusions CSV: ",
         file.path(out_dir, make_name("exclusions_raw", "csv")))
 
 # ---------- 3) SELECTION FLOW (counts table + PNG) ----------
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
     "ABG or continuous SpO₂ in 24h",
     "Meets ARF criteria in 24h"
   ),
   remaining = c(
     cand_n, nrow(step1), nrow(step2), nrow(step3), nrow(step4), nrow(step5), nrow(step6)
   )
 ) %>%
   mutate(excluded_at_step = dplyr::lag(remaining, default = remaining[1]) - remaining,
          step_f = factor(step, levels = step))
 
 # Save counts CSV
 write_csv(flow_df, file.path(out_dir, make_name("selection_flow_counts", "csv")))
 
 # Save bar-style flow PNG (kept for quick QA)
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
 
 ggsave(file.path(out_dir, make_name("selection_flow", "png")),
        plot = p_flow, width = 8, height = 5, dpi = 300)
 
 message("Saved selection flowchart PNG: ",
         file.path(out_dir, make_name("selection_flow", "png")))
 
 # ---------- 4) CONSORT-style flowchart ----------

 # Ensure excluded_at_step exists
 if (!"excluded_at_step" %in% names(flow_df)) {
   flow_df <- flow_df |>
     dplyr::mutate(excluded_at_step = dplyr::lag(remaining, default = remaining[1]) - remaining)
 }
 
 get_excl_lab_vec <- function(steps) {
   vapply(steps, function(s) {
     s <- as.character(s)
     if (grepl("^ICU .*candidates", s)) return("")                           # no side branch at step 1
     if (grepl("18", s))                           return("Under 18")
     if (grepl("Demograph", s, ignore.case=TRUE))  return("Missing demographics")
     if (grepl("ICU.*24", s, ignore.case=TRUE))    return("ICU stay < 24h")
     if (grepl("Geograph|ZIP|tract|block", s, ignore.case=TRUE))
       return("Missing geography")
     if (grepl("ABG|SpO", s, ignore.case=TRUE))    return("No ABG or continuous SpO2")
     if (grepl("ARF", s, ignore.case=TRUE))        return("No ARF criteria in +/-24h")
     "Excluded"
   }, character(1))
 }
 
 excl_labels_vec <- get_excl_lab_vec(flow_df$step)
 
 
 # Exclusion labels for each step (edit if you renamed steps)
 excl_labels <- setNames(
   c(
     "",                                # step 1 (candidates)
     "Under 18",                        # step 2
     "Missing demographics",            # step 3
     "ICU stay < 24h",                  # step 4
     "Missing geography",               # step 5
     "No ABG or continuous SpO₂",       # step 6  
     "No ARF criteria in ±24h"          # step 7
   ),
   flow_df$step
 )
 
 # Escape helper for DOT labels
 escape_label <- function(s) {
   s <- as.character(s)
   
   # Normalize common Unicode we use
   s <- gsub("\u00B1", "+/-", s, perl = TRUE)        # ±  -> +/-
   s <- gsub("\u2082", "2",    s, perl = TRUE)       # ₂  -> 2
   s <- gsub("\u2212", "-",    s, perl = TRUE)       # −  -> -
   s <- gsub("\u2013|\u2014", "-", s, perl = TRUE)   # –— -> -
   s <- gsub("[\u200B-\u200D\uFEFF]", "", s, perl = TRUE)  # zero-width chars
   
   # Transliterate but NEVER return NA
   s2 <- suppressWarnings(iconv(s, to = "ASCII//TRANSLIT", sub = ""))
   s2[is.na(s2)] <- s[is.na(s2)]
   
   # Escape for DOT
   s2 <- gsub("\\\\", "\\\\\\\\", s2)   # backslashes
   s2 <- gsub("\"", "\\\\\"", s2)       # quotes
   s2
 }
 
 # Keep nodes
 keep_nodes <- vapply(seq_len(nrow(flow_df)), function(i) {
   lbl <- sprintf("%s\nn = %s",
                  escape_label(flow_df$step[i]),
                  scales::comma(flow_df$remaining[i]))
   sprintf('K%d [label="%s", fillcolor="#f6f7f9", color="#2b2f36"];', i, lbl)
 }, character(1))
 
 # Edges between keeps
 edges <- if (nrow(flow_df) > 1) sprintf("K%d -> K%d;", 1:(nrow(flow_df)-1), 2:nrow(flow_df)) else character(0)
 
 # Exclusion side branches
 drop_nodes <- character(0); drop_edges <- character(0)
 for (i in 2:nrow(flow_df)) {
   excl_n <- flow_df$excluded_at_step[i]
   if (is.na(excl_n) || excl_n <= 0) next
   lab <- excl_labels_vec[i]                 # robust, index-based
   if (!nzchar(lab)) next
   lbl <- sprintf("%s\nn = %s", escape_label(lab), scales::comma(excl_n))
   drop_nodes <- c(drop_nodes, sprintf('X%d [label="%s", fillcolor="#fdecea", color="#c0392b"];', i, lbl))
   drop_edges <- c(drop_edges, sprintf("K%d -> X%d;", i-1, i))
 }
 
 dot <- paste(
   'digraph G {',
   '  rankdir=TB;',
   '  fontname="Helvetica";',
   '  node [shape=box, style="rounded,filled", fontname="Helvetica"];',
   '  edge [fontname="Helvetica"];',
   paste(keep_nodes, collapse = "\n"),
   paste(drop_nodes, collapse = "\n"),
   paste(edges, collapse = "\n"),
   paste(drop_edges, collapse = "\n"),
   '}',
   sep = "\n"
 )
 
 dot_path <- file.path(out_dir, make_name("cohort_flow", "dot"))
 svg_path <- file.path(out_dir, make_name("cohort_flow", "svg"))
 png_path <- file.path(out_dir, make_name("cohort_flow", "png"))
 
 writeLines(dot, dot_path)
 gr2 <- DiagrammeR::grViz(dot)
 svg_txt <- DiagrammeRsvg::export_svg(gr2)
 writeLines(svg_txt, svg_path)
 rsvg::rsvg_png(svg_path, png_path, width = 1600, height = 1100)
 
 message("Saved selection flowchart DOT/SVG/PNG: ",
         file.path(out_dir, make_name("cohort_flow", "png")))

 # ===============================
 # Perioperative ARF Control Cohort
 # ===============================
 
 # 1) Load minimal diagnosis table
 hospital_dx <- get_min("hospital_diagnosis",
                        c("hospitalization_id", "diagnosis_code")
 ) %>%
   mutate(diagnosis_code = toupper(trimws(diagnosis_code)))
 
 # 2) Flag perioperative control codes: J95.82-.84
 dx_periop <- hospital_dx %>%
  mutate(periop_ctrl = grepl("^J95\\.(82|83|84)$", diagnosis_code)) %>%
  group_by(hospitalization_id) %>%
  summarize(dx_periop_ctrl = any(periop_ctrl, na.rm = TRUE), .groups = "drop")

# 3) Join to your existing `flags` table and select controls
# (Assumes you already built `flags` in your inclusion script)
stopifnot(exists("flags"))

flags_periop <- flags %>%
  left_join(dx_periop, by = "hospitalization_id") %>%
  mutate(
    is_periop_ctrl = coalesce(dx_periop_ctrl, FALSE),
    # Same base filters as your main cohort
    include_periop =
      coalesce(adult, FALSE) &
      coalesce(icu_24h, FALSE) &
      coalesce(has_demo, FALSE) &
      coalesce(has_geo, FALSE) &
      is_periop_ctrl &
      !coalesce(arf_criterion_met, FALSE)  # explicitly exclude ARF
  )

# 4) Build a cohort with the SAME columns as `cohort_min`
# If WINDOW_H isn't in scope (should be), default to 24h
if (!exists("WINDOW_H")) WINDOW_H <- 24

cohort_min_periop <- flags_periop %>%
  filter(include_periop) %>%
  transmute(
    patient_id, hospitalization_id,
    admission_dttm, discharge_dttm,
    first_icu_in, last_icu_out,
    icu_los_hours,
    age_years, sex_category, race_category, ethnicity_category,
    census_tract, county_code,
    # ARF subtype flags mirror the main cohort schema; should all be FALSE by design
    hypoxemic_arf   = coalesce(any_hypox, FALSE) & FALSE,
    hypercapnic_arf = coalesce(any_hypercap, FALSE) & FALSE,
    mixed_arf       = FALSE,
    # window columns identical to main cohort
    data_window_start = first_icu_in - lubridate::dhours(WINDOW_H),
    data_window_end   = first_icu_in + lubridate::dhours(WINDOW_H),
    # add an explicit control marker (extra column; remove if you want exact parity)
    is_periop_control = TRUE
  )

# 5) Quick sanity print
cat("Perioperative control cohort size: ", nrow(cohort_min_periop), "\n", sep = "")

# 6) Save with your naming convention
# (Assumes make_name(), SITE_NAME, SYSTEM_TIME, out_dir already exist; if not, define them)
if (!exists("make_name")) {
  sanitize_tag <- function(x) {
    x <- if (is.null(x)) "SITE" else as.character(x)
    x <- iconv(x, to = "ASCII//TRANSLIT")
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) "SITE" else x
  }
  SITE_NAME   <- sanitize_tag(if (exists("config")) config$site_name else Sys.getenv("SITE_NAME", "SITE"))
  SYSTEM_TIME <- format(Sys.time(), "%Y%m%dT%H%M%S")
  make_name   <- function(result_name, ext = NULL) {
    base <- paste(result_name, SITE_NAME, SYSTEM_TIME, sep = "_")
    if (is.null(ext)) base else paste0(base, ".", ext)
  }
  out_dir <- file.path("outputs", paste0("run_", SITE_NAME, "_", SYSTEM_TIME))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
}

# readr::write_csv(
#   cohort_min_periop,
#   file.path(out_dir, make_name("cohort_inclusion_periop", "csv"))
# )
# message("Saved perioperative control cohort CSV: ",
#         file.path(out_dir, make_name("cohort_inclusion_periop", "csv")))

# =========================
# Wrap up
# =========================

keep_vars <- c("clif_tables", "cohort_min", "cohort_min_periop", "repo")
rm(list = setdiff(ls(envir = .GlobalEnv), keep_vars), envir = .GlobalEnv)

message("\n🎯 Cohort identification complete!")
message("📂 Next steps: Please immediately run `02_REFER_linkage_analysis.R`.\n")









