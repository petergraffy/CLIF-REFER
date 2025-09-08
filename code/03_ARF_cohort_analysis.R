# ================================================================================================
# ICU REspiratory Failure Environmental Risk (REFER) Index | PI: Peter Graffy (graffy@uchicago.edu)
# 
# Build cohorts
# # Read ARF & periop controls; attach their index hospitalization (admission/discharge).
# # Derive encounter year and patient ZIP at index.
#
# Build medical history (lookback window)
# # Pull prior diagnoses/procedures/meds/vitals from CLIF up to index admission start (e.g., 365 or 730 days).
# # Engineer comorbidity (Charlson/Elixhauser), prior ICU stays, vent history, baseline labs, meds (e.g., ACEi/ARB, diuretics), and utilization (ED/hospitalizations).
# 
# Define outcomes for index encounter
# # Primary: ICU LOS (days), hospital LOS, in-hospital mortality, 30-day mortality (if you have death_dttm), mechanical ventilation days, vasopressor use, AKI, RRT, readmission (7/30 days), discharge disposition.
# # Create binary and/or time/continuous versions as needed.
# 
# Link exposome
# # Join county-year exposome: SVI, PM2.5, NO2, Daymet summaries
# # Build an exposome index (z-score composite) and/or analyze components separately.
# 
# Model
# # LOS: neg-bin or log-normal regression; mortality/adverse events: logistic; vent days: neg-bin or zero-inflated.
# # Adjust for age/sex/race/ethnicity, comorbidity, admission dx/service/ICU type, seasonality/year, and hospital factors (if multi-site later).
# # Sensitivity: ARF vs periop controls, alternative lookback windows, component-wise exposome, DAG-guided covariate sets, IPTW/matching.
# =================================================================================================

## --- 0) Keep a minimal workspace as requested -------------------------------
keep_vars <- c("clif_tables", "cohort_min", "cohort_min_periop")
rm(list = setdiff(ls(envir = .GlobalEnv), keep_vars), envir = .GlobalEnv)

library(tidyverse)
library(lubridate)
library(janitor)
library(MASS)   # glm.nb
library(broom)
library(patchwork)
library(ggplot2)
library(ggeffects)

## --- 1) Helpers --------------------------------------------------------------
# Replace your helper with this:
get_tbl <- function(nm) {
  # find clif_tables up the environment chain (incl. .GlobalEnv)
  ct <- get0("clif_tables", inherits = TRUE)
  if (is.null(ct)) {
    stop("Couldn't find an object named 'clif_tables'. Make sure it exists in your workspace (exact name).")
  }
  # allow case-insensitive key match
  key <- nm
  if (!key %in% names(ct)) {
    ci_match <- names(ct)[tolower(names(ct)) == tolower(nm)]
    if (length(ci_match) == 1) key <- ci_match
  }
  if (!key %in% names(ct)) {
    stop(sprintf(
      "Table '%s' not found in clif_tables. Available: %s",
      nm, paste(names(ct), collapse = ", ")
    ))
  }
  janitor::clean_names(ct[[key]])
}

stopifnot(exists("clif_tables", inherits = TRUE))
print(names(clif_tables))  # confirm the exact table keys (e.g., "patient" vs "Patient")

# Core CLIF tables
patient         <- get_tbl("clif_patient")
hospitalization <- get_tbl("clif_hospitalization")
diagnosis       <- get_tbl("clif_hospital_diagnosis")
support         <- get_tbl("clif_respiratory_support")
med_admin       <- get_tbl("clif_medication_admin_continuous")     
icu_stay        <- get_tbl("clif_adt")
vitals          <- get_tbl("clif_vitals")
labs            <- get_tbl("clif_labs")

# ---------- helpers you need ----------
pick_col <- function(df, candidates, required = TRUE) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  if (required) stop(sprintf("None of these columns found: %s", paste(candidates, collapse = ", ")))
  rep(NA, nrow(df))
}

# helper stays the same
coalesce_any <- function(data, candidates) {
  cols <- dplyr::select(data, dplyr::any_of(candidates))
  if (ncol(cols) == 0) return(rep(NA_character_, nrow(data)))
  dplyr::coalesce(!!!cols)
}

safe_ts <- function(x, tz = "America/Chicago") {
  if (inherits(x, "POSIXt")) return(x)
  if (is.numeric(x)) {
    x2 <- ifelse(x > 1e12, x/1000, x)  # ms -> s
    return(lubridate::as_datetime(x2, tz = tz))
  }
  suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c("ymd_HMS","ymd_HM","ymd","ymdTz","ymdT",
               "mdy_HMS","mdy_HM","mdy",
               "dmy_HMS","dmy_HM","dmy","HMS"),
    tz = tz, quiet = TRUE
  ))
}

add_index_fields <- function(df) {
  # ensure admit/discharge exist (pull from hospitalization if missing)
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
  
  # compute raw values *outside* mutate (avoid using '.' with base pipe)
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


## 3) Index encounter fields for cohorts --------------------------------------

arf_idx    <- cohort_min        |> janitor::clean_names() |> add_index_fields() |> dplyr::mutate(cohort = "ARF")
periop_idx <- cohort_min_periop |> janitor::clean_names() |> add_index_fields() |> dplyr::mutate(cohort = "PERIOP")
cohort_all <- dplyr::bind_rows(arf_idx, periop_idx) |> dplyr::filter(!is.na(index_admit))

## 4) Medical history (unchanged logic; uses med_admin) -----------------------
lookback_days <- 365
cohort_lb <- cohort_all |>
  dplyr::transmute(patient_id, hospitalization_id, index_admit, index_year, cohort,
                   lb_start = index_admit - lubridate::days(lookback_days))

# ICU-only segments from ADT -----------------------------------------------
adt_tmp <- icu_stay |>
  dplyr::left_join(
    hospitalization |> dplyr::select(hospitalization_id, patient_id),
    by = "hospitalization_id"
  )

icu_segs <- adt_tmp |>
  dplyr::mutate(
    # pull raw strings/numbers, then parse with safe_ts()
    in_raw  = pick_col(adt_tmp, c("in_dttm","adt_in_dttm","icu_in","arrival_dttm","in_time")),
    out_raw = pick_col(adt_tmp, c("out_dttm","adt_out_dttm","icu_out","departure_dttm","out_time")),
    in_ts   = safe_ts(in_raw,  tz = "America/Chicago"),
    out_ts  = safe_ts(out_raw, tz = "America/Chicago"),
    loccat  = tolower(pick_col(adt_tmp, c("location_category","loc_category","unit_category")))
  ) |>
  # keep ICU-like locations only
  dplyr::filter(stringr::str_detect(loccat, "\\bicu\\b|\\bccu\\b|\\bmicu\\b|\\bsicu\\b|\\bcicu\\b|\\bnicu\\b|\\bpicu\\b")) |>
  # valid rows
  dplyr::filter(!is.na(patient_id), !is.na(in_ts), !is.na(out_ts)) |>
  # drop obvious bad intervals
  dplyr::filter(out_ts > in_ts)

# 2) count prior ICU stays in the lookback window
icu_hist <- icu_segs |>
  dplyr::inner_join(cohort_lb,
                    by = dplyr::join_by(patient_id),
                    relationship = "many-to-many") |>
  dplyr::filter(out_ts < index_admit, out_ts >= lb_start) |>
  dplyr::filter(hospitalization_id.x != hospitalization_id.y) |>
  dplyr::group_by(patient_id, hospitalization_id.y) |>
  dplyr::summarize(
    prior_icu_stays = dplyr::n_distinct(hospitalization_id.x),
    .groups = "drop"
  ) |>
  dplyr::rename(hospitalization_id = hospitalization_id.y)

# Baseline cardio-relevant meds in lookback (from clif_medication_admin_continuous)
# Attach patient_id to med_admin first, then build baseline meds
med_tmp <- med_admin |>
  dplyr::left_join(
    hospitalization |> dplyr::select(hospitalization_id, patient_id),
    by = "hospitalization_id"
  )

baseline_meds <- med_tmp |>
  dplyr::mutate(
    # parse start time robustly (handles POSIXct, character, or numeric epoch)
    admin_ts = safe_ts(pick_col(med_tmp, c("admin_dttm","start_dttm","start_time","start_ts"), required = TRUE)),
    # pick a med name/text field and normalize to lower
    med_low  = tolower(pick_col(med_tmp, c("medication_name","medication","drug_name","med_name"), required = TRUE))
  ) |>
  # join to cohort window by patient_id
  dplyr::inner_join(cohort_lb, by = "patient_id", relationship = "many-to-many") |>
  # lookback window: strictly prior to index admit
  dplyr::filter(!is.na(admin_ts), admin_ts >= lb_start, admin_ts < index_admit) |>
  dplyr::mutate(
    is_acei_arb = stringr::str_detect(med_low, "lisinopril|losartan|valsartan|enalapril|olmesartan|candesartan"),
    is_diuretic = stringr::str_detect(med_low, "furosemide|bumetanide|torsemide|hydrochlorothiazide"),
    is_bb       = stringr::str_detect(med_low, "metoprolol|carvedilol|atenolol|propranolol")
  ) |>
  dplyr::group_by(patient_id, hospitalization_id.y) |>
  dplyr::summarize(
    any_acei_arb = as.integer(any(is_acei_arb, na.rm = TRUE)),
    any_diuretic = as.integer(any(is_diuretic, na.rm = TRUE)),
    any_bb       = as.integer(any(is_bb, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::rename(hospitalization_id = hospitalization_id.y)

history_features <- cohort_all |>
  dplyr::select(patient_id, hospitalization_id) |>
  dplyr::distinct() |>
  dplyr::left_join(icu_hist,      by = c("patient_id","hospitalization_id")) |>
  dplyr::left_join(baseline_meds, by = c("patient_id","hospitalization_id")) |>
  dplyr::mutate(dplyr::across(c(prior_icu_stays, any_acei_arb, any_diuretic, any_bb), ~tidyr::replace_na(., 0)))

## 5) Outcomes (ICU LOS from clif_adt; vent from clif_respiratory_support) ----
# ICU LOS
icu_los <- icu_segs |>
  dplyr::semi_join(cohort_all, by = "hospitalization_id") |>
  dplyr::mutate(seg_days = as.numeric(difftime(out_ts, in_ts, units = "days"))) |>
  dplyr::group_by(hospitalization_id) |>
  dplyr::summarize(icu_los_days = sum(seg_days, na.rm = TRUE), .groups = "drop")

# Hospital LOS from hospitalization
# 1) LOS straight from cohort_all (already parsed)
hosp_los_core <- cohort_all %>%
  dplyr::transmute(
    hospitalization_id,
    hosp_los_days = as.numeric(difftime(index_discharge, index_admit, units = "days"))
  )

# 2) Optional meta from hospitalization (no timestamp parsing)
hosp_meta <- hospitalization %>%
  dplyr::select(hospitalization_id, discharge_category, county_code)

# 3) Final LOS + meta
hosp_los <- hosp_los_core %>%
  dplyr::left_join(hosp_meta, by = "hospitalization_id")

# In-hospital & 30-day mortality (if death_dttm exists in clif_patient)
mortality_instay <- cohort_all |>
  dplyr::left_join(patient |> dplyr::select(patient_id, death_dttm), by = "patient_id") |>
  dplyr::mutate(
    death_ts      = safe_ts(death_dttm),
    in_hosp_death = as.integer(!is.na(death_ts) & death_ts >= index_admit & death_ts <= index_discharge),
    death_30d     = as.integer(!is.na(death_ts) & death_ts <= (index_admit + lubridate::days(30)))
  ) |>
  dplyr::select(hospitalization_id, in_hosp_death, death_30d)

# Ventilation flag from clif_respiratory_support
# (looks for invasive ventilation or similar in mode/type fields)
vent_flag <- support %>%
  dplyr::mutate(
    dev_low = tolower(coalesce(device_name, ""))
  ) %>%
  # include anything that has "vent" but exclude the noninvasive terms
  dplyr::filter(
    stringr::str_detect(dev_low, "\\bvent\\b") | stringr::str_detect(dev_low, "vent;") |
      stringr::str_detect(dev_low, "vent ")
  ) %>%
  dplyr::filter(
    !stringr::str_detect(dev_low, "bipap|cpap|high flow|nasal cannula|nc|trach collar|oxytrach|room air")
  ) %>%
  dplyr::semi_join(cohort_all, by = "hospitalization_id") %>%
  dplyr::distinct(hospitalization_id) %>%
  dplyr::mutate(vent_proc_flag = 1L)

# --- 1) add patient_id to support and parse recorded_time --------------------
support_tmp <- support %>%
  dplyr::left_join(
    hospitalization %>% dplyr::select(hospitalization_id, patient_id),
    by = "hospitalization_id"
  ) %>%
  dplyr::mutate(
    rec_time = safe_ts(dplyr::coalesce(
      # put your actual time column(s) here in priority order
      .data$recorded_time, .data$recorded_dttm
    )),
    dev_low = tolower(dplyr::coalesce(.data$device_name, ""))  # normalize device text
  ) %>%
  dplyr::filter(!is.na(rec_time)) %>%               # keep rows we can time-order
  dplyr::semi_join(cohort_all, by = "hospitalization_id")  # only index cohort stays

# --- 2) classify invasive vent vs NIV ---------------------------------------
# invasive vent: strings containing 'vent' but excluding common noninvasive terms
support_class <- support_tmp %>%
  dplyr::mutate(
    is_niv = stringr::str_detect(
      dev_low, "bipap|cpap|high flow|hf vent|nasal cannula|nc|venturi|face mask|face tent|trach collar|oxytrach|room air|t-piece|ram cannula|aerosol mask|o2 hood"
    ),
    has_vent_token = stringr::str_detect(dev_low, "(^|[ ;])vent([ ;]|$)"),
    is_invasive_vent = has_vent_token & !is_niv
  )

# --- 3) sum consecutive durations (hours) with a gap threshold ---------------
# We add the time between consecutive rows when BOTH are invasive vent and
# the gap is not huge (e.g., <= 6 hours). Tweak threshold as needed.
gap_hours <- 6

vent_durations <- support_class %>%
  dplyr::arrange(hospitalization_id, rec_time) %>%
  dplyr::group_by(hospitalization_id) %>%
  dplyr::mutate(
    next_time   = dplyr::lead(rec_time),
    next_invas  = dplyr::lead(is_invasive_vent),
    gap_hr      = as.numeric(difftime(next_time, rec_time, units = "hours")),
    add_hours   = dplyr::if_else(
      is_invasive_vent & next_invas & !is.na(gap_hr) & gap_hr > 0 & gap_hr <= gap_hours,
      gap_hr, 0
    ),
    next_niv    = dplyr::lead(is_niv),
    add_niv_hrs = dplyr::if_else(
      is_niv & next_niv & !is.na(gap_hr) & gap_hr > 0 & gap_hr <= gap_hours,
      gap_hr, 0
    )
  ) %>%
  dplyr::summarize(
    vent_hours = sum(add_hours, na.rm = TRUE),
    niv_hours  = sum(add_niv_hrs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    vent_proc_flag = as.integer(vent_hours > 0)
  )

# Vasopressors from med_admin_continuous
vaso_flag <- med_admin |>
  dplyr::mutate(med_low = tolower(pick_col(med_admin, c("medication_name","medication","drug_name","med_name")))) |>
  dplyr::filter(stringr::str_detect(med_low, "norepinephrine|epinephrine|phenylephrine|vasopressin|dopamine")) |>
  dplyr::semi_join(cohort_all, by="hospitalization_id") |>
  dplyr::distinct(hospitalization_id) |>
  dplyr::mutate(vaso_flag = 1L)

# AKI (rough, creatinine swing) from clif_labs
aki_flag <- labs |>
  dplyr::mutate(name_low = tolower(pick_col(labs, c("lab_name","test_name","component","loinc_name")))) |>
  dplyr::filter(stringr::str_detect(name_low, "creatinine")) |>
  dplyr::semi_join(cohort_all, by="hospitalization_id") |>
  dplyr::group_by(hospitalization_id) |>
  dplyr::summarize(
    aki_flag = as.integer((max(lab_value_numeric, na.rm=TRUE) - min(lab_value_numeric, na.rm=TRUE)) >= 0.3),
    .groups = "drop"
  )

outcomes <- cohort_all |>
  dplyr::select(patient_id, hospitalization_id, cohort, index_admit, index_discharge, index_year, county_code) |>
  dplyr::left_join(icu_los,  by="hospitalization_id") |>
  dplyr::left_join(hosp_los, by="hospitalization_id") |>
  dplyr::left_join(mortality_instay, by="hospitalization_id") |>
  dplyr::left_join(vent_flag, by="hospitalization_id") |>
  dplyr::left_join(vaso_flag, by="hospitalization_id") |>
  dplyr::left_join(aki_flag,  by="hospitalization_id") |>
  dplyr::mutate(dplyr::across(c(vent_proc_flag, vaso_flag, aki_flag, in_hosp_death, death_30d), ~tidyr::replace_na(., 0L)))

# --- 4) join into outcomes ---------------------------------------------------
outcomes <- outcomes %>%
  dplyr::left_join(vent_durations, by = "hospitalization_id") %>%
  dplyr::mutate(
    vent_hours     = dplyr::coalesce(vent_hours, 0),
    niv_hours      = dplyr::coalesce(niv_hours, 0)  # drop this line if you don't want NIV
  )

outcomes <- outcomes %>%
  # coalesce duplicated fields, then drop the extras
  dplyr::mutate(
    county_code    = dplyr::coalesce(county_code.x, county_code.y),
    vent_proc_flag = dplyr::coalesce(vent_proc_flag.x, vent_proc_flag.y)
  ) %>%
  dplyr::select(
    patient_id, hospitalization_id, cohort, index_admit, index_discharge,
    index_year, county_code, icu_los_days, hosp_los_days, discharge_category,
    in_hosp_death, death_30d, vent_proc_flag, vaso_flag, aki_flag,
    vent_hours, niv_hours
  )


# --- 5) Join in exposome features -------------------------------------------

outcomes <- outcomes %>%
  dplyr::mutate(
    fips_county = stringr::str_pad(as.character(county_code), width = 5, pad = "0")
  )

svi   <- readr::read_csv("SVI_county_year.csv")
pm25  <- readr::read_csv("pm25_county_year.csv")
no2   <- readr::read_csv("no2_county_year.csv")
daymet<- readr::read_csv("daymet_county_year_allvars.csv")

exposome <- svi %>%
  dplyr::left_join(pm25,  by = c("GEOID","year")) %>%
  dplyr::left_join(no2,   by = c("GEOID","year")) %>%
  dplyr::left_join(daymet,by = c("GEOID","year"))

outcomes$GEOID <- outcomes$county_code

outcomes_exp <- outcomes %>%
  dplyr::left_join(exposome, by = c("GEOID","index_year" = "year"))

arf_summary <- outcomes_exp %>%
  dplyr::filter(cohort == "ARF") %>%
  dplyr::summarize(
    n            = dplyr::n(),
    mort_inhosp  = mean(in_hosp_death, na.rm = TRUE),
    mort_30d     = mean(death_30d, na.rm = TRUE),
    mean_ICU_los = mean(icu_los_days, na.rm = TRUE),
    mean_PM25    = mean(pm25_mean, na.rm = TRUE),
    mean_NO2     = mean(no2_mean, na.rm = TRUE),
    mean_SVI     = mean(svi_overall, na.rm = TRUE)
  )

arf_strat <- outcomes_exp %>%
  filter(cohort == "ARF") %>%
  mutate(
    svi_tertile = ntile(svi_overall, 3),
    pm25_quint  = ntile(pm25_mean, 5)
  ) %>%
  group_by(svi_tertile) %>%
  summarise(
    n = n(),
    mort_inhosp = mean(in_hosp_death, na.rm = TRUE),
    mort_30d    = mean(death_30d, na.rm = TRUE),
    mean_ICU_los= mean(icu_los_days, na.rm = TRUE),
    .groups = "drop"
  )

####### OUTCOME MODELING

arf_exp <- outcomes_exp %>%
  dplyr::filter(cohort == "ARF")

###### adjusting for confounders

arf_exp <- arf_exp %>%
  dplyr::left_join(
    patient %>%
      dplyr::select(patient_id, race_name, ethnicity_category, sex_category, birth_date),
    by = "patient_id"
  )

arf_exp <- arf_exp %>%
  dplyr::mutate(
    # combine race + ethnicity
    race_ethnicity = dplyr::case_when(
      stringr::str_to_lower(ethnicity_category) %in% c("hispanic", "latino", "latinx") ~
        paste0("Hispanic ", race_name),
      TRUE ~ paste0("Non-Hispanic ", race_name)
    ),
    # calculate age at index admission
    age = as.numeric(difftime(index_admit, birth_date, units = "days")) / 365.25
  )

arf_exp <- arf_exp %>%
  dplyr::mutate(
    re_low = stringr::str_to_lower(race_ethnicity),
    is_nonhisp = stringr::str_detect(re_low, "\\bnon[- ]?hispanic\\b"),
    is_hisp    = stringr::str_detect(re_low, "\\bhispanic\\b") & !is_nonhisp,
    is_white   = stringr::str_detect(re_low, "\\bwhite\\b"),
    is_black   = stringr::str_detect(re_low, "black"),
    is_asian_any = stringr::str_detect(
      re_low, "asian|mideast|filipino|chinese|korean|vietnamese|pacific islander|samoan"
    ),
    race_ethnicity_simple = dplyr::case_when(
      is_white & is_hisp    ~ "Hispanic White",
      is_white & is_nonhisp ~ "Non-Hispanic White",
      is_black & is_hisp    ~ "Hispanic Black",
      is_black & is_nonhisp ~ "Non-Hispanic Black",
      is_asian_any          ~ "Asian",
      TRUE                  ~ "Other"
    )
  ) %>%
  dplyr::select(-re_low, -is_nonhisp, -is_hisp, -is_white, -is_black, -is_asian_any)

table(arf_exp$race_ethnicity_simple, useNA = "ifany")

# Set reference levels (optional)
arf_exp <- arf_exp %>%
  dplyr::mutate(
    sex_category = factor(sex_category),
    race_ethnicity_simple = factor(
      race_ethnicity_simple,
      levels = c("Non-Hispanic White","Hispanic White","Non-Hispanic Black","Hispanic Black","Asian","Other")
    )
  )


# logistic regreesion for mortality

fit_mort <- glm(
  in_hosp_death ~ pm25_mean + no2_mean,
  data = arf_exp,
  family = binomial()
)

summary(fit_mort)
exp(cbind(OR = coef(fit_mort), confint(fit_mort)))  # odds ratios

# Logistic regression for 30-day mortality
fit_mort30 <- glm(
  death_30d ~ pm25_mean + no2_mean,
  data = arf_exp,
  family = binomial()
)

exp(cbind(OR = coef(fit_mort30), confint(fit_mort30)))

# Length of stay (overdispersed → use negative binomial)
fit_los <- glm.nb(
  icu_los_days ~ pm25_mean + no2_mean,
  data = arf_exp
)

summary(fit_los)
exp(coef(fit_los))  # incidence rate ratios

# Vent hours (can do NB or log-transformed linear regression)
fit_vent <- glm.nb(
  vent_hours ~ pm25_mean + no2_mean,
  data = arf_exp
)
summary(fit_vent)

# adjusted in hospital mortality

fit_mort_adj <- glm(
  in_hosp_death ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp,
  family = binomial()
)

summary(fit_mort_adj)
exp(cbind(OR = coef(fit_mort_adj), confint(fit_mort_adj)))  # odds ratios

# Logistic regression for 30-day mortality
fit_mort30_adj <- glm(
  death_30d ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp,
  family = binomial()
)

summary(fit_mort30_adj)
exp(cbind(OR = coef(fit_mort30_adj), confint(fit_mort30_adj)))

# Length of stay (overdispersed → use negative binomial)
fit_los_adj <- glm.nb(
  icu_los_days ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp
)

summary(fit_los_adj)
exp(coef(fit_los_adj))  # incidence rate ratios

# Vent hours (can do NB or log-transformed linear regression)
fit_vent_adj <- glm.nb(
  vent_hours ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp
)
summary(fit_vent_adj)
exp(coef(fit_vent_adj))








#################

fit_mort_adj <- glm(
  in_hosp_death ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp,
  family = binomial()
)

summary(fit_mort_adj)
# Odds ratios with 95% CI
OR_mort <- broom::tidy(fit_mort_adj, exponentiate = TRUE, conf.int = TRUE)
OR_mort

#################

fit_mort_adj30 <- glm(
  death_30d ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp,
  family = binomial()
)

summary(fit_mort_adj30)
# Odds ratios with 95% CI
OR_mort <- broom::tidy(fit_mort_adj30, exponentiate = TRUE, conf.int = TRUE)
OR_mort

#### ICU LOS

fit_icu_nb <- glm.nb(
  icu_los_days ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp
)

IRR_icu <- broom::tidy(fit_icu_nb, exponentiate = TRUE, conf.int = TRUE)  # incidence rate ratios
IRR_icu

####################### VENT HOURS

fit_vent_nb <- glm.nb(
  vent_hours ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp
)

IRR_vent <- broom::tidy(fit_vent_nb, exponentiate = TRUE, conf.int = TRUE)
IRR_vent


######################### AKI AND VASO

fit_aki <- glm(
  aki_flag ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp, family = binomial()
)

fit_vaso <- glm(
  vaso_flag ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp, family = binomial()
)

OR_aki  <- broom::tidy(fit_aki,  exponentiate = TRUE, conf.int = TRUE)
OR_vaso <- broom::tidy(fit_vaso, exponentiate = TRUE, conf.int = TRUE)

OR_aki
OR_vaso


############## PLOTTING FOR NO2

# Make a prediction grid for NO2 across observed range
newdat <- arf_exp %>%
  dplyr::summarise(min_no2 = min(no2_mean, na.rm=TRUE),
                   max_no2 = max(no2_mean, na.rm=TRUE)) %>%
  tidyr::expand_grid(
    no2_mean = seq(min_no2, max_no2, length.out=100),
    pm25_mean = mean(arf_exp$pm25_mean, na.rm=TRUE),
    age = mean(arf_exp$age, na.rm=TRUE),
    sex_category = "Female", # choose reference
    race_ethnicity_simple = "Non-Hispanic White",
    svi_overall = mean(arf_exp$svi_overall, na.rm=TRUE)
  )

# Fit logistic for in-hosp death
fit_mort_no2 <- glm(
  in_hosp_death ~ pm25_mean + no2_mean + age + sex_category +
    race_ethnicity_simple + svi_overall,
  data = arf_exp,
  family = binomial()
)

newdat$pred <- predict(fit_mort_no2, newdata = newdat, type = "response")

ggplot(newdat, aes(x=no2_mean, y=pred)) +
  geom_line(size=1.2, color="firebrick") +
  labs(x="NO2 (ppb)", y="Predicted probability of in-hospital death",
       title="Higher NO2 associated with greater in-hospital mortality in ARF Patients") +
  theme_classic(base_size=14)


fit_mort30_no2 <- glm(
  death_30d ~ pm25_mean + no2_mean + age + sex_category +
    race_ethnicity_simple + svi_overall,
  data = arf_exp,
  family = binomial()
)

newdat$pred30 <- predict(fit_mort30_no2, newdata = newdat, type = "response")

ggplot(newdat, aes(x=no2_mean, y=pred30)) +
  geom_line(size=1.2, color="navy") +
  labs(x="NO2 (ppb)", y="Predicted probability of 30-day death",
       title="Higher NO2 associated with greater 30-day mortality in ARF") +
  theme_classic(base_size=14)

fit_vent_no2 <- glm.nb(
  vent_hours ~ pm25_mean + no2_mean + age + sex_category +
    race_ethnicity_simple + svi_overall,
  data = arf_exp
)

newdat$pred_vent <- predict(fit_vent_no2, newdata = newdat, type = "response")

ggplot(newdat, aes(x=no2_mean, y=pred_vent)) +
  geom_line(size=1.2, color="darkgreen") +
  labs(x="NO2 (ppb)", y="Predicted ventilation hours",
       title="Higher NO2 associated with greater invasive ventilation use") +
  theme_classic(base_size=14)

############# COMBO FIGURE

# Pick reference levels from the most common categories -----------------
ref_sex <- arf_exp %>% count(sex_category, sort = TRUE) %>% slice(1) %>% pull(sex_category)
ref_re  <- arf_exp %>% count(race_ethnicity_simple, sort = TRUE) %>% slice(1) %>% pull(race_ethnicity_simple)

# Helper to build a prediction grid across observed NO2 -----------------
make_grid <- function(df, n = 100) {
  rng <- df %>% summarize(no2_min = min(no2_mean, na.rm = TRUE),
                          no2_max = max(no2_mean, na.rm = TRUE))
  tibble(
    no2_mean = seq(rng$no2_min, rng$no2_max, length.out = n),
    pm25_mean = mean(df$pm25_mean, na.rm = TRUE),
    age = mean(df$age, na.rm = TRUE),
    sex_category = ref_sex,
    race_ethnicity_simple = ref_re,
    svi_overall = mean(df$svi_overall, na.rm = TRUE)
  ) %>% mutate(no2_10 = no2_mean / 10)  # scale to per 10-ppb
}

grid <- make_grid(arf_exp)

# MODELS ----------------------------------------------------------------
fit_mort_inhosp <- glm(
  in_hosp_death ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp, family = binomial()
)

fit_mort_30d <- glm(
  death_30d ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp, family = binomial()
)

fit_vent_nb <- glm.nb(
  vent_hours ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp
)

# PREDICTIONS + 95% CIs -------------------------------------------------
# Logistic: get link + se, then plogis back-transform
pred_ci_logistic <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>%
    mutate(
      link = pr$fit, se = pr$se.fit,
      link_lo = link - 1.96 * se,
      link_hi = link + 1.96 * se,
      pred   = plogis(link),
      lo     = plogis(link_lo),
      hi     = plogis(link_hi)
    )
}

# NegBin: link is log(mu); back-transform with exp()
pred_ci_negbin <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>%
    mutate(
      link = pr$fit, se = pr$se.fit,
      pred = exp(link),
      lo   = exp(link - 1.96 * se),
      hi   = exp(link + 1.96 * se)
    )
}

df_mort_inhosp <- pred_ci_logistic(fit_mort_inhosp, grid)
df_mort_30d    <- pred_ci_logistic(fit_mort_30d,    grid)
df_vent        <- pred_ci_negbin (fit_vent_nb,      grid)

# PLOTS -----------------------------------------------------------------
p1 <- ggplot(df_mort_inhosp, aes(x = no2_10, y = pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(size = 1.2) +
  labs(
    x = "NO2 (per 10 ppb)",
    y = "Predicted probability",
    title = "In-hospital mortality vs NO2",
    subtitle = paste("Adjusted for PM2.5, age, sex, \nrace/ethnicity, SVI — \nrefs:",
                     ref_sex, "/", ref_re)
  ) +
  theme_classic(base_size = 14)

p2 <- ggplot(df_mort_30d, aes(x = no2_10, y = pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(size = 1.2) +
  labs(
    x = "NO2 (per 10 ppb)",
    y = "Predicted probability",
    title = "30-day mortality vs NO2",
    subtitle = paste("Adjusted for PM2.5, age, sex, \nrace/ethnicity, SVI — \nrefs:",
                     ref_sex, "/", ref_re)
  ) +
  theme_classic(base_size = 14)

p3 <- ggplot(df_vent, aes(x = no2_10, y = pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(size = 1.2) +
  labs(
    x = "NO2 (per 10 ppb)",
    y = "Predicted mean ventilation hours",
    title = "Ventilation hours vs NO2",
    subtitle = paste("Negative binomial; adjusted for PM2.5, age, sex, \nrace/ethnicity, SVI — \nrefs:",
                     ref_sex, "/", ref_re)
  ) +
  theme_classic(base_size = 14)

# COMBINE & SAVE --------------------------------------------------------
combo <- (p1 | p2) / p3
combo

ggsave("no2_arfi_outcomes_combined.png", combo, width = 11, height = 9, dpi = 300)

############### Hypercapnic vs hypoxemic vs mixed arf differences

arf_exp <- arf_exp %>%
  left_join(
    cohort_min %>%
      dplyr::select(patient_id, hospitalization_id,
             hypoxemic_arf, hypercapnic_arf, mixed_arf),
    by = c("patient_id","hospitalization_id")
  )

arf_exp <- arf_exp %>%
  mutate(
    arf_subtype = case_when(
      mixed_arf == 1 ~ "Mixed",
      hypoxemic_arf == 1 ~ "Hypoxemic",
      hypercapnic_arf == 1 ~ "Hypercapnic",
      TRUE ~ "Other"   # catch-all, if some records are missing flags
    ),
    arf_subtype = factor(arf_subtype,
                         levels = c("Hypoxemic","Hypercapnic","Mixed","Other"))
  )

fit_mort_sub <- glm(
  in_hosp_death ~ pm25_mean + no2_mean * arf_subtype + age +
    sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp,
  family = binomial()
)

summary(fit_mort_sub)

# Extract odds ratios per subtype
broom::tidy(fit_mort_sub, exponentiate = TRUE, conf.int = TRUE)


fit_vent_sub <- glm.nb(
  vent_hours ~ pm25_mean + no2_mean * arf_subtype + age +
    sex_category + race_ethnicity_simple + svi_overall,
  data = arf_exp
)

broom::tidy(fit_vent_sub, exponentiate = TRUE, conf.int = TRUE)

# Predicted curves for NO2 × subtype
gge_mort <- ggpredict(fit_mort_sub, terms = c("no2_mean [all]", "arf_subtype"))
gge_vent <- ggpredict(fit_vent_sub, terms = c("no2_mean [all]", "arf_subtype"))

plot(gge_mort) + labs(title="In-hospital mortality by NO₂ across ARF subtypes")
plot(gge_vent) + labs(title="Ventilation hours by NO₂ across ARF subtypes")


arf_exp <- arf_exp %>%
  mutate(
    no2_10   = no2_mean / 10,
    pm25_mean = as.numeric(pm25_mean),
    no2_mean  = as.numeric(no2_mean),
    age       = as.numeric(age)
  )

# Make sure subtype and covariates are factors
arf_exp <- arf_exp %>%
  mutate(
    arf_subtype = droplevels(factor(arf_subtype,
                                    levels = c("Hypoxemic","Hypercapnic","Mixed","Other"))),
    sex_category = droplevels(factor(sex_category)),
    race_ethnicity_simple = droplevels(factor(race_ethnicity_simple))
  )

# Choose reference levels (or set explicitly if you prefer)
ref_sex <- arf_exp %>% count(sex_category, sort = TRUE) %>% slice(1) %>% pull(sex_category)
ref_re  <- arf_exp %>% count(race_ethnicity_simple, sort = TRUE) %>% slice(1) %>% pull(race_ethnicity_simple)

# -------------------- Models (NO2 × subtype) --------------------
fit_mort_inhosp <- glm(
  in_hosp_death ~ pm25_mean + no2_10 * arf_subtype + age + sex_category +
    race_ethnicity_simple + svi_overall,
  data = arf_exp, family = binomial()
)

fit_mort_30d <- glm(
  death_30d ~ pm25_mean + no2_10 * arf_subtype + age + sex_category +
    race_ethnicity_simple + svi_overall,
  data = arf_exp, family = binomial()
)

fit_vent_nb <- glm.nb(
  vent_hours ~ pm25_mean + no2_10 * arf_subtype + age + sex_category +
    race_ethnicity_simple + svi_overall,
  data = arf_exp
)

# -------------------- Prediction grid --------------------
make_grid <- function(df, n = 150) {
  rng <- df %>% summarize(lo = quantile(no2_10, 0.01, na.rm = TRUE),
                          hi = quantile(no2_10, 0.99, na.rm = TRUE))
  expand.grid(
    no2_10 = seq(rng$lo, rng$hi, length.out = n),
    arf_subtype = levels(df$arf_subtype)
  ) %>%
    as_tibble() %>%
    mutate(
      pm25_mean = mean(df$pm25_mean, na.rm = TRUE),
      age = mean(df$age, na.rm = TRUE),
      sex_category = ref_sex,
      race_ethnicity_simple = ref_re,
      svi_overall = mean(df$svi_overall, na.rm = TRUE)
    )
}
grid <- make_grid(arf_exp)

# -------------------- Predict with 95% CIs --------------------
pred_ci_logistic <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>%
    mutate(
      link = pr$fit, se = pr$se.fit,
      pred = plogis(link),
      lo   = plogis(link - 1.96 * se),
      hi   = plogis(link + 1.96 * se)
    )
}

pred_ci_negbin <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  newdata %>%
    mutate(
      link = pr$fit, se = pr$se.fit,
      pred = exp(link),
      lo   = exp(link - 1.96 * se),
      hi   = exp(link + 1.96 * se)
    )
}

df_mort_inhosp <- pred_ci_logistic(fit_mort_inhosp, grid)
df_mort_30d    <- pred_ci_logistic(fit_mort_30d,    grid)
df_vent        <- pred_ci_negbin (fit_vent_nb,      grid)

# -------------------- Theme & palette --------------------
theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

pal <- scale_color_brewer(palette = "Dark2")
fill_pal <- scale_fill_brewer(palette = "Dark2")

# -------------------- Plots --------------------
p1 <- ggplot(df_mort_inhosp,
             aes(x = no2_10, y = pred, color = arf_subtype, fill = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, color = NA) +
  geom_line(size = 1.2) +
  pal + fill_pal +
  labs(
    title = "In-hospital mortality vs NO2 (per 10 ppb) \nby ARF subtype",
    x = "NO2 (per 10 ppb increase)",
    y = "Predicted probability",
    color = "ARF subtype", fill = "ARF subtype",
    subtitle = paste("Adjusted for PM2.5, age, sex, race/ethnicity, SVI; \nrefs:",
                     ref_sex, "/", ref_re)
  ) +
  theme_pub

p2 <- ggplot(df_mort_30d,
             aes(x = no2_10, y = pred, color = arf_subtype, fill = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, color = NA) +
  geom_line(size = 1.2) +
  pal + fill_pal +
  labs(
    title = "30-day mortality vs NO2 (per 10 ppb) \nby ARF subtype",
    x = "NO2 (per 10 ppb increase)",
    y = "Predicted probability",
    color = "ARF subtype", fill = "ARF subtype",
    subtitle = paste("Adjusted for PM2.5, age, sex, race/ethnicity, SVI; \nrefs:",
                     ref_sex, "/", ref_re)
  ) +
  theme_pub

p3 <- ggplot(df_vent,
             aes(x = no2_10, y = pred, color = arf_subtype, fill = arf_subtype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, color = NA) +
  geom_line(size = 1.2) +
  pal + fill_pal +
  labs(
    title = "Invasive ventilation hours vs NO2 (per 10 ppb) \nby ARF subtype",
    x = "NO2 (per 10 ppb increase)",
    y = "Predicted mean hours",
    color = "ARF subtype", fill = "ARF subtype",
    subtitle = paste("Negative binomial; adjusted for PM2.5, age, sex, race/ethnicity, SVI; \nrefs:",
                     ref_sex, "/", ref_re)
  ) +
  theme_pub

# -------------------- Combine & export --------------------
combined <- (p1 | p2) / p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined

 ggsave("fig_no2_by_subtype_combined.png", combined, width = 12, height = 9, dpi = 300)
 ggsave("fig_no2_by_subtype_combined.pdf", combined, width = 12, height = 9)











