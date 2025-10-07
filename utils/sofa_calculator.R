# ===============================================================================================
# SOFA Score Calculator for CLIF-ARFVI
# Purpose: Calculate SOFA-97 scores for ICU patients
# Author: Kaveri Chhikara
# ===============================================================================================

library(tidyverse)
library(lubridate)

# Helper function to calculate PaO2 from SpO2 (for imputation when PaO2 is missing)
calc_pao2 <- function(s) {
  s <- s / 100
  a <- (11700) / ((1 / s) - 1)
  b <- sqrt((50^3) + (a^2))
  pao2 <- ((b + a)^(1/3)) - ((b - a)^(1/3))
  return(pao2)
}

# Define medications and their unit conversion information
med_unit_info <- list(
  norepinephrine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr","mcg/kg/hour", "mg/kg/hr", 
                         "mg/kg/hour", "mcg/min", "mg/hr","mg/hour")
  ),
  epinephrine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", 
                         "mg/hr", "mg/hour", "mg/kg/hour", "mcg/kg/hour")
  ),
  phenylephrine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mcg/kg/hour", "mg/kg/hr",
                         "mg/kg/hour","mcg/min", "mg/hr","mg/hour")
  ),
  vasopressin = list(
    required_unit = "units/min",
    acceptable_units = c("units/min", "units/hr","units/hour", "milliunits/min", 
                         "milliunits/hr", "milliunits/hour")
  ),
  dopamine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr","mcg/kg/hour", 
                         "mg/kg/hr", "mg/kg/hour","mcg/min", "mg/hr", "mg/hour")
  ),
  angiotensin = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("ng/kg/min", "ng/kg/hr","ng/kg/hour")
  ),
  dobutamine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", 
                         "mg/hr", "mg/hour", "mcg/kg/hour", "mg/kg/hour")
  ),
  milrinone = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", 
                         "mg/hr", "mg/hour", "mcg/kg/hour", "mg/kg/hour")
  )
)

# Helper function to get medication dose conversion factor
get_conversion_factor <- function(med_category, med_dose_unit, weight_kg) {
  med_info <- med_unit_info[[med_category]]
  if (is.null(med_info)) return(NA_real_)

  med_dose_unit <- tolower(med_dose_unit)
  if (!(med_dose_unit %in% med_info$acceptable_units)) return(NA_real_)

  factor <- NA_real_

  if (med_category %in% c("norepinephrine", "epinephrine", "phenylephrine",
                          "dopamine", "milrinone", "dobutamine")) {
    if (med_dose_unit == "mcg/kg/min") {
      factor <- 1
    } else if (med_dose_unit %in% c("mcg/kg/hr", "mcg/kg/hour")) {
      factor <- 1 / 60
    } else if (med_dose_unit  %in% c("mg/kg/hr","mg/kg/hour") ) {
      factor <- 1000 / 60
    } else if (med_dose_unit == "mcg/min") {
      factor <- 1 / weight_kg
    } else if (med_dose_unit %in% c("mg/hr", "mg/hour")) {
      factor <- (1000 / 60) / weight_kg
    }
  } else if (med_category == "angiotensin") {
    if (med_dose_unit == "ng/kg/min") {
      factor <- 1 / 1000
    } else if (med_dose_unit %in% c("ng/kg/hr","ng/kg/hour")) {
      factor <- 1 / 1000 / 60
    }
  } else if (med_category == "vasopressin") {
    if (med_dose_unit == "units/min") {
      factor <- 1
    } else if (med_dose_unit %in% c("units/hr", "units/hour")) {
      factor <- 1 / 60
    } else if (med_dose_unit == "milliunits/min") {
      factor <- 1 / 1000
    } else if (med_dose_unit %in% c("milliunits/hr", "milliunits/hour")) {
      factor <- 1 / 1000 / 60
    }
  }

  return(factor)
}

# Main optimized SOFA calculation function
calculate_sofa <- function(cohort_data,
                          vitals_df,
                          labs_df,
                          support_df,
                          med_admin_df,
                          scores_df,
                          window_hours = 24,
                          safe_ts = NULL) {

  # Use default safe_ts if not provided
  if (is.null(safe_ts)) {
    safe_ts <- function(x) {
      if (inherits(x, "POSIXt")) return(x)
      lubridate::as_datetime(x, tz = "UTC")
    }
  }

  # Pre-filter to cohort IDs to reduce data volume
  cohort_ids <- unique(cohort_data$hospitalization_id)

  # Prepare time windows
  cohort_windows <- cohort_data %>%
    mutate(
      icu_admit_time = safe_ts(icu_admit_time),
      window_end = icu_admit_time + lubridate::hours(window_hours)
    ) %>%
    dplyr::select(hospitalization_id, icu_admit_time, window_end)

  # ---- 1. VITALS (MAP, SpO2, Weight) - Pre-filtered ----
  cat("Processing vitals...\n")

  # Get weights separately (usually recorded less frequently)
  weights <- vitals_df %>%
    filter(hospitalization_id %in% cohort_ids,
           tolower(vital_category) == "weight_kg") %>%
    mutate(weight = as.numeric(vital_value)) %>%
    filter(!is.na(weight), weight > 10, weight < 500) %>%  # reasonable weight bounds
    group_by(hospitalization_id) %>%
    summarise(weight_kg = median(weight, na.rm = TRUE), .groups = "drop")

  # Get MAP and SpO2 for the window
  vitals_window <- vitals_df %>%
    filter(hospitalization_id %in% cohort_ids) %>%
    mutate(
      vital_cat = tolower(vital_category),
      vital_val = as.numeric(vital_value),
      recorded_ts = safe_ts(recorded_dttm)
    ) %>%
    filter(vital_cat %in% c("map", "spo2"),
           !is.na(vital_val)) %>%
    # Apply reasonable bounds
    filter(
      (vital_cat == "map" & vital_val >= 20 & vital_val <= 250) |
      (vital_cat == "spo2" & vital_val >= 50 & vital_val <= 100)
    ) %>%
    inner_join(cohort_windows, by = "hospitalization_id") %>%
    filter(recorded_ts >= icu_admit_time & recorded_ts <= window_end) %>%
    group_by(hospitalization_id, vital_cat) %>%
    summarise(worst_val = min(vital_val, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = vital_cat,
      values_from = worst_val,
      values_fill = NA_real_
    ) %>%
    mutate(
      pao2_imputed = case_when(
        !is.na(spo2) & spo2 < 97 ~ calc_pao2(spo2),
        TRUE ~ NA_real_
      )
    )

  # ---- 2. RESPIRATORY SUPPORT - Pre-filtered ----
  cat("Processing respiratory support...\n")

  resp_window <- support_df %>%
    filter(hospitalization_id %in% cohort_ids) %>%
    mutate(
      recorded_ts = safe_ts(recorded_dttm),
      fio2_num = as.numeric(fio2_set),
      # Standardize FiO2 to decimal
      fio2_std = case_when(
        is.na(fio2_num) ~ NA_real_,
        fio2_num > 1 & fio2_num <= 100 ~ fio2_num / 100,
        fio2_num >= 0.21 & fio2_num <= 1 ~ fio2_num,
        TRUE ~ NA_real_
      ),
      device_cat = tolower(device_category)
    ) %>%
    inner_join(cohort_windows, by = "hospitalization_id") %>%
    filter(recorded_ts >= icu_admit_time & recorded_ts <= window_end) %>%
    group_by(hospitalization_id) %>%
    summarise(
      # Get max FiO2 (worst oxygenation need)
      fio2_max = suppressWarnings(max(fio2_std, na.rm = TRUE)),
      # Get highest level of support
      has_imv = any(grepl("imv|vent", device_cat), na.rm = TRUE),
      has_nippv = any(grepl("nippv", device_cat), na.rm = TRUE),
      has_cpap = any(grepl("cpap", device_cat), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      fio2_max = ifelse(is.finite(fio2_max), fio2_max, NA_real_),
      resp_support_max = case_when(
        has_imv ~ "Vent",
        has_nippv ~ "NIPPV",
        has_cpap ~ "CPAP",
        TRUE ~ "Other"
      )
    )

  # ---- 3. LABS - Pre-filtered and optimized ----
  cat("Processing labs...\n")

  labs_window <- labs_df %>%
    filter(hospitalization_id %in% cohort_ids) %>%
    mutate(
      lab_cat = tolower(lab_category),
      lab_val = as.numeric(lab_value_numeric),
      result_ts = safe_ts(lab_result_dttm)
    ) %>%
    filter(
      lab_cat %in% c("creatinine", "bilirubin_total", "po2_arterial", "platelet_count"),
      !is.na(lab_val)
    ) %>%
    # Apply reasonable bounds
    filter(
      (lab_cat == "creatinine" & lab_val >= 0 & lab_val <= 25) |
      (lab_cat == "bilirubin_total" & lab_val >= 0 & lab_val <= 100) |
      (lab_cat == "po2_arterial" & lab_val >= 20 & lab_val <= 800) |
      (lab_cat == "platelet_count" & lab_val >= 0 & lab_val <= 3000)
    ) %>%
    inner_join(cohort_windows, by = "hospitalization_id") %>%
    filter(result_ts >= icu_admit_time & result_ts <= window_end)

  # Process each lab type separately to avoid issues with duplicates
  labs_creatinine <- labs_window %>%
    filter(lab_cat == "creatinine") %>%
    group_by(hospitalization_id) %>%
    summarise(creatinine = max(lab_val, na.rm = TRUE), .groups = "drop") %>%
    mutate(creatinine = ifelse(is.finite(creatinine), creatinine, NA_real_))

  labs_bilirubin <- labs_window %>%
    filter(lab_cat == "bilirubin_total") %>%
    group_by(hospitalization_id) %>%
    summarise(bilirubin_total = max(lab_val, na.rm = TRUE), .groups = "drop") %>%
    mutate(bilirubin_total = ifelse(is.finite(bilirubin_total), bilirubin_total, NA_real_))

  labs_po2 <- labs_window %>%
    filter(lab_cat == "po2_arterial") %>%
    group_by(hospitalization_id) %>%
    summarise(po2_arterial = min(lab_val, na.rm = TRUE), .groups = "drop") %>%
    mutate(po2_arterial = ifelse(is.finite(po2_arterial), po2_arterial, NA_real_))

  labs_platelets <- labs_window %>%
    filter(lab_cat == "platelet_count") %>%
    group_by(hospitalization_id) %>%
    summarise(platelet_count = min(lab_val, na.rm = TRUE), .groups = "drop") %>%
    mutate(platelet_count = ifelse(is.finite(platelet_count), platelet_count, NA_real_))

  # Combine all labs
  labs_combined <- cohort_windows %>%
    dplyr::select(hospitalization_id) %>%
    left_join(labs_creatinine, by = "hospitalization_id") %>%
    left_join(labs_bilirubin, by = "hospitalization_id") %>%
    left_join(labs_po2, by = "hospitalization_id") %>%
    left_join(labs_platelets, by = "hospitalization_id")

  # ---- 4. GCS SCORES - Pre-filtered ----
  cat("Processing GCS scores...\n")

  gcs_window <- scores_df %>%
    filter(hospitalization_id %in% cohort_ids) %>%
    mutate(
      assessment_cat = tolower(assessment_category),
      gcs_val = as.numeric(numerical_value),
      recorded_ts = safe_ts(recorded_dttm)
    ) %>%
    filter(assessment_cat == "gcs_total",
           !is.na(gcs_val),
           gcs_val >= 3,
           gcs_val <= 15) %>%
    inner_join(cohort_windows, by = "hospitalization_id") %>%
    filter(recorded_ts >= icu_admit_time & recorded_ts <= window_end) %>%
    group_by(hospitalization_id) %>%
    summarise(
      min_gcs_score = min(gcs_val, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(min_gcs_score = ifelse(is.finite(min_gcs_score), min_gcs_score, NA_real_))

  # ---- 5. MEDICATIONS (Vasopressors) - Pre-filtered ----
  cat("Processing medications...\n")

  required_meds <- c("norepinephrine", "epinephrine", "phenylephrine",
                    "vasopressin", "dopamine", "angiotensin",
                    "dobutamine", "milrinone")

  meds_window <- med_admin_df %>%
    filter(hospitalization_id %in% cohort_ids) %>%
    mutate(
      med_cat = tolower(med_category),
      admin_ts = safe_ts(admin_dttm),
      dose_num = as.numeric(med_dose)
    ) %>%
    filter(med_cat %in% required_meds,
           !is.na(dose_num),
           dose_num > 0) %>%
    inner_join(cohort_windows, by = "hospitalization_id") %>%
    filter(admin_ts >= icu_admit_time & admin_ts <= window_end) %>%
    left_join(weights, by = "hospitalization_id") %>%
    mutate(
      # Use median weight of 80kg if missing
      weight_kg = coalesce(weight_kg, 80),
      conversion_factor = mapply(
        get_conversion_factor,
        med_cat,
        med_dose_unit,
        weight_kg,
        USE.NAMES = FALSE
      ),
      dose_converted = dose_num * conversion_factor
    ) %>%
    filter(!is.na(dose_converted)) %>%
    group_by(hospitalization_id, med_cat) %>%
    summarise(
      max_dose = max(dose_converted, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(max_dose = ifelse(is.finite(max_dose), max_dose, NA_real_)) %>%
    pivot_wider(
      names_from = med_cat,
      values_from = max_dose,
      values_fill = NA_real_
    )

  # ---- 6. COMBINE ALL COMPONENTS ----
  cat("Combining all components...\n")

  sofa_data <- cohort_windows %>%
    dplyr::select(hospitalization_id) %>%
    left_join(vitals_window, by = "hospitalization_id") %>%
    left_join(resp_window, by = "hospitalization_id") %>%
    left_join(labs_combined, by = "hospitalization_id") %>%
    left_join(gcs_window, by = "hospitalization_id") %>%
    left_join(meds_window, by = "hospitalization_id")

  # Calculate P/F ratios (ensure numeric types)
  sofa_data <- sofa_data %>%
    mutate(
      # Ensure numeric for division
      po2_arterial = as.numeric(po2_arterial),
      fio2_max = as.numeric(fio2_max),
      pao2_imputed = as.numeric(pao2_imputed),

      # Calculate ratios only when both values exist and FiO2 > 0
      p_f = case_when(
        !is.na(po2_arterial) & !is.na(fio2_max) & fio2_max > 0 ~ po2_arterial / fio2_max,
        TRUE ~ NA_real_
      ),
      p_f_imputed = case_when(
        !is.na(pao2_imputed) & !is.na(fio2_max) & fio2_max > 0 ~ pao2_imputed / fio2_max,
        TRUE ~ NA_real_
      )
    )

  # ---- 7. CALCULATE SOFA SCORES ----
  cat("Calculating SOFA scores...\n")

  sofa_scores <- sofa_data %>%
    mutate(
      # Cardiovascular (MAP and vasopressors)
      sofa_cv = case_when(
        dopamine > 15 | epinephrine > 0.1 | norepinephrine > 0.1 ~ 4,
        dopamine > 5 | (epinephrine > 0 & epinephrine <= 0.1) |
          (norepinephrine > 0 & norepinephrine <= 0.1) ~ 3,
        (dopamine > 0 & dopamine <= 5) | dobutamine > 0 ~ 2,
        map < 70 ~ 1,
        TRUE ~ 0
      ),

      # Coagulation (Platelets)
      sofa_coag = case_when(
        platelet_count < 20 ~ 4,
        platelet_count < 50 ~ 3,
        platelet_count < 100 ~ 2,
        platelet_count < 150 ~ 1,
        TRUE ~ 0
      ),

      # Liver (Bilirubin)
      sofa_liver = case_when(
        bilirubin_total >= 12.0 ~ 4,
        bilirubin_total >= 6.0 ~ 3,
        bilirubin_total >= 2.0 ~ 2,
        bilirubin_total >= 1.2 ~ 1,
        TRUE ~ 0
      ),

      # Renal (Creatinine)
      sofa_renal = case_when(
        creatinine >= 5.0 ~ 4,
        creatinine >= 3.5 ~ 3,
        creatinine >= 2.0 ~ 2,
        creatinine >= 1.2 ~ 1,
        TRUE ~ 0
      ),

      # Respiratory (P/F ratio with ventilation status)
      sofa_resp = case_when(
        # Use actual PaO2/FiO2 if available
        !is.na(p_f) & p_f < 100 & resp_support_max %in% c("Vent", "NIPPV", "CPAP") ~ 4,
        !is.na(p_f) & p_f < 200 & resp_support_max %in% c("Vent", "NIPPV", "CPAP") ~ 3,
        !is.na(p_f) & p_f < 300 ~ 2,
        !is.na(p_f) & p_f < 400 ~ 1,
        # Use imputed if actual not available
        !is.na(p_f_imputed) & p_f_imputed < 100 & resp_support_max %in% c("Vent", "NIPPV", "CPAP") ~ 4,
        !is.na(p_f_imputed) & p_f_imputed < 200 & resp_support_max %in% c("Vent", "NIPPV", "CPAP") ~ 3,
        !is.na(p_f_imputed) & p_f_imputed < 300 ~ 2,
        !is.na(p_f_imputed) & p_f_imputed < 400 ~ 1,
        TRUE ~ 0
      ),

      # CNS (GCS)
      sofa_cns = case_when(
        min_gcs_score < 6 ~ 4,
        min_gcs_score >= 6 & min_gcs_score <= 9 ~ 3,
        min_gcs_score >= 10 & min_gcs_score <= 12 ~ 2,
        min_gcs_score >= 13 & min_gcs_score <= 14 ~ 1,
        TRUE ~ 0
      ),

      # Total SOFA score
      sofa_total = sofa_cv + sofa_coag + sofa_liver + sofa_renal + sofa_resp + sofa_cns
    )

  # Return scores with key variables
  sofa_scores %>%
    dplyr::select(
      hospitalization_id,
      # Component scores
      sofa_cv, sofa_coag, sofa_liver, sofa_renal, sofa_resp, sofa_cns,
      # Total score
      sofa_total,
      # Optional: Include key values for verification (comment out if not needed)
      # map, platelet_count, bilirubin_total, creatinine,
      # p_f, p_f_imputed, min_gcs_score
    )
}