# ===============================================================================================
# ICU Respiratory Failure Environmental Risk (REFER) Index
# PI: Peter Graffy (graffy@uchicago.edu)
# Script: 01_build_refer.R
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
  library(yaml)
})

# ---- 0.1 Load config (YAML or fallback defaults) ----------------------------------------------
keep_vars <- c("clif_tables", "cohort_min", "cohort_min_periop")
rm(list = setdiff(ls(envir = .GlobalEnv), keep_vars), envir = .GlobalEnv)

#setwd("~/Desktop/Peter/Postdoc/CLIF-ARFVI") <------ set your working directory here if needed

# ------------------------------- 0) Config & Setup --------------------------------
# Expect: utils/config.R defines a list named `config`
# Example fields (customize in utils/config.R as needed):
#   site_name, tables_path, file_type,
#   project_code, version, output_dir, figures_dir

source("utils/config.R")

# safety checks
if (!exists("config") || !is.list(config)) {
  stop("`config` was not created by utils/config.R. Make sure it defines a list named `config`.")
}

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
  output_dir   = config$output_dir   %||% "output",
  figures_dir  = config$figures_dir  %||% "figures",
  site_name    = config$site_name    %||% "site",
  date_stamp   = format(Sys.Date(), "%Y%m%d"),
  run_id       = format(Sys.time(), "%Y%m%d_%H%M%S")
)

# prefix like: refer_v0_1_20250908
cfg$prefix <- paste(cfg$project_code, cfg$site_name, cfg$date_stamp, sep = "_")

# ensure output folders exist
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$output_dir, cfg$figures_dir), recursive = TRUE, showWarnings = FALSE)

# ---- Save helpers (use cfg/output_dir/figures_dir) ----------------------------
.path_csv <- function(name)
  file.path(cfg$output_dir, paste0(cfg$prefix, "_", name, "_", cfg$run_id, ".csv"))

.path_rds <- function(name)
  file.path(cfg$output_dir, paste0(cfg$prefix, "_", name, "_", cfg$run_id, ".rds"))

.path_png <- function(name)
  file.path(cfg$output_dir, cfg$figures_dir, paste0(cfg$prefix, "_", name, "_", cfg$run_id, ".png"))

.path_pdf <- function(name)
  file.path(cfg$output_dir, cfg$figures_dir, paste0(cfg$prefix, "_", name, "_", cfg$run_id, ".pdf"))

save_tbl   <- function(x, name) readr::write_csv(x, .path_csv(name))
save_model <- function(fit, name) saveRDS(fit, .path_rds(name))
save_plot  <- function(p, name, w = 11, h = 9, dpi = 300) {
  ggsave(.path_png(name), p, width = w, height = h, dpi = dpi)
  ggsave(.path_pdf(name), p, width = w, height = h)
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
  dplyr::select(patient_id, hospitalization_id, cohort, index_admit, index_discharge, index_year, county_code) |>
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

svi    <- readr::read_csv("SVI_county_year.csv")
pm25   <- readr::read_csv("pm25_county_year.csv")
no2    <- readr::read_csv("no2_county_year.csv")
daymet <- readr::read_csv("daymet_county_year_allvars.csv")

exposome <- svi |>
  left_join(pm25,   by = c("GEOID","year")) |>
  left_join(no2,    by = c("GEOID","year")) |>
  left_join(daymet, by = c("GEOID","year"))

outcomes_exp <- outcomes |>
  left_join(exposome, by = c("GEOID","index_year" = "year"))

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

save_tbl(arf_exp, "arf_exp")

# ------------------------------------ 7) Models (Adjusted) --------------------------------------
# Helper to tidy-save any model
tidy_and_save <- function(fit, name, exponentiate = FALSE) {
  tt <- broom::tidy(fit, exponentiate = exponentiate, conf.int = TRUE)
  save_tbl(tt, paste0("model_tidy_", name))
  tt
}

# Logistic: in-hospital death
fit_mort_adj <- glm(in_hosp_death ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
                    data = arf_exp, family = binomial())
tt_mort <- tidy_and_save(fit_mort_adj, "inhosp_death_adj", exponentiate = TRUE)

# Logistic: 30-day death
fit_mort30_adj <- glm(death_30d ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
                      data = arf_exp, family = binomial())
tt_mort30 <- tidy_and_save(fit_mort30_adj, "death30d_adj", exponentiate = TRUE)

# NegBin: ICU LOS
fit_icu_nb <- glm.nb(icu_los_days ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
                     data = arf_exp)
tt_icu <- tidy_and_save(fit_icu_nb, "icu_los_adj", exponentiate = TRUE)

# NegBin: Vent hours
fit_vent_nb <- glm.nb(vent_hours ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
                      data = arf_exp)
tt_vent <- tidy_and_save(fit_vent_nb, "vent_hours_adj", exponentiate = TRUE)

# Optional: AKI & vaso (logistic)
fit_aki  <- glm(aki_flag  ~ pm25_mean + no2_mean + age + sex_category + race_ethnicity_simple + svi_overall,
                data = arf_exp, family = binomial())

tidy_and_save(fit_aki,  "aki_adj",  exponentiate = TRUE)


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
    svi_overall = mean(df$svi_overall, na.rm = TRUE)
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

fit_mort_inhosp_sub <- glm(in_hosp_death ~ pm25_mean + no2_10 * arf_subtype + age + sex_category + race_ethnicity_simple + svi_overall,
                           data = arf_exp, family = binomial())
fit_mort_30d_sub    <- glm(death_30d     ~ pm25_mean + no2_10 * arf_subtype + age + sex_category + race_ethnicity_simple + svi_overall,
                           data = arf_exp, family = binomial())
fit_vent_nb_sub     <- glm.nb(vent_hours  ~ pm25_mean + no2_10 * arf_subtype + age + sex_category + race_ethnicity_simple + svi_overall,
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
      svi_overall = mean(df$svi_overall, na.rm = TRUE)
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
    race_ethnicity_simple + svi_overall,
  data = arf_exp, family = binomial()
)

fit_mort_30d_sub <- glm(
  death_30d ~ pm25_mean + no2_10 * arf_subtype + age + sex_category +
    race_ethnicity_simple + svi_overall,
  data = arf_exp, family = binomial()
)

fit_vent_nb_sub <- MASS::glm.nb(
  vent_hours ~ pm25_mean + no2_10 * arf_subtype + age + sex_category +
    race_ethnicity_simple + svi_overall,
  data = arf_exp
)

# --- Build grid from the model's levels and coerce factors before predict ---
arf_lvls <- fit_mort_inhosp_sub$xlevels$arf_subtype  # training levels used by model

make_grid_sub_safe <- function(df, arf_lvls, n = 150) {
  rng <- df %>%
    dplyr::summarize(lo = stats::quantile(no2_10, 0.01, na.rm = TRUE),
                     hi = stats::quantile(no2_10, 0.99, na.rm = TRUE))
  expand.grid(
    no2_10 = seq(rng$lo, rng$hi, length.out = n),
    arf_subtype = arf_lvls,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      pm25_mean = mean(df$pm25_mean, na.rm = TRUE),
      age = mean(df$age, na.rm = TRUE),
      sex_category = names(sort(table(df$sex_category), decreasing = TRUE))[1],
      race_ethnicity_simple = names(sort(table(df$race_ethnicity_simple), decreasing = TRUE))[1],
      svi_overall = mean(df$svi_overall, na.rm = TRUE),
      # coerce to the model's factor levels
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

message("All done. Outputs written to: ", normalizePath(cfg$output_dir))




















