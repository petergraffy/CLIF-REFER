# ============================================================
# MASTER SUPPLEMENT TABLE (Part 1)
#   - Pooled CIF-SHR results (REML + I²)  [from your CIF SHR meta-analysis]
#   - Pooled regression model meta-analysis (random DL) [from /models aggregation]
#   - 30-day pooled CIF values (risk-set weighted) at Day 30 by exposure bin
#
# OUTPUTS:
#   1) master_supplement_table_long.csv
#   2) master_supplement_table_wide.csv
#   3) master_supplement_table.docx  (optional, requires officer/flextable)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(stringr)
  library(glue)
  library(readr)
  library(fs)
})

# ------------------------------------------------------------
# 0) PATHS (edit if needed)
# ------------------------------------------------------------
base_dir   <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis"

# CIF-SHR pooled tables (from your earlier block writing these files)
shr_dir    <- base_dir
shr_file   <- file.path(shr_dir, "CIF_SHR_PM25_NO2_overall_REML_I2.csv")

# Models meta-analysis pooled table (from "Model aggregation" section)
models_dir <- file.path(base_dir, "models", "aggregated_outputs")
models_file <- file.path(models_dir, "pooled_meta_by_pollutant_outcome.csv")

# CIF pooled curves master file (from your CIF section)
# If you already wrote combined_cif_all_sites.csv:
cif_out_dir <- file.path(base_dir, "cif", "fig_out")
cif_master_file <- file.path(cif_out_dir, "combined_cif_all_sites.csv")

# Where to save the master table
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Standardize pollutant strings
norm_pollutant <- function(x) {
  x2 <- str_to_upper(as.character(x))
  case_when(
    str_detect(x2, "PM2\\.5|PM25|PM_?2\\.?5") ~ "PM2.5",
    str_detect(x2, "NO2|NO_?2")              ~ "NO2",
    TRUE ~ as.character(x)
  )
}

# Standardize outcome labels across sections
norm_outcome <- function(x) {
  xl <- str_to_lower(as.character(x))
  case_when(
    str_detect(xl, "death") & str_detect(xl, "30")                 ~ "30-day death",
    str_detect(xl, "inhosp") | str_detect(xl, "in[-_ ]?hosp")      ~ "In-hospital death",
    str_detect(xl, "icu") & str_detect(xl, "los")                  ~ "ICU length of stay",
    str_detect(xl, "vent")                                         ~ "Ventilation hours",
    str_detect(xl, "successful") & str_detect(xl, "extub")         ~ "Successful extubation",
    str_detect(xl, "persistent") | str_detect(xl, "\\bprf\\b")     ~ "Persistent respiratory failure",
    str_detect(xl, "^death$")                                      ~ "Death",
    TRUE ~ as.character(x)
  )
}

# Extract pooled SHR and I2 if you have string columns
fmt_ci <- function(est, lo, hi, digits = 2) {
  ifelse(is.na(est) | is.na(lo) | is.na(hi), NA_character_,
         sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"), est, lo, hi))
}

# Find the nearest (or exact) day=30 row robustly
pick_day_value <- function(df, day_target = 30) {
  # expects grouped df with 'day' and a metric col called 'cif'
  if (!nrow(df)) return(tibble(day = NA_real_, cif = NA_real_))
  if (any(df$day == day_target, na.rm = TRUE)) {
    df %>% filter(day == day_target) %>% slice(1) %>% select(day, cif)
  } else {
    # nearest day
    df %>% mutate(d = abs(day - day_target)) %>% arrange(d) %>% slice(1) %>% select(day, cif)
  }
}

# ------------------------------------------------------------
# 2) READ: CIF-SHR pooled (REML + I²)
# ------------------------------------------------------------
stopifnot(file.exists(shr_file))

shr_pooled <- read_csv(shr_file, show_col_types = FALSE) %>%
  clean_names()

# Expect columns like: outcome, exposure, sites, pooled_shr_95_ci, i_2
# We’ll be permissive.
shr_pooled2 <- shr_pooled %>%
  transmute(
    section   = "CIF SHR meta-analysis",
    pollutant = norm_pollutant(Exposure %||% exposure),
    outcome   = norm_outcome(Outcome %||% outcome),
    method    = "Random-effects (REML)",
    k         = as.numeric(Sites %||% sites),
    estimate  = NA_real_,
    lcl       = NA_real_,
    ucl       = NA_real_,
    effect_ci = `Pooled SHR (95% CI)` %||% pooled_shr_95_ci %||% pooled_shr_95__ci,
    i2        = as.numeric(`I^2 (%)` %||% i_2____ %||% i_2___ %||% i2 %||% i_2)
  )

# If your file didn’t carry numeric est/ci, leave effect_ci as the formatted string.
# (You already formatted it, so this is fine.)

# ------------------------------------------------------------
# 3) READ: Models pooled meta-analysis (random DL)
# ------------------------------------------------------------
stopifnot(file.exists(models_file))

models_pooled <- read_csv(models_file, show_col_types = FALSE) %>%
  clean_names()

# Expect columns: pollutant, outcome, method, pooled_est, pooled_lcl, pooled_ucl, k, Q, Q_p, tau2
models_random <- models_pooled %>%
  mutate(
    pollutant = norm_pollutant(pollutant),
    outcome   = norm_outcome(outcome),
    method    = case_when(
      str_detect(method, "random") ~ "Random-effects (DL)",
      str_detect(method, "fixed")  ~ "Fixed-effect",
      TRUE ~ method
    )
  ) %>%
  filter(method == "Random-effects (DL)") %>%
  transmute(
    section   = "Regression model meta-analysis",
    pollutant,
    outcome,
    method,
    k         = as.numeric(k),
    estimate  = as.numeric(pooled_est),
    lcl       = as.numeric(pooled_lcl),
    ucl       = as.numeric(pooled_ucl),
    effect_ci = fmt_ci(estimate, lcl, ucl, digits = 2),
    i2        = NA_real_  # not computed in your DL block
  )

# ------------------------------------------------------------
# 4) READ: CIF pooled Day-30 values by bin
#     We compute pooled CIF curves from the per-site CIF master
#     and then extract day=30 pooled CIF by pollutant/outcome/bin.
# ------------------------------------------------------------
stopifnot(file.exists(cif_master_file))

cif_all <- read_csv(cif_master_file, show_col_types = FALSE) %>%
  clean_names()

# We need at least: pollutant, cause_label, day, cif, risk_set
required_cols <- c("pollutant","cause_label","day","cif")
missing <- setdiff(required_cols, names(cif_all))
if (length(missing) > 0) stop("combined_cif_all_sites.csv is missing required columns: ", paste(missing, collapse = ", "))

# risk_set may be missing in some exports; fallback weight = 1
if (!("risk_set" %in% names(cif_all))) {
  message("NOTE: 'risk_set' not found in CIF master file. Using equal weights across rows.")
  cif_all <- cif_all %>% mutate(risk_set = 1)
}

# Choose which bin column to use (priority: US bins if present)
bin_col <- case_when(
  "exposure_bin_us_q4"  %in% names(cif_all) ~ "exposure_bin_us_q4",
  "exposure_bin_us_med" %in% names(cif_all) ~ "exposure_bin_us_med",
  "exposure_bin_us"     %in% names(cif_all) ~ "exposure_bin_us",
  "exposure_bin"        %in% names(cif_all) ~ "exposure_bin",
  TRUE ~ NA_character_
)

if (is.na(bin_col)) stop("No exposure bin column found in CIF master file (looked for exposure_bin_us_q4/us_med/us/exposure_bin).")

# Standardize bins to a single column called 'bin'
cif_all2 <- cif_all %>%
  mutate(
    pollutant = norm_pollutant(pollutant),
    outcome   = norm_outcome(cause_label),
    bin       = as.character(.data[[bin_col]])
  ) %>%
  filter(!is.na(bin), !is.na(outcome), !is.na(pollutant), !is.na(day), !is.na(cif))

# Risk-set weighted pooled curve at each day
cif_pooled <- cif_all2 %>%
  group_by(pollutant, outcome, bin, day) %>%
  summarize(
    n_sites       = n_distinct(site_name %||% site %||% site_id, na.rm = TRUE),
    total_at_risk = sum(pmax(risk_set, 0), na.rm = TRUE),
    cif           = weighted.mean(cif, w = pmax(risk_set, 0), na.rm = TRUE),
    .groups = "drop"
  )

# Day-30 (or nearest) per pollutant/outcome/bin
cif_day30 <- cif_pooled %>%
  group_by(pollutant, outcome, bin) %>%
  group_modify(~ pick_day_value(.x, day_target = 30)) %>%
  ungroup() %>%
  transmute(
    section   = glue("Pooled CIF (risk-set weighted) — {bin_col}"),
    pollutant,
    outcome,
    method    = "Risk-set weighted",
    k         = NA_real_,
    estimate  = cif,
    lcl       = NA_real_,
    ucl       = NA_real_,
    effect_ci = ifelse(is.na(estimate), NA_character_, sprintf("%.3f", estimate)),
    i2        = NA_real_,
    day_used  = day
  )

# ------------------------------------------------------------
# 5) COMBINE: master long + master wide
# ------------------------------------------------------------
master_long <- bind_rows(
  shr_pooled2,
  models_random,
  cif_day30 %>% select(-day_used)
) %>%
  arrange(section, pollutant, outcome, method)

write_csv(master_long, file.path(out_dir, "master_supplement_table_long.csv"))

# Build a wide version: one row per outcome; columns by section/pollutant
master_wide <- master_long %>%
  mutate(
    col = case_when(
      section == "CIF SHR meta-analysis"        ~ paste0("CIF_SHR_REML_", pollutant),
      section == "Regression model meta-analysis" ~ paste0("MODEL_META_DL_", pollutant),
      str_detect(section, "^Pooled CIF")        ~ paste0("CIF_DAY30_", pollutant, "_", str_replace_all(section, ".*—\\s*", ""), "_", "BIN_", if_else(is.na(method), "", method)),
      TRUE ~ paste0(section, "_", pollutant)
    )
  ) %>%
  select(outcome, col, effect_ci, i2, k) %>%
  # For CIF rows, effect_ci is the day-30 CIF; keep it as effect_ci
  group_by(outcome, col) %>%
  summarize(
    value = dplyr::first(effect_ci),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = col, values_from = value) %>%
  arrange(factor(outcome, levels = c(
    "Death","Successful extubation","Persistent respiratory failure",
    "Ventilation hours","ICU length of stay","In-hospital death","30-day death"
  )))

write_csv(master_wide, file.path(out_dir, "master_supplement_table_wide.csv"))

message("Wrote:\n  - ", file.path(out_dir, "master_supplement_table_long.csv"),
        "\n  - ", file.path(out_dir, "master_supplement_table_wide.csv"))

# ------------------------------------------------------------
# 6) OPTIONAL: Write a Word docx with a clean table
# ------------------------------------------------------------
# If you want a Word table, uncomment this block.
# (Requires officer + flextable)
#
# suppressPackageStartupMessages({
#   library(officer)
#   library(flextable)
# })
#
# ft <- flextable(master_wide)
# ft <- autofit(ft)
# ft <- theme_booktabs(ft)
# ft <- fontsize(ft, size = 10, part = "all")
# ft <- bold(ft, part = "header")
#
# doc <- read_docx()
# doc <- body_add_par(doc, "Master Supplement Table", style = "heading 1")
# doc <- body_add_par(doc, "Pooled effects and 30-day pooled CIF values.", style = "Normal")
# doc <- body_add_flextable(doc, ft)
#
# print(doc, target = file.path(out_dir, "master_supplement_table.docx"))
# message("Wrote: ", file.path(out_dir, "master_supplement_table.docx"))
