# ================================================================================================
# Manifest-driven pooled analysis
#
# Reads only analysis_ready/manifest.csv active rows. This replaces ad hoc globbing
# over sites/analysis and avoids duplicated Finder-copy files.
# ================================================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(fs)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(cmd_file)) normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = TRUE) else NA_character_
}
repo <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE) else normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(repo, "analysis_ready", "manifest.csv"))) repo <- normalizePath(getwd(), mustWork = TRUE)

ready_dir <- file.path(repo, "analysis_ready")
manifest_path <- file.path(ready_dir, "manifest.csv")
tables_dir <- file.path(ready_dir, "pooled", "tables")
figures_dir <- file.path(ready_dir, "pooled", "figures")
dir_create(tables_dir)
dir_create(figures_dir)

manifest <- readr::read_csv(manifest_path, show_col_types = FALSE) %>%
  mutate(active = as.logical(active), exists = as.logical(exists))

active_files <- function(family) {
  out <- manifest %>%
    filter(.data$family == !!family, active, exists) %>%
    mutate(path = file.path(repo, curated_relpath))
  if (!nrow(out)) stop("No active manifest rows for family: ", family)
  out
}

norm_pollutant <- function(x) {
  xu <- str_to_upper(as.character(x))
  case_when(
    str_detect(xu, "PM2\\.5|PM25|PM_?2") ~ "PM2.5",
    str_detect(xu, "NO2|NO_?2") ~ "NO2",
    TRUE ~ as.character(x)
  )
}

norm_outcome <- function(x) {
  xl <- str_to_lower(as.character(x))
  case_when(
    str_detect(xl, "death30d|30") & str_detect(xl, "death|mort") ~ "30-day mortality",
    str_detect(xl, "inhosp|in_hosp|in-hospital") ~ "In-hospital mortality",
    str_detect(xl, "vent") ~ "Ventilation hours",
    str_detect(xl, "icu") & str_detect(xl, "los") ~ "ICU LOS",
    str_detect(xl, "successful") ~ "Successful extubation",
    str_detect(xl, "persistent|prf") ~ "Persistent respiratory failure",
    str_detect(xl, "^death$") ~ "Death",
    TRUE ~ as.character(x)
  )
}

yi_from_exp_ci <- function(df) {
  if (!("std.error" %in% names(df))) df$std.error <- NA_real_
  df %>%
    mutate(
      estimate = as.numeric(estimate),
      conf.low = as.numeric(conf.low),
      conf.high = as.numeric(conf.high),
      yi = log(estimate),
      sei = case_when(
        is.finite(conf.low) & is.finite(conf.high) & conf.low > 0 & conf.high > 0 ~
          (log(conf.high) - log(conf.low)) / (2 * 1.96),
        is.finite(std.error) & std.error > 0 ~ as.numeric(std.error),
        TRUE ~ NA_real_
      )
    )
}

meta_dl <- function(yi, sei) {
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[keep]
  sei <- sei[keep]
  k <- length(yi)
  if (k == 0) return(tibble(k = 0, beta = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_,
                            tau2 = NA_real_, q = NA_real_, q_p = NA_real_, i2 = NA_real_))
  vi <- sei^2
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - fixed)^2)
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (k > 1 && c_val > 0) max(0, (q - (k - 1)) / c_val) else 0
  wr <- 1 / (vi + tau2)
  beta <- sum(wr * yi) / sum(wr)
  se <- sqrt(1 / sum(wr))
  tibble(
    k = k,
    beta = beta,
    se = se,
    lo = beta - 1.96 * se,
    hi = beta + 1.96 * se,
    tau2 = tau2,
    q = q,
    q_p = if (k > 1) pchisq(q, df = k - 1, lower.tail = FALSE) else NA_real_,
    i2 = if (k > 1 && q > 0) max(0, (q - (k - 1)) / q) * 100 else 0
  )
}

# ------------------------------------------------------------------------------------------------
# 1) CIF master from manifest-approved site files
# ------------------------------------------------------------------------------------------------

read_one_cif <- function(row) {
  df <- readr::read_csv(row$path, show_col_types = FALSE, progress = FALSE)
  names(df) <- names(df) %>% str_to_lower() %>% str_replace_all("[^a-z0-9]+", "_") %>% str_replace_all("^_|_$", "")
  source_pollutant <- if ("pollutant" %in% names(df)) df$pollutant else row$pollutant
  source_cause <- if ("cause_label" %in% names(df)) df$cause_label else if ("cause" %in% names(df)) df$cause else NA_character_
  source_cause <- case_when(
    as.character(source_cause) == "1" ~ "Successful extubation",
    as.character(source_cause) == "2" ~ "Death",
    as.character(source_cause) == "3" ~ "Persistent respiratory failure",
    TRUE ~ as.character(source_cause)
  )
  df %>%
    mutate(
      manifest_site = row$site,
      manifest_pollutant = row$pollutant,
      manifest_file = row$curated_relpath,
      pollutant = norm_pollutant(source_pollutant),
      cause_label = norm_outcome(source_cause),
      site_name = row$site
    )
}

cif_manifest <- active_files("cif_plotdf")
cif_all <- pmap_dfr(cif_manifest, function(...) read_one_cif(list(...))) %>%
  distinct()

readr::write_csv(cif_all, file.path(tables_dir, "combined_cif_all_sites_manifest.csv"))

bin_col <- case_when(
  "exposure_bin_us_q4" %in% names(cif_all) ~ "exposure_bin_us_q4",
  "exposure_bin_us_med" %in% names(cif_all) ~ "exposure_bin_us_med",
  "exposure_bin_us" %in% names(cif_all) ~ "exposure_bin_us",
  "exposure_bin" %in% names(cif_all) ~ "exposure_bin",
  "exposure_group" %in% names(cif_all) ~ "exposure_group",
  "legend_label" %in% names(cif_all) ~ "legend_label",
  TRUE ~ NA_character_
)
if (is.na(bin_col)) stop("No exposure-bin column found in manifest CIF inputs.")
if (!("risk_set" %in% names(cif_all))) cif_all$risk_set <- 1

cif_day30 <- cif_all %>%
  mutate(
    bin = as.character(.data[[bin_col]]),
    day = as.numeric(day),
    cif = as.numeric(cif),
    risk_set = pmax(as.numeric(risk_set), 0)
  ) %>%
  filter(is.finite(day), is.finite(cif), !is.na(bin), !is.na(pollutant), !is.na(cause_label)) %>%
  group_by(pollutant, outcome = cause_label, bin, day) %>%
  summarise(
    n_sites = n_distinct(site_name),
    total_at_risk = sum(risk_set, na.rm = TRUE),
    cif = weighted.mean(cif, w = if_else(risk_set > 0, risk_set, 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(pollutant, outcome, bin) %>%
  slice_min(abs(day - 30), n = 1, with_ties = FALSE) %>%
  ungroup()

readr::write_csv(cif_day30, file.path(tables_dir, "pooled_cif_day30_by_bin_manifest.csv"))

# ------------------------------------------------------------------------------------------------
# 2) Cumulative pollutant model meta-analysis
# ------------------------------------------------------------------------------------------------

read_model_term <- function(row) {
  df <- readr::read_csv(row$path, show_col_types = FALSE, progress = FALSE)
  term_target <- if (row$pollutant == "NO2") {
    c("no2_mean_cummean_2018toYr", "no2_mean", "no2_10")
  } else {
    c("pm25_mean_cummean_2018toYr", "pm25_mean")
  }
  df %>%
    filter(term %in% term_target) %>%
    slice(1) %>%
    yi_from_exp_ci() %>%
    transmute(site = row$site, pollutant = row$pollutant, outcome = row$outcome,
              term, estimate, conf.low, conf.high, yi, sei)
}

model_effects <- active_files("model_cumulative") %>%
  pmap_dfr(function(...) read_model_term(list(...))) %>%
  filter(is.finite(yi), is.finite(sei), sei > 0)

readr::write_csv(model_effects, file.path(tables_dir, "per_site_cumulative_model_effects_manifest.csv"))

model_pooled <- model_effects %>%
  group_by(pollutant, outcome) %>%
  summarise(meta = list(meta_dl(yi, sei)), .groups = "drop") %>%
  unnest(meta) %>%
  mutate(
    beta = as.numeric(beta),
    lo = as.numeric(lo),
    hi = as.numeric(hi),
    outcome = norm_outcome(outcome),
    effect = exp(beta),
    lo = exp(lo),
    hi = exp(hi),
    method = "Random-effects (DL)"
  ) %>%
  select(pollutant, outcome, k, method, effect, lo, hi, tau2, q, q_p, i2)

readr::write_csv(model_pooled, file.path(tables_dir, "pooled_cumulative_model_effects_manifest.csv"))

# ------------------------------------------------------------------------------------------------
# 3) NO2 by ARF subtype model meta-analysis
#
# Model tidy CSVs were saved exponentiated. Convert OR/IRR and CIs back to log scale before pooling.
# For subtype-specific effects, covariance between main NO2 and interaction terms is unavailable in
# tidy CSVs, so the SE is conservative only if covariance is non-positive and approximate otherwise.
# ------------------------------------------------------------------------------------------------

clean_subtype <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace("^arf subtype", "") %>%
    str_trim() %>%
    str_to_title()
}

read_subtype_effects <- function(row) {
  df <- readr::read_csv(row$path, show_col_types = FALSE, progress = FALSE) %>%
    yi_from_exp_ci()
  main <- df %>% filter(term == "no2_10") %>% slice(1)
  if (!nrow(main)) return(tibble())
  out <- tibble(
    site = row$site,
    outcome = row$outcome,
    subtype = "Hypoxemic",
    yi = main$yi[[1]],
    sei = main$sei[[1]],
    se_method = "main_effect"
  )
  interactions <- df %>%
    filter(str_detect(term, "^no2_10:arf_subtype")) %>%
    mutate(subtype = clean_subtype(str_replace(term, "^no2_10:arf_subtype", "")))
  if (nrow(interactions)) {
    out <- bind_rows(
      out,
      interactions %>%
        transmute(
          site = row$site,
          outcome = row$outcome,
          subtype,
          yi = main$yi[[1]] + yi,
          sei = sqrt(main$sei[[1]]^2 + sei^2),
          se_method = "main_plus_interaction_no_covariance"
        )
    )
  }
  out
}

subtype_effects <- active_files("model_subtype") %>%
  pmap_dfr(function(...) read_subtype_effects(list(...))) %>%
  filter(is.finite(yi), is.finite(sei), sei > 0)

readr::write_csv(subtype_effects, file.path(tables_dir, "per_site_no2_subtype_effects_manifest.csv"))

subtype_pooled <- subtype_effects %>%
  group_by(outcome, subtype) %>%
  summarise(meta = list(meta_dl(yi, sei)), .groups = "drop") %>%
  unnest(meta) %>%
  mutate(
    outcome = norm_outcome(outcome),
    measure = if_else(outcome == "Ventilation hours", "IRR per 10 ppb NO2", "OR per 10 ppb NO2"),
    effect = exp(beta),
    lo = exp(lo),
    hi = exp(hi),
    method = "Random-effects (DL)"
  ) %>%
  select(outcome, subtype, k, method, measure, effect, lo, hi, tau2, q, q_p, i2)

readr::write_csv(subtype_pooled, file.path(tables_dir, "pooled_no2_by_subtype_manifest.csv"))

plot_subtype <- subtype_pooled %>%
  mutate(label = paste(outcome, subtype, sep = " | "))

p <- ggplot(plot_subtype, aes(x = effect, y = reorder(label, effect), xmin = lo, xmax = hi)) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_errorbar(orientation = "y", width = 0.2) +
  geom_point(size = 2.6) +
  scale_x_log10() +
  labs(x = "Pooled NO2 effect", y = NULL) +
  theme_classic(base_size = 12)

ggsave(file.path(figures_dir, "pooled_no2_by_subtype_manifest.png"), p, width = 10, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "pooled_no2_by_subtype_manifest.pdf"), p, width = 10, height = 6)

message("Manifest-driven pooled analysis complete.")
message("Tables:  ", tables_dir)
message("Figures: ", figures_dir)
