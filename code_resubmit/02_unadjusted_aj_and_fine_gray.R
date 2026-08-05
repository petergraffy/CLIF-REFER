#!/usr/bin/env Rscript

# Unadjusted Aalen-Johansen CIFs and adjusted Fine-Gray models
#
# Input:
#   A reviewer-optimized analysis dataset with:
#     ftime_days, event_code, pm25_12m_zcta/no2_12m_zcta, age_10, sex,
#     race_ethnicity, charlson_score, index_year_f, and ACS covariates.
#
# Event coding:
#   0 = censored
#   1 = successful extubation
#   2 = death
#   3 = persistent respiratory failure

suppressPackageStartupMessages({
  library(broom)
  library(cmprsk)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(glue)
  library(janitor)
  library(prodlim)
  library(readr)
  library(riskRegression)
  library(scales)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !file.exists(args[[1]])) {
  stop("Usage: Rscript code_resubmit/02_unadjusted_aj_and_fine_gray.R <analysis_dataset.csv> [output_dir]")
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
repo <- normalizePath(file.path(dirname(input_path), "..", ".."), mustWork = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- args[[2]] %||% file.path(repo, "output", "code_resubmit", run_id)
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading analysis dataset: ", input_path)
message("Output: ", out_dir)

analysis_df <- readr::read_csv(input_path, show_col_types = FALSE, progress = FALSE) %>%
  janitor::clean_names() %>%
  mutate(
    site = site %||% "site",
    site = as.character(site),
    ftime_days = as.numeric(ftime_days),
    event_code = as.integer(event_code),
    sex = factor(sex),
    race_ethnicity = forcats::fct_lump_min(factor(race_ethnicity), min = 100, other_level = "Other/Unknown"),
    index_year_f = factor(index_year_f %||% index_year),
    age_10 = age_10 %||% (age / 10),
    pm25_per_5 = pm25_per_5 %||% (pm25_12m_zcta / 5),
    no2_per_10 = no2_per_10 %||% (no2_12m_zcta / 10),
    acs_median_household_income_10k = acs_median_household_income_10k %||% (acs_median_household_income / 10000)
  ) %>%
  filter(is.finite(ftime_days), ftime_days > 0, event_code %in% 0:3)

site_name <- dplyr::first(stats::na.omit(analysis_df$site)) %||% "site"

event_labels <- c(
  `1` = "Successful extubation",
  `2` = "Death",
  `3` = "Persistent respiratory failure"
)

pollutants <- tibble::tribble(
  ~pollutant, ~raw_col, ~term, ~unit_label,
  "PM2.5", "pm25_12m_zcta", "pm25_per_5", "per 5 ug/m3",
  "NO2", "no2_12m_zcta", "no2_per_10", "per 10 ppb"
)

make_quartile <- function(x) {
  q <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)
  q <- unique(q)
  if (length(q) < 5) return(factor(rep(NA_character_, length(x)), levels = paste0("Q", 1:4)))
  cut(x, breaks = q, include.lowest = TRUE, labels = paste0("Q", 1:4))
}

tidy_cuminc <- function(ci, pollutant) {
  rows <- lapply(names(ci), function(nm) {
    x <- ci[[nm]]
    if (!is.list(x) || is.null(x$time) || is.null(x$est)) return(NULL)
    toks <- strsplit(nm, "\\s+")[[1]]
    cause <- suppressWarnings(as.integer(toks[length(toks)]))
    group <- paste(toks[-length(toks)], collapse = " ")
    tibble(
      pollutant = pollutant,
      exposure_quartile = group,
      event_code = cause,
      outcome = unname(event_labels[as.character(cause)]),
      time = x$time,
      cif = x$est,
      var = x$var %||% NA_real_
    )
  })
  bind_rows(rows) %>%
    filter(event_code %in% 1:3, is.finite(time), is.finite(cif)) %>%
    mutate(
      outcome = factor(outcome, levels = unname(event_labels)),
      exposure_quartile = factor(exposure_quartile, levels = paste0("Q", 1:4))
    )
}

aj_df <- bind_rows(lapply(seq_len(nrow(pollutants)), function(i) {
  spec <- pollutants[i, ]
  dat <- analysis_df %>%
    mutate(exposure_quartile = make_quartile(.data[[spec$raw_col]])) %>%
    filter(!is.na(exposure_quartile))
  ci <- cmprsk::cuminc(dat$ftime_days, dat$event_code, group = dat$exposure_quartile, cencode = 0)
  tidy_cuminc(ci, spec$pollutant)
}))

readr::write_csv(aj_df, file.path(out_dir, "unadjusted_aalen_johansen_cif_by_exposure_quartile.csv"))

quartile_counts <- bind_rows(lapply(seq_len(nrow(pollutants)), function(i) {
  spec <- pollutants[i, ]
  analysis_df %>%
    mutate(exposure_quartile = make_quartile(.data[[spec$raw_col]])) %>%
    filter(!is.na(exposure_quartile)) %>%
    count(pollutant = spec$pollutant, exposure_quartile, event_code, name = "n") %>%
    mutate(outcome = case_when(
      event_code == 0L ~ "Censored",
      event_code %in% 1:3 ~ unname(event_labels[as.character(event_code)]),
      TRUE ~ NA_character_
    ))
}))
readr::write_csv(quartile_counts, file.path(out_dir, "unadjusted_aalen_johansen_event_counts_by_quartile.csv"))

quartile_colors <- c(
  "Q1" = "#2166AC",
  "Q2" = "#67A9CF",
  "Q3" = "#F4A582",
  "Q4" = "#B2182B"
)

p_aj <- ggplot(aj_df, aes(time, cif, color = exposure_quartile)) +
  geom_step(linewidth = 0.9) +
  facet_grid(outcome ~ pollutant, scales = "free_y") +
  scale_color_manual(values = quartile_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 7), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = "Unadjusted Aalen-Johansen Cumulative Incidence by Exposure Quartile",
    x = "Days since ARF onset",
    y = "Cumulative incidence",
    color = "Exposure quartile"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "unadjusted_aalen_johansen_cif_by_exposure_quartile.png"), p_aj, width = 11, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "unadjusted_aalen_johansen_cif_by_exposure_quartile.pdf"), p_aj, width = 11, height = 10)

primary_social_covars <- c(
  "acs_pct_poverty",
  "acs_pct_unemployed",
  "acs_pct_no_vehicle",
  "acs_pct_nonwhite",
  "acs_median_household_income_10k",
  "acs_pct_bachelor_plus"
)

adjustment_covars <- c(
  "age_10",
  "sex",
  "race_ethnicity",
  "charlson_score",
  "index_year_f",
  primary_social_covars
)
adjustment_covars <- intersect(adjustment_covars, names(analysis_df))

drop_uninformative_covars <- function(model_df, covars, exposure_terms) {
  purrr::keep(covars, function(covar) {
    if (covar %in% exposure_terms) return(TRUE)
    x <- model_df[[covar]]
    if (is.factor(x) || is.character(x)) return(dplyr::n_distinct(x, na.rm = TRUE) >= 2)
    isTRUE(stats::var(as.numeric(x), na.rm = TRUE) > 0)
  }) %>% unname()
}

fit_fgr <- function(df, cause_code, exposure_terms, model_label) {
  outcome_label <- unname(event_labels[as.character(cause_code)])
  message("Fine-Gray: ", outcome_label, " / ", model_label)
  covars <- c(exposure_terms, adjustment_covars)
  model_df <- df %>%
    dplyr::select(ftime_days, event_code, all_of(covars)) %>%
    tidyr::drop_na() %>%
    mutate(
      event = factor(
        event_code,
        levels = c(0, 1, 2, 3),
        labels = c("censor", unname(event_labels))
      )
    )
  covars <- drop_uninformative_covars(model_df, covars, exposure_terms)
  model_df <- model_df %>%
    dplyr::select(ftime_days, event, all_of(covars))

  fg_df <- survival::finegray(
    as.formula(paste0("Surv(ftime_days, event) ~ ", paste(covars, collapse = " + "))),
    data = model_df,
    etype = outcome_label
  )

  fit <- survival::coxph(
    as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", paste(covars, collapse = " + "))),
    data = fg_df,
    weights = fgwt,
    ties = "efron"
  )

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
    site = site_name,
    outcome = outcome_label,
    model = model_label,
      term,
    n = nrow(model_df),
      events = sum(model_df$event == outcome_label),
      subdistribution_hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      adjustment_set = paste(setdiff(covars, exposure_terms), collapse = " + ")
    )
}

fg_specs <- list(
  list(exposure_terms = "pm25_per_5", model_label = "PM25 single-pollutant"),
  list(exposure_terms = "no2_per_10", model_label = "NO2 single-pollutant"),
  list(exposure_terms = c("pm25_per_5", "no2_per_10"), model_label = "PM25 + NO2")
)

fine_gray_results <- bind_rows(lapply(1:3, function(cause_code) {
  bind_rows(lapply(fg_specs, function(spec) {
    fit_fgr(analysis_df, cause_code, spec$exposure_terms, spec$model_label)
  }))
}))

readr::write_csv(fine_gray_results, file.path(out_dir, "fine_gray_results_same_covariates_as_cox.csv"))

message("Wrote outputs:")
message(" - ", file.path(out_dir, "unadjusted_aalen_johansen_cif_by_exposure_quartile.csv"))
message(" - ", file.path(out_dir, "fine_gray_results_same_covariates_as_cox.csv"))
message(" - ", file.path(fig_dir, "unadjusted_aalen_johansen_cif_by_exposure_quartile.png"))
