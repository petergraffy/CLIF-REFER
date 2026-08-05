#!/usr/bin/env Rscript

# Exposure-response curves for the reviewer-optimized primary models.
#
# Curves show model-based predictions on the outcome scale, matching the
# original REFER exposure-response style:
#   - Cox mortality model: predicted mortality probability at a fixed day
#   - Quasi-Poisson VFD model: predicted mean ventilator-free days
#
# The script produces whole-cohort and ARF-subtype versions for PM2.5 and NO2.

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(glue)
  library(janitor)
  library(readr)
  library(scales)
  library(splines)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !file.exists(args[[1]])) {
  stop("Usage: Rscript code_resubmit/03_exposure_response_primary_models.R <analysis_dataset.csv> [output_dir]")
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
repo <- normalizePath(file.path(dirname(input_path), "..", ".."), mustWork = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- args[[2]] %||% file.path(repo, "output", "code_resubmit", run_id)
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
mortality_prediction_day <- as.numeric(Sys.getenv("REFER_MORTALITY_PREDICTION_DAY", "28"))

message("Reading analysis dataset: ", input_path)
message("Output: ", out_dir)

safe_mean <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

top_level <- function(x) {
  x <- droplevels(factor(x))
  tab <- sort(table(x), decreasing = TRUE)
  if (!length(tab)) NA_character_ else names(tab)[[1]]
}

analysis_df <- readr::read_csv(input_path, show_col_types = FALSE, progress = FALSE) %>%
  janitor::clean_names() %>%
  mutate(
    site = site %||% "site",
    mortality_ftime_days = as.numeric(mortality_ftime_days),
    mortality_event = as.integer(mortality_event),
    ventilator_free_days = as.numeric(ventilator_free_days),
    age_10 = age_10 %||% (age / 10),
    pm25_per_5 = pm25_per_5 %||% (pm25_12m_zcta / 5),
    no2_per_10 = no2_per_10 %||% (no2_12m_zcta / 10),
    acs_median_household_income_10k = acs_median_household_income_10k %||% (acs_median_household_income / 10000),
    sex = factor(sex),
    race_ethnicity = forcats::fct_lump_min(factor(race_ethnicity), min = 100, other_level = "Other/Unknown"),
    index_year_f = factor(index_year_f %||% index_year),
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed"))
  ) %>%
  filter(
    is.finite(mortality_ftime_days),
    mortality_ftime_days > 0,
    !is.na(mortality_event),
    !is.na(ventilator_free_days)
  )

site_name <- dplyr::first(stats::na.omit(analysis_df$site)) %||% "site"
spline_df <- 3L

social_covars <- c(
  "acs_pct_poverty",
  "acs_pct_unemployed",
  "acs_pct_no_vehicle",
  "acs_pct_nonwhite",
  "acs_median_household_income_10k",
  "acs_pct_bachelor_plus"
)

base_covars <- c(
  "age_10",
  "sex",
  "race_ethnicity",
  "charlson_score",
  "index_year_f",
  intersect(social_covars, names(analysis_df))
)

pollutant_specs <- tibble::tribble(
  ~pollutant, ~term, ~raw_col, ~other_term, ~display_label, ~x_label,
  "NO2", "no2_per_10", "no2_12m_zcta", "pm25_per_5", "Nitrogen~dioxide~(NO[2])", expression(NO[2]~"(ppb)"),
  "PM2.5", "pm25_per_5", "pm25_12m_zcta", "no2_per_10", "Fine~particulate~matter~(PM[2.5])", expression(PM[2.5]~"("*mu*"g/m"^3*")")
)

outcome_specs <- tibble::tribble(
  ~outcome, ~endpoint, ~y_label,
  glue("Mortality by day {mortality_prediction_day}"), "mortality", "Predicted probability",
  "Ventilator-free days", "vfd", "Predicted mean VFD"
)

drop_uninformative_covars <- function(model_df, covars, exposure_terms) {
  purrr::keep(covars, function(covar) {
    if (covar %in% exposure_terms) return(TRUE)
    x <- model_df[[covar]]
    if (is.factor(x) || is.character(x)) return(dplyr::n_distinct(x, na.rm = TRUE) >= 2)
    isTRUE(stats::var(as.numeric(x), na.rm = TRUE) > 0)
  }) %>% unname()
}

make_prediction_grid <- function(data, term, raw_col, other_term, by_subtype) {
  exposure_q <- quantile(data[[raw_col]], probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
  exposure_grid <- seq(exposure_q[[1]], exposure_q[[3]], length.out = 160)

  if (by_subtype) {
    grid <- tidyr::expand_grid(
      exposure_raw = exposure_grid,
      arf_subtype = levels(droplevels(data$arf_subtype))
    )
  } else {
    grid <- tibble(exposure_raw = exposure_grid)
  }

  grid[[term]] <- if (term == "pm25_per_5") grid$exposure_raw / 5 else grid$exposure_raw / 10

  if (other_term %in% names(data)) {
    grid[[other_term]] <- median(data[[other_term]], na.rm = TRUE)
  }

  grid %>%
    mutate(
      age_10 = safe_mean(data$age_10),
      sex = top_level(data$sex),
      race_ethnicity = top_level(data$race_ethnicity),
      charlson_score = safe_mean(data$charlson_score),
      index_year_f = top_level(data$index_year_f),
      acs_pct_poverty = safe_mean(data$acs_pct_poverty),
      acs_pct_unemployed = safe_mean(data$acs_pct_unemployed),
      acs_pct_no_vehicle = safe_mean(data$acs_pct_no_vehicle),
      acs_pct_nonwhite = safe_mean(data$acs_pct_nonwhite),
      acs_median_household_income_10k = safe_mean(data$acs_median_household_income_10k),
      acs_pct_bachelor_plus = safe_mean(data$acs_pct_bachelor_plus),
      reference_exposure = exposure_q[[2]],
      exposure_p05 = exposure_q[[1]],
      exposure_p95 = exposure_q[[3]]
    ) %>%
    mutate(
      across(any_of(c("sex", "race_ethnicity", "index_year_f", "arf_subtype")), factor)
    )
}

align_factor_levels <- function(grid, data) {
  for (nm in intersect(c("sex", "race_ethnicity", "index_year_f", "arf_subtype"), names(grid))) {
    grid[[nm]] <- factor(as.character(grid[[nm]]), levels = levels(droplevels(data[[nm]])))
  }
  grid
}

predict_outcome_scale <- function(fit, pred_grid, endpoint) {
  if (endpoint == "mortality") {
    pr <- predict(fit, newdata = pred_grid, type = "lp", se.fit = TRUE, reference = "zero")
    bh <- survival::basehaz(fit, centered = FALSE)
    idx <- which(bh$time <= mortality_prediction_day)
    h0 <- if (length(idx)) bh$hazard[max(idx)] else 0
    transform_lp <- function(lp) 1 - exp(-h0 * exp(lp))
    return(list(
      estimate = transform_lp(pr$fit),
      conf_low = transform_lp(pr$fit - 1.96 * pr$se.fit),
      conf_high = transform_lp(pr$fit + 1.96 * pr$se.fit)
    ))
  }

  pr <- predict(fit, newdata = pred_grid, type = "link", se.fit = TRUE)
  list(
    estimate = exp(pr$fit),
    conf_low = exp(pr$fit - 1.96 * pr$se.fit),
    conf_high = exp(pr$fit + 1.96 * pr$se.fit)
  )
}

fit_curve <- function(outcome, endpoint, pollutant, term, raw_col, other_term, display_label, by_subtype) {
  covars <- c(term, other_term, base_covars, raw_col)
  if (by_subtype) covars <- c(covars, "arf_subtype")

  model_df <- analysis_df %>%
    dplyr::select(mortality_ftime_days, mortality_event, ventilator_free_days, all_of(covars)) %>%
    tidyr::drop_na() %>%
    mutate(
      sex = droplevels(sex),
      race_ethnicity = droplevels(race_ethnicity),
      index_year_f = droplevels(index_year_f),
      across(any_of("arf_subtype"), droplevels)
    )

  exposure_terms <- c(term, other_term)
  rhs_covars <- drop_uninformative_covars(model_df, c(other_term, base_covars), exposure_terms = c(term, other_term))
  spline_term <- paste0("splines::ns(", term, ", df = ", spline_df, ")")
  if (by_subtype) spline_term <- paste0(spline_term, " * arf_subtype")
  rhs <- paste(c(spline_term, rhs_covars), collapse = " + ")

  fml <- if (endpoint == "mortality") {
    as.formula(paste0("Surv(mortality_ftime_days, mortality_event) ~ ", rhs))
  } else {
    as.formula(paste0("ventilator_free_days ~ ", rhs))
  }

  fit <- if (endpoint == "mortality") {
    survival::coxph(fml, data = model_df, ties = "efron")
  } else {
    stats::glm(fml, data = model_df, family = quasipoisson(link = "log"))
  }

  pred_grid <- make_prediction_grid(model_df, term, raw_col, other_term, by_subtype) %>%
    align_factor_levels(model_df)
  pred <- predict_outcome_scale(fit, pred_grid, endpoint)

  pred_grid %>%
    mutate(
      site = site_name,
      curve_type = if_else(by_subtype, "By ARF subtype", "Whole cohort"),
      outcome = outcome,
      endpoint = endpoint,
      pollutant = pollutant,
      pollutant_label = display_label,
      estimate = pred$estimate,
      conf_low = pred$conf_low,
      conf_high = pred$conf_high,
      n = nrow(model_df),
      events = if_else(endpoint == "mortality", sum(model_df$mortality_event == 1L), NA_integer_),
      adjustment_set = paste(rhs_covars, collapse = " + "),
      mortality_prediction_day = if_else(endpoint == "mortality", mortality_prediction_day, NA_real_)
    )
}

prediction_df <- tidyr::crossing(outcome_specs, pollutant_specs) %>%
  rowwise() %>%
  do(bind_rows(
    fit_curve(
      .$outcome, .$endpoint, .$pollutant, .$term, .$raw_col, .$other_term,
      .$display_label, by_subtype = FALSE
    ),
    fit_curve(
      .$outcome, .$endpoint, .$pollutant, .$term, .$raw_col, .$other_term,
      .$display_label, by_subtype = TRUE
    )
  )) %>%
  ungroup() %>%
  mutate(
    outcome = factor(outcome, levels = outcome_specs$outcome),
    curve_type = factor(curve_type, levels = c("Whole cohort", "By ARF subtype")),
    pollutant = factor(pollutant, levels = pollutant_specs$pollutant),
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$display_label),
    arf_subtype = forcats::fct_na_value_to_level(
      factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed")),
      level = "Overall"
    )
  )

readr::write_csv(prediction_df, file.path(out_dir, "primary_exposure_response_predictions.csv"))

rug_df <- bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
  spec <- pollutant_specs[i, ]
  ranges <- prediction_df %>%
    filter(pollutant == spec$pollutant) %>%
    summarise(exposure_p05 = first(exposure_p05), exposure_p95 = first(exposure_p95), .groups = "drop")

  analysis_df %>%
    transmute(
      pollutant = spec$pollutant,
      pollutant_label = spec$display_label,
      exposure_raw = .data[[spec$raw_col]],
      arf_subtype
    ) %>%
    filter(
      !is.na(exposure_raw),
      exposure_raw >= ranges$exposure_p05,
      exposure_raw <= ranges$exposure_p95
    )
})) %>%
  crossing(outcome = outcome_specs$outcome, curve_type = levels(prediction_df$curve_type)) %>%
  mutate(
    outcome = factor(outcome, levels = levels(prediction_df$outcome)),
    curve_type = factor(curve_type, levels = levels(prediction_df$curve_type)),
    pollutant = factor(pollutant, levels = levels(prediction_df$pollutant)),
    pollutant_label = factor(pollutant_label, levels = levels(prediction_df$pollutant_label)),
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed")),
    rug_group = if_else(curve_type == "Whole cohort", "Overall", as.character(arf_subtype))
  )

curve_colors <- c(
  "Overall" = "black",
  "Hypoxemic" = "#0072B2",
  "Hypercapnic" = "#009E73",
  "Mixed" = "#D55E00"
)

theme_response <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

plot_curves <- function(curve_type_value, filename_stub) {
  df <- prediction_df %>% filter(curve_type == curve_type_value)
  rugs <- rug_df %>% filter(curve_type == curve_type_value)

  p <- ggplot(df, aes(exposure_raw, estimate, color = arf_subtype, fill = arf_subtype)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.14, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.05) +
    geom_rug(
      data = rugs,
      aes(x = exposure_raw, color = factor(rug_group, levels = names(curve_colors))),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.12,
      linewidth = 0.2,
      show.legend = FALSE
    ) +
    facet_grid(outcome ~ pollutant_label, scales = "free", labeller = labeller(pollutant_label = label_parsed)) +
    scale_color_manual(values = curve_colors, drop = TRUE) +
    scale_fill_manual(values = curve_colors, drop = TRUE) +
    scale_y_continuous(labels = label_number(accuracy = 0.01, trim = TRUE)) +
    labs(
      title = glue("Primary Exposure-Response Curves: {curve_type_value}"),
      x = expression("Exposure concentration (NO"[2]*" ppb; PM"[2.5]*" "*mu*"g/m"^3*")"),
      y = "Model-predicted outcome",
      color = NULL,
      fill = NULL
    ) +
    theme_response()

  ggsave(file.path(fig_dir, paste0(filename_stub, ".png")), p, width = 12.5, height = 8.5, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename_stub, ".pdf")), p, width = 12.5, height = 8.5)
  p
}

p_overall <- plot_curves("Whole cohort", "primary_exposure_response_predicted_outcome_whole_cohort")
p_subtype <- plot_curves("By ARF subtype", "primary_exposure_response_predicted_outcome_by_arf_subtype")

message("Wrote:")
message(" - ", file.path(out_dir, "primary_exposure_response_predictions.csv"))
message(" - ", file.path(fig_dir, "primary_exposure_response_predicted_outcome_whole_cohort.png"))
message(" - ", file.path(fig_dir, "primary_exposure_response_predicted_outcome_by_arf_subtype.png"))
