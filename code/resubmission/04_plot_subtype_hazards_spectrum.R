#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(survival)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(cmd_file)) normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = TRUE) else NA_character_
}
repo <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE) else getwd()

args <- commandArgs(trailingOnly = TRUE)
out_root <- args[[1]] %||% NA_character_
if (is.na(out_root)) {
  candidates <- list.dirs(file.path(repo, "output", "resubmission"), recursive = FALSE, full.names = TRUE)
  candidates <- candidates[file.exists(file.path(candidates, "resubmission_analysis_dataset.csv"))]
  if (!length(candidates)) stop("No resubmission output directory with analysis dataset found.")
  out_root <- candidates[which.max(file.info(candidates)$mtime)]
}
out_root <- normalizePath(out_root, mustWork = TRUE)
fig_dir <- file.path(out_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

subtype_colors <- c(
  "Hypoxemic" = "#0072B2",
  "Hypercapnic" = "#009E73",
  "Mixed" = "#D55E00"
)

theme_spectrum <- function(base_size = 17) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
      panel.border = element_rect(fill = NA, color = "grey55", linewidth = 0.55),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      plot.title = element_text(face = "bold", size = rel(1.25), hjust = 0.02, margin = margin(b = 7)),
      plot.subtitle = element_text(color = "grey25", hjust = 0.02, margin = margin(b = 10)),
      strip.text = element_text(face = "bold", color = "grey10", margin = margin(5, 5, 5, 5)),
      strip.background = element_rect(fill = "grey94", color = "grey68", linewidth = 0.55),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.98)),
      legend.key.width = unit(32, "pt"),
      panel.spacing.x = unit(12, "pt"),
      panel.spacing.y = unit(12, "pt"),
      plot.margin = margin(8, 16, 8, 8)
    )
}

top_level <- function(x) {
  x <- droplevels(factor(x))
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[[1]]
}

safe_mean <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

outcome_specs <- tibble::tribble(
  ~cause_code, ~outcome,
  1L, "Extubation",
  2L, "Death",
  3L, "Persistent RF"
)

pollutant_specs <- tibble::tribble(
  ~pollutant, ~term, ~adjustment_term, ~raw_col, ~raw_label,
  "Fine particulate matter (PM2.5)", "pm25_per_5", "no2_per_10", "pm25_12m_zcta", "PM2.5, 12-month mean (ug/m3)",
  "Nitrogen dioxide (NO2)", "no2_per_10", "pm25_per_5", "no2_12m_zcta", "NO2, 12-month mean/prior-year fallback (ppb)"
)

analysis_df <- read_csv(
  file.path(out_root, "resubmission_analysis_dataset.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed")),
    sex = factor(sex),
    race_ethnicity = factor(race_ethnicity),
    index_year_f = factor(index_year_f),
    index_year_f = droplevels(index_year_f),
    acs_median_household_income_10k = acs_median_household_income_10k %||% (acs_median_household_income / 10000)
  ) %>%
  filter(!is.na(arf_subtype))

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
  social_covars
)

make_prediction_rows <- function(data, fit, pollutant_term, adjustment_term, raw_col, outcome_label, pollutant_label) {
  exposure_q <- quantile(data[[raw_col]], probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
  exposure_grid <- seq(exposure_q[[1]], exposure_q[[3]], length.out = 120)

  subtype_levels <- levels(droplevels(data$arf_subtype))
  pred_grid <- tidyr::expand_grid(
    exposure_raw = exposure_grid,
    arf_subtype = subtype_levels
  )
  pred_grid[[pollutant_term]] <- if (pollutant_term == "pm25_per_5") {
    pred_grid$exposure_raw / 5
  } else if (pollutant_term == "no2_per_10") {
    pred_grid$exposure_raw / 10
  } else {
    pred_grid$exposure_raw
  }

  pred_grid <- pred_grid %>%
    mutate(
      !!adjustment_term := median(data[[adjustment_term]], na.rm = TRUE),
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
      arf_subtype = factor(arf_subtype, levels = levels(data$arf_subtype)),
      sex = factor(sex, levels = levels(data$sex)),
      race_ethnicity = factor(race_ethnicity, levels = levels(data$race_ethnicity)),
      index_year_f = factor(index_year_f, levels = levels(data$index_year_f))
    )

  ref_grid <- pred_grid %>%
    mutate(exposure_raw = exposure_q[[2]])
  ref_grid[[pollutant_term]] <- if (pollutant_term == "pm25_per_5") {
    exposure_q[[2]] / 5
  } else if (pollutant_term == "no2_per_10") {
    exposure_q[[2]] / 10
  } else {
    exposure_q[[2]]
  }

  x_pred <- model.matrix(delete.response(terms(fit)), pred_grid)
  x_ref <- model.matrix(delete.response(terms(fit)), ref_grid)
  beta <- coef(fit)
  vc <- vcov(fit)

  keep <- names(beta)
  x_diff <- x_pred[, keep, drop = FALSE] - x_ref[, keep, drop = FALSE]
  log_hr <- as.vector(x_diff %*% beta)
  se_log_hr <- sqrt(rowSums((x_diff %*% vc) * x_diff))

  pred_grid %>%
    transmute(
      outcome = outcome_label,
      pollutant = pollutant_label,
      exposure_raw,
      arf_subtype,
      hazard_ratio = exp(log_hr),
      conf_low = exp(log_hr - 1.96 * se_log_hr),
      conf_high = exp(log_hr + 1.96 * se_log_hr),
      reference_exposure = exposure_q[[2]],
      exposure_p05 = exposure_q[[1]],
      exposure_p95 = exposure_q[[3]]
    )
}

fit_one <- function(cause_code, outcome_label, pollutant_term, adjustment_term, raw_col, pollutant_label) {
  covars <- c(pollutant_term, adjustment_term, "arf_subtype", base_covars)
  model_df <- analysis_df %>%
    select(ftime_days, event_code, all_of(covars), all_of(raw_col)) %>%
    drop_na() %>%
    mutate(
      arf_subtype = droplevels(arf_subtype),
      sex = droplevels(sex),
      race_ethnicity = droplevels(race_ethnicity),
      index_year_f = droplevels(index_year_f)
    )

  if (!nrow(model_df) || sum(model_df$event_code == cause_code) < 10) return(tibble())

  fml <- as.formula(paste0(
    "Surv(ftime_days, event_code == ", cause_code, ") ~ ",
    pollutant_term,
    " * arf_subtype + ",
    paste(c(adjustment_term, base_covars), collapse = " + ")
  ))
  fit <- survival::coxph(fml, data = model_df, ties = "efron")
  make_prediction_rows(model_df, fit, pollutant_term, adjustment_term, raw_col, outcome_label, pollutant_label) %>%
    mutate(
      n = nrow(model_df),
      events = sum(model_df$event_code == cause_code)
    )
}

prediction_df <- tidyr::crossing(outcome_specs, pollutant_specs) %>%
  rowwise() %>%
  do(fit_one(
    .$cause_code,
    .$outcome,
    .$term,
    .$adjustment_term,
    .$raw_col,
    .$pollutant
  )) %>%
  ungroup() %>%
  mutate(
    outcome = factor(outcome, levels = c("Extubation", "Death", "Persistent RF")),
    pollutant = factor(
      pollutant,
      levels = c("Nitrogen dioxide (NO2)", "Fine particulate matter (PM2.5)")
    ),
    pollutant_label = factor(
      case_when(
        pollutant == "Nitrogen dioxide (NO2)" ~ "Nitrogen~dioxide~(NO[2])",
        pollutant == "Fine particulate matter (PM2.5)" ~ "Fine~particulate~matter~(PM[2.5])",
        TRUE ~ as.character(pollutant)
      ),
      levels = c("Nitrogen~dioxide~(NO[2])", "Fine~particulate~matter~(PM[2.5])")
    ),
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed"))
  )

readr::write_csv(prediction_df, file.path(out_root, "resubmission_subtype_hazard_spectrum_predictions.csv"))

p <- ggplot(prediction_df, aes(x = exposure_raw, y = hazard_ratio, color = arf_subtype, fill = arf_subtype)) +
  geom_hline(yintercept = 1, linewidth = 0.4, color = "grey35") +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.14, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.1) +
  facet_grid(
    outcome ~ pollutant_label,
    scales = "free_x",
    labeller = labeller(pollutant_label = label_parsed)
  ) +
  scale_color_manual(values = subtype_colors, drop = FALSE) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = label_number(accuracy = 0.01, trim = TRUE)
  ) +
  labs(
    title = "Subtype-Specific Cause-Specific Hazards Across Air Pollution Exposure",
    subtitle = "Curves are adjusted hazard ratios relative to the median exposure, with the alternate pollutant and covariates held constant",
    x = expression("Exposure concentration (NO"[2]*" ppb; PM"[2.5]*" "*mu*"g/m"^3*")"),
    y = "Relative hazard, log scale",
    color = "ARF subtype",
    fill = "ARF subtype"
  ) +
  theme_spectrum()

png_path <- file.path(fig_dir, "resubmission_subtype_hazard_spectrum.png")
pdf_path <- file.path(fig_dir, "resubmission_subtype_hazard_spectrum.pdf")
ggsave(png_path, p, width = 16, height = 12, dpi = 300)
ggsave(pdf_path, p, width = 16, height = 12)

message("Wrote subtype hazard spectrum figure:")
message(" - ", png_path)
message(" - ", pdf_path)
