#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(splines)
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

spline_df <- 3L

theme_response <- function(base_size = 18) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      plot.title = element_text(face = "bold", size = rel(1.15), hjust = 0.02, margin = margin(b = 7)),
      plot.subtitle = element_text(color = "grey25", hjust = 0.02, margin = margin(b = 10)),
      strip.text = element_text(face = "bold", color = "grey10", margin = margin(5, 5, 5, 5)),
      strip.background = element_rect(fill = "grey94", color = "grey68", linewidth = 0.55),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.98)),
      legend.key.width = unit(32, "pt"),
      panel.spacing.x = unit(14, "pt"),
      panel.spacing.y = unit(14, "pt"),
      plot.margin = margin(8, 16, 8, 8)
    )
}

safe_mean <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

top_level <- function(x) {
  x <- droplevels(factor(x))
  names(sort(table(x), decreasing = TRUE))[[1]]
}

analysis_df <- read_csv(
  file.path(out_root, "resubmission_analysis_dataset.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed")),
    sex = factor(sex),
    race_ethnicity = factor(race_ethnicity),
    index_year_f = factor(index_year_f),
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

pollutant_specs <- tibble::tribble(
  ~pollutant, ~term, ~adjustment_term, ~raw_col, ~display_label,
  "Nitrogen dioxide (NO2)", "no2_per_10", "pm25_per_5", "no2_12m_zcta", "Nitrogen~dioxide~(NO[2])",
  "Fine particulate matter (PM2.5)", "pm25_per_5", "no2_per_10", "pm25_12m_zcta", "Fine~particulate~matter~(PM[2.5])"
)

outcome_specs <- tibble::tribble(
  ~outcome, ~endpoint, ~event_col, ~time_col,
  "In-hospital mortality", "mortality", "mortality_event", "mortality_ftime_days",
  "Ventilator-free days", "vfd", NA_character_, NA_character_
)

make_prediction_grid <- function(data, pollutant_term, adjustment_term, raw_col) {
  exposure_q <- quantile(data[[raw_col]], probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
  exposure_grid <- seq(exposure_q[[1]], exposure_q[[3]], length.out = 140)
  subtype_levels <- levels(droplevels(data$arf_subtype))
  
  pred_grid <- tidyr::expand_grid(
    exposure_raw = exposure_grid,
    arf_subtype = subtype_levels
  )
  
  pred_grid[[pollutant_term]] <- if (pollutant_term == "pm25_per_5") {
    pred_grid$exposure_raw / 5
  } else {
    pred_grid$exposure_raw / 10
  }
  
  pred_grid %>%
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
      index_year_f = factor(index_year_f, levels = levels(data$index_year_f)),
      reference_exposure = exposure_q[[2]],
      exposure_p05 = exposure_q[[1]],
      exposure_p95 = exposure_q[[3]]
    )
}

reference_grid <- function(pred_grid, pollutant_term) {
  ref <- pred_grid %>% mutate(exposure_raw = reference_exposure)
  ref[[pollutant_term]] <- if (pollutant_term == "pm25_per_5") {
    ref$reference_exposure / 5
  } else {
    ref$reference_exposure / 10
  }
  ref
}

predict_ratio <- function(fit, pred_grid, ref_grid) {
  x_pred <- model.matrix(delete.response(terms(fit)), pred_grid)
  x_ref <- model.matrix(delete.response(terms(fit)), ref_grid)
  beta <- coef(fit)
  vc <- vcov(fit)
  keep <- names(beta)
  x_diff <- x_pred[, keep, drop = FALSE] - x_ref[, keep, drop = FALSE]
  eta <- as.vector(x_diff %*% beta)
  se_eta <- sqrt(rowSums((x_diff %*% vc) * x_diff))
  list(
    ratio = exp(eta),
    conf_low = exp(eta - 1.96 * se_eta),
    conf_high = exp(eta + 1.96 * se_eta)
  )
}

fit_one_curve <- function(outcome, endpoint, event_col, time_col, pollutant, term, adjustment_term, raw_col, display_label) {
  outcome <- outcome[[1]]
  endpoint <- endpoint[[1]]
  event_col <- event_col[[1]]
  time_col <- time_col[[1]]
  pollutant <- pollutant[[1]]
  term <- term[[1]]
  adjustment_term <- adjustment_term[[1]]
  raw_col <- raw_col[[1]]
  display_label <- display_label[[1]]
  
  outcome_cols <- if (endpoint == "mortality") c(event_col, time_col) else character()
  covars <- c(term, adjustment_term, "arf_subtype", base_covars, raw_col)
  model_df <- analysis_df %>%
    select(all_of(covars), all_of(outcome_cols), any_of("ventilator_free_days")) %>%
    drop_na() %>%
    mutate(
      arf_subtype = droplevels(arf_subtype),
      sex = droplevels(sex),
      race_ethnicity = droplevels(race_ethnicity),
      index_year_f = droplevels(index_year_f)
    )
  
  if (!nrow(model_df)) return(tibble())
  
  rhs <- paste0(
    "splines::ns(", term, ", df = ", spline_df, ") * arf_subtype + ",
    paste(c(adjustment_term, base_covars), collapse = " + ")
  )
  fml <- if (endpoint == "mortality") {
    as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", rhs))
  } else {
    as.formula(paste0("ventilator_free_days ~ ", rhs))
  }
  
  fit <- if (endpoint == "mortality") {
    survival::coxph(fml, data = model_df, ties = "efron")
  } else {
    stats::glm(fml, data = model_df, family = quasipoisson(link = "log"))
  }
  
  pred_grid <- make_prediction_grid(model_df, term, adjustment_term, raw_col)
  ref_grid <- reference_grid(pred_grid, term)
  pred <- predict_ratio(fit, pred_grid, ref_grid)
  event_n <- if (endpoint == "mortality") sum(model_df[[event_col]] == 1, na.rm = TRUE) else NA_integer_
  
  pred_grid %>%
    transmute(
      outcome,
      endpoint,
      pollutant,
      pollutant_label = display_label,
      exposure_raw,
      arf_subtype,
      estimate = pred$ratio,
      conf_low = pred$conf_low,
      conf_high = pred$conf_high,
      reference_exposure,
      exposure_p05,
      exposure_p95,
      n = nrow(model_df),
      events = event_n
    )
}

prediction_df <- tidyr::crossing(outcome_specs, pollutant_specs) %>%
  rowwise() %>%
  do(fit_one_curve(
    .$outcome,
    .$endpoint,
    .$event_col,
    .$time_col,
    .$pollutant,
    .$term,
    .$adjustment_term,
    .$raw_col,
    .$display_label
  )) %>%
  ungroup() %>%
  mutate(
    outcome = factor(outcome, levels = c("In-hospital mortality", "Ventilator-free days")),
    pollutant = factor(pollutant, levels = pollutant_specs$pollutant),
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$display_label),
    arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed"))
  )

readr::write_csv(
  prediction_df,
  file.path(out_root, "resubmission_primary_exposure_response_by_arf_subtype_predictions.csv")
)

panel_rug_anchor <- prediction_df %>%
  group_by(outcome, pollutant, pollutant_label) %>%
  summarise(
    rug_anchor = min(conf_low, estimate, na.rm = TRUE),
    .groups = "drop"
  )

rug_df <- bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
  spec <- pollutant_specs[i, ]
  exposure_range <- prediction_df %>%
    filter(pollutant == spec$pollutant) %>%
    summarise(
      exposure_p05 = first(exposure_p05),
      exposure_p95 = first(exposure_p95),
      .groups = "drop"
    )
  
  analysis_df %>%
    transmute(
      pollutant = spec$pollutant,
      pollutant_label = spec$display_label,
      exposure_raw = .data[[spec$raw_col]],
      arf_subtype
    ) %>%
    filter(
      !is.na(exposure_raw),
      !is.na(arf_subtype),
      exposure_raw >= exposure_range$exposure_p05,
      exposure_raw <= exposure_range$exposure_p95
    )
})) %>%
  tidyr::crossing(outcome = levels(prediction_df$outcome)) %>%
  left_join(panel_rug_anchor, by = c("outcome", "pollutant", "pollutant_label")) %>%
  mutate(
    outcome = factor(outcome, levels = levels(prediction_df$outcome)),
    pollutant = factor(pollutant, levels = levels(prediction_df$pollutant)),
    pollutant_label = factor(pollutant_label, levels = levels(prediction_df$pollutant_label)),
    arf_subtype = factor(arf_subtype, levels = levels(prediction_df$arf_subtype)),
    rug_row = as.integer(arf_subtype),
    rug_y = if_else(
      outcome == "Ventilator-free days",
      c(0.72, 0.76, 0.80)[rug_row],
      c(0.58, 0.62, 0.66)[rug_row]
    ),
    rug_yend = rug_y + 0.012
  )

p <- ggplot(prediction_df, aes(x = exposure_raw, y = estimate, color = arf_subtype, fill = arf_subtype)) +
  geom_hline(yintercept = 1, linewidth = 0.4, color = "grey35") +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.13, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.1) +
  geom_segment(
    data = rug_df,
    aes(x = exposure_raw, xend = exposure_raw, y = rug_y, yend = rug_yend, color = arf_subtype),
    inherit.aes = FALSE,
    alpha = 0.18,
    linewidth = 0.18,
    show.legend = FALSE
  ) +
  facet_grid(
    outcome ~ pollutant_label,
    scales = "free",
    labeller = labeller(pollutant_label = label_parsed)
  ) +
  scale_color_manual(values = subtype_colors, drop = FALSE) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_y_log10(
    breaks = c(0.6, 0.75, 0.9, 1, 1.1, 1.25, 1.5, 2),
    labels = label_number(accuracy = 0.01, trim = TRUE),
    limits = function(x) {
      c(min(x, na.rm = TRUE) * 0.97, max(x, na.rm = TRUE) * 1.03)
    }
  ) +
  labs(
    title = "Exposure-Response by ARF Subtype",
    x = expression("Exposure concentration (NO"[2]*" ppb; PM"[2.5]*" "*mu*"g/m"^3*")"),
    y = "Ratio relative to median exposure, log scale",
    color = "ARF subtype",
    fill = "ARF subtype"
  ) +
  theme_response()

png_path <- file.path(fig_dir, "resubmission_primary_exposure_response_by_arf_subtype.png")
pdf_path <- file.path(fig_dir, "resubmission_primary_exposure_response_by_arf_subtype.pdf")
ggsave(png_path, p, width = 17, height = 10.5, dpi = 300)
ggsave(pdf_path, p, width = 17, height = 10.5)

message("Wrote primary exposure-response by subtype figure:")
message(" - ", png_path)
message(" - ", pdf_path)
