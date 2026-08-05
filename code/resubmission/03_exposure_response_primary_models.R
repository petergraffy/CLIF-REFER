#!/usr/bin/env Rscript

# Exposure-response curves for the reviewer-optimized primary models.
#
# Curves show model-based predictions on the outcome scale, matching the
# original REFER exposure-response style:
#   - Logistic mortality model: predicted mortality probability by day 28
#   - Quasi-Poisson VFD model: predicted mean ventilator-free days
#
# The script produces whole-cohort, ARF-subtype, sex, and race/ethnicity
# versions for PM2.5 and NO2.

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(glue)
  library(janitor)
  library(patchwork)
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
  stop("Usage: Rscript code/resubmission/03_exposure_response_primary_models.R <analysis_dataset.csv> [output_dir]")
}
arg_or <- function(i, default = NA_character_) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
repo <- normalizePath(file.path(dirname(input_path), "..", ".."), mustWork = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- arg_or(2, file.path(repo, "output", "resubmission", run_id))
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
    mortality_day28_event = as.integer(
      if ("mortality_day28_event" %in% names(.)) {
        mortality_day28_event
      } else {
        mortality_event == 1L & mortality_ftime_days <= mortality_prediction_day
      }
    ),
    ventilator_free_days = as.numeric(ventilator_free_days),
    age_10 = age_10 %||% (age / 10),
    pm25_per_5 = pm25_per_5 %||% (pm25_12m_zcta / 5),
    no2_per_10 = no2_per_10 %||% (no2_12m_zcta / 10),
    acs_median_household_income_10k = acs_median_household_income_10k %||% (acs_median_household_income / 10000),
    sex = factor(sex),
    race_ethnicity = forcats::fct_collapse(
      factor(race_ethnicity),
      "Other/Unknown" = c("Hispanic Other/Unknown", "Other/Unknown")
    ),
    race_ethnicity = forcats::fct_lump_min(race_ethnicity, min = 100, other_level = "Other/Unknown"),
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
  glue("Mortality by day {mortality_prediction_day}"), "mortality", glue("Predicted mortality probability by day {mortality_prediction_day}"),
  "Ventilator-free days", "vfd", "Predicted ventilator-free days"
)

outcome_colors <- setNames(
  c("#2166AC", "#B2182B"),
  c(glue("Predicted mortality probability by day {mortality_prediction_day}"), "Predicted ventilator-free days")
)

pollutant_x_labels <- list(
  "NO2" = expression(NO[2]~"(ppb)"),
  "PM2.5" = expression(PM[2.5]~"("*mu*"g/m"^3*")")
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
    pr <- predict(fit, newdata = pred_grid, type = "link", se.fit = TRUE)
    return(list(
      estimate = plogis(pr$fit),
      conf_low = plogis(pr$fit - 1.96 * pr$se.fit),
      conf_high = plogis(pr$fit + 1.96 * pr$se.fit)
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
    dplyr::select(mortality_day28_event, ventilator_free_days, all_of(covars)) %>%
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
    as.formula(paste0("mortality_day28_event ~ ", rhs))
  } else {
    as.formula(paste0("ventilator_free_days ~ ", rhs))
  }

  fit <- if (endpoint == "mortality") {
    stats::glm(fml, data = model_df, family = binomial(link = "logit"))
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
      events = if_else(endpoint == "mortality", sum(model_df$mortality_day28_event == 1L), NA_integer_),
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
  left_join(outcome_specs %>% select(outcome, y_label), by = "outcome") %>%
  mutate(
    outcome = factor(outcome, levels = outcome_specs$outcome),
    y_label = factor(y_label, levels = outcome_specs$y_label),
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
  left_join(outcome_specs %>% select(outcome, y_label), by = "outcome") %>%
  mutate(
    outcome = factor(outcome, levels = levels(prediction_df$outcome)),
    y_label = factor(y_label, levels = levels(prediction_df$y_label)),
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
subtype_colors <- curve_colors[names(curve_colors) != "Overall"]

theme_response <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(face = "bold", size = 8.5, margin = margin(t = 8, r = 8, b = 8, l = 8)),
      strip.placement = "outside",
      legend.position = "bottom",
      plot.margin = margin(5.5, 5.5, 5.5, 14)
    )
}

wrap_pollutant_panels <- function(plot_fun, heights = NULL) {
  plots <- lapply(levels(prediction_df$pollutant), function(pollutant_value) plot_fun(pollutant_value))
  patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom")
}

make_separate_rug_panel <- function(rugs, group_col, palette_values, pollutant_value) {
  group_col <- rlang::ensym(group_col)
  rugs <- rugs %>%
    mutate(rug_level = factor(as.character(!!group_col), levels = names(palette_values))) %>%
    filter(!is.na(rug_level))

  ggplot(rugs, aes(x = exposure_raw, y = rug_level, color = rug_level)) +
    geom_point(shape = "|", size = 1.8, alpha = 0.32, show.legend = FALSE) +
    scale_color_manual(values = palette_values, drop = TRUE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = pollutant_x_labels[[as.character(pollutant_value)]], y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10, color = "grey20"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.25),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 13),
      plot.margin = margin(0, 5.5, 5.5, 5.5),
      legend.position = "none"
    )
}

apply_outcome_y_windows <- function(df) {
  df %>%
    mutate(
      y_window_low = if_else(endpoint == "mortality", 0, 15),
      y_window_high = if_else(endpoint == "mortality", 0.40, 28),
      estimate_plot = pmin(pmax(estimate, y_window_low), y_window_high),
      conf_low_plot = pmin(pmax(conf_low, y_window_low), y_window_high),
      conf_high_plot = pmin(pmax(conf_high, y_window_low), y_window_high)
    )
}

make_y_window_df <- function(plot_df) {
  plot_df %>%
    distinct(y_label, endpoint) %>%
    mutate(
      exposure_raw = min(plot_df$exposure_raw, na.rm = TRUE),
      y_window_low = if_else(endpoint == "mortality", 0, 15),
      y_window_high = if_else(endpoint == "mortality", 0.40, 28)
    ) %>%
    pivot_longer(c(y_window_low, y_window_high), names_to = "window_bound", values_to = "y_value")
}

label_outcome_axis <- function(x) {
  if_else(
    abs(x) >= 1,
    label_number(accuracy = 1, trim = TRUE)(x),
    label_number(accuracy = 0.01, trim = TRUE)(x)
  )
}

plot_curves <- function(curve_type_value, filename_stub) {
  df <- prediction_df %>%
    filter(curve_type == curve_type_value) %>%
    apply_outcome_y_windows()
  rugs <- rug_df %>% filter(curve_type == curve_type_value)
  is_overall <- identical(curve_type_value, "Whole cohort")
  palette_values <- if (is_overall) outcome_colors else subtype_colors
  df <- df %>%
    mutate(curve_group = if (is_overall) as.character(y_label) else as.character(arf_subtype))
  rugs <- rugs %>%
    mutate(rug_group = if (is_overall) as.character(y_label) else as.character(rug_group))

  p <- wrap_pollutant_panels(function(pollutant_value) {
    plot_df <- df %>% filter(pollutant == pollutant_value)
    plot_rugs <- rugs %>% filter(pollutant == pollutant_value)

    curve_plot <- ggplot(plot_df, aes(exposure_raw, estimate_plot, color = curve_group, fill = curve_group)) +
      geom_blank(data = make_y_window_df(plot_df), aes(x = exposure_raw, y = y_value), inherit.aes = FALSE) +
      geom_ribbon(aes(ymin = conf_low_plot, ymax = conf_high_plot), alpha = 0.14, color = NA, show.legend = FALSE) +
      geom_line(linewidth = 1.05) +
      facet_grid(y_label ~ ., scales = "free_y", switch = "y") +
      scale_color_manual(values = palette_values, drop = TRUE) +
      scale_fill_manual(values = palette_values, drop = TRUE) +
      scale_y_continuous(labels = label_outcome_axis, expand = expansion(mult = c(0, 0))) +
      labs(
        title = rlang::parse_expr(as.character(plot_df$pollutant_label[[1]])),
        x = if (is_overall) pollutant_x_labels[[as.character(pollutant_value)]] else NULL,
        y = NULL,
        color = NULL,
        fill = NULL
      ) +
      theme_response()

    if (is_overall) {
      curve_plot +
        geom_rug(
          data = plot_rugs,
          aes(x = exposure_raw, color = factor(rug_group, levels = names(palette_values))),
          inherit.aes = FALSE,
          sides = "b",
          alpha = 0.16,
          linewidth = 0.24,
          show.legend = FALSE
        )
    } else {
      curve_plot / make_separate_rug_panel(plot_rugs, rug_group, palette_values, pollutant_value) +
        plot_layout(heights = c(1, 0.18))
    }
  })

  ggsave(file.path(fig_dir, paste0(filename_stub, ".png")), p, width = 12.5, height = 8.5, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename_stub, ".pdf")), p, width = 12.5, height = 8.5)
  p
}

p_overall <- plot_curves("Whole cohort", "primary_exposure_response_predicted_outcome_whole_cohort")
p_subtype <- plot_curves("By ARF subtype", "primary_exposure_response_predicted_outcome_by_arf_subtype")

plot_subtype_mortality_only <- function() {
  df <- prediction_df %>%
    filter(curve_type == "By ARF subtype", endpoint == "mortality") %>%
    apply_outcome_y_windows()
  rugs <- rug_df %>%
    filter(curve_type == "By ARF subtype", outcome == glue("Mortality by day {mortality_prediction_day}"))

  p <- wrap_pollutant_panels(function(pollutant_value) {
    plot_df <- df %>% filter(pollutant == pollutant_value)
    plot_rugs <- rugs %>% filter(pollutant == pollutant_value)

    curve_plot <- ggplot(plot_df, aes(exposure_raw, estimate_plot, color = arf_subtype, fill = arf_subtype)) +
      geom_ribbon(aes(ymin = conf_low_plot, ymax = conf_high_plot), alpha = 0.14, color = NA, show.legend = FALSE) +
      geom_line(linewidth = 1.05) +
      scale_color_manual(values = subtype_colors, drop = TRUE) +
      scale_fill_manual(values = subtype_colors, drop = TRUE) +
      scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 0.40), expand = expansion(mult = c(0, 0))) +
      labs(
        title = rlang::parse_expr(as.character(plot_df$pollutant_label[[1]])),
        x = NULL,
        y = glue("Predicted mortality probability by day {mortality_prediction_day}"),
        color = NULL,
        fill = NULL
      ) +
      theme_response()

    curve_plot / make_separate_rug_panel(plot_rugs, rug_group, subtype_colors, pollutant_value) +
      plot_layout(heights = c(1, 0.22))
  })

  ggsave(
    file.path(fig_dir, "primary_mortality_day28_logistic_exposure_response_by_arf_subtype.png"),
    p,
    width = 12.5,
    height = 5.5,
    dpi = 300
  )
  ggsave(
    file.path(fig_dir, "primary_mortality_day28_logistic_exposure_response_by_arf_subtype.pdf"),
    p,
    width = 12.5,
    height = 5.5
  )
  p
}

p_subtype_mortality <- plot_subtype_mortality_only()

make_group_prediction_grid <- function(data, term, raw_col, other_term, group_var) {
  exposure_q <- quantile(data[[raw_col]], probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
  exposure_grid <- seq(exposure_q[[1]], exposure_q[[3]], length.out = 160)
  group_levels <- levels(droplevels(data[[group_var]]))

  grid <- tidyr::expand_grid(
    exposure_raw = exposure_grid,
    group_level = group_levels
  )
  grid[[group_var]] <- factor(grid$group_level, levels = levels(droplevels(data[[group_var]])))
  grid[[term]] <- if (term == "pm25_per_5") grid$exposure_raw / 5 else grid$exposure_raw / 10

  if (other_term %in% names(data)) {
    grid[[other_term]] <- median(data[[other_term]], na.rm = TRUE)
  }

  grid %>%
    mutate(
      age_10 = safe_mean(data$age_10),
      sex = if (group_var == "sex") sex else top_level(data$sex),
      race_ethnicity = if (group_var == "race_ethnicity") race_ethnicity else top_level(data$race_ethnicity),
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
      across(any_of(c("sex", "race_ethnicity", "index_year_f")), factor)
    ) %>%
    align_factor_levels(data)
}

fit_mortality_group_curve <- function(group_var, group_label, pollutant, term, raw_col, other_term, display_label) {
  covars <- c(term, other_term, base_covars, raw_col, group_var)

  model_df <- analysis_df %>%
    dplyr::select(mortality_day28_event, all_of(covars)) %>%
    tidyr::drop_na() %>%
    mutate(
      sex = droplevels(sex),
      race_ethnicity = droplevels(race_ethnicity),
      index_year_f = droplevels(index_year_f)
    )

  exposure_terms <- c(term, other_term)
  rhs_covars <- setdiff(drop_uninformative_covars(model_df, c(other_term, base_covars), exposure_terms), group_var)
  spline_term <- paste0("splines::ns(", term, ", df = ", spline_df, ") * ", group_var)
  rhs <- paste(c(spline_term, rhs_covars), collapse = " + ")
  fml <- as.formula(paste0("mortality_day28_event ~ ", rhs))
  fit <- stats::glm(fml, data = model_df, family = binomial(link = "logit"))

  pred_grid <- make_group_prediction_grid(model_df, term, raw_col, other_term, group_var)
  pred <- predict_outcome_scale(fit, pred_grid, "mortality")

  pred_grid %>%
    mutate(
      site = site_name,
      curve_type = group_label,
      outcome = glue("Mortality by day {mortality_prediction_day}"),
      endpoint = "mortality",
      pollutant = pollutant,
      pollutant_label = display_label,
      group_var = group_var,
      group_label = group_label,
      estimate = pred$estimate,
      conf_low = pred$conf_low,
      conf_high = pred$conf_high,
      n = nrow(model_df),
      events = sum(model_df$mortality_day28_event == 1L),
      adjustment_set = paste(rhs_covars, collapse = " + "),
      mortality_prediction_day = mortality_prediction_day
    )
}

fit_vfd_group_curve <- function(group_var, group_label, pollutant, term, raw_col, other_term, display_label) {
  covars <- c(term, other_term, base_covars, raw_col, group_var)

  model_df <- analysis_df %>%
    dplyr::select(ventilator_free_days, all_of(covars)) %>%
    tidyr::drop_na() %>%
    mutate(
      sex = droplevels(sex),
      race_ethnicity = droplevels(race_ethnicity),
      index_year_f = droplevels(index_year_f)
    )

  exposure_terms <- c(term, other_term)
  rhs_covars <- setdiff(drop_uninformative_covars(model_df, c(other_term, base_covars), exposure_terms), group_var)
  spline_term <- paste0("splines::ns(", term, ", df = ", spline_df, ") * ", group_var)
  rhs <- paste(c(spline_term, rhs_covars), collapse = " + ")
  fml <- as.formula(paste0("ventilator_free_days ~ ", rhs))
  fit <- stats::glm(fml, data = model_df, family = quasipoisson(link = "log"))

  pred_grid <- make_group_prediction_grid(model_df, term, raw_col, other_term, group_var)
  pred <- predict_outcome_scale(fit, pred_grid, "vfd")

  pred_grid %>%
    mutate(
      site = site_name,
      curve_type = group_label,
      outcome = "Ventilator-free days",
      endpoint = "vfd",
      pollutant = pollutant,
      pollutant_label = display_label,
      group_var = group_var,
      group_label = group_label,
      estimate = pred$estimate,
      conf_low = pred$conf_low,
      conf_high = pred$conf_high,
      n = nrow(model_df),
      events = NA_integer_,
      adjustment_set = paste(rhs_covars, collapse = " + "),
      mortality_prediction_day = NA_real_
    )
}

mortality_group_predictions <- bind_rows(
  tidyr::crossing(group_var = c("sex", "race_ethnicity"), pollutant_specs) %>%
    rowwise() %>%
    do(fit_mortality_group_curve(
      .$group_var,
      if_else(.$group_var == "sex", "By sex", "By race/ethnicity"),
      .$pollutant,
      .$term,
      .$raw_col,
      .$other_term,
      .$display_label
    )) %>%
    ungroup()
) %>%
  mutate(
    curve_type = factor(curve_type, levels = c("By sex", "By race/ethnicity")),
    pollutant = factor(pollutant, levels = pollutant_specs$pollutant),
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$display_label)
  )

vfd_group_predictions <- bind_rows(
  tidyr::crossing(group_var = c("sex", "race_ethnicity"), pollutant_specs) %>%
    rowwise() %>%
    do(fit_vfd_group_curve(
      .$group_var,
      if_else(.$group_var == "sex", "By sex", "By race/ethnicity"),
      .$pollutant,
      .$term,
      .$raw_col,
      .$other_term,
      .$display_label
    )) %>%
    ungroup()
)

group_predictions <- bind_rows(mortality_group_predictions, vfd_group_predictions) %>%
  left_join(outcome_specs %>% select(outcome, y_label), by = "outcome") %>%
  mutate(
    outcome = factor(outcome, levels = outcome_specs$outcome),
    y_label = factor(y_label, levels = outcome_specs$y_label),
    curve_type = factor(curve_type, levels = c("By sex", "By race/ethnicity")),
    pollutant = factor(pollutant, levels = pollutant_specs$pollutant),
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$display_label)
  )

readr::write_csv(
  mortality_group_predictions,
  file.path(out_dir, "primary_mortality_exposure_response_by_sex_race_predictions.csv")
)
readr::write_csv(
  group_predictions,
  file.path(out_dir, "primary_exposure_response_by_sex_race_predictions.csv")
)

make_group_rug_df <- function(group_var) {
  bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
    spec <- pollutant_specs[i, ]
    ranges <- group_predictions %>%
      filter(pollutant == spec$pollutant, group_var == !!group_var) %>%
      summarise(exposure_p05 = first(exposure_p05), exposure_p95 = first(exposure_p95), .groups = "drop")

    analysis_df %>%
      transmute(
        pollutant = spec$pollutant,
        pollutant_label = spec$display_label,
        exposure_raw = .data[[spec$raw_col]],
        group_level = as.character(.data[[group_var]])
      ) %>%
      filter(
        !is.na(exposure_raw),
        !is.na(group_level),
        exposure_raw >= ranges$exposure_p05,
        exposure_raw <= ranges$exposure_p95
      )
  })) %>%
    crossing(outcome = outcome_specs$outcome) %>%
    left_join(outcome_specs %>% select(outcome, y_label), by = "outcome") %>%
    mutate(
      outcome = factor(outcome, levels = outcome_specs$outcome),
      y_label = factor(y_label, levels = outcome_specs$y_label),
      pollutant = factor(pollutant, levels = levels(group_predictions$pollutant)),
      pollutant_label = factor(pollutant_label, levels = levels(group_predictions$pollutant_label))
    )
}

plot_group_combined_curves <- function(group_var, group_label, filename_stub, palette_values) {
  df <- group_predictions %>%
    filter(.data$group_var == !!group_var) %>%
    apply_outcome_y_windows() %>%
    mutate(group_level = factor(group_level, levels = names(palette_values)))
  rugs <- make_group_rug_df(group_var) %>%
    mutate(group_level = factor(group_level, levels = names(palette_values)))
  ribbon_alpha <- if (length(palette_values) > 3) 0.06 else 0.14

  p <- wrap_pollutant_panels(function(pollutant_value) {
    plot_df <- df %>% filter(pollutant == pollutant_value)
    plot_rugs <- rugs %>% filter(pollutant == pollutant_value)

    curve_plot <- ggplot(plot_df, aes(exposure_raw, estimate_plot, color = group_level, fill = group_level)) +
      geom_blank(data = make_y_window_df(plot_df), aes(x = exposure_raw, y = y_value), inherit.aes = FALSE) +
      geom_ribbon(aes(ymin = conf_low_plot, ymax = conf_high_plot), alpha = ribbon_alpha, color = NA, show.legend = FALSE) +
      geom_line(linewidth = 1.05) +
      facet_grid(y_label ~ ., scales = "free_y", switch = "y") +
      scale_color_manual(values = palette_values, drop = TRUE) +
      scale_fill_manual(values = palette_values, drop = TRUE) +
      scale_y_continuous(labels = label_outcome_axis, expand = expansion(mult = c(0, 0))) +
      labs(
        title = rlang::parse_expr(as.character(plot_df$pollutant_label[[1]])),
        x = NULL,
        y = NULL,
        color = NULL,
        fill = NULL
      ) +
      theme_response()

    curve_plot / make_separate_rug_panel(plot_rugs, group_level, palette_values, pollutant_value) +
      plot_layout(heights = c(1, 0.18))
  })

  plot_width <- if (identical(group_var, "race_ethnicity")) 16 else 12.5
  ggsave(file.path(fig_dir, paste0(filename_stub, ".png")), p, width = plot_width, height = 8.5, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename_stub, ".pdf")), p, width = plot_width, height = 8.5)
  p
}

plot_mortality_group_curves <- function(group_var, group_label, filename_stub, palette_values) {
  df <- mortality_group_predictions %>%
    filter(.data$group_var == !!group_var) %>%
    apply_outcome_y_windows() %>%
    mutate(group_level = factor(group_level, levels = names(palette_values)))
  rugs <- make_group_rug_df(group_var) %>%
    mutate(group_level = factor(group_level, levels = names(palette_values)))
  ribbon_alpha <- if (length(palette_values) > 3) 0.06 else 0.14

  p <- wrap_pollutant_panels(function(pollutant_value) {
    plot_df <- df %>% filter(pollutant == pollutant_value)
    plot_rugs <- rugs %>% filter(pollutant == pollutant_value)

    curve_plot <- ggplot(plot_df, aes(exposure_raw, estimate_plot, color = group_level, fill = group_level)) +
      geom_ribbon(aes(ymin = conf_low_plot, ymax = conf_high_plot), alpha = ribbon_alpha, color = NA, show.legend = FALSE) +
      geom_line(linewidth = 1.05) +
      scale_color_manual(values = palette_values, drop = TRUE) +
      scale_fill_manual(values = palette_values, drop = TRUE) +
      scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 0.40), expand = expansion(mult = c(0, 0))) +
      labs(
        title = rlang::parse_expr(as.character(plot_df$pollutant_label[[1]])),
        x = NULL,
        y = glue("Predicted mortality probability by day {mortality_prediction_day}"),
        color = NULL,
        fill = NULL
      ) +
      theme_response()

    curve_plot / make_separate_rug_panel(plot_rugs, group_level, palette_values, pollutant_value) +
      plot_layout(heights = c(1, 0.22))
  })

  plot_width <- if (identical(group_var, "race_ethnicity")) 16 else 12.5
  ggsave(file.path(fig_dir, paste0(filename_stub, ".png")), p, width = plot_width, height = 5.5, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename_stub, ".pdf")), p, width = plot_width, height = 5.5)
  p
}

sex_palette <- c("Female" = "#009E73", "Male" = "#D55E00")
race_levels <- levels(droplevels(analysis_df$race_ethnicity))
race_palette <- setNames(
  c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")[seq_along(race_levels)],
  race_levels
)

p_sex <- plot_mortality_group_curves("sex", "By sex", "primary_mortality_exposure_response_predicted_by_sex", sex_palette)
p_race <- plot_mortality_group_curves("race_ethnicity", "By race/ethnicity", "primary_mortality_exposure_response_predicted_by_race_ethnicity", race_palette)
p_sex_combined <- plot_group_combined_curves("sex", "By sex", "primary_exposure_response_predicted_outcome_by_sex", sex_palette)
p_race_combined <- plot_group_combined_curves("race_ethnicity", "By race/ethnicity", "primary_exposure_response_predicted_outcome_by_race_ethnicity", race_palette)

message("Wrote:")
message(" - ", file.path(out_dir, "primary_exposure_response_predictions.csv"))
message(" - ", file.path(out_dir, "primary_mortality_exposure_response_by_sex_race_predictions.csv"))
message(" - ", file.path(out_dir, "primary_exposure_response_by_sex_race_predictions.csv"))
message(" - ", file.path(fig_dir, "primary_exposure_response_predicted_outcome_whole_cohort.png"))
message(" - ", file.path(fig_dir, "primary_exposure_response_predicted_outcome_by_arf_subtype.png"))
message(" - ", file.path(fig_dir, "primary_exposure_response_predicted_outcome_by_sex.png"))
message(" - ", file.path(fig_dir, "primary_exposure_response_predicted_outcome_by_race_ethnicity.png"))
message(" - ", file.path(fig_dir, "primary_mortality_day28_logistic_exposure_response_by_arf_subtype.png"))
message(" - ", file.path(fig_dir, "primary_mortality_exposure_response_predicted_by_sex.png"))
message(" - ", file.path(fig_dir, "primary_mortality_exposure_response_predicted_by_race_ethnicity.png"))
