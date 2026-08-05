#!/usr/bin/env Rscript

# Unadjusted Aalen-Johansen CIFs and adjusted Fine-Gray models
#
# Input:
#   A reviewer-optimized analysis dataset with:
#     ftime_days, event_code, pm25_12m_zcta/no2_12m_zcta, age_10, sex,
#     race_ethnicity, charlson_score, index_year_f, and ACS covariates.
#
# Event coding among IMV patients:
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
  library(patchwork)
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

competing_risk_diagnostics <- analysis_df %>%
  summarise(
    n_input = n(),
    n_with_imv_after_arf = if ("has_imv_after_arf" %in% names(.)) sum(has_imv_after_arf %in% TRUE, na.rm = TRUE) else NA_integer_,
    n_event_code_0_before_exhaustive_recode = sum(event_code == 0L, na.rm = TRUE),
    n_imv_event_code_0_before_exhaustive_recode = if ("has_imv_after_arf" %in% names(.)) {
      sum(has_imv_after_arf %in% TRUE & event_code == 0L, na.rm = TRUE)
    } else {
      NA_integer_
    }
  )

if ("has_imv_after_arf" %in% names(analysis_df)) {
  analysis_df <- analysis_df %>%
    filter(has_imv_after_arf %in% TRUE) %>%
    mutate(event_name = if ("event_name" %in% names(.)) as.character(event_name) else NA_character_) %>%
    mutate(
      event_code = if_else(event_code == 0L, 3L, event_code),
      event_name = case_when(
        event_code == 1L ~ "extubation",
        event_code == 2L ~ "death",
        event_code == 3L ~ "persistent_rf",
        TRUE ~ event_name
      )
    )
} else {
  warning("Input has no `has_imv_after_arf` column; assuming all rows are eligible for the IMV competing-risk analysis.")
  analysis_df <- analysis_df %>%
    mutate(event_name = if ("event_name" %in% names(.)) as.character(event_name) else NA_character_) %>%
    mutate(
      event_code = if_else(event_code == 0L, 3L, event_code),
      event_name = case_when(
        event_code == 1L ~ "extubation",
        event_code == 2L ~ "death",
        event_code == 3L ~ "persistent_rf",
        TRUE ~ event_name
      )
    )
}

competing_risk_diagnostics <- competing_risk_diagnostics %>%
  mutate(
    n_secondary_competing_risk = nrow(analysis_df),
    n_event_code_0_after_exhaustive_recode = sum(analysis_df$event_code == 0L, na.rm = TRUE),
    n_extubation = sum(analysis_df$event_code == 1L, na.rm = TRUE),
    n_death = sum(analysis_df$event_code == 2L, na.rm = TRUE),
    n_persistent_rf = sum(analysis_df$event_code == 3L, na.rm = TRUE)
  )

readr::write_csv(
  competing_risk_diagnostics,
  file.path(out_dir, "competing_risk_exhaustive_imv_assignment_diagnostics.csv")
)

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
})) %>%
  filter(time <= 28)

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
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 7, 14, 21, 28), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.06))) +
  labs(
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

fit_fgr_object <- function(df, cause_code, exposure_terms, model_label) {
  outcome_label <- unname(event_labels[as.character(cause_code)])
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

  list(
    fit = fit,
    model_df = model_df,
    covars = covars,
    cause_code = cause_code,
    outcome_label = outcome_label,
    model_label = model_label
  )
}

typical_value <- function(x) {
  if (is.factor(x)) {
    tab <- sort(table(x), decreasing = TRUE)
    return(factor(names(tab)[[1]], levels = levels(x)))
  }
  if (is.character(x)) {
    tab <- sort(table(x), decreasing = TRUE)
    return(names(tab)[[1]])
  }
  mean(as.numeric(x), na.rm = TRUE)
}

make_prediction_row <- function(model_df, covars) {
  out <- lapply(covars, function(covar) typical_value(model_df[[covar]]))
  names(out) <- covars
  as_tibble(out)
}

predict_fg_cif <- function(fit_obj, source_df, pollutant, raw_col, exposure_term) {
  time_grid <- seq(0, 28, by = 0.25)
  quartile_df <- source_df %>%
    mutate(exposure_quartile = make_quartile(.data[[raw_col]])) %>%
    filter(!is.na(exposure_quartile), !is.na(.data[[raw_col]])) %>%
    group_by(exposure_quartile) %>%
    summarise(exposure_median = median(.data[[raw_col]], na.rm = TRUE), .groups = "drop")

  pred_rows <- bind_rows(lapply(seq_len(nrow(quartile_df)), function(i) {
    row <- make_prediction_row(fit_obj$model_df, fit_obj$covars)
    row[[exposure_term]] <- if (exposure_term == "pm25_per_5") {
      quartile_df$exposure_median[[i]] / 5
    } else {
      quartile_df$exposure_median[[i]] / 10
    }
    row %>%
      mutate(
        pollutant = pollutant,
        exposure_quartile = quartile_df$exposure_quartile[[i]],
        exposure_median = quartile_df$exposure_median[[i]]
      )
  }))

  sf <- survival::survfit(fit_obj$fit, newdata = pred_rows, se.fit = FALSE)
  surv_summary <- summary(sf, times = time_grid, extend = TRUE)

  if (is.matrix(surv_summary$surv)) {
    strata_id <- rep(seq_len(ncol(surv_summary$surv)), each = length(surv_summary$time))
    time <- rep(surv_summary$time, times = ncol(surv_summary$surv))
    survival_prob <- as.vector(surv_summary$surv)
  } else {
    strata_id <- if (is.null(surv_summary$strata)) {
      rep(1L, length(surv_summary$time))
    } else {
      as.integer(gsub("^.*=", "", surv_summary$strata))
    }
    time <- surv_summary$time
    survival_prob <- surv_summary$surv
  }

  tibble(
    site = site_name,
    pollutant = pollutant,
    exposure_quartile = pred_rows$exposure_quartile[strata_id],
    exposure_median = pred_rows$exposure_median[strata_id],
    outcome = fit_obj$outcome_label,
    event_code = fit_obj$cause_code,
    time = time,
    cif = pmax(0, pmin(1, 1 - survival_prob)),
    model = fit_obj$model_label,
    adjustment_set = paste(setdiff(fit_obj$covars, exposure_term), collapse = " + ")
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

fg_cif_predictions <- bind_rows(lapply(seq_len(nrow(pollutants)), function(i) {
  spec <- pollutants[i, ]
  bind_rows(lapply(1:3, function(cause_code) {
    fit_obj <- fit_fgr_object(analysis_df, cause_code, spec$term, paste0(spec$pollutant, " single-pollutant"))
    predict_fg_cif(fit_obj, analysis_df, spec$pollutant, spec$raw_col, spec$term)
  }))
})) %>%
  mutate(
    pollutant = factor(pollutant, levels = c("NO2", "PM2.5")),
    outcome = recode(
      as.character(outcome),
      "Successful extubation" = "Extubation",
      "Persistent respiratory failure" = "Persistent RF"
    ),
    outcome = factor(outcome, levels = c("Extubation", "Death", "Persistent RF")),
    exposure_quartile = factor(exposure_quartile, levels = paste0("Q", 1:4))
  )

readr::write_csv(
  fg_cif_predictions,
  file.path(out_dir, "fine_gray_adjusted_cif_predictions_by_exposure_quartile.csv")
)

fg_quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
fg_quartile_colors <- setNames(quartile_colors, fg_quartile_labels)

fg_plot_data <- fg_cif_predictions %>%
  mutate(
    exposure_quartile = recode(
      as.character(exposure_quartile),
      "Q1" = "Q1 lowest",
      "Q2" = "Q2",
      "Q3" = "Q3",
      "Q4" = "Q4 highest"
    ),
    exposure_quartile = factor(exposure_quartile, levels = fg_quartile_labels),
    outcome = factor(outcome, levels = c("Extubation", "Death", "Persistent RF"))
  )

fg_panel_scales <- fg_plot_data %>%
  group_by(pollutant, outcome) %>%
  summarise(y_max = max(cif, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_limit = y_max * 1.04)

theme_fg_cif <- function() {
  theme_minimal(base_size = 18) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      legend.key.width = unit(30, "pt"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16, color = "grey20"),
      panel.spacing.x = unit(28, "pt"),
      panel.spacing.y = unit(22, "pt"),
      plot.margin = margin(4, 8, 4, 8)
    )
}

make_fg_curve_panel <- function(pollutant_value, outcome_value, show_y_axis = TRUE, show_x_axis = FALSE) {
  dat <- fg_plot_data %>% filter(pollutant == pollutant_value, outcome == outcome_value)
  y_limit <- fg_panel_scales %>%
    filter(pollutant == pollutant_value, outcome == outcome_value) %>%
    pull(y_limit)

  ggplot(dat, aes(time, cif, color = exposure_quartile)) +
    geom_line(linewidth = 1.05) +
    scale_color_manual(values = fg_quartile_colors, drop = FALSE) +
    scale_x_continuous(
      limits = c(0, 28),
      breaks = c(0, 7, 14, 21, 28),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, y_limit),
      expand = expansion(mult = c(0.015, 0.035))
    ) +
    labs(
      x = if (show_x_axis) "Days since ARF onset" else NULL,
      y = if (show_y_axis) "Adjusted cumulative incidence" else NULL,
      color = "Quartile"
    ) +
    theme_fg_cif() +
    theme(
      axis.title.y = if (show_y_axis) element_text(size = 19) else element_blank(),
      axis.text.y = if (show_y_axis) element_text(size = 17, color = "grey20") else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line(color = "black", linewidth = 0.35) else element_blank(),
      axis.title.x = if (show_x_axis) element_text(size = 18, margin = margin(t = 7)) else element_blank(),
      legend.position = "bottom"
    )
}

make_fg_column_header <- function(label_expr) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "grey45", linewidth = 0.55) +
    annotate("text", x = 0.5, y = 0.42, label = label_expr, parse = TRUE, hjust = 0.5, vjust = 0.5, fontface = "bold", size = 6.4) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(12, 8, -18, 8))
}

make_fg_row_header <- function(label_text) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "grey45", linewidth = 0.55) +
    annotate("text", x = 0.5, y = 0.5, label = label_text, angle = -90, hjust = 0.5, vjust = 0.5, fontface = "bold", size = 6.8) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(4, 8, 2, -2))
}

fg_header_row <- wrap_plots(
  make_fg_column_header("'Nitrogen dioxide'~(NO[2])"),
  make_fg_column_header("'Fine particulate matter'~(PM[2.5])"),
  plot_spacer(),
  nrow = 1,
  widths = c(1, 1, 0.075)
)

fg_outcome_levels <- c("Extubation", "Death", "Persistent RF")
fg_pollutant_levels <- c("NO2", "PM2.5")

fg_row_plots <- lapply(seq_along(fg_outcome_levels), function(i) {
  outcome_value <- fg_outcome_levels[[i]]
  wrap_plots(
    make_fg_curve_panel(fg_pollutant_levels[[1]], outcome_value, show_y_axis = TRUE, show_x_axis = i == length(fg_outcome_levels)),
    make_fg_curve_panel(fg_pollutant_levels[[2]], outcome_value, show_y_axis = FALSE, show_x_axis = i == length(fg_outcome_levels)),
    make_fg_row_header(outcome_value),
    nrow = 1,
    widths = c(1, 1, 0.075)
  )
})

p_fg_cif <- (fg_header_row / wrap_plots(fg_row_plots, ncol = 1)) +
  plot_layout(heights = c(0.075, 3), guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = "bottom",
      legend.text = element_text(size = 19),
      plot.margin = margin(8, 14, 8, 14)
    )
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fine_gray_adjusted_cif_by_exposure_quartile.png"), p_fg_cif, width = 21, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "fine_gray_adjusted_cif_by_exposure_quartile.pdf"), p_fg_cif, width = 21, height = 14)

message("Wrote outputs:")
message(" - ", file.path(out_dir, "unadjusted_aalen_johansen_cif_by_exposure_quartile.csv"))
message(" - ", file.path(out_dir, "fine_gray_results_same_covariates_as_cox.csv"))
message(" - ", file.path(out_dir, "fine_gray_adjusted_cif_predictions_by_exposure_quartile.csv"))
message(" - ", file.path(fig_dir, "unadjusted_aalen_johansen_cif_by_exposure_quartile.png"))
message(" - ", file.path(fig_dir, "fine_gray_adjusted_cif_by_exposure_quartile.png"))
