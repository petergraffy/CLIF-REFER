#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
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
  candidates <- candidates[file.exists(file.path(candidates, "resubmission_aalen_johansen_cif_plot_data.csv"))]
  if (!length(candidates)) stop("No resubmission output directory with CIF plot data found.")
  out_root <- candidates[which.max(file.info(candidates)$mtime)]
}
out_root <- normalizePath(out_root, mustWork = TRUE)
fig_dir <- file.path(out_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
quartile_colors <- c(
  "Q1 lowest" = "#2166AC",
  "Q2" = "#67A9CF",
  "Q3" = "#F4A582",
  "Q4 highest" = "#B2182B"
)

theme_cif_transplant <- function() {
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
      plot.title = element_text(face = "bold", size = 21),
      plot.caption = element_text(size = 13, color = "grey25", hjust = 0.5, margin = margin(t = 8)),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16, color = "grey20"),
      strip.text = element_text(face = "bold", size = 18, color = "grey10"),
      panel.spacing.x = unit(28, "pt"),
      panel.spacing.y = unit(22, "pt"),
      plot.margin = margin(8, 32, 8, 8)
    )
}

theme_table_box <- function() {
  theme_minimal(base_size = 16) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.title = element_text(size = 17, color = "grey10", margin = margin(t = 7)),
      axis.text.x = element_text(size = 16, color = "grey20", margin = margin(t = 2)),
      axis.text.y = element_text(size = 17, face = "bold"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.25),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 17, face = "bold", color = "grey20", margin = margin(b = 3)),
      plot.margin = margin(3, 8, 4, 8)
    )
}

cif_plot_data <- readr::read_csv(
  file.path(out_root, "resubmission_aalen_johansen_cif_plot_data.csv"),
  show_col_types = FALSE
) %>%
  filter(!state %in% c("censor", "(s0)")) %>%
  mutate(
    exposure = recode(exposure, "PM2.5" = "PM2.5", "NO2" = "NO2"),
    state = factor(state, levels = c("Successful extubation", "Death", "Persistent respiratory failure")),
    quartile = recode(
      as.character(quartile),
      "Q1" = "Q1 lowest",
      "Q2" = "Q2",
      "Q3" = "Q3",
    "Q4" = "Q4 highest"
    ),
    quartile = factor(quartile, levels = quartile_labels)
  ) %>%
  mutate(
    state = recode(
      as.character(state),
      "Successful extubation" = "Extubation",
      "Persistent respiratory failure" = "Persistent RF"
    )
  ) %>%
  mutate(state = factor(state, levels = c("Extubation", "Death", "Persistent RF"))) %>%
  filter(!is.na(quartile), !is.na(state))

analysis_data <- readr::read_csv(
  file.path(out_root, "resubmission_analysis_dataset.csv"),
  show_col_types = FALSE
)

tick_times <- c(0, 7, 14, 21, 28)
x_limits <- c(-2.5, 30)
event_code_lookup <- tibble(
  state = factor(c("Extubation", "Death", "Persistent RF"), levels = c("Extubation", "Death", "Persistent RF")),
  event_code = c(1L, 2L, 3L)
)

make_risk_event_table <- function(data, exposure_var, exposure_label) {
  base <- data %>%
    filter(!is.na(.data[[exposure_var]])) %>%
    mutate(
      quartile = recode(
        as.character(.data[[exposure_var]]),
        "Q1" = "Q1 lowest",
        "Q2" = "Q2",
        "Q3" = "Q3",
        "Q4" = "Q4 highest"
      ),
      quartile = factor(quartile, levels = quartile_labels)
    ) %>%
    filter(!is.na(quartile))

  tidyr::expand_grid(
    exposure = exposure_label,
    state = event_code_lookup$state,
    quartile = factor(quartile_labels, levels = quartile_labels),
    time = tick_times
  ) %>%
    left_join(event_code_lookup, by = "state") %>%
    rowwise() %>%
    mutate(
      n_risk = sum(base$quartile == quartile & base$ftime_days >= time, na.rm = TRUE),
      n_events = sum(base$quartile == quartile & base$event_code == event_code & base$ftime_days <= time, na.rm = TRUE),
      label = paste0(n_events, "/", n_risk)
    ) %>%
    ungroup()
}

risk_event_table <- bind_rows(
  make_risk_event_table(analysis_data, "no2_q", "NO2"),
  make_risk_event_table(analysis_data, "pm25_q", "PM2.5")
)

panel_scales <- cif_plot_data %>%
  group_by(exposure, state) %>%
  summarise(y_max = max(cif, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_limit = y_max * 1.04)

make_curve_panel <- function(exposure_value, state_value, show_y_axis = TRUE) {
  dat <- cif_plot_data %>% filter(exposure == exposure_value, state == state_value)
  y_limit <- panel_scales %>%
    filter(exposure == exposure_value, state == state_value) %>%
    pull(y_limit)

  ggplot(dat, aes(x = time, y = cif, color = quartile)) +
    geom_step(linewidth = 1.05) +
    scale_color_manual(values = quartile_colors, drop = FALSE) +
    scale_x_continuous(
      breaks = tick_times,
      labels = label_number(accuracy = 1, trim = TRUE),
      limits = x_limits,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      labels = label_percent(accuracy = 1),
      limits = c(0, y_limit),
      expand = expansion(mult = c(0.015, 0.035))
    ) +
    labs(
      title = state_value,
      x = NULL,
      y = if (show_y_axis) "Cumulative incidence" else NULL,
      color = "Quartile"
    ) +
    theme_cif_transplant() +
    theme(
      axis.title.y = if (show_y_axis) element_text(size = 19) else element_blank(),
      axis.text.y = if (show_y_axis) element_text(size = 17, color = "grey20") else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line(color = "black", linewidth = 0.35) else element_blank(),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 4), size = 22),
      legend.position = "bottom",
      plot.margin = margin(4, 8, 1, 8)
    )
}

make_table_panel <- function(exposure_value, state_value, show_x_axis = FALSE) {
  dat <- risk_event_table %>%
    filter(exposure == exposure_value, state == state_value) %>%
    mutate(quartile = factor(quartile, levels = rev(quartile_labels)))

  ggplot(dat, aes(x = time, y = quartile, label = label, color = quartile)) +
    geom_text(size = 5.2, fontface = "bold") +
    scale_color_manual(values = quartile_colors, drop = FALSE, guide = "none") +
    scale_x_continuous(
      breaks = tick_times,
      labels = label_number(accuracy = 1),
      limits = x_limits,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_discrete(labels = function(x) str_replace(x, " lowest| highest", "")) +
    labs(title = "Events / at risk", x = if (show_x_axis) "Days since ARF onset" else NULL, y = NULL) +
    theme_table_box() +
    theme(
      axis.text.x = if (show_x_axis) element_text(size = 16, color = "grey20", margin = margin(t = 2)) else element_blank(),
      axis.ticks.x = if (show_x_axis) element_line(color = "black", linewidth = 0.25) else element_blank(),
      legend.position = "none",
      plot.margin = margin(1, 8, if (show_x_axis) 6 else 2, 8)
    )
}

make_cell <- function(exposure_value, state_value, show_y_axis = TRUE, show_x_axis = FALSE) {
  make_curve_panel(exposure_value, state_value, show_y_axis = show_y_axis) /
    make_table_panel(exposure_value, state_value, show_x_axis = show_x_axis) +
    plot_layout(heights = c(4.2, 1.28))
}

state_levels <- c("Extubation", "Death", "Persistent RF")
exposure_levels <- c("NO2", "PM2.5")

make_column_header <- function(label_expr) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "grey45", linewidth = 0.55) +
    annotate("text", x = 0.5, y = 0.5, label = label_expr, parse = TRUE, hjust = 0.5, vjust = 0.5, fontface = "bold", size = 6.4) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(6, 8, -14, 8))
}

header_row <- wrap_plots(
  make_column_header("'Nitrogen dioxide'~(NO[2])"),
  make_column_header("'Fine particulate matter'~(PM[2.5])"),
  nrow = 1
)

row_plots <- lapply(seq_along(state_levels), function(i) {
  state_value <- state_levels[[i]]
  wrap_plots(
    make_cell(exposure_levels[[1]], state_value, show_y_axis = TRUE, show_x_axis = i == length(state_levels)),
    make_cell(exposure_levels[[2]], state_value, show_y_axis = FALSE, show_x_axis = i == length(state_levels)),
    nrow = 1
  )
})

p_cif <- (header_row / wrap_plots(row_plots, ncol = 1)) +
  plot_layout(heights = c(0.075, 3), guides = "collect") +
  plot_annotation(
    title = "Aalen-Johansen Cumulative Incidence by Exposure Quartile",
    theme = theme(
      plot.title = element_text(face = "bold", size = 30, hjust = 0.02, margin = margin(b = 8)),
      legend.position = "bottom",
      legend.text = element_text(size = 19),
      plot.margin = margin(8, 14, 8, 14)
    )
  ) &
  theme(legend.position = "bottom")

png_path <- file.path(fig_dir, "resubmission_aalen_johansen_cif_quartiles_transplant_style.png")
pdf_path <- file.path(fig_dir, "resubmission_aalen_johansen_cif_quartiles_transplant_style.pdf")
ggsave(png_path, p_cif, width = 20, height = 19, dpi = 300)
ggsave(pdf_path, p_cif, width = 20, height = 19)

message("Wrote transplant-style AJ curves:")
message(" - ", png_path)
message(" - ", pdf_path)
