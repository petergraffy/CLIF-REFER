#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
pooled_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path("output", "resubmission_pooled", "20260819_all_sites_with_michigan")
}

pooled_dir <- normalizePath(pooled_dir, mustWork = TRUE)
site_path <- file.path(pooled_dir, "site_fine_gray_effect_estimates.csv")
pooled_path <- file.path(pooled_dir, "pooled_fine_gray_effect_estimates.csv")
if (!file.exists(site_path)) stop("Missing site Fine-Gray estimates: ", site_path)
if (!file.exists(pooled_path)) stop("Missing pooled Fine-Gray estimates: ", pooled_path)

fig_dir <- file.path(pooled_dir, "figures")
journal_dir <- file.path(pooled_dir, "journal_editable_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(journal_dir, recursive = TRUE, showWarnings = FALSE)

site_labels <- c(
  Emory = "Emory",
  hopkins = "Hopkins",
  Michigan = "Michigan",
  NU = "Northwestern",
  OHSU = "OHSU",
  penn = "Penn",
  RUSH = "RUSH",
  UCMC = "UCMC",
  UCSF = "UCSF",
  UMN = "Minnesota"
)

site_effects <- read_csv(site_path, show_col_types = FALSE) %>%
  transmute(
    site = as.character(site),
    site_label = coalesce(as.character(site_label), recode(as.character(site), !!!site_labels, .default = as.character(site))),
    row_type = "Site-specific",
    analysis_family,
    outcome,
    model,
    term,
    pollutant,
    effect_measure,
    n_total = as.numeric(n),
    events_total = as.numeric(events),
    estimate = as.numeric(estimate),
    conf_low = as.numeric(conf_low),
    conf_high = as.numeric(conf_high),
    p_value = as.numeric(p_value)
  )

pooled_effects <- read_csv(pooled_path, show_col_types = FALSE) %>%
  transmute(
    site = "Pooled",
    site_label = "Pooled",
    row_type = "Pooled",
    analysis_family,
    outcome,
    model,
    term,
    pollutant,
    effect_measure,
    n_total = as.numeric(n_total),
    events_total = as.numeric(events_total),
    estimate = as.numeric(estimate),
    conf_low = as.numeric(conf_low),
    conf_high = as.numeric(conf_high),
    p_value = as.numeric(p_value)
  )

site_order_top <- c("Emory", "Hopkins", "Michigan", "Minnesota", "Northwestern", "OHSU", "Penn", "RUSH", "UCMC", "UCSF", "Pooled")
site_display_map <- stats::setNames(
  c(paste("Site", 1:10), "Pooled"),
  site_order_top
)
estimand_order <- c(
  "NO2 single-pollutant",
  "PM2.5 single-pollutant",
  "NO2 adjusted for PM2.5",
  "PM2.5 adjusted for NO2"
)
offsets <- c(
  "NO2 single-pollutant" = 0.24,
  "PM2.5 single-pollutant" = 0.08,
  "NO2 adjusted for PM2.5" = -0.08,
  "PM2.5 adjusted for NO2" = -0.24
)

plot_df <- bind_rows(site_effects, pooled_effects) %>%
  filter(analysis_family == "Adjusted Fine-Gray") %>%
  filter(model %in% c("NO2 single-pollutant", "PM25 single-pollutant", "PM25 + NO2")) %>%
  mutate(
    site_label = recode(site_label, "hopkins" = "Hopkins", "penn" = "Penn", "NU" = "Northwestern", "UMN" = "Minnesota", .default = site_label),
    outcome_panel = case_when(
      outcome == "Successful extubation" ~ "Successful extubation",
      outcome == "Death" ~ "Death",
      outcome == "Persistent respiratory failure" ~ "Persistent RF",
      TRUE ~ outcome
    ),
    x_title = "Adjusted subdistribution hazard ratio (95% CI)",
    estimand = case_when(
      model == "NO2 single-pollutant" & term == "no2_per_10" ~ "NO2 single-pollutant",
      model == "PM25 single-pollutant" & term == "pm25_per_5" ~ "PM2.5 single-pollutant",
      model == "PM25 + NO2" & term == "no2_per_10" ~ "NO2 adjusted for PM2.5",
      model == "PM25 + NO2" & term == "pm25_per_5" ~ "PM2.5 adjusted for NO2",
      TRUE ~ NA_character_
    ),
    site_label = factor(site_label, levels = site_order_top),
    site_display = factor(site_display_map[as.character(site_label)], levels = unname(site_display_map)),
    y_base = length(site_order_top) - as.integer(site_label) + 1,
    estimand = factor(estimand, levels = estimand_order),
    y_plot = y_base + unname(offsets[as.character(estimand)])
  ) %>%
  filter(!is.na(estimand), is.finite(estimate), is.finite(conf_low), is.finite(conf_high))

readr::write_csv(plot_df, file.path(pooled_dir, "pooled_fine_gray_effect_forest_publication_ready_data.csv"))

site_cols <- c(
  "Emory" = "#E69F00",
  "Hopkins" = "#56B4E9",
  "Michigan" = "#009E73",
  "Minnesota" = "#CC79A7",
  "Northwestern" = "#0072B2",
  "OHSU" = "#D55E00",
  "Penn" = "#7E57C2",
  "RUSH" = "#5E5E5E",
  "UCMC" = "#00A6D6",
  "UCSF" = "#8DA0CB",
  "Pooled" = "#000000"
)

shape_vals <- c(
  "NO2 single-pollutant" = 21,
  "PM2.5 single-pollutant" = 22,
  "NO2 adjusted for PM2.5" = 24,
  "PM2.5 adjusted for NO2" = 23
)

shape_labs <- c(
  "NO2 single-pollutant" = expression(NO[2]~"single-pollutant"),
  "PM2.5 single-pollutant" = expression(PM[2.5]~"single-pollutant"),
  "NO2 adjusted for PM2.5" = expression(NO[2]~"adjusted for"~PM[2.5]),
  "PM2.5 adjusted for NO2" = expression(PM[2.5]~"adjusted for"~NO[2])
)

theme_forest <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1, margin = margin(b = 8)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.key.width = unit(12, "pt"),
      legend.box = "vertical",
      plot.margin = margin(8, 10, 8, 8),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

panel_limits <- function(df) {
  lo <- min(df$conf_low, na.rm = TRUE)
  hi <- max(df$conf_high, na.rm = TRUE)
  c(max(0.08, lo * 0.82), hi * 1.12)
}

panel_breaks <- function(limits) {
  candidates <- if (limits[[2]] > 10) {
    c(0.10, 0.50, 1.00, 5.00, 20.00)
  } else if (limits[[2]] > 4) {
    c(0.25, 0.50, 1.00, 2.00, 5.00)
  } else {
    c(0.50, 0.75, 1.00, 1.50, 2.00, 3.00)
  }
  candidates[candidates >= limits[[1]] & candidates <= limits[[2]]]
}

make_panel <- function(df, panel_title, show_y = TRUE) {
  x_limits <- panel_limits(df)
  x_breaks <- panel_breaks(x_limits)
  ggplot(df, aes(colour = site_label, fill = site_label, shape = estimand)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.45) +
    geom_segment(
      aes(x = conf_low, xend = conf_high, y = y_plot, yend = y_plot, linewidth = row_type),
      lineend = "butt",
      show.legend = FALSE
    ) +
    geom_point(aes(x = estimate, y = y_plot, size = row_type), stroke = 0.45, colour = "black") +
    scale_colour_manual(values = site_cols, drop = FALSE) +
    scale_fill_manual(values = site_cols, drop = FALSE) +
    scale_shape_manual(values = shape_vals, labels = shape_labs, drop = FALSE) +
    scale_size_manual(values = c("Site-specific" = 2.25, "Pooled" = 3.05), guide = "none") +
    scale_linewidth_manual(values = c("Site-specific" = 0.46, "Pooled" = 0.85), guide = "none") +
    scale_y_continuous(
      limits = c(0.45, length(site_order_top) + 0.55),
      breaks = seq_along(rev(site_order_top)),
      labels = rev(unname(site_display_map)),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_x_log10(
      limits = x_limits,
      breaks = x_breaks,
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      title = panel_title,
      x = "Adjusted subdistribution hazard ratio (95% CI)",
      y = NULL,
      colour = "Site",
      fill = "Site",
      shape = "Pollutant model"
    ) +
    theme_forest(11) +
    theme(
      axis.text.y = if (show_y) element_text(size = 10) else element_blank(),
      axis.ticks.y = if (show_y) element_line(color = "black", linewidth = 0.35) else element_blank(),
      axis.line.y = if (show_y) element_line(color = "black", linewidth = 0.35) else element_blank(),
      axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
      legend.text = element_text(size = 9.2),
      legend.title = element_text(size = 9.6)
    ) +
    guides(
      colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(shape = 21, size = 2.8, linewidth = 0)),
      fill = "none",
      shape = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(fill = "grey80", colour = "black", size = 2.8))
    )
}

p_extubation <- make_panel(
  plot_df %>% filter(outcome_panel == "Successful extubation"),
  "Successful extubation",
  show_y = TRUE
)
p_death <- make_panel(
  plot_df %>% filter(outcome_panel == "Death"),
  "Death",
  show_y = FALSE
)
p_persistent <- make_panel(
  plot_df %>% filter(outcome_panel == "Persistent RF"),
  "Persistent RF",
  show_y = FALSE
)

p <- (p_extubation | p_death | p_persistent) +
  plot_layout(widths = c(1.12, 1, 1), guides = "collect") &
  theme(legend.position = "bottom")

save_plot_set <- function(plot, stem, target_dir) {
  ggsave(file.path(target_dir, paste0(stem, ".pdf")), plot, width = 14.2, height = 7.4, units = "in", device = cairo_pdf)
  ggsave(file.path(target_dir, paste0(stem, ".png")), plot, width = 14.2, height = 7.4, units = "in", dpi = 450)
}

save_plot_set(p, "pooled_fine_gray_effect_forest_publication_ready", fig_dir)
save_plot_set(p, "pooled_fine_gray_effect_forest", fig_dir)
save_plot_set(p, "Supp_Figure_fine_gray_effect_forest", journal_dir)

message("Saved publication-ready Fine-Gray SHR forest to: ", file.path(fig_dir, "pooled_fine_gray_effect_forest_publication_ready.png"))
