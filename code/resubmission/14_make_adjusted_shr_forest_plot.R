#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
pooled_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path("output", "resubmission_pooled", "20260819_all_sites_with_michigan")
}

pooled_dir <- normalizePath(pooled_dir, mustWork = TRUE)
effect_path <- file.path(pooled_dir, "pooled_fine_gray_effect_estimates.csv")
if (!file.exists(effect_path)) {
  stop("Could not find pooled Fine-Gray estimate file: ", effect_path)
}

fig_dir <- file.path(pooled_dir, "figures")
journal_dir <- file.path(pooled_dir, "journal_editable_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(journal_dir, recursive = TRUE, showWarnings = FALSE)

theme_refer_forest <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.width = unit(22, "pt"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 24, 10, 12)
    )
}

effects <- read.csv(effect_path, stringsAsFactors = FALSE, check.names = FALSE)

plot_df <- effects[
  effects$model %in% c("NO2 single-pollutant", "PM25 single-pollutant") &
    effects$pollutant %in% c("NO2", "PM2.5") &
    effects$outcome %in% c("Successful extubation", "Death", "Persistent respiratory failure"),
  ,
  drop = FALSE
]

if (nrow(plot_df) == 0) {
  stop("No single-pollutant adjusted Fine-Gray estimates found in: ", effect_path)
}

plot_df$outcome_label <- ifelse(
  plot_df$outcome == "Persistent respiratory failure",
  "Persistent RF",
  plot_df$outcome
)
plot_df$y_base <- c(
  "Persistent RF" = 1,
  "Death" = 2,
  "Successful extubation" = 3
)[plot_df$outcome_label]
plot_df$pollutant_label <- factor(plot_df$pollutant, levels = c("NO2", "PM2.5"))
plot_df$y_plot <- plot_df$y_base + ifelse(plot_df$pollutant == "NO2", 0.08, -0.08)
plot_df$estimate_label <- sprintf("%.2f (%.2f-%.2f)", plot_df$estimate, plot_df$conf_low, plot_df$conf_high)
plot_df$label_x <- plot_df$conf_high * 1.018

pollutant_cols <- c("NO2" = "#0072B2", "PM2.5" = "#D55E00")
pollutant_labs <- c("NO2" = expression(NO[2]~"(per 10 ppb)"), "PM2.5" = expression(PM[2.5]~"(per 5 "*mu*"g/m"^3*")"))

x_min <- max(0.8, min(plot_df$conf_low, na.rm = TRUE) * 0.94)
x_max <- max(plot_df$conf_high, na.rm = TRUE) * 1.42
break_candidates <- c(0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50)
x_breaks <- break_candidates[break_candidates >= x_min & break_candidates <= x_max]
if (!1 %in% x_breaks) {
  x_breaks <- sort(unique(c(x_breaks, 1)))
}

p <- ggplot(plot_df, aes(colour = pollutant_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.45) +
  geom_segment(aes(x = conf_low, xend = conf_high, y = y_plot, yend = y_plot), linewidth = 0.72) +
  geom_segment(aes(x = conf_low, xend = conf_low, y = y_plot - 0.045, yend = y_plot + 0.045), linewidth = 0.72) +
  geom_segment(aes(x = conf_high, xend = conf_high, y = y_plot - 0.045, yend = y_plot + 0.045), linewidth = 0.72) +
  geom_point(aes(x = estimate, y = y_plot), size = 2.8) +
  geom_text(
    aes(x = label_x, y = y_plot, label = estimate_label),
    hjust = 0,
    size = 3.25,
    colour = "black",
    show.legend = FALSE
  ) +
  scale_colour_manual(values = pollutant_cols, labels = pollutant_labs, drop = FALSE) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = x_breaks,
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.01, 0.035))
  ) +
  scale_y_continuous(
    limits = c(0.62, 3.38),
    breaks = c(1, 2, 3),
    labels = c("Persistent RF", "Death", "Successful extubation"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Adjusted subdistribution hazard ratio (95% CI)",
    y = NULL,
    colour = NULL
  ) +
  theme_refer_forest(11) +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.text.x = element_text(size = 10.5),
    axis.text.y = element_text(size = 11.5),
    legend.text = element_text(size = 10.5),
    legend.margin = margin(t = 4)
  )

for (target_dir in c(fig_dir, journal_dir)) {
  ggsave(
    filename = file.path(target_dir, "pooled_adjusted_shr_forest_single.pdf"),
    plot = p,
    width = 8.8,
    height = 3.6,
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(target_dir, "pooled_adjusted_shr_forest_single.png"),
    plot = p,
    width = 8.8,
    height = 3.6,
    dpi = 450
  )
}

message("Saved adjusted SHR forest plot to: ", file.path(fig_dir, "pooled_adjusted_shr_forest_single.png"))
