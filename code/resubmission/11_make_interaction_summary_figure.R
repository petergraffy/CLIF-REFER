#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
pooled_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path("output", "resubmission_pooled", "20260819_all_sites_with_michigan")
}

pooled_dir <- normalizePath(pooled_dir, mustWork = TRUE)
interaction_path <- file.path(pooled_dir, "pooled_subgroup_interaction_tests_fisher.csv")
if (!file.exists(interaction_path)) {
  stop("Could not find pooled interaction test file: ", interaction_path)
}

fig_dir <- file.path(pooled_dir, "figures")
journal_dir <- file.path(pooled_dir, "journal_editable_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(journal_dir, recursive = TRUE, showWarnings = FALSE)

fmt_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", sprintf("%.3f", p))
  )
}

interaction_df <- readr::read_csv(interaction_path, show_col_types = FALSE) %>%
  mutate(
    interaction = recode(
      subgroup,
      arf_subtype = "ARF subtype",
      race_ethnicity = "Race and ethnicity",
      sex = "Sex",
      .default = str_replace_all(subgroup, "_", " ")
    ),
    outcome_label = recode(
      outcome,
      "Mortality by day 28" = "Mortality by day 28",
      "Ventilator-free days" = "Ventilator-free days",
      .default = outcome
    ),
    pollutant_label = recode(
      pollutant,
      "NO2" = "NO\u2082",
      "PM2.5" = "PM\u2082\u2024\u2085",
      .default = pollutant
    ),
    p_label = fmt_p(p_fisher),
    p_group = case_when(
      p_fisher < 0.001 ~ "p < 0.001",
      p_fisher < 0.01 ~ "0.001 <= p < 0.01",
      p_fisher < 0.05 ~ "0.01 <= p < 0.05",
      TRUE ~ "p >= 0.05"
    ),
    significant = p_fisher < 0.05,
    cell_label = p_label,
    text_colour = if_else(p_fisher < 0.001, "white", "black"),
    interaction = factor(interaction, levels = c("ARF subtype", "Race and ethnicity", "Sex")),
    outcome_label = factor(outcome_label, levels = c("Mortality by day 28", "Ventilator-free days")),
    pollutant_label = factor(pollutant_label, levels = c("NO\u2082", "PM\u2082.\u2085")),
    p_group = factor(
      p_group,
      levels = c("p < 0.001", "0.001 <= p < 0.01", "0.01 <= p < 0.05", "p >= 0.05")
    )
  )

column_labels <- c(
  "Mortality by day 28.NO\u2082" = "Mortality by day 28\nNO\u2082",
  "Mortality by day 28.PM\u2082.\u2085" = "Mortality by day 28\nPM\u2082.\u2085",
  "Ventilator-free days.NO\u2082" = "Ventilator-free days\nNO\u2082",
  "Ventilator-free days.PM\u2082.\u2085" = "Ventilator-free days\nPM\u2082.\u2085"
)

plot_df <- interaction_df %>%
  mutate(
    column = paste(outcome_label, pollutant_label, sep = "."),
    column = factor(column, levels = names(column_labels), labels = unname(column_labels))
  )

interaction_plot <- ggplot(plot_df, aes(x = column, y = interaction, fill = p_group)) +
  geom_tile(color = "black", linewidth = 0.45, width = 0.96, height = 0.82) +
  geom_text(aes(label = cell_label, color = text_colour), size = 4.2, lineheight = 0.95) +
  scale_fill_manual(
    values = c(
      "p < 0.001" = "#2166AC",
      "0.001 <= p < 0.01" = "#67A9CF",
      "0.01 <= p < 0.05" = "#D1E5F0",
      "p >= 0.05" = "#F7F7F7"
    ),
    name = "Interaction p-value"
  ) +
  scale_color_identity() +
  guides(fill = guide_legend(override.aes = list(color = "black", linewidth = 0.45))) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(plot_df$interaction))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 11.5, lineheight = 0.95),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10.5),
    legend.text = element_text(size = 10),
    legend.key.width = unit(24, "pt"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(14, 16, 10, 12)
  )

ggsave(
  file.path(fig_dir, "pooled_primary_interaction_summary.png"),
  interaction_plot,
  width = 8.4,
  height = 3.9,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(fig_dir, "pooled_primary_interaction_summary.pdf"),
  interaction_plot,
  width = 8.4,
  height = 3.9,
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(journal_dir, "Supp_Figure_primary_interaction_summary.png"),
  interaction_plot,
  width = 8.4,
  height = 3.9,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(journal_dir, "Supp_Figure_primary_interaction_summary.pdf"),
  interaction_plot,
  width = 8.4,
  height = 3.9,
  device = cairo_pdf,
  bg = "white"
)

message("Saved interaction summary figure:")
message(" - ", file.path(fig_dir, "pooled_primary_interaction_summary.png"))
message(" - ", file.path(fig_dir, "pooled_primary_interaction_summary.pdf"))
message(" - ", file.path(journal_dir, "Supp_Figure_primary_interaction_summary.png"))
message(" - ", file.path(journal_dir, "Supp_Figure_primary_interaction_summary.pdf"))
