#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
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
effect_path <- file.path(pooled_dir, "pooled_subgroup_effect_estimates.csv")
interaction_path <- file.path(pooled_dir, "pooled_subgroup_interaction_tests_fisher.csv")
if (!file.exists(effect_path)) {
  stop("Could not find pooled subgroup estimate file: ", effect_path)
}
if (!file.exists(interaction_path)) {
  stop("Could not find pooled interaction test file: ", interaction_path)
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
      strip.background = element_rect(fill = "#F2F2F2", color = "#B8B8B8", linewidth = 0.45),
      strip.text = element_text(color = "black", face = "bold", size = base_size - 0.5),
      panel.spacing.x = unit(1.8, "lines"),
      panel.spacing.y = unit(1.0, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.width = unit(20, "pt"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 34, 10, 12)
    )
}

fmt_p <- function(p) {
  vapply(
    p,
    function(x) {
      if (is.na(x)) {
        ""
      } else if (x < 0.0001) {
        "<0.0001"
      } else {
        sprintf("%.2g", x)
      }
    },
    character(1)
  )
}

p_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

p_stars_spaced <- function(p) {
  stars <- p_stars(p)
  if_else(stars == "", "", paste0(" ", stars))
}

subgroup_labels <- c(
  arf_subtype = "ARF subtype",
  race_ethnicity = "Race and ethnicity",
  sex = "Sex"
)

row_map <- tibble(
  subgroup = c(
    "arf_subtype", "arf_subtype", "arf_subtype",
    "race_ethnicity", "race_ethnicity", "race_ethnicity", "race_ethnicity", "race_ethnicity",
    "sex", "sex"
  ),
  subgroup_level = c(
    "Hypercapnic", "Hypoxemic", "Mixed",
    "Asian", "Hispanic White", "Non-Hispanic Black", "Non-Hispanic White", "Other/Unknown",
    "Female", "Male"
  ),
  y = c(12, 11, 10, 7.5, 6.5, 5.5, 4.5, 3.5, 1.5, 0.5)
)

header_map <- tibble(
  subgroup = c("arf_subtype", "race_ethnicity", "sex"),
  subgroup_type = c("ARF subtype", "Race and ethnicity", "Sex"),
  y = c(13.65, 9.15, 3.15)
)

pollutant_labels <- c(
  NO2 = "NO\u2082",
  `PM2.5` = "PM2.5"
)

effect_df <- readr::read_csv(effect_path, show_col_types = FALSE) %>%
  filter(
    outcome %in% c("Mortality by day 28", "Ventilator-free days"),
    pollutant %in% c("NO2", "PM2.5"),
    subgroup %in% names(subgroup_labels)
  ) %>%
  mutate(
    subgroup_type = recode(subgroup, !!!subgroup_labels),
    subgroup_type = factor(subgroup_type, levels = unname(subgroup_labels)),
    subgroup_level = str_replace(subgroup_level, "^Hispanic Other/Unknown$", "Other/Unknown"),
    outcome_label = recode(
      outcome,
      "Mortality by day 28" = "Mortality by day 28",
      "Ventilator-free days" = "Ventilator-free days"
    ),
    pollutant_label = recode(pollutant, !!!pollutant_labels),
    panel = factor(outcome_label, levels = c("Mortality by day 28", "Ventilator-free days")),
    point_stars = p_stars(p_value),
    estimate_ci = sprintf("%.2f (%.2f-%.2f)%s", estimate, conf_low, conf_high, point_stars),
    n_label = paste0("n=", format(n_total, big.mark = ",")),
    star_x = if_else(
      outcome == "Ventilator-free days",
      pmax(conf_low / 1.015, 0.865),
      pmin(conf_high * 1.08, 2.28)
    ),
    star_hjust = if_else(outcome == "Ventilator-free days", 1, 0)
  ) %>%
  left_join(row_map, by = c("subgroup", "subgroup_level")) %>%
  mutate(y_plot = y + if_else(pollutant == "NO2", 0.12, -0.12))

interaction_df <- readr::read_csv(interaction_path, show_col_types = FALSE) %>%
  filter(
    outcome %in% c("Mortality by day 28", "Ventilator-free days"),
    pollutant %in% c("NO2", "PM2.5"),
    subgroup %in% names(subgroup_labels)
  ) %>%
  mutate(
    subgroup_type = recode(subgroup, !!!subgroup_labels),
    subgroup_type = factor(subgroup_type, levels = unname(subgroup_labels)),
    outcome_label = recode(
      outcome,
      "Mortality by day 28" = "Mortality by day 28",
      "Ventilator-free days" = "Ventilator-free days"
    ),
    pollutant_label = recode(pollutant, !!!pollutant_labels),
    panel = factor(outcome_label, levels = levels(effect_df$panel)),
    p_text = fmt_p(p_fisher),
    interaction_label = if_else(
      str_starts(p_text, "<"),
      paste0(pollutant_label, " ", p_text, p_stars(p_fisher)),
      paste0(pollutant_label, " p=", p_text, p_stars(p_fisher))
    ),
    interaction_p_label = if_else(
      str_starts(p_text, "<"),
      paste0("p", p_text, p_stars_spaced(p_fisher)),
      paste0("p=", p_text, p_stars_spaced(p_fisher))
    )
  ) %>%
  left_join(header_map, by = c("subgroup", "subgroup_type")) %>%
  mutate(
    interaction_y = y - if_else(pollutant == "NO2", 0.78, 1.14)
  ) %>%
  arrange(panel, subgroup, pollutant)

interaction_headers <- interaction_df %>%
  distinct(panel, subgroup, subgroup_type, y)

axis_breaks <- row_map$y
axis_labels <- row_map$subgroup_level

colour_values <- c(
  NO2 = "#0072B2",
  `PM2.5` = "#D55E00"
)

make_panel <- function(outcome_name, x_limits, x_breaks, show_y_axis = TRUE) {
  panel_effect <- effect_df %>%
    filter(panel == outcome_name)

  if (outcome_name == "Ventilator-free days") {
    panel_effect <- panel_effect %>%
      mutate(
        star_x = pmax(conf_low / 1.01, x_limits[[1]] * 1.01),
        star_hjust = 1
      )
  } else {
    panel_effect <- panel_effect %>%
      mutate(
        star_x = pmin(conf_high * 1.035, x_limits[[2]] * 0.985),
        star_hjust = 0
      )
  }

  panel_interaction <- interaction_df %>%
    filter(panel == outcome_name)

  panel_interaction_headers <- interaction_headers %>%
    filter(panel == outcome_name)

  sig_legend_df <- tibble(
    x = x_limits[[1]],
    y = rep(0.05, 3),
    p_value_star = factor(
      c("*", "**", "***"),
      levels = c("*", "**", "***"),
      labels = c("* p<0.05", "** p<0.01", "*** p<0.001")
    )
  )

  p <- ggplot(panel_effect, aes(x = estimate, y = y_plot, colour = pollutant)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.14, linewidth = 0.6) +
    geom_point(size = 2.1) +
    geom_text(
      data = panel_interaction_headers,
      aes(x = x_limits[[1]], y = y, label = subgroup_type),
      hjust = 0,
      vjust = 0.5,
      size = 3.4,
      fontface = "bold",
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = panel_interaction_headers,
      aes(x = x_limits[[1]], y = y - 0.42, label = "Interaction p:"),
      hjust = 0,
      vjust = 0.5,
      size = 2.35,
      lineheight = 0.9,
      fontface = "italic",
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = panel_interaction,
      aes(x = x_limits[[1]], y = interaction_y, label = interaction_p_label, colour = pollutant),
      hjust = 0,
      vjust = 0.5,
      size = 2.35,
      lineheight = 0.9,
      fontface = "italic",
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = panel_effect %>% filter(point_stars != ""),
      aes(x = star_x, y = y_plot, label = point_stars, hjust = star_hjust),
      vjust = 0.55,
      size = 4.0,
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_point(
      data = sig_legend_df,
      aes(x = x, y = y, shape = p_value_star),
      alpha = 0,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      name = "Pollutant",
      values = colour_values,
      breaks = c("NO2", "PM2.5"),
      labels = c(expression(NO[2]), expression(PM[2.5]))
    ) +
    scale_shape_manual(
      name = "P value",
      values = c("* p<0.05" = 32, "** p<0.01" = 32, "*** p<0.001" = 32)
    ) +
    guides(
      colour = guide_legend(order = 1),
      shape = guide_legend(
        order = 2,
        override.aes = list(alpha = 0, colour = "black", size = 0)
      )
    ) +
    scale_x_log10(
      limits = x_limits,
      breaks = x_breaks,
      labels = function(x) sprintf("%.2f", x)
    ) +
    scale_y_continuous(
      breaks = axis_breaks,
      labels = if (show_y_axis) axis_labels else rep("", length(axis_labels)),
      limits = c(-0.1, 14.25),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Pooled ratio estimate with 95% CI",
      y = NULL,
      title = outcome_name
    ) +
    theme_refer_forest(10.5) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9.5),
      axis.text.x = element_text(size = 9.5),
      axis.title.x = element_text(size = 11, margin = margin(t = 8)),
      legend.position = "bottom",
      legend.title = element_text(size = 9.5, face = "bold"),
      legend.text = element_text(size = 10)
    )

  if (!show_y_axis) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  p
}

mortality_plot <- make_panel(
  "Mortality by day 28",
  x_limits = c(0.40, 2.35),
  x_breaks = c(0.5, 0.75, 1.0, 1.5, 2.0),
  show_y_axis = TRUE
)

vfd_plot <- make_panel(
  "Ventilator-free days",
  x_limits = c(0.86, 1.15),
  x_breaks = c(0.90, 0.95, 1.00, 1.05, 1.10, 1.15),
  show_y_axis = FALSE
)

forest_plot <- mortality_plot + vfd_plot +
  plot_layout(nrow = 1, guides = "collect", widths = c(1.08, 1)) &
  theme(
    legend.position = "bottom"
  )

out_base <- "pooled_subgroup_interaction_effect_forest"
journal_base <- "Supp_Figure_subgroup_interaction_effect_forest"

ggsave(file.path(fig_dir, paste0(out_base, ".png")), forest_plot, width = 10.8, height = 8.6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, paste0(out_base, ".pdf")), forest_plot, width = 10.8, height = 8.6, device = cairo_pdf, bg = "white")
ggsave(file.path(journal_dir, paste0(journal_base, ".png")), forest_plot, width = 10.8, height = 8.6, dpi = 300, bg = "white")
ggsave(file.path(journal_dir, paste0(journal_base, ".pdf")), forest_plot, width = 10.8, height = 8.6, device = cairo_pdf, bg = "white")

readr::write_csv(
  effect_df %>%
    select(outcome, pollutant, subgroup, subgroup_level, n_total, events_total, effect_measure, estimate, conf_low, conf_high, p_value, point_stars) %>%
    left_join(
      interaction_df %>%
        transmute(outcome, pollutant, subgroup, interaction_p_value = p_fisher, interaction_stars = p_stars(p_fisher)),
      by = c("outcome", "pollutant", "subgroup")
    ),
  file.path(pooled_dir, "pooled_subgroup_interaction_forest_values.csv")
)

message("Saved subgroup interaction forest plot:")
message(" - ", file.path(fig_dir, paste0(out_base, ".png")))
message(" - ", file.path(fig_dir, paste0(out_base, ".pdf")))
message(" - ", file.path(journal_dir, paste0(journal_base, ".png")))
message(" - ", file.path(journal_dir, paste0(journal_base, ".pdf")))
