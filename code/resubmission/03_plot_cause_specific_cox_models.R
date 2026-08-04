#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
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
  candidates <- candidates[file.exists(file.path(candidates, "resubmission_cause_specific_cox_results.csv"))]
  if (!length(candidates)) stop("No resubmission output directory with Cox results found.")
  out_root <- candidates[which.max(file.info(candidates)$mtime)]
}
out_root <- normalizePath(out_root, mustWork = TRUE)
fig_dir <- file.path(out_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pollutant_colors <- c(
  "Fine particulate matter (PM2.5)" = "#0072B2",
  "Nitrogen dioxide (NO2)" = "#D55E00"
)

model_row_labels <- c(
  "PM2.5 only" = "PM[2.5]*' only'",
  "NO2 only" = "NO[2]*' only'",
  "PM2.5 in PM2.5 + NO2" = "PM[2.5]*' in PM'[2.5]*' + NO'[2]",
  "NO2 in PM2.5 + NO2" = "NO[2]*' in PM'[2.5]*' + NO'[2]"
)

theme_cox <- function(base_size = 18) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid.major.y = element_line(linewidth = 0.18, color = "grey91"),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey86"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.55),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      plot.title = element_text(face = "bold", size = rel(1.25), hjust = 0.02, margin = margin(b = 7)),
      plot.subtitle = element_text(color = "grey25", hjust = 0.02, margin = margin(b = 10)),
      strip.text = element_text(face = "bold", color = "grey10", size = rel(1.0), margin = margin(5, 5, 5, 5)),
      strip.background = element_rect(fill = "grey94", color = "grey68", linewidth = 0.55),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.98)),
      panel.spacing.x = unit(12, "pt"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

cox_results <- read_csv(
  file.path(out_root, "resubmission_cause_specific_cox_results.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    outcome = recode(
      cause,
      "Successful extubation" = "Extubation",
      "Death" = "Death",
      "Persistent respiratory failure" = "Persistent RF"
    ),
    outcome = factor(outcome, levels = c("Extubation", "Death", "Persistent RF")),
    pollutant = case_when(
      term == "pm25_per_5" ~ "Fine particulate matter (PM2.5)",
      term == "no2_per_10" ~ "Nitrogen dioxide (NO2)",
      TRUE ~ term
    ),
    pollutant = factor(pollutant, levels = names(pollutant_colors)),
    scale_label = case_when(
      term == "pm25_per_5" ~ "per 5 ug/m3",
      term == "no2_per_10" ~ "per 10 ppb",
      TRUE ~ ""
    ),
    model_row = case_when(
      model == "PM25 single-pollutant" ~ "PM2.5 only",
      model == "NO2 single-pollutant" ~ "NO2 only",
      model == "PM25 + NO2" & term == "pm25_per_5" ~ "PM2.5 in PM2.5 + NO2",
      model == "PM25 + NO2" & term == "no2_per_10" ~ "NO2 in PM2.5 + NO2",
      TRUE ~ model
    ),
    model_row = factor(
      model_row,
      levels = rev(c("PM2.5 only", "NO2 only", "PM2.5 in PM2.5 + NO2", "NO2 in PM2.5 + NO2"))
    ),
    estimate_label = sprintf("%.2f (%.2f, %.2f)", hazard_ratio, conf_low, conf_high),
    text_x = conf_high * 1.08
  ) %>%
  group_by(outcome) %>%
  mutate(label_limit_x = max(conf_high, na.rm = TRUE) * 1.65) %>%
  ungroup()

p_cox <- ggplot(cox_results, aes(x = hazard_ratio, y = model_row, color = pollutant)) +
  geom_blank(aes(x = label_limit_x)) +
  geom_vline(xintercept = 1, linewidth = 0.45, color = "grey35") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.16, linewidth = 0.95) +
  geom_point(size = 3.2) +
  geom_text(
    aes(x = text_x, label = estimate_label),
    hjust = 0,
    size = 4.2,
    color = "grey12",
    show.legend = FALSE
  ) +
  facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
  scale_color_manual(
    values = pollutant_colors,
    labels = c(
      "Fine particulate matter (PM2.5)" = expression("Fine particulate matter (PM"[2.5]*")"),
      "Nitrogen dioxide (NO2)" = expression("Nitrogen dioxide (NO"[2]*")")
    ),
    drop = FALSE
  ) +
  scale_x_log10(
    breaks = scales::breaks_log(n = 6),
    labels = label_number(accuracy = 0.01, trim = TRUE),
    expand = expansion(mult = c(0.04, 0.20))
  ) +
  scale_y_discrete(labels = function(x) parse(text = model_row_labels[x])) +
  labs(
    title = "Cause-Specific Cox Models",
    subtitle = "Hazard ratios adjusted for demographics, Charlson score, ARF subtype, index year, and ZCTA ACS social vulnerability covariates",
    x = "Hazard ratio, log scale",
    y = NULL,
    color = "Pollutant"
  ) +
  theme_cox()

png_path <- file.path(fig_dir, "resubmission_cause_specific_cox_models.png")
pdf_path <- file.path(fig_dir, "resubmission_cause_specific_cox_models.pdf")
ggsave(png_path, p_cox, width = 20, height = 9.5, dpi = 300)
ggsave(pdf_path, p_cox, width = 20, height = 9.5)

message("Wrote Cox model figure:")
message(" - ", png_path)
message(" - ", pdf_path)
