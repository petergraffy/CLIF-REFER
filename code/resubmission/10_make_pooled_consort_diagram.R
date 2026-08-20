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
flow_path <- file.path(pooled_dir, "pooled_inclusion_flow_counts.csv")
if (!file.exists(flow_path)) {
  stop("Could not find pooled inclusion flow file: ", flow_path)
}

fig_dir <- file.path(pooled_dir, "figures")
journal_dir <- file.path(pooled_dir, "journal_editable_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(journal_dir, recursive = TRUE, showWarnings = FALSE)

flow <- readr::read_csv(flow_path, show_col_types = FALSE)

first_draft_icu_admissions <- 665074

fmt_n <- function(x) {
  ifelse(is.na(x), "", format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE))
}

fmt_pct <- function(x, denom) {
  ifelse(is.na(x) | is.na(denom) | denom == 0, "", sprintf("%.1f%%", 100 * x / denom))
}

step_value <- function(step_order, col = "n_hospitalizations_or_encounters") {
  out <- flow %>% filter(.data$step_order == !!step_order) %>% slice(1)
  if (!nrow(out) || !col %in% names(out)) return(NA_real_)
  out[[col]]
}

exclusion_label <- function(title, n_excluded, denominator) {
  paste0(title, "\nn = ", fmt_n(n_excluded), " (", fmt_pct(n_excluded, denominator), ")")
}

adult_n <- 665042
demographics_n <- 663876
icu24_n <- 532905
icu24_excl_n <- 130971
geography_n <- 522090
geography_excl_n <- 10815
physiology_n <- 522004
physiology_excl_n <- 86
arf_n <- step_value(4)
no_arf_n <- physiology_n - arf_n
imv_n <- step_value(11)
no_imv_n <- step_value(11, "n_excluded_from_prior")

main_label <- function(title, step_order) {
  n_hosp <- step_value(step_order, "n_hospitalizations_or_encounters")
  n_pat <- step_value(step_order, "n_patients")
  if (!is.na(n_pat)) {
    paste0(title, "\nN=", fmt_n(n_hosp), " encounters\n", fmt_n(n_pat), " patients")
  } else {
    paste0(title, "\nN=", fmt_n(n_hosp), " encounters")
  }
}

small_label <- function(title, step_order, col = "n_hospitalizations_or_encounters") {
  paste0(title, "\nN=", fmt_n(step_value(step_order, col)), " encounters")
}

boxes <- tribble(
  ~id, ~x, ~y, ~label, ~type,
  "icu", 0.0, 8.2,
  paste0("ICU 2018-2024 (candidates)\nn = ", fmt_n(first_draft_icu_admissions)), "included",
  "adult", 0.0, 7.25,
  ">=18 years\nn = 665,042", "included",
  "demographics", 0.0, 6.25,
  "Demographics present\nn = 663,876", "included",
  "icu24", -0.95, 5.1,
  paste0("ICU stay >=24h\nn = ", fmt_n(icu24_n)), "included",
  "icu24_excl", 1.35, 5.1,
  exclusion_label("ICU stay <24h", icu24_excl_n, demographics_n), "excluded",
  "geography", -1.9, 3.95,
  paste0("Geography present (ZCTA)\nn = ", fmt_n(geography_n)), "included",
  "geography_excl", 1.15, 3.95,
  exclusion_label("Missing geography", geography_excl_n, icu24_n), "excluded",
  "physiology", -2.75, 2.8,
  paste0("ABG or continuous SpO2 in 24h\nn = ", fmt_n(physiology_n)), "included",
  "physiology_excl", 0.1, 2.8,
  exclusion_label("No ABG or continuous SpO2", physiology_excl_n, geography_n), "excluded",
  "arf", -3.55, 1.65,
  paste0("Meets ARF criteria in 24h\nn = ", fmt_n(arf_n)), "included",
  "no_arf", -1.6, 1.65,
  exclusion_label("No ARF criteria in +/-24h", no_arf_n, physiology_n), "excluded",
  "imv", -4.45, 0.45,
  paste0("IMV after ARF onset\nn = ", fmt_n(imv_n)), "included",
  "no_imv", -2.25, 0.45,
  exclusion_label("No IMV after ARF onset", no_imv_n, arf_n), "excluded"
)

segments <- tribble(
  ~x, ~y, ~xend, ~yend,
  0.0, 7.9, 0.0, 7.55,
  0.0, 6.95, 0.0, 6.55,
  0.0, 5.95, -0.75, 5.4,
  0.0, 5.95, 1.15, 5.4,
  -0.95, 4.8, -1.7, 4.25,
  -0.95, 4.8, 0.95, 4.25,
  -1.9, 3.65, -2.55, 3.1,
  -1.9, 3.65, -0.1, 3.1,
  -2.75, 2.5, -3.35, 1.95,
  -2.75, 2.5, -1.8, 1.95,
  -3.55, 1.35, -4.25, 0.75,
  -3.55, 1.35, -2.45, 0.75
)

consort_plot <- ggplot() +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.55,
    arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
    color = "black"
  ) +
  geom_label(
    data = boxes,
    aes(x = x, y = y, label = label, color = type),
    fill = "#F2F2F2",
    size = 3.9,
    linewidth = 0.55,
    label.padding = unit(0.22, "lines"),
    label.r = unit(0.35, "lines"),
    lineheight = 0.95,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      included = "black",
      excluded = "#E32222"
    )
  ) +
  coord_cartesian(xlim = c(-5.05, 2.45), ylim = c(-0.25, 8.65), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(18, 18, 18, 18)
  )

out_files <- file.path(
  c(fig_dir, journal_dir),
  c("pooled_inclusion_flow_consort.png", "Figure_S1_pooled_inclusion_flow_consort.pdf")
)

ggsave(file.path(fig_dir, "pooled_inclusion_flow_consort.png"), consort_plot, width = 8.5, height = 8.4, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "pooled_inclusion_flow_consort.pdf"), consort_plot, width = 8.5, height = 8.4, device = cairo_pdf, bg = "white")
ggsave(file.path(journal_dir, "Figure_S1_pooled_inclusion_flow_consort.png"), consort_plot, width = 8.5, height = 8.4, dpi = 300, bg = "white")
ggsave(file.path(journal_dir, "Figure_S1_pooled_inclusion_flow_consort.pdf"), consort_plot, width = 8.5, height = 8.4, device = cairo_pdf, bg = "white")

message("Saved pooled inclusion-flow diagram:")
message(" - ", file.path(fig_dir, "pooled_inclusion_flow_consort.png"))
message(" - ", file.path(fig_dir, "pooled_inclusion_flow_consort.pdf"))
message(" - ", file.path(journal_dir, "Figure_S1_pooled_inclusion_flow_consort.png"))
message(" - ", file.path(journal_dir, "Figure_S1_pooled_inclusion_flow_consort.pdf"))
