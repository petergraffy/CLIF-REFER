#!/usr/bin/env Rscript

# Backfill non-PHI exposure distribution bins from a restored site analysis dataset.
#
# Usage:
#   Rscript code/resubmission/99_backfill_exposure_distribution_bins.R <site_output_dir>
#
# The script writes only primary_exposure_distribution_bins.csv and does not
# modify or delete the source row-level files.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !nzchar(args[[1]])) {
  stop("Usage: Rscript code/resubmission/99_backfill_exposure_distribution_bins.R <site_output_dir>")
}

site_dir <- normalizePath(args[[1]], mustWork = TRUE)
dataset_path <- file.path(site_dir, "analysis_dataset_reviewer_optimized.csv")
if (!file.exists(dataset_path)) {
  dataset_path <- file.path(site_dir, "resubmission_analysis_dataset.csv")
}
if (!file.exists(dataset_path)) {
  stop("No restored analysis dataset found in: ", site_dir)
}

analysis_df <- readr::read_csv(dataset_path, show_col_types = FALSE, progress = FALSE)

required_cols <- c("no2_12m_zcta", "pm25_12m_zcta", "arf_subtype", "sex", "race_ethnicity")
missing_cols <- setdiff(required_cols, names(analysis_df))
if (length(missing_cols)) {
  stop("Dataset is missing required columns for exposure bins: ", paste(missing_cols, collapse = ", "))
}

site_name <- if ("site" %in% names(analysis_df)) {
  first(stats::na.omit(as.character(analysis_df$site))) %||% basename(dirname(site_dir))
} else {
  basename(dirname(site_dir))
}

pollutant_specs <- tibble::tribble(
  ~pollutant, ~raw_col, ~display_label,
  "NO2", "no2_12m_zcta", "Nitrogen~dioxide~(NO[2])",
  "PM2.5", "pm25_12m_zcta", "Fine~particulate~matter~(PM[2.5])"
)

pollutant_bin_widths <- c(
  "NO2" = as.numeric(Sys.getenv("REFER_RUG_BIN_WIDTH_NO2", "0.025")),
  "PM2.5" = as.numeric(Sys.getenv("REFER_RUG_BIN_WIDTH_PM25", "0.01"))
)

make_exposure_distribution_bins <- function(df, group_var) {
  group_level <- if (identical(group_var, "overall")) {
    rep("Overall", nrow(df))
  } else {
    as.character(df[[group_var]])
  }

  bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
    spec <- pollutant_specs[i, ]
    width <- pollutant_bin_widths[[spec$pollutant]]
    if (!is.finite(width) || width <= 0) {
      stop("Invalid exposure rug bin width for ", spec$pollutant)
    }

    df %>%
      transmute(
        site = site_name,
        group_var = group_var,
        group_level = group_level,
        pollutant = spec$pollutant,
        pollutant_label = spec$display_label,
        bin_width = width,
        exposure_raw = as.numeric(.data[[spec$raw_col]])
      ) %>%
      filter(!is.na(group_level), is.finite(exposure_raw)) %>%
      mutate(
        bin_lower = floor(exposure_raw / bin_width) * bin_width,
        bin_upper = bin_lower + bin_width,
        bin_mid = bin_lower + bin_width / 2
      ) %>%
      count(site, group_var, group_level, pollutant, pollutant_label, bin_width, bin_lower, bin_upper, bin_mid, name = "n")
  }))
}

exposure_distribution_bins <- bind_rows(
  make_exposure_distribution_bins(analysis_df, "overall"),
  make_exposure_distribution_bins(analysis_df, "arf_subtype"),
  make_exposure_distribution_bins(analysis_df, "sex"),
  make_exposure_distribution_bins(analysis_df, "race_ethnicity")
)

out_path <- file.path(site_dir, "primary_exposure_distribution_bins.csv")
readr::write_csv(exposure_distribution_bins, out_path)

message("Wrote: ", out_path)
message("Rows: ", nrow(exposure_distribution_bins))
