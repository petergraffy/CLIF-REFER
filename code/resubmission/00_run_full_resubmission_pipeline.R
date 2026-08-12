#!/usr/bin/env Rscript

# One-command site runner for the REFER resubmission analysis.
#
# Usage:
#   Rscript code/resubmission/00_run_full_resubmission_pipeline.R
#   Rscript code/resubmission/00_run_full_resubmission_pipeline.R <primary_dataset.csv> [no_icu_los_dataset.csv] [output_dir]
#
# With no primary dataset, this runs the raw CLIF/ZCTA cohort builder first.
# With a primary dataset, it reuses that dataset and only runs downstream
# primary, secondary, sensitivity, figure, and aggregate-export steps.

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
arg_or <- function(i, default = NA_character_) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(cmd_file)) normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = TRUE) else NA_character_
}
repo <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
primary_dataset_arg <- arg_or(1)
no_icu_dataset_arg <- arg_or(2)
resubmission_dir <- file.path(repo, "code", "resubmission")
out_dir <- arg_or(3, file.path(repo, "output", "resubmission", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
private_work_dir <- file.path(tempdir(), paste0("refer_resubmission_private_", stamp))
dir.create(private_work_dir, recursive = TRUE, showWarnings = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")

run_step <- function(step, script, script_args = character(), env = character()) {
  message("\n== ", step, " ==")
  message("Running: ", paste(c(rscript, "--vanilla", script, script_args), collapse = " "))
  start_time <- Sys.time()
  status <- system2(
    rscript,
    args = c("--vanilla", script, script_args),
    env = env,
    stdout = "",
    stderr = ""
  )
  end_time <- Sys.time()
  if (!identical(status, 0L)) {
    stop("Step failed: ", step, " (exit status ", status, ")")
  }
  tibble(
    step = step,
    script = script,
    status = "completed",
    started_at = as.character(start_time),
    ended_at = as.character(end_time),
    elapsed_seconds = as.numeric(difftime(end_time, start_time, units = "secs"))
  )
}

find_existing <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (length(hits)) normalizePath(hits[[1]], mustWork = TRUE) else NA_character_
}

message("Repository: ", repo)
message("Unified output directory: ", out_dir)
message("Private row-level working directory: ", private_work_dir)

step_log <- list()

if (is.na(primary_dataset_arg) || !file.exists(primary_dataset_arg)) {
  cohort_builder <- file.path(repo, "code", "resubmission", "01_zcta_arf_onset_cause_specific_cox.R")
  step_log[[length(step_log) + 1]] <- run_step(
    "Build ARF cohorts from raw CLIF/ZCTA data",
    cohort_builder,
    env = c(paste0("REFER_RESUBMISSION_OUTPUT_DIR=", normalizePath(private_work_dir, mustWork = FALSE)))
  )
  primary_dataset <- file.path(private_work_dir, "resubmission_analysis_dataset.csv")
  no_icu_dataset <- file.path(private_work_dir, "resubmission_analysis_dataset_no_icu_los_restriction.csv")

  aggregate_builder_outputs <- c(
    "resubmission_cohort_summary.csv",
    "resubmission_cohort_summary_no_icu_los_restriction.csv",
    "resubmission_cohort_summary_no_peak_covid.csv"
  )
  for (filename in aggregate_builder_outputs) {
    src <- file.path(private_work_dir, filename)
    if (file.exists(src)) {
      file.copy(src, file.path(out_dir, filename), overwrite = TRUE)
    }
  }
} else {
  primary_dataset <- normalizePath(primary_dataset_arg, mustWork = TRUE)
  no_icu_dataset <- if (!is.na(no_icu_dataset_arg) && file.exists(no_icu_dataset_arg)) {
    normalizePath(no_icu_dataset_arg, mustWork = TRUE)
  } else {
    find_existing(c(
      file.path(dirname(primary_dataset), "resubmission_analysis_dataset_no_icu_los_restriction.csv"),
      file.path(dirname(primary_dataset), "analysis_dataset_no_icu_los_restriction.csv"),
      file.path(dirname(dirname(primary_dataset)), "resubmission_analysis_dataset_no_icu_los_restriction.csv")
    ))
  }
}

if (!file.exists(primary_dataset)) {
  stop("Primary analysis dataset was not found: ", primary_dataset)
}

primary_script <- file.path(resubmission_dir, "01_primary_reviewer_optimized_models.R")
reviewer_dataset <- file.path(private_work_dir, "analysis_dataset_reviewer_optimized.csv")
step_log[[length(step_log) + 1]] <- run_step(
  "Primary day-28 mortality, Cox sensitivity, VFD, Table 1, COVID sensitivity",
  primary_script,
  c(primary_dataset, out_dir),
  env = c(
    paste0("REFER_INTERNAL_ANALYSIS_DATASET=", normalizePath(reviewer_dataset, mustWork = FALSE)),
    "REFER_EXPORT_ROW_LEVEL_DATASETS=false"
  )
)

if (!file.exists(reviewer_dataset)) {
  stop("Internal reviewer-optimized dataset was not created: ", reviewer_dataset)
}

analysis_steps <- list(
  list(
    step = "Secondary competing-risk AJ and Fine-Gray analyses",
    script = file.path(resubmission_dir, "02_unadjusted_aj_and_fine_gray.R"),
    args = c(reviewer_dataset, out_dir)
  ),
  list(
    step = "Primary exposure-response figures and prediction exports",
    script = file.path(resubmission_dir, "03_exposure_response_primary_models.R"),
    args = c(reviewer_dataset, out_dir)
  ),
  list(
    step = "Subgroup estimates and interaction tests",
    script = file.path(resubmission_dir, "04_subgroup_interaction_estimates.R"),
    args = c(reviewer_dataset, out_dir)
  ),
  list(
    step = "All-CLIF raw mortality aggregate supplement",
    script = file.path(resubmission_dir, "05_clif_site_mortality_distribution_export.R"),
    args = c(repo, out_dir)
  )
)

for (x in analysis_steps) {
  step_log[[length(step_log) + 1]] <- run_step(x$step, x$script, x$args)
}

sensitivity_args <- c(reviewer_dataset, out_dir)
if (!is.na(no_icu_dataset) && file.exists(no_icu_dataset)) {
  sensitivity_args <- c(sensitivity_args, no_icu_dataset)
} else {
  warning("No no-ICU-LOS-restriction dataset found; that sensitivity will be marked skipped.")
}

step_log[[length(step_log) + 1]] <- run_step(
  "Primary sensitivity model suite",
  file.path(resubmission_dir, "06_primary_sensitivity_models.R"),
  sensitivity_args
)

step_log[[length(step_log) + 1]] <- run_step(
  "Site inclusion-flow count export",
  file.path(resubmission_dir, "07_export_inclusion_flow_counts.R"),
  c(out_dir)
)

pipeline_manifest <- tibble(
  repo = repo,
  output_dir = normalizePath(out_dir, mustWork = FALSE),
  primary_input_dataset = if (!is.na(primary_dataset_arg) && file.exists(primary_dataset_arg)) primary_dataset else "temporary_private_working_file_not_exported",
  reviewer_optimized_dataset = "temporary_private_working_file_not_exported",
  no_icu_los_dataset = if (!is.na(primary_dataset_arg) && !is.na(no_icu_dataset) && file.exists(no_icu_dataset)) {
    no_icu_dataset
  } else if (!is.na(no_icu_dataset) && file.exists(no_icu_dataset)) {
    "temporary_private_working_file_not_exported"
  } else {
    NA_character_
  },
  private_row_level_outputs = "not_exported_to_site_output",
  generated_at = as.character(Sys.time())
)

readr::write_csv(bind_rows(step_log), file.path(out_dir, "resubmission_pipeline_step_log.csv"))
readr::write_csv(pipeline_manifest, file.path(out_dir, "resubmission_pipeline_manifest.csv"))

message("\nPipeline complete.")
message("Unified output directory: ", normalizePath(out_dir, mustWork = FALSE))
message("Manifest: ", file.path(out_dir, "resubmission_pipeline_manifest.csv"))
message("Step log: ", file.path(out_dir, "resubmission_pipeline_step_log.csv"))
