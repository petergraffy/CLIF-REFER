# REFER resubmission all-in-one RStudio runner.
#
# Open this file in RStudio and click Source. This script runs the complete
# resubmission pipeline in one R session, without requiring PowerShell,
# Command Prompt, Terminal, or child Rscript processes.
#
# Configuration comes from code/resubmission/resubmission_config.json if present,
# otherwise config/config.json. Row-level intermediate datasets are written only
# to a private temporary working folder so they are not included in the site
# output folder for pooling.

find_this_file <- function() {
  path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
  if (!is.na(path)) return(path)

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NA_character_)
    if (!is.na(path) && nzchar(path) && file.exists(path)) return(normalizePath(path, mustWork = TRUE))
  }

  NA_character_
}

find_repo_root <- function() {
  this_file <- find_this_file()
  candidates <- c(getwd())
  if (!is.na(this_file)) {
    candidates <- c(file.path(dirname(this_file), "..", ".."), candidates)
  }

  for (candidate in candidates) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (file.exists(file.path(candidate, "renv.lock")) &&
        file.exists(file.path(candidate, "code", "resubmission", "01_zcta_arf_onset_cause_specific_cox.R"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    "Could not find the repository root. Open this file from the cloned ",
    "CLIF-REFER repository, or set the RStudio working directory to the repository root."
  )
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

repo <- find_repo_root()
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  message("Package `jsonlite` is not available yet; running setup before reading config.")
  setup_script <- file.path(repo, "code", "resubmission", "00_setup_renv.R")
  setup_env <- new.env(parent = globalenv())
  setup_env$.refer_script_path <- normalizePath(setup_script, mustWork = TRUE)
  setup_env$commandArgs <- function(trailingOnly = FALSE) {
    if (isTRUE(trailingOnly)) character() else c(base::commandArgs(FALSE), paste0("--file=", setup_env$.refer_script_path))
  }
  source(setup_script, local = setup_env, echo = FALSE, chdir = FALSE)
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package `jsonlite` is required but was not available after setup.")
}

config_path <- file.path(repo, "code", "resubmission", "resubmission_config.json")
if (!file.exists(config_path)) {
  config_path <- file.path(repo, "config", "config.json")
}
if (!file.exists(config_path)) {
  stop("Could not find code/resubmission/resubmission_config.json or config/config.json.")
}

config <- jsonlite::fromJSON(config_path)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_base <- config$output_dir %||% file.path("output", "resubmission")
out_dir <- if (grepl("^(/|[A-Za-z]:[/\\\\]|~)", output_base)) {
  file.path(path.expand(output_base), stamp)
} else {
  file.path(repo, output_base, stamp)
}
out_dir <- normalizePath(out_dir, mustWork = FALSE)
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")
private_work_dir <- file.path(tempdir(), paste0("refer_resubmission_private_", stamp))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_work_dir, recursive = TRUE, showWarnings = FALSE)

message("Repository: ", repo)
message("Config: ", config_path)
message("Unified output directory: ", out_dir)
message("Private row-level working directory: ", private_work_dir)

pipeline_env <- new.env(parent = globalenv())
pipeline_env$commandArgs <- function(trailingOnly = FALSE) {
  if (isTRUE(trailingOnly)) {
    return(get(".refer_script_args", envir = pipeline_env, inherits = FALSE))
  }
  c(base::commandArgs(FALSE), paste0("--file=", get(".refer_script_path", envir = pipeline_env, inherits = FALSE)))
}

step_log_name <- function(step) {
  x <- tolower(step)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  paste0(format(Sys.time(), "%H%M%S"), "_", x, ".log")
}

run_sourced_step <- function(step, script, args = character()) {
  script <- normalizePath(script, mustWork = TRUE)
  log_file <- file.path(log_dir, step_log_name(step))

  message("\n== ", step, " ==")
  message("Script: ", script)
  message("Step log: ", log_file)

  writeLines(
    c(
      paste0("Step: ", step),
      paste0("Started at: ", Sys.time()),
      paste0("Repository: ", repo),
      paste0("Working directory: ", getwd()),
      paste0("Script: ", script),
      paste0("Arguments: ", if (length(args)) paste(args, collapse = " | ") else "none"),
      ""
    ),
    log_file,
    useBytes = TRUE
  )

  pipeline_env$.refer_script_args <- args
  pipeline_env$.refer_script_path <- script

  start_time <- Sys.time()
  ok <- FALSE
  sink(log_file, append = TRUE, split = TRUE)
  sink(log_file, append = TRUE, type = "message")
  on.exit({
    while (sink.number(type = "message") > 0) sink(type = "message")
    while (sink.number() > 0) sink()
  }, add = TRUE)

  tryCatch(
    {
      source(script, local = pipeline_env, echo = FALSE, chdir = FALSE)
      ok <- TRUE
    },
    error = function(e) {
      message("\nERROR: ", conditionMessage(e))
      pipeline_env$.refer_last_error <- e
    }
  )

  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number() > 0) sink()

  write(
    c(
      "",
      paste0("Ended at: ", Sys.time()),
      paste0("Elapsed seconds: ", as.numeric(difftime(Sys.time(), start_time, units = "secs"))),
      paste0("Status: ", if (ok) "completed" else "failed")
    ),
    file = log_file,
    append = TRUE
  )

  if (!ok) {
    lines <- readLines(log_file, warn = FALSE)
    message("\n--- Last lines from step log ---")
    message(paste(utils::tail(lines, 80), collapse = "\n"))
    stop("Step failed: ", step, ". See step log: ", log_file)
  }

  data.frame(
    step = step,
    script = script,
    status = "completed",
    started_at = as.character(start_time),
    ended_at = as.character(Sys.time()),
    elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    step_log = normalizePath(log_file, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

resubmission_dir <- file.path(repo, "code", "resubmission")
raw_dataset <- file.path(private_work_dir, "resubmission_analysis_dataset.csv")
no_icu_dataset <- file.path(private_work_dir, "resubmission_analysis_dataset_no_icu_los_restriction.csv")
reviewer_dataset <- file.path(private_work_dir, "analysis_dataset_reviewer_optimized.csv")

step_log <- list()
step_log[[length(step_log) + 1]] <- run_sourced_step(
  "Restore R package environment",
  file.path(resubmission_dir, "00_setup_renv.R")
)
step_log[[length(step_log) + 1]] <- run_sourced_step(
  "Build ARF cohorts from raw CLIF/ZCTA data",
  file.path(resubmission_dir, "01_zcta_arf_onset_cause_specific_cox.R"),
  private_work_dir
)

for (filename in c(
  "resubmission_cohort_summary.csv",
  "resubmission_cohort_summary_no_icu_los_restriction.csv",
  "resubmission_cohort_summary_no_peak_covid.csv"
)) {
  src <- file.path(private_work_dir, filename)
  if (file.exists(src)) {
    file.copy(src, file.path(out_dir, filename), overwrite = TRUE)
  }
}

if (!file.exists(raw_dataset)) {
  stop("Primary analysis dataset was not created: ", raw_dataset)
}

step_log[[length(step_log) + 1]] <- run_sourced_step(
  "Primary day-28 mortality, Cox sensitivity, VFD, Table 1, COVID sensitivity",
  file.path(resubmission_dir, "01_primary_reviewer_optimized_models.R"),
  c(raw_dataset, out_dir, reviewer_dataset)
)
if (!file.exists(reviewer_dataset)) {
  stop("Internal reviewer-optimized dataset was not created: ", reviewer_dataset)
}

step_log[[length(step_log) + 1]] <- run_sourced_step(
  "Secondary competing-risk AJ and Fine-Gray analyses",
  file.path(resubmission_dir, "02_unadjusted_aj_and_fine_gray.R"),
  c(reviewer_dataset, out_dir)
)
step_log[[length(step_log) + 1]] <- run_sourced_step(
  "Primary exposure-response figures and prediction exports",
  file.path(resubmission_dir, "03_exposure_response_primary_models.R"),
  c(reviewer_dataset, out_dir)
)
step_log[[length(step_log) + 1]] <- run_sourced_step(
  "Subgroup estimates and interaction tests",
  file.path(resubmission_dir, "04_subgroup_interaction_estimates.R"),
  c(reviewer_dataset, out_dir)
)
step_log[[length(step_log) + 1]] <- run_sourced_step(
  "All-CLIF raw mortality aggregate supplement",
  file.path(resubmission_dir, "05_clif_site_mortality_distribution_export.R"),
  c(repo, out_dir)
)
step_log[[length(step_log) + 1]] <- run_sourced_step(
  "Primary sensitivity model suite",
  file.path(resubmission_dir, "06_primary_sensitivity_models.R"),
  c(reviewer_dataset, out_dir, no_icu_dataset)
)
step_log[[length(step_log) + 1]] <- run_sourced_step(
  "Site inclusion-flow count export",
  file.path(resubmission_dir, "07_export_inclusion_flow_counts.R"),
  out_dir
)

pipeline_manifest <- data.frame(
  repo = repo,
  config_path = normalizePath(config_path, mustWork = FALSE),
  output_dir = normalizePath(out_dir, mustWork = FALSE),
  primary_input_dataset = "temporary_private_working_file_not_exported",
  reviewer_optimized_dataset = "temporary_private_working_file_not_exported",
  no_icu_los_dataset = "temporary_private_working_file_not_exported",
  private_row_level_outputs = "not_exported_to_site_output",
  generated_at = as.character(Sys.time()),
  runner = "00_run_all_in_one_from_rstudio.R",
  stringsAsFactors = FALSE
)

readr::write_csv(dplyr::bind_rows(step_log), file.path(out_dir, "resubmission_pipeline_step_log.csv"))
readr::write_csv(pipeline_manifest, file.path(out_dir, "resubmission_pipeline_manifest.csv"))

message("\nPipeline complete.")
message("Unified output directory: ", normalizePath(out_dir, mustWork = FALSE))
message("Manifest: ", file.path(out_dir, "resubmission_pipeline_manifest.csv"))
message("Step log: ", file.path(out_dir, "resubmission_pipeline_step_log.csv"))
