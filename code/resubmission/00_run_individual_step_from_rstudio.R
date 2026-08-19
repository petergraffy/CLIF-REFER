# Open this file in RStudio, edit STEP_TO_RUN if needed, and click Source.
# This lets sites run one resubmission step at a time without using
# PowerShell, Command Prompt, or Terminal.

# Choose one:
#   "00_setup"
#   "01_build_cohort"
#   "01_primary_models"
#   "02_competing_risks"
#   "03_exposure_response"
#   "04_subgroup_interactions"
#   "05_all_clif_mortality"
#   "06_sensitivity_models"
#   "07_inclusion_flow"
#   "all_manual"
STEP_TO_RUN <- "all_manual"

# Reuse this same folder across steps. For a fresh manual run, change the final
# folder name, then Source this file again.
RUN_DIR <- file.path("output", "resubmission", "manual_run")

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
        file.exists(file.path(candidate, "code", "resubmission", "00_run_full_resubmission_pipeline.R"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    "Could not find the repository root. Open this file from the cloned ",
    "CLIF-REFER repository, or set the RStudio working directory to the repository root."
  )
}

rscript_path <- function() {
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  candidates <- c(
    file.path(R.home("bin"), exe),
    file.path(R.home("bin"), "x64", exe)
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits)) return(normalizePath(hits[[1]], mustWork = TRUE))
  stop("Could not find Rscript for this R installation. R.home() is: ", R.home())
}

run_child_script <- function(label, script, args = character()) {
  rscript <- rscript_path()
  command <- paste(c(rscript, "--vanilla", script, args), collapse = " ")

  message("\n== ", label, " ==")
  message("Running: ", command)

  status <- system2(
    rscript,
    args = c("--vanilla", script, args),
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L)) {
    stop(label, " failed with exit status ", status, ". See the console output above.")
  }

  invisible(status)
}

repo <- find_repo_root()
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo)

run_dir <- normalizePath(RUN_DIR, mustWork = FALSE)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

raw_dataset <- file.path(run_dir, "resubmission_analysis_dataset.csv")
no_icu_dataset <- file.path(run_dir, "resubmission_analysis_dataset_no_icu_los_restriction.csv")
analysis_dataset <- file.path(run_dir, "analysis_dataset_reviewer_optimized.csv")

scripts <- list(
  setup = file.path(repo, "code", "resubmission", "00_setup_renv.R"),
  build_cohort = file.path(repo, "code", "resubmission", "01_zcta_arf_onset_cause_specific_cox.R"),
  primary = file.path(repo, "code", "resubmission", "01_primary_reviewer_optimized_models.R"),
  competing = file.path(repo, "code", "resubmission", "02_unadjusted_aj_and_fine_gray.R"),
  exposure = file.path(repo, "code", "resubmission", "03_exposure_response_primary_models.R"),
  subgroups = file.path(repo, "code", "resubmission", "04_subgroup_interaction_estimates.R"),
  all_clif = file.path(repo, "code", "resubmission", "05_clif_site_mortality_distribution_export.R"),
  sensitivities = file.path(repo, "code", "resubmission", "06_primary_sensitivity_models.R"),
  inclusion_flow = file.path(repo, "code", "resubmission", "07_export_inclusion_flow_counts.R")
)

run_step <- function(step_name) {
  switch(
    step_name,
    "00_setup" = run_child_script("Restore R package environment", scripts$setup),
    "01_build_cohort" = run_child_script("Build ARF cohorts", scripts$build_cohort, run_dir),
    "01_primary_models" = run_child_script("Primary models and Table 1", scripts$primary, c(raw_dataset, run_dir, analysis_dataset)),
    "02_competing_risks" = run_child_script("AJ and Fine-Gray competing risks", scripts$competing, c(analysis_dataset, run_dir)),
    "03_exposure_response" = run_child_script("Exposure-response figures", scripts$exposure, c(analysis_dataset, run_dir)),
    "04_subgroup_interactions" = run_child_script("Subgroup estimates and interactions", scripts$subgroups, c(analysis_dataset, run_dir)),
    "05_all_clif_mortality" = run_child_script("All-CLIF raw mortality supplement", scripts$all_clif, c(repo, run_dir)),
    "06_sensitivity_models" = run_child_script("Sensitivity models", scripts$sensitivities, c(analysis_dataset, run_dir, no_icu_dataset)),
    "07_inclusion_flow" = run_child_script("Inclusion-flow count export", scripts$inclusion_flow, run_dir),
    stop("Unknown STEP_TO_RUN: ", step_name)
  )
}

message("Repository: ", repo)
message("Manual run folder: ", run_dir)

if (identical(STEP_TO_RUN, "all_manual")) {
  for (step_name in c(
    "00_setup",
    "01_build_cohort",
    "01_primary_models",
    "02_competing_risks",
    "03_exposure_response",
    "04_subgroup_interactions",
    "05_all_clif_mortality",
    "06_sensitivity_models",
    "07_inclusion_flow"
  )) {
    run_step(step_name)
  }
} else {
  run_step(STEP_TO_RUN)
}

message("\nRequested resubmission step(s) finished.")
