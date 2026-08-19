# Open this file in RStudio and click Source to run the full REFER
# resubmission pipeline without using PowerShell or Command Prompt.

# Sites usually should not edit these defaults.
RUN_SETUP <- TRUE
PRIMARY_DATASET <- ""
NO_ICU_LOS_DATASET <- ""
OUTPUT_DIR <- ""

find_this_file <- function() {
  path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
  if (!is.na(path)) return(path)

  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(cmd_file)) {
    path <- tryCatch(normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = TRUE), error = function(e) NA_character_)
    if (!is.na(path)) return(path)
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NA_character_)
    if (!is.na(path) && nzchar(path) && file.exists(path)) return(normalizePath(path, mustWork = TRUE))
  }

  NA_character_
}

find_repo_root <- function() {
  this_file <- find_this_file()
  candidates <- character()
  if (!is.na(this_file)) {
    candidates <- c(candidates, file.path(dirname(this_file), "..", ".."))
  }
  candidates <- c(candidates, getwd())

  for (candidate in candidates) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (file.exists(file.path(candidate, "renv.lock")) &&
        file.exists(file.path(candidate, "code", "resubmission", "00_run_full_resubmission_pipeline.R"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    "Could not find the repository root. In RStudio, open this file from the cloned ",
    "CLIF-REFER repository, or set the RStudio working directory to the repository root."
  )
}

rscript_path <- function() {
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  path <- file.path(R.home("bin"), exe)
  if (file.exists(path)) return(normalizePath(path, mustWork = TRUE))

  path <- file.path(R.home("bin"), "x64", exe)
  if (file.exists(path)) return(normalizePath(path, mustWork = TRUE))

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

message("Repository: ", repo)
message("R executable: ", R.home())

setup_script <- file.path(repo, "code", "resubmission", "00_setup_renv.R")
pipeline_script <- file.path(repo, "code", "resubmission", "00_run_full_resubmission_pipeline.R")

if (isTRUE(RUN_SETUP)) {
  run_child_script("Restore R package environment", setup_script)
}

pipeline_args <- character()
if (nzchar(PRIMARY_DATASET)) {
  pipeline_args <- c(pipeline_args, normalizePath(PRIMARY_DATASET, mustWork = TRUE))
  if (nzchar(NO_ICU_LOS_DATASET)) {
    pipeline_args <- c(pipeline_args, normalizePath(NO_ICU_LOS_DATASET, mustWork = TRUE))
  }
  if (nzchar(OUTPUT_DIR)) {
    pipeline_args <- c(pipeline_args, normalizePath(OUTPUT_DIR, mustWork = FALSE))
  }
} else if (nzchar(OUTPUT_DIR)) {
  stop("OUTPUT_DIR can only be used here when PRIMARY_DATASET is also supplied.")
}

run_child_script("Run full REFER resubmission pipeline", pipeline_script, pipeline_args)

message("\nREFER resubmission pipeline finished.")
