#!/usr/bin/env Rscript

# Site-facing dependency setup for the REFER resubmission analysis.
#
# Recommended use from the repository root:
#   Rscript --vanilla code/resubmission/00_setup_renv.R

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(cmd_file)) normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = TRUE) else NA_character_
}

repo <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

lockfile <- file.path(repo, "renv.lock")
if (!file.exists(lockfile)) {
  stop("Could not find renv.lock at repository root: ", lockfile)
}

if (is.null(getOption("repos")) || identical(getOption("repos")[["CRAN"]], "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

project_library <- renv::paths$library(project = repo)
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_library, .libPaths()))

message("Restoring R package environment from: ", lockfile)
message("Project library: ", project_library)
renv::restore(project = repo, library = project_library, lockfile = lockfile, prompt = FALSE)

message("renv restore complete.")
