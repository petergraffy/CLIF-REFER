#!/usr/bin/env Rscript

# Site-facing dependency setup for the REFER resubmission analysis.
#
# Recommended use from the repository root:
#   Rscript --vanilla code/resubmission/00_setup_renv.R

default_user_library <- function() {
  r_version <- paste(R.version$major, sub("\\..*", "", R.version$minor), sep = ".")
  if (.Platform$OS.type == "windows") {
    local_appdata <- Sys.getenv("LOCALAPPDATA", unset = "")
    if (nzchar(local_appdata)) {
      return(file.path(local_appdata, "R", "win-library", r_version))
    }
  }
  Sys.getenv("R_LIBS_USER", unset = file.path(path.expand("~"), "R", paste0("library-", r_version)))
}

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

user_library <- default_user_library()
dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(user_library, .libPaths())))

cran_repo <- Sys.getenv(
  "REFER_RENV_CRAN_REPO",
  unset = "https://packagemanager.posit.co/cran/latest"
)
options(repos = c(CRAN = cran_repo))
options(renv.config.ppm.enabled = TRUE)

if (.Platform$OS.type == "windows") {
  options(pkgType = "binary")
  options(install.packages.compile.from.source = "never")
}

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", lib = user_library)
}
stopifnot(requireNamespace("renv", quietly = TRUE))

project_library <- renv::paths$library(project = repo)
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_library, .libPaths()))

message("Restoring R package environment from: ", lockfile)
message("Bootstrap user library: ", user_library)
message("Project library: ", project_library)
message("CRAN repository: ", getOption("repos")[["CRAN"]])
exclude_packages <- strsplit(Sys.getenv("REFER_RENV_EXCLUDE_PACKAGES", "duckdb"), "[,; ]+")[[1]]
exclude_packages <- exclude_packages[nzchar(exclude_packages)]
if (length(exclude_packages)) {
  message("Skipping packages not required by the resubmission pipeline: ", paste(exclude_packages, collapse = ", "))
}
renv::restore(
  project = repo,
  library = project_library,
  lockfile = lockfile,
  exclude = exclude_packages,
  prompt = FALSE
)

message("renv restore complete.")
