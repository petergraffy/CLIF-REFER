# 00_renv_restore.R
# Ensure renv is installed and load the project environment

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Restore the project's package environment
renv::activate()
renv::restore()
