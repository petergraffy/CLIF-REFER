# Setup R environment using renv
# Install renv if not already installed:
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# Initialize renv for the project:
renv::init()
# Install required packages:
renv::install(c("knitr", "here", "tidyverse", "arrow", "gtsummary"))
# Save the project's package state:
renv::snapshot()