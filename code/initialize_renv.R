# Setup R environment using renv
# Install renv if not already installed:
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# Initialize renv for the project:
renv::init()

# All required packages (deduped)
pkgs <- c(
  "knitr","fst","here","tidyverse","arrow","gtsummary",
  "dplyr","tidyr","stringr","lubridate","purrr","fuzzyjoin","janitor",
  "MASS","broom","patchwork","ggplot2","ggeffects","gt",
  "rlang","data.table","readr","glue","scales",
  "DiagrammeR","DiagrammeRsvg","rsvg", "tidycensus", "fixest", 
  "marginaleffects", "pscl", "glmmTMB", "digest", "pROC", 
  "tibble", "forcats", "grid", 'jsonlite', 'cmprsk'
)
pkgs <- unique(pkgs)

# Install and lock
renv::install(pkgs)
renv::settings$snapshot.type("all")
renv::snapshot(prompt = FALSE) # Then run snapshot again

