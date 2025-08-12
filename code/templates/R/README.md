## Instructions for initializing your project environment with R

1. Copy `initialize_renv_template.R` and rename to `initialize_renv.R`
2. Run the script to install `renv` and create a project snapshot which will be stored in an `renv` folder in the project
3. run `renv::snapshot()` at the end of code development prior to distributing the project to be run across the consortium to ensure the most up to date packages are included in the environment
4. DELETE this file from your project repo
