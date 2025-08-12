# Step 1: Define the repository URL and file paths
repo_url <- "https://github.com/08wparker/CLIF_cohort_identifier/archive/refs/heads/main.zip"
download_zip <- "CLIF_cohort_identifier.zip"
unzip_dir <- "CLIF_cohort_identifier-main"  # Adjust based on the folder structure after unzipping

# Step 2: Remove old zip file and unzipped directory if they exist
if (file.exists(download_zip)) {
  message("Deleting the old ZIP file...")
  file.remove(download_zip)
}

if (dir.exists(unzip_dir)) {
  message("Deleting the old unzipped directory...")
  unlink(unzip_dir, recursive = TRUE)
}

# Step 3: Download the repository ZIP (after deleting the old one)
message("Downloading the GitHub repository as a ZIP file...")
download.file(repo_url, destfile = download_zip)
message("Download complete.")

# Step 4: Unzip the repository
message("Unzipping the repository...")
unzip(download_zip)
message("Unzipping complete.")

# Step 5: Parse the `renv.lock` file to detect required packages from `CLIF_cohort_identifier`
renv_lock_file <- file.path(unzip_dir, "renv.lock")

if (!file.exists(renv_lock_file)) {
  stop("renv.lock file not found. Unable to detect required packages.")
}

message("Detecting required packages from renv.lock file...")

# Step 6: Extract packages from renv.lock file
renv_content <- jsonlite::fromJSON(renv_lock_file)
required_packages <- names(renv_content$Packages)

# Step 7: Check for missing packages in your current `renv` environment
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing_packages) > 0) {
  message("Installing missing packages into your local renv environment: ", paste(missing_packages, collapse = ", "))
  
  # Use `renv::install()` to install the missing packages into the current renv environment
  renv::install(missing_packages)
} else {
  message("All required packages are already installed.")
}

# Step 8: Load the required packages into the session from the current renv environment
lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Loading package: ", pkg)
    library(pkg, character.only = TRUE)
  }
})

# Step 9: Run the Shiny app from the unzipped folder
app_dir <- file.path(unzip_dir, "CLIF_cohort_ID_application")  # Update to the path of your app
if (dir.exists(app_dir)) {
  message("Running the Shiny app...")
  shiny::runApp(app_dir)
} else {
  stop("App directory not found. Please check the path.")
}
