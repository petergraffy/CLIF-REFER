 ## Code directory

This directory contains scripts for the project workflow. The general workflow consists of two main steps: cohort identification, and linkage/analysis.

### General Workflow

0. In the exposome folder, unzip the acs_estimates.csv.zip folder and save it in the exposome folder with the other csv files. Then initialize your R environment using `00_renv_restore`.

1. Run the `01_REFER_cohort_identification.R` script
   This script should:
   - Apply inclusion and exclusion criteria
   - Select required fields from each table
   - Filter tables to include only required observations
   
Then, **without making any changes to your R environment**, you can run the next step immediately after.

2. Run the `02_REFER_linkage_analysis.R` script
   This script should:
   - Perform project-specific quality control checks on the filtered cohort data.
   - Perform geospatial linkage to exposome features in the accompanied csv files.
   - Perform modeling and return outputs and figures.

If you accidentally cleared your environment between steps, you can run step 1 again to get the objects back in your environment, or you can re-read the saved csv output from file one in the output folder.

**Once all code has been run, be sure to delete cohort_inclusion_[site].csv, cohort_inclusion_periop_[site].csv, and exclusions_raw_[site] from the run_[site] folder in output as this has line-level patient data.**



