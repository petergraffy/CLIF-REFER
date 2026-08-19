# REFER Resubmission Primary Analysis

This folder is a cleaner, reviewer-focused replacement for the original first-run model block in `code/02_REFER_linkage_analysis.R`.

The goal is not to expand the analysis in every direction. It implements only the main changes implied by the reviewer feedback:

1. Use mortality by day 28 after ARF onset as the primary binary mortality estimand.
2. Keep Cox proportional hazards models for mortality after ARF onset as a sensitivity analysis.
3. Use ventilator-free days through day 28 as the primary ventilation-duration estimand, while retaining raw invasive mechanical ventilation duration as a supplemental site export.
4. Use the stronger covariate set now available: age, sex, race/ethnicity, calendar year, Charlson score, and available ACS/social vulnerability covariates.
5. Add a COVID-era sensitivity analysis excluding the 12-month period most affected by early pandemic care disruptions.

## Main Script

`00_setup_renv.R`

Sites should restore the R package environment before running the analysis:

```sh
Rscript --vanilla code/resubmission/00_setup_renv.R
```

On Windows, this script bootstraps `renv` into a writable user library, uses
Posit Package Manager for CRAN binaries, and disables silent source-build
fallbacks. If a site must use an institutional CRAN mirror, set
`REFER_RENV_CRAN_REPO` before running the setup script.

By default, `00_setup_renv.R` skips `duckdb` during `renv::restore()` because
the resubmission pipeline does not use it and several Windows systems try to
compile it from source. Sites should pull the latest repository version if they
see a `duckdb` compilation error. The skip list can be changed with
`REFER_RENV_EXCLUDE_PACKAGES`, but sites should not need to do this for the
REFER resubmission run.

## Windows Command-Line Instructions

Use PowerShell or Command Prompt from the cloned repository root. Do not run the
files line-by-line from the R console, and do not double-click the scripts. The
working directory must be the folder that contains `code`, `config`, `exposome`,
and `renv.lock`.

First, open PowerShell and move to the repository folder:

```powershell
cd "C:\Users\<username>\Projects\CLIF-REFER"
```

Confirm that PowerShell is in the correct place:

```powershell
Test-Path "code\resubmission\00_setup_renv.R"
Test-Path "code\resubmission\00_run_full_resubmission_pipeline.R"
Test-Path "renv.lock"
Test-Path "exposome\zcta\zcta_acs_community_covariates_2005_2023.csv.gz"
```

Each command should print `True`. If any command prints `False`, the current
folder is not the cloned repository root or the repository is incomplete.

If `Rscript` works from the command line, run:

```powershell
Rscript --vanilla "code\resubmission\00_setup_renv.R"
Rscript --vanilla "code\resubmission\00_run_full_resubmission_pipeline.R"
```

If Windows says `Rscript is not recognized`, R is installed but not on the
system `PATH`. Use the full path to `Rscript.exe` instead. Common locations are:

```powershell
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" --vanilla "code\resubmission\00_setup_renv.R"
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" --vanilla "code\resubmission\00_run_full_resubmission_pipeline.R"
```

or:

```powershell
& "C:\Users\<username>\AppData\Local\Programs\R\R-4.4.2\bin\x64\Rscript.exe" --vanilla "code\resubmission\00_setup_renv.R"
& "C:\Users\<username>\AppData\Local\Programs\R\R-4.4.2\bin\x64\Rscript.exe" --vanilla "code\resubmission\00_run_full_resubmission_pipeline.R"
```

Replace `R-4.4.2` with the R version installed at the site, if needed. To find
the installed `Rscript.exe`, run:

```powershell
Get-ChildItem -Path "C:\Program Files\R" -Filter Rscript.exe -Recurse -ErrorAction SilentlyContinue
Get-ChildItem -Path "$env:LOCALAPPDATA\Programs\R" -Filter Rscript.exe -Recurse -ErrorAction SilentlyContinue
```

The setup step should end with `renv restore complete.` The full pipeline should
end with output written under:

```text
output\resubmission\<timestamp>\
```

During the full pipeline, each step writes stdout and stderr logs under:

```text
output\resubmission\<timestamp>\logs\
```

If a step fails, the runner prints the last lines from those logs and reports
the log folder. This is usually the fastest way to see whether the issue is a
missing CLIF table, an unexpected column name, a missing exposome file, or an
incorrect `config\config.json` path.

The two most common Windows problems are:

- `Rscript is not recognized`: use the full path to `Rscript.exe` as shown
  above, or add R to the Windows `PATH`.
- Scripts use the wrong directories: reopen PowerShell, `cd` to the cloned
  repository root, confirm the `Test-Path` checks above, and rerun the two
  commands from there.

If the setup script fails while installing a package, send the complete console
log beginning with the `00_setup_renv.R` command. If the pipeline fails after
setup succeeds, send the complete console log beginning with the
`00_run_full_resubmission_pipeline.R` command plus the two log files for the
failed step from `output\resubmission\<timestamp>\logs\`.

`00_run_full_resubmission_pipeline.R`

This is the site-facing one-command runner. With no arguments, it first builds
the ARF cohorts from raw CLIF/ZCTA data using the resubmission cohort builder,
then runs every primary, secondary, sensitivity, figure, diagnostic, and
aggregate export into one timestamped folder:

```sh
Rscript code/resubmission/00_run_full_resubmission_pipeline.R
```

For local development only, reruns can start from an existing row-level cohort
export. These files are PHI-bearing working files and should not be uploaded
for pooling:

```sh
Rscript code/resubmission/00_run_full_resubmission_pipeline.R \
  output/resubmission/<run>/resubmission_analysis_dataset.csv \
  output/resubmission/<run>/resubmission_analysis_dataset_no_icu_los_restriction.csv
```

## Required Exposome Files

The resubmission cohort builder uses repo-local ZCTA exposure files by default.
Sites should receive these files in `exposome/zcta/`:

- `air_pollution_zcta_pm25_monthly_2005_2023.parquet`
- `air_pollution_zcta_o3_monthly_2005_2023.parquet`
- `air_pollution_zcta_no2_annual_2005_2025.parquet`
- `zcta_acs_community_covariates_2005_2023.csv.gz`
- `zcta_acs_community_covariates_dictionary.csv`
- `zcta_acs_community_covariates_qc.csv`
- `no2_monthly/no2_zcta_monthly_2019.parquet` through `no2_monthly/no2_zcta_monthly_2025.parquet`

If a site needs to override these paths, copy
`resubmission_config_template.json` to `resubmission_config.json` and edit the
corresponding `zcta_*` fields.

The ACS file is pre-downloaded and analysis-ready. It includes the ZCTA-level
social covariates used in the models, so sites do not need Census API access to
run the resubmission analysis.

`01_primary_reviewer_optimized_models.R`

This script assumes the original cohort/linkage workflow has already constructed `arf_exp`. If available, it also uses:

- `support_class` to compute ventilator-free days from invasive mechanical ventilation records
- `diagnosis` to calculate Charlson score when `charlson_score` is not already present
- `patient` to recover death time when `arf_exp` does not already include one

The script writes outputs to:

`output/resubmission/<timestamp>/`

`06_primary_sensitivity_models.R`

This script reads the reviewer-optimized analysis dataset and exports poolable
site-level sensitivity model results for: adding first-24-hour SOFA total,
removing the ICU length-of-stay restriction when a companion dataset is
available, adding O3 as a copollutant, using a 36-month exposure window, and
using a 14-day follow-up window.

The 14-day VFD sensitivity requires exact `ventilator_free_days_day14` or
`imv_days_through_day14` fields. The primary dataset builder now creates these
when raw respiratory support data are available; otherwise the sensitivity
status file marks the VFD component as partial rather than imputing it.

## Primary Models

Primary mortality model:

```r
mortality_day28_event ~
  exposure +
  age_10 + sex + race_ethnicity + charlson_score + index_year_f +
  social_covariates
```

Cox sensitivity model:

```r
Surv(mortality_ftime_days, mortality_event) ~
  exposure +
  age_10 + sex + race_ethnicity + charlson_score + index_year_f +
  social_covariates
```

VFD model:

```r
ventilator_free_days ~
  exposure +
  age_10 + sex + race_ethnicity + charlson_score + index_year_f +
  social_covariates
```

Exposure models are run as:

- PM2.5 single-pollutant, scaled per 5 ug/m3
- NO2 single-pollutant, scaled per 10 ppb
- PM2.5 + NO2 multipollutant model

## COVID Sensitivity

The primary script also reruns the primary model set after excluding ARF onset dates from March 1, 2020 through February 28, 2021. Sites can override the window with environment variables:

```sh
REFER_COVID_EXCLUDE_START=2020-03-01
REFER_COVID_EXCLUDE_END=2021-02-28
```

This sensitivity exports site-level, poolable result tables only; no patient-level data are required for cross-site pooling.

## Output Files

- `cohort_summary_reviewer_optimized.csv`
- `primary_mortality_day28_logistic_results.csv`
- `primary_mortality_cox_results.csv`
- `primary_mortality_cox_ph_diagnostics.csv`
- `primary_vfd_quasipoisson_results.csv`
- `primary_vfd_model_diagnostics.csv`
- `primary_vfd_poisson_vs_quasipoisson_diagnostics.csv`
- `primary_vfd_calibration_by_fitted_decile.csv`
- `primary_imv_duration_quasipoisson_results.csv`
- `primary_imv_duration_model_diagnostics.csv`
- `primary_imv_duration_poisson_vs_quasipoisson_diagnostics.csv`
- `primary_imv_duration_calibration_by_fitted_decile.csv`
- `primary_exposure_distribution_bins.csv`
- `sensitivity_exclude_peak_covid_12m_cohort_summary.csv`
- `sensitivity_exclude_peak_covid_12m_mortality_day28_logistic_results.csv`
- `sensitivity_exclude_peak_covid_12m_mortality_cox_results.csv`
- `sensitivity_exclude_peak_covid_12m_vfd_quasipoisson_results.csv`
- `sensitivity_exclude_peak_covid_12m_imv_duration_quasipoisson_results.csv`
- `sensitivity_exclude_peak_covid_12m_pooling_table.csv`
- `primary_vs_exclude_peak_covid_12m_pooling_table.csv`
- `primary_sensitivity_models_pooling_table.csv`
- `primary_sensitivity_models_run_status.csv`
- `primary_sensitivity_models_cohort_summaries.csv`
- `primary_sensitivity_3y_exposure_diagnostics.csv`
- `site_inclusion_flow_counts.csv`
- `site_inclusion_flow_counts_wide.csv`
- `site_inclusion_flow_counts_dictionary.csv`

The pipeline uses row-level analysis datasets internally while it is running,
but it does not export them to the site output folder by default. Files such as
`analysis_dataset_reviewer_optimized.csv`, `resubmission_analysis_dataset.csv`,
`resubmission_analysis_dataset_no_icu_los_restriction.csv`, and
`readmission_analysis_dataset.csv` should be treated as PHI-bearing working
files and should not be uploaded for pooling.

`primary_exposure_distribution_bins.csv` is an aggregate, non-row-level export
used to draw pooled rug-style exposure distributions. By default, it uses
0.025 ppb NO2 bins and 0.01 ug/m3 PM2.5 bins; sites can override these with
`REFER_RUG_BIN_WIDTH_NO2` and `REFER_RUG_BIN_WIDTH_PM25`.

The VFD and supplemental IMV-duration diagnostic files are generated
automatically at each site when the primary script runs. They provide the
model-level dispersion diagnostics, Poisson versus quasi-Poisson exposure
inference, and fitted-value decile calibration summaries needed to justify the
quasi-Poisson mean model without sharing row-level data.

The site inclusion-flow files are generated as the final pipeline step. They
combine the all-CLIF denominator, ARF cohort construction, exposure/covariate
availability, primary complete-case counts, IMV competing-risk assignment, and
sensitivity cohort sizes into aggregate-only exports that can be pooled across
sites.
