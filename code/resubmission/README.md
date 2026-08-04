# Resubmission Analysis

This folder contains the reviewer-responsive analysis framework for the REFER
air-pollution resubmission.

## Goal

The analysis changes the estimand from the original county-level cumulative
exposure/outcome regression workflow to:

- fixed pre-ARF exposure interval
- ZCTA-level PM2.5 and NO2 exposure assignment
- ARF onset as time zero
- primary cohort requiring ICU length of stay of at least 24 hours
- cause-specific Cox proportional hazards models
- baseline comorbidity adjustment with Charlson score
- ZCTA-level ACS social vulnerability adjustment
- unadjusted Aalen-Johansen cumulative incidence curves for death,
  successful extubation, and persistent respiratory failure

## Inputs

The script expects the same local CLIF tables used by the main project plus
ZCTA air-pollution parquet files.

Copy `resubmission_config_template.json` to `resubmission_config.json` and edit
the paths. The ZCTA exposure files can come from the transplant/pollution release:

- `air_pollution_zcta_pm25_monthly_2005_2023.parquet`
- `no2_zcta_monthly_YYYY.parquet` files from
  <https://github.com/petergraffy/environment_transplant_survival/releases/tag/no2-zcta-monthly-v1>
- `zcta_acs_community_covariates_2005_2023.csv.gz` from the transplant
  project community covariate build

PM2.5 and NO2 are averaged across the 12 complete months before ARF onset when
monthly data are available. The monthly NO2 release starts in 2019; encounters
without a complete 12-month pre-onset monthly NO2 window use the annual
prior-year ZCTA NO2 file as a fallback. If monthly NO2 files are not found, all
NO2 analyses use annual prior-year ZCTA NO2.

Primary cause-specific Cox models adjust for age, sex, race/ethnicity,
Charlson score, index year, and ZCTA ACS variables for poverty, unemployment,
no vehicle access, non-White population, median household income, and
bachelor-or-higher educational attainment. A separate sensitivity model also
adjusts for ARF subtype, and another sensitivity model adjusts for total SOFA
score calculated during the first 24 hours of ICU admission.

To address possible bias from excluding short ICU stays, the script also writes
a no-ICU-length-of-stay-restriction sensitivity cohort and repeats the primary
cause-specific Cox models in that secondary cohort.

## Run

From the repository root:

```sh
Rscript code/resubmission/01_zcta_arf_onset_cause_specific_cox.R
```

Outputs are written to `output/resubmission/<timestamp>/` by default:

- `resubmission_analysis_dataset.csv`
- `resubmission_analysis_dataset_no_icu_los_restriction.csv`
- `resubmission_cause_specific_cox_results.csv`
- `resubmission_cause_specific_cox_results_arf_subtype_sensitivity.csv`
- `resubmission_cause_specific_cox_results_sofa_sensitivity.csv`
- `resubmission_cause_specific_cox_results_no_icu_los_restriction.csv`
- `resubmission_aalen_johansen_cif_plot_data.csv`
- `resubmission_aalen_johansen_cif_quartiles.png`
- `resubmission_cohort_summary.csv`
