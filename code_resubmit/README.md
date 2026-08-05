# REFER Resubmission Primary Analysis

This folder is a cleaner, reviewer-focused replacement for the original first-run model block in `code/02_REFER_linkage_analysis.R`.

The goal is not to expand the analysis in every direction. It implements only the main changes implied by the reviewer feedback:

1. Use mortality by day 28 after ARF onset as the primary binary mortality estimand.
2. Keep Cox proportional hazards models for mortality after ARF onset as a sensitivity analysis.
3. Replace raw invasive ventilation duration with ventilator-free days through day 28.
4. Use the stronger covariate set now available: age, sex, race/ethnicity, calendar year, Charlson score, and available ACS/social vulnerability covariates.

## Main Script

`01_primary_reviewer_optimized_models.R`

This script assumes the original cohort/linkage workflow has already constructed `arf_exp`. If available, it also uses:

- `support_class` to compute ventilator-free days from invasive mechanical ventilation records
- `diagnosis` to calculate Charlson score when `charlson_score` is not already present
- `patient` to recover death time when `arf_exp` does not already include one

The script writes outputs to:

`output/code_resubmit/<timestamp>/`

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

## Output Files

- `analysis_dataset_reviewer_optimized.csv`
- `cohort_summary_reviewer_optimized.csv`
- `primary_mortality_day28_logistic_results.csv`
- `primary_mortality_cox_results.csv`
- `primary_mortality_cox_ph_diagnostics.csv`
- `primary_vfd_quasipoisson_results.csv`
