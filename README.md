# *CLIF Project Title*

## CLIF VERSION 

[major].[minor]

## Objective

*Describe the project objective*

## Required CLIF tables and fields

Please refer to the online [CLIF data dictionary](https://clif-consortium.github.io/website/data-dictionary.html), [ETL tools](https://github.com/clif-consortium/CLIF/tree/main/etl-to-clif-resources), and [specific table contacts](https://github.com/clif-consortium/CLIF?tab=readme-ov-file#relational-clif) for more information on constructing the required tables and fields. 

*List all required tables for the project here, and provide a brief rationale for why they are required.*

Example:
The following tables are required:
1. **patient**: `patient_id`, `race_category`, `ethnicity_category`, `sex_category`
2. **hospitalization**: `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`
3. **vitals**: `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value`
   - `vital_category` = 'heart_rate', 'resp_rate', 'sbp', 'dbp', 'map', 'resp_rate', 'spo2'
4. **labs**: `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value`
   - `lab_category` = 'lactate'
5. **medication_admin_continuous**: `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit`
   - `med_category` = "norepinephrine", "epinephrine", "phenylephrine", "vasopressin", "dopamine", "angiotensin", "nicardipine", "nitroprusside", "clevidipine", "cisatracurium"
6. **respiratory_support**: `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `lpm_set`, `resp_rate_set`, `peep_set`, `resp_rate_obs`

## Cohort identification
*Describe study cohort inclusion and exclusion criteria here*

## Expected Results

*Describe the output of the analysis. The final project results should be saved in the [`output/final`](output/README.md) directory.*

## Detailed Instructions for running the project

## 1. Update `config/config.json`
Follow instructions in the [config/README.md](config/README.md) file for detailed configuration steps.

**Note: if using the `01_run_cohort_id_app.R` file, this step is not necessary as the app will create the config file for the user**

## 2. Set up the project environment

*Describe the steps to setup the project environment.*

Example for R:
Run `00_renv_restore.R` in the [code](code/templates/R) to set up the project environment

Example for Python:
Open your terminal and run the following commands:
```
python3 -m venv .mobilization
source .mobilization/bin/activate
pip install -r requirements.txt 
```

## 3. Run code

Detailed instructions on the code workflow are provided in the [code directory](code/README.md)

## Example Repositories
* [CLIF Adult Sepsis Events](https://github.com/08wparker/CLIF_sepsis) for R
* [CLIF Eligibility for mobilization](https://github.com/kaveriC/CLIF-eligibility-for-mobilization) for Python
* [CLIF Variation in Ventilation](https://github.com/ingra107/clif_vent_variation)
---


