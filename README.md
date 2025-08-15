# Environmental and Clinical Determinants of Acute Respiratory Failure Trajectories in ICU Patients

## CLIF VERSION
2.0.0

## Objective
To determine clinical and non-clinical factors, including environmental exposures such as chronic air pollution and heat vulnerability, that are associated with the onset and outcomes of acute respiratory failure in ICU patients in Chicago.

## Required CLIF tables and fields

**Demographics**
- **patient**: `patient_id`, `birth_date`, `race_category`, `ethnicity_category`, `sex_category`, `zip_code`, `preferred_language`

**Hospitalization & ICU stay**
- **hospitalization**: `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `admitting_service`, `discharge_service`

**Clinical trajectories (type-specific)**
- **vitals**: `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value`  
  - `vital_category` = 'heart_rate', 'resp_rate', 'sbp', 'dbp', 'map', 'spo2', 'temperature'
  - Rationale: Key for trending oxygen saturation and respiratory rate.

- **labs**: `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value`  
  - For hypoxemic ARF: `lab_category` = 'pao2', 'fio2', 'pao2_fio2_ratio', 'spo2'  
  - For hypercapnic ARF: `lab_category` = 'paco2', 'ph', 'bicarbonate'  

**Therapeutics**
- **medication_admin_continuous**: `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit`
  - `med_category` = "norepinephrine", "epinephrine", "phenylephrine", "vasopressin", "dopamine", "angiotensin", "nicardipine", "nitroprusside", "clevidipine", "cisatracurium"

**Respiratory support**
- **respiratory_support**: `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `fio2_set`, `peep_set`, `resp_rate_set`, `tidal_volume_set`, `plateau_pressure`, `pao2_fio2_ratio`  
  - Rationale: Tracks escalation or titration of ventilatory support for each ARF subtype.

**Diagnosis & outcomes**
- **diagnosis**: `hospitalization_id`, `diagnosis_code`, `diagnosis_category`, `diagnosis_type`  
  - Hypoxemic ARF: ICD-10 `J96.0x`  
  - Hypercapnic ARF: ICD-10 `J96.1x`  
  - Acute on chronic respiratory failure: `J96.2x` (specify hypoxemic vs hypercapnic if coded)
- **icu_outcomes**: `icu_mortality`, `hospital_mortality`, `icu_length_of_stay`, `hospital_length_of_stay`

---

## Cohort identification

**Inclusion criteria**
1. Adult patients (≥18 years) admitted to ICU between 2015–2025.
2. At least one documented ICD-10 diagnosis of ARF:
   - Hypoxemic: J96.0x
   - Hypercapnic: J96.1x
   - Acute on chronic (specify type if possible): J96.2x
3. Available ABG and/or continuous pulse oximetry data within ±24h of ICU admission.
4. Residential ZIP code for environmental data linkage.

**Exclusion criteria**
- Missing key demographic data (age, sex, race).
- Hospitalizations <24 hours in ICU.
- Repeat ICU stays within same hospitalization (only first considered for primary analysis).

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


