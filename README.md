# Environmental and Clinical Determinants of Acute Respiratory Failure Trajectories in ICU Patients

## CLIF VERSION

2.1.0

## Objective

To determine clinical and non-clinical factors, including environmental
exposures such as chronic air pollution and heat vulnerability, that are
associated with the onset and outcomes of acute respiratory failure in
ICU patients in Chicago.

## Required CLIF tables and fields

**Demographics**

- **patient**: `patient_id`, `birth_date`, `race_category`, `ethnicity_category`, `sex_category`, `preferred_language`
  - *Note:* `zip_code` and all other geocoding fields are captured on the **hospitalization** record (see below).

**Hospitalization & ICU stay**
- **hospitalization**: `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `admitting_service`, `discharge_service`, `zip_code`, `death_dttm`
  - (Depending on availability) geocoded linkage such as `latitude`, `longitude`, `census_tract`, `census_block_group`, `county`, etc.

**Clinical trajectories (type-specific)**
- **vitals**: `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value`  
  - Include **all** available `vital_category` values. At minimum: `'heart_rate'`, `'resp_rate'`, `'sbp'`, `'dbp'`, `'map'`, `'spo2'`, `'temperature'`, `'height'`, `'weight'`.

- **labs**: `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value`
  - For hypoxemic ARF: `lab_category` ∈ `'pao2'`, `'fio2'`, `'pao2_fio2_ratio'`, `'spo2'`
  - For hypercapnic ARF: `lab_category` ∈ `'paco2'`, `'ph'`, `'bicarbonate'`

**Therapeutics**
- **medication_admin_continuous**: `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit`  
  - Vasopressors/vasoactives: `"norepinephrine"`, `"epinephrine"`, `"phenylephrine"`, `"vasopressin"`, `"dopamine"`, `"angiotensin"`  
  - Antihypertensives (continuous): `"nicardipine"`, `"nitroprusside"`, `"nitroglycerin"` 
  - Neuromuscular blockade: `"cisatracurium"`, `"vecuronium"`, `"rocuronium"`  
  - Respiratory/airway: `"naloxone"` (narcan), `"albuterol_continuous"` (inhaled continuous albuterol)  
  - Sedation/analgesia (continuous): `"propofol"`, `"midazolam"`, `"dexmedetomidine"`, `"fentanyl"`

**Respiratory support**
- **respiratory_support**: `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `fio2_set`, `peep_set`, `resp_rate_set`, `tidal_volume_set`, `plateau_pressure`, `pao2_fio2_ratio`

**Diagnosis & outcomes**
- **hospital_diagnosis**: `hospitalization_id`, `diagnosis_code`, `diagnosis_category`, `diagnosis_type`  
  - Hypoxemic ARF: ICD-10 `J96.0x`  
  - Hypercapnic ARF: ICD-10 `J96.1x`  
  - Acute on chronic respiratory failure: `J96.2x` (specify hypoxemic vs hypercapnic if coded)

- **icu_outcomes**: `icu_mortality`, `hospital_mortality`, `icu_length_of_stay`, `hospital_length_of_stay`

**Control cohort (perioperative respiratory failure)**
- ICD-10 `J95.82`: Acute pulmonary insufficiency following thoracic surgery  
- ICD-10 `J95.83`: Acute pulmonary insufficiency following nonthoracic surgery  
- ICD-10 `J95.84`: Acute and chronic respiratory failure following surgery


------------------------------------------------------------------------

## Cohort identification

**Inclusion criteria**:

1.  Adult patients (≥18 years) admitted to ICU between 2018–2024.

2.  At least one of the following criteria of acute respiratory failure
    is met:

    -   Acute hypoxemic respiratory failure (any one of the following)

        -   SpO2 less than 90% on room air

        -   PaO2 of 60 mm Hg or less on room air

        -   PaO2–FiO2 ratio of 300 or less (on any amount of FiO2)

    -   Acute hypercapnic respiratory failure (both of the following)

        -   PaCO2 of 45 mm Hg or more AND

        -   Arterial pH \< 7.35

3.  Available ABG and/or continuous pulse oximetry data within ±24h of
    ICU admission.

4.  Residential ZIP code for environmental data linkage.

***-\> Note that mixed hypoxic and hypercapnic respiratory failure is
common and should be accounted for.***

***-\> Also note that for SpO2 and PaO2, these numbers are directly
affected by supplemental oxygen (FiO2) via whatever delivery mechanism.
The definitions for ARF we choose for these values will be on room air
(21% FiO2). P/F ratio can define ARF even on supplemental oxygen.***

**Exclusion criteria** 
- Missing key demographic data (age, sex, race). 
- Hospitalizations \<24 hours in ICU. 
- Repeat ICU stays within same hospitalization (only first considered for primary analysis).

**Bibliography for definitions of ARF**

1\. Lagina, M. & Valley, T. S. Diagnosis and Management of Acute
Respiratory Failure. Critical Care Clinics 40, 235–253 (2024).

2\. Baldomero, A. K. et al. Effectiveness and Harms of High-Flow Nasal
Oxygen for Acute Respiratory Failure: An Evidence Report for a Clinical
Guideline From the American College of Physicians. Ann Intern Med 174,
952–966 (2021).

3\. RENOVATE Investigators and the BRICNet Authors et al. High-Flow
Nasal Oxygen vs Noninvasive Ventilation in Patients With Acute
Respiratory Failure: The RENOVATE Randomized Clinical Trial. JAMA 333,
875 (2025).

4\. Mirabile, V. S., Shebl, E., Sankari, A. & Burns, B. Respiratory
Failure in Adults. in StatPearls (StatPearls Publishing, Treasure Island
(FL), 2025).

## Detailed Instructions for running the project

## 1. Update `config/config.json`

Follow instructions in the [config/README.md](config/README.md) file for detailed configuration steps.

**Note: if using the `01_run_cohort_id_app.R` file, this step is not necessary as the app will create the config file for the user**

## 2. Set up the project environment

*Describe the steps to setup the project environment.*

Example for R: Run `00_renv_restore.R` in the [code](code/templates/R) to set up the project environment

Example for Python: Open your terminal and run the following commands:

```         
python3 -m venv .mobilization
source .mobilization/bin/activate
pip install -r requirements.txt 
```

## 3. Run code

Detailed instructions on the code workflow are provided in the [code directory](code/README.md)

## Example Repositories

-   [CLIF Adult Sepsis Events](https://github.com/08wparker/CLIF_sepsis) for R

-   [CLIF Eligibility for mobilization](https://github.com/kaveriC/CLIF-eligibility-for-mobilization) for Python

-   

    ## [CLIF Variation in Ventilation](https://github.com/ingra107/clif_vent_variation)