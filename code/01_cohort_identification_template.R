# Cohort identification script for inpatient admissions - ARF Study

# Load required libraries
library(knitr)
library(here)
library(tidyverse)
library(arrow)
library(gtsummary)

# Objective: identify a cohort of hospitalizations from CLIF tables
# for acute respiratory failure analysis

# Specify inpatient cohort parameters
start_date <- "2018-01-01"
end_date <- "2024-12-31"
include_pediatric <- FALSE
include_er_deaths <- TRUE

# Specify required CLIF tables (updated for CLIF v2.1 and this project)
tables <- c("patient", "hospitalization", "vitals", "labs", 
            "medication_admin_continuous", "adt", 
            "respiratory_support", "hospital_diagnosis", 
            "microbiology_culture")

true_tables <- c("patient", "hospitalization", "adt", "vitals", "labs", 
                 "medication_admin_continuous", "respiratory_support", 
                 "hospital_diagnosis")

table_flags <- setNames(tables %in% true_tables, tables)

# Load configuration utility
source("utils/config.R")
site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type

print(paste("Site Name:", site_name))
print(paste("Tables Path:", tables_path))
print(paste("File Type:", file_type))

# Load required CLIF files
clif_table_filenames <- list.files(path = tables_path, 
                                   pattern = paste0("^clif_.*\\.", file_type, "$"), 
                                   full.names = TRUE)
clif_table_basenames <- basename(clif_table_filenames) %>%
  str_remove(paste0("\\.", file_type, "$"))

required_files <- paste0("clif_", names(table_flags)[table_flags])
missing_tables <- setdiff(required_files, clif_table_basenames)
if (length(missing_tables) > 0) {
  stop(paste("Missing required tables:", 
             paste(missing_tables, collapse = ", ")))
}
required_filenames <- clif_table_filenames[clif_table_basenames %in% required_files] 

# Read CLIF tables
if (file_type == "parquet") {
  data_list <- lapply(required_filenames, function(file) {
    open_dataset(file) %>% collect()
  })
} else if (file_type == "csv") {
  data_list <- lapply(required_filenames, read_csv)
} else if (file_type == "fst") {
  data_list <- lapply(required_filenames, read.fst)
} else {
  stop("Unsupported file format")
}

# Assign tables
for (i in seq_along(required_filenames)) {
  object_name <- str_remove(basename(required_filenames[i]),
                            paste0("\\.", file_type, "$"))
  object_name <- make.names(object_name)
  assign(object_name, data_list[[i]])
}
remove(data_list)

# Filter by admission date and age
clif_hospitalization_filtered <- clif_hospitalization %>%
  filter(admission_dttm >= start_date & admission_dttm <= end_date)

if (!include_pediatric) {
  clif_hospitalization_filtered <- clif_hospitalization_filtered %>%
    filter(age_at_admission >= 18)
}

# Identify hospitalizations with ICU or ward stays
inpatient_hospitalization_ids <- clif_adt %>%
  filter(location_category %in% c("ward", "icu")) %>%
  select(hospitalization_id) %>%
  distinct() %>%
  pull(hospitalization_id)

cohort_hospitalization_ids <- clif_hospitalization_filtered %>%
  filter(hospitalization_id %in% inpatient_hospitalization_ids) %>%
  pull(hospitalization_id)

# Include ER-only deaths if specified
if (include_er_deaths) {
  ER_only_hospitalization_ids <- clif_adt %>%
    group_by(hospitalization_id) %>%
    filter(all(location_category == "er")) %>%
    pull(hospitalization_id)
  
  ER_death_ids <- clif_hospitalization %>%
    filter(hospitalization_id %in% ER_only_hospitalization_ids) %>%
    filter(discharge_category == "Expired") %>%
    pull(hospitalization_id)
  
  cohort_hospitalization_ids <- union(cohort_hospitalization_ids, ER_death_ids)
}

# # Optional: Filter to ARF-related diagnoses (J96.xx)
# resp_dx_ids <- clif_hospital_diagnosis %>%
#   filter(str_detect(as.character(diagnostic_code), "^J96")) %>%
#   pull(hospitalization_id) %>%
#   unique()
# 
# cohort_hospitalization_ids <- intersect(cohort_hospitalization_ids, resp_dx_ids)

# Export list of hospitalization IDs
save(cohort_hospitalization_ids, 
     file = here("output/intermediate/01_cohort_hospitalization_ids.RData"))

# Utility: filter tables by hospitalization_id
filter_clif_table <- function(table, filter_col, cohort_ids, select_cols = NULL) {
  filtered_table <- table %>% filter(!!sym(filter_col) %in% cohort_ids)
  if (!is.null(select_cols)) {
    filtered_table <- filtered_table %>% select(all_of(select_cols))
  }
  return(filtered_table)
}

# Remove patient from general loop
table_flags["patient"] <- FALSE

# Filter all CLIF tables by cohort
for (table_name in names(table_flags)[table_flags]) {
  full_table_name <- paste0("clif_", table_name)
  assign(
    paste0(full_table_name, "_cohort"), 
    filter_clif_table(get(full_table_name), "hospitalization_id", cohort_hospitalization_ids)
  )
}

# Filter patient table separately
cohort_patient_ids <- clif_hospitalization_filtered %>%
  filter(hospitalization_id %in% cohort_hospitalization_ids) %>%
  pull(patient_id) %>%
  unique()

clif_patient_cohort <- filter_clif_table(clif_patient, "patient_id", cohort_patient_ids)

# Save filtered tables
save(list = ls(pattern = "clif_.*_cohort"), 
     file = here("output/intermediate/01_clif_cohort_tables.RData"))

# Table 1: patient-level
admits_per_patient <- clif_hospitalization_filtered %>%
  group_by(patient_id) %>%
  summarise(n_hospitalizations = n())

table_one_patient <- clif_patient_cohort %>%
  left_join(admits_per_patient, by = "patient_id") %>%
  left_join(clif_hospitalization_filtered %>% select(patient_id, admission_dttm), by = "patient_id") %>%
  mutate(age = as.numeric(difftime(admission_dttm, birth_date, units = "days")) / 365.25) %>%
  select(age, sex_category, race_category, ethnicity_category, 
         n_hospitalizations, language_name) %>%
  tbl_summary()

print(table_one_patient)

# Table 1: hospitalization-level
ever_icu <- clif_adt %>%
  filter(hospitalization_id %in% cohort_hospitalization_ids) %>%
  filter(location_category == "icu") %>%
  select(hospitalization_id) %>%
  mutate(ever_icu = 1) %>%
  distinct()

table_one_hospitalization <- clif_hospitalization_filtered %>%
  mutate(length_of_stay = as.numeric(as.Date(discharge_dttm) - 
                                       as.Date(admission_dttm), units = "days")) %>%
  select(patient_id, hospitalization_id, age_at_admission, discharge_category, 
         admission_type_name, length_of_stay) %>%
  left_join(clif_patient_cohort %>% 
              select(patient_id, race_category, sex_category, 
                     ethnicity_category, language_name)) %>% 
  left_join(ever_icu, by = "hospitalization_id") %>%
  mutate(ever_icu = ifelse(is.na(ever_icu), 0, 1)) %>%
  select(-patient_id, -hospitalization_id) %>%
  tbl_summary(by = ever_icu)

# table_one_hospitalization %>% 
#   as_gt() %>% 
#   gt::gtsave(filename = here(paste0("output/intermediate/Table_One_", Sys.Date(), "_", 
#                                     config$site_name, ".pdf")))
