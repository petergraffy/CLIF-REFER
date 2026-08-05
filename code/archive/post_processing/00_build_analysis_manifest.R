# ================================================================================================
# Build curated analysis manifest
#
# Non-destructively copies canonical pooled-analysis inputs into analysis_ready/
# and writes analysis_ready/manifest.csv. Downstream pooled code should only use
# active rows from this manifest.
# ================================================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(fs)
  library(tools)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(cmd_file)) normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = TRUE) else NA_character_
}
repo <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE) else normalizePath(getwd(), mustWork = TRUE)
if (!dir_exists(file.path(repo, "sites"))) repo <- normalizePath(getwd(), mustWork = TRUE)

analysis_dir <- file.path(repo, "sites", "analysis")
ready_dir <- file.path(repo, "analysis_ready")
input_dir <- file.path(ready_dir, "site_inputs")

dir_create(ready_dir)
dir_create(input_dir)
dir_create(file.path(ready_dir, "pooled", "tables"))
dir_create(file.path(ready_dir, "pooled", "figures"))

norm_site <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "^JHU$", "Hopkins")
  x <- str_replace_all(x, "^hopkins$", "Hopkins")
  x <- str_replace_all(x, "^UPenn$", "UPenn")
  x
}

site_from_path <- function(path) {
  bn <- basename(path)
  m <- str_match(bn, "refer_([^_]+)_")
  if (!is.na(m[, 2])) return(norm_site(m[, 2]))
  m <- str_match(bn, "^([^_]+)_\\d{8}_model")
  if (!is.na(m[, 2])) return(norm_site(m[, 2]))
  parts <- str_split(path, .Platform$file.sep, simplify = TRUE)
  candidates <- c("UCMC", "Emory", "EU", "UMN", "umn", "UCSF", "NU", "OHSU",
                  "RUSH", "Penn", "UPenn", "Michigan", "Hopkins", "JHU")
  hit <- intersect(candidates, parts)
  if (length(hit)) return(norm_site(hit[[length(hit)]]))
  NA_character_
}

read_first_row <- function(path) {
  tryCatch(readr::read_csv(path, n_max = 1, show_col_types = FALSE, progress = FALSE),
           error = function(e) tibble())
}

first_value <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% names(df) && nrow(df) > 0) return(as.character(df[[nm]][[1]]))
  }
  NA_character_
}

classify_model <- function(path) {
  bn <- basename(path)
  tibble(
    outcome = case_when(
      str_detect(bn, "death30d") ~ "death30d",
      str_detect(bn, "inhosp") ~ "inhosp_death",
      str_detect(bn, "vent_hours|venthrs") ~ "vent_hours",
      str_detect(bn, "iculos|icu_los") ~ "icu_los",
      TRUE ~ NA_character_
    ),
    pollutant = case_when(
      str_detect(bn, "no2") ~ "NO2",
      str_detect(bn, "pm25") ~ "PM2.5",
      TRUE ~ NA_character_
    ),
    model_type = case_when(
      str_detect(bn, "x_subtype") ~ "x_subtype",
      str_detect(bn, "cumulative") ~ "cumulative",
      str_detect(bn, "adj_plus_covid") ~ "adj_plus_covid",
      TRUE ~ "other"
    )
  )
}

copy_curated <- function(df) {
  if (!nrow(df)) return(df)
  df %>%
    mutate(
      curated_dir = file.path(input_dir, family, if_else(is.na(site) | site == "", "unknown_site", site)),
      curated_path = file.path(curated_dir, basename(source_path))
    ) %>%
    rowwise() %>%
    mutate(
      copied = {
        dir_create(curated_dir)
        file.copy(source_path, curated_path, overwrite = TRUE)
      },
      md5 = unname(tools::md5sum(curated_path))
    ) %>%
    ungroup() %>%
    select(-curated_dir)
}

latest_by <- function(df, keys) {
  if (!nrow(df)) return(df)
  df %>%
    mutate(
      run_key = str_extract(basename(source_path), "\\d{8}_\\d{6}|\\d{8}"),
      run_key = if_else(is.na(run_key), "00000000", run_key)
    ) %>%
    arrange(across(all_of(keys)), desc(run_key), desc(file_info(source_path)$modification_time), source_path) %>%
    group_by(across(all_of(keys))) %>%
    slice(1) %>%
    ungroup()
}

make_record <- function(paths, family) {
  if (!length(paths)) return(tibble(source_path = character(), family = character()))
  tibble(source_path = normalizePath(paths, mustWork = TRUE), family = family)
}

records <- list()

# Site summaries and QA
records$site_summary <- dir_ls(file.path(repo, c("sites", "output_10082025")), recurse = TRUE, type = "file",
                               regexp = "arf_summary_.*\\.csv$") %>%
  make_record("site_summary") %>%
  mutate(site = map_chr(source_path, site_from_path)) %>%
  latest_by(c("site", "family"))

records$missingness <- dir_ls(file.path(repo, c("sites", "output_10082025")), recurse = TRUE, type = "file",
                              regexp = "missing_pct_.*\\.csv$") %>%
  make_record("missingness") %>%
  mutate(site = map_chr(source_path, site_from_path)) %>%
  latest_by(c("site", "family"))

# County counts are pooled spatial inputs. Keep one canonical aggregate if present.
county_paths <- dir_ls(analysis_dir, type = "file", regexp = "^combined_arf_counts_by_county_year\\.csv$")
if (!length(county_paths)) county_paths <- dir_ls(analysis_dir, type = "file", regexp = "^arf_counts_by_county_year\\.csv$")
records$county_counts <- make_record(county_paths, "county_counts") %>%
  mutate(site = "pooled")

# CIF files in sites/analysis/cif are the complete current collection. Read one row to identify site.
cif_dir <- file.path(analysis_dir, "cif")
cif_paths <- dir_ls(cif_dir, type = "file", regexp = "site_cif_plotdf__.*_cumulative(?: copy \\d+)?\\.csv$")
records$cif_plotdf <- tibble(source_path = normalizePath(cif_paths, mustWork = TRUE)) %>%
  mutate(
    first = map(source_path, read_first_row),
    site = map_chr(first, ~ norm_site(first_value(.x, c("site_name", "site", "site_id")))),
    family = "cif_plotdf",
    pollutant = case_when(str_detect(basename(source_path), "NO2") ~ "NO2",
                          str_detect(basename(source_path), "PM25") ~ "PM2.5",
                          TRUE ~ NA_character_),
    outcome = NA_character_
  ) %>%
  select(-first) %>%
  latest_by(c("site", "family", "pollutant"))

bin_paths <- dir_ls(cif_dir, type = "file", regexp = "site_exposure_bins__.*_cumulative(?: copy \\d+)?\\.csv$")
records$cif_bins <- tibble(source_path = normalizePath(bin_paths, mustWork = TRUE)) %>%
  mutate(
    first = map(source_path, read_first_row),
    site = map_chr(first, ~ norm_site(first_value(.x, c("site_name", "site", "site_id")))),
    family = "cif_bins",
    pollutant = case_when(str_detect(basename(source_path), "NO2") ~ "NO2",
                          str_detect(basename(source_path), "PM25") ~ "PM2.5",
                          TRUE ~ NA_character_),
    outcome = NA_character_
  ) %>%
  select(-first) %>%
  latest_by(c("site", "family", "pollutant"))

# Model inputs for pooled analyses
model_paths <- dir_ls(file.path(analysis_dir, "models"), recurse = FALSE, type = "file", regexp = "\\.csv$")
records$model_cumulative <- tibble(source_path = normalizePath(model_paths, mustWork = TRUE)) %>%
  bind_cols(map_dfr(.$source_path, classify_model)) %>%
  mutate(site = map_chr(source_path, site_from_path), family = "model_cumulative") %>%
  filter(model_type == "cumulative", !is.na(site), !is.na(outcome), !is.na(pollutant)) %>%
  latest_by(c("site", "family", "outcome", "pollutant"))

subtype_paths <- dir_ls(file.path(analysis_dir, "subtype"), recurse = FALSE, type = "file",
                        regexp = "model_tidy_.*x_subtype.*\\.csv$")
records$model_subtype <- tibble(source_path = normalizePath(subtype_paths, mustWork = TRUE)) %>%
  bind_cols(map_dfr(.$source_path, classify_model)) %>%
  mutate(site = map_chr(source_path, site_from_path), family = "model_subtype", pollutant = "NO2") %>%
  filter(!is.na(site), outcome %in% c("death30d", "inhosp_death", "vent_hours")) %>%
  latest_by(c("site", "family", "outcome"))

# Fine-Gray SHR site exports
shr_paths <- dir_ls(analysis_dir, type = "file", regexp = "shr_all_outcomes_cumulative_only_NOLOOP.*\\.csv$")
records$finegray_shr <- make_record(shr_paths, "finegray_shr") %>%
  mutate(site = map_chr(source_path, site_from_path)) %>%
  latest_by(c("site", "family"))

# Table 1 artifacts
table_paths <- dir_ls(file.path(analysis_dir, "tab1"), type = "file", regexp = "\\.(csv|rtf)$")
records$table1 <- make_record(table_paths, "table1") %>%
  mutate(site = map_chr(source_path, site_from_path),
         site = if_else(is.na(site), "pooled", site)) %>%
  latest_by(c("site", "family", "source_path"))

manifest <- bind_rows(records) %>%
  mutate(
    site = norm_site(site),
    outcome = outcome %||% NA_character_,
    pollutant = pollutant %||% NA_character_,
    source_relpath = path_rel(source_path, start = repo),
    active = TRUE,
    exists = file_exists(source_path),
    notes = ""
  ) %>%
  arrange(family, site, pollutant, outcome, source_relpath)

manifest <- copy_curated(manifest) %>%
  mutate(curated_relpath = path_rel(curated_path, start = repo)) %>%
  select(site, family, pollutant, outcome, active, exists, source_relpath,
         curated_relpath, md5, notes)

manifest_path <- file.path(ready_dir, "manifest.csv")
readr::write_csv(manifest, manifest_path)

summary_path <- file.path(ready_dir, "manifest_summary.csv")
manifest %>%
  count(family, site, name = "n_files") %>%
  arrange(family, site) %>%
  readr::write_csv(summary_path)

message("Wrote manifest: ", manifest_path)
message("Wrote summary:  ", summary_path)
message("Curated files:  ", nrow(manifest))
