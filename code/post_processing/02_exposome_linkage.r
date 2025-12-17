

# packages
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("janitor", quietly = TRUE)) install.packages("janitor")
library(tidyverse); library(janitor)

# ---- helper to standardize column names for both schemas ----
standardize_pollutant <- function(df, pollutant = c("no2","pm25")) {
  pollutant <- match.arg(pollutant)
  df <- janitor::clean_names(df)
  # find the value column based on common header variants
  val_col <- switch(
    pollutant,
    "no2"  = names(df)[names(df) %in% c("no2_mean","mean_no2","no2","mean_no_2","avg_no2")][1],
    "pm25" = names(df)[names(df) %in% c("pm25_mean","mean_pm25","pm2_5_mean","pm25","mean_pm2_5")][1]
  )
  if (is.na(val_col)) stop("Couldn't find a value column for ", pollutant, " in: ", paste(names(df), collapse = ", "))
  
  df %>%
    transmute(
      GEOID = as.character(.data[[names(df)[names(df) %in% c("geoid","county_fips","fips")][1]]]),
      year  = as.integer(.data[[names(df)[names(df) %in% c("year","yr")][1]]]),
      value = as.numeric(.data[[val_col]])
    )
}

# ---- file paths (edit if your filenames differ) ----
no2_files  <- c("conus_county_no2_2005_2020.csv", "no2_county_year.csv")
pm25_files <- c("conus_county_pm25_1998_2023.csv", "pm25_county_year.csv")

# optional: a lookup for STATEFP/NAME using any conus_county_* file you have
county_lookup <- tryCatch(
  readr::read_csv("conus_county_pm25_1998_2023.csv", show_col_types = FALSE) |>
    clean_names() |>
    transmute(GEOID = as.character(geoid),
              STATEFP = as.character(statefp),
              NAME = name) |>
    distinct(),
  error = function(e) NULL
)

# ---- NO2 ----
no2_list <- lapply(no2_files[file.exists(no2_files)], function(f) {
  readr::read_csv(f, show_col_types = FALSE) |>
    standardize_pollutant("no2") |>
    mutate(source = basename(f),
           priority = ifelse(grepl("county_year", f), 1L, 2L))  # prefer *_county_year
})

no2_all <- bind_rows(no2_list) |>
  filter(dplyr::between(year, 2005, 2024)) |>
  arrange(GEOID, year, priority) |>
  distinct(GEOID, year, .keep_all = TRUE) |>
  select(GEOID, year, no2_mean = value)

if (!is.null(county_lookup)) {
  no2_all <- left_join(county_lookup, no2_all, by = "GEOID") |>
    select(GEOID, STATEFP, NAME, year, no2_mean)
}

no2_all <- arrange(no2_all, GEOID, year)

readr::write_csv(no2_all, "conus_county_no2_2005_2024.csv")

# ---- PM2.5 ----
pm25_list <- lapply(pm25_files[file.exists(pm25_files)], function(f) {
  readr::read_csv(f, show_col_types = FALSE) |>
    standardize_pollutant("pm25") |>
    mutate(source = basename(f),
           priority = ifelse(grepl("county_year", f), 1L, 2L))
})

pm25_all <- bind_rows(pm25_list) |>
  filter(dplyr::between(year, 2005, 2024)) |>
  arrange(GEOID, year, priority) |>
  distinct(GEOID, year, .keep_all = TRUE) |>
  select(GEOID, year, pm25_mean = value) 

if (!is.null(county_lookup)) {
  pm25_all <- left_join(county_lookup, pm25_all, by = "GEOID") |>
    select(GEOID, STATEFP, NAME, year, pm25_mean)
}

pm25_all <- arrange(pm25_all, GEOID, year)

readr::write_csv(pm25_all, "conus_county_pm25_2005_2024.csv")

# quick checks
no2_all %>% count(year) %>% print(n = 30)
pm25_all %>% count(year) %>% print(n = 30)


# needs zoo for interpolation
if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
library(zoo)

# ensure we have one row per county-year across 2005–2024, then linearly interpolate
pm25_all <- pm25_all |>
  # keep only needed columns
  select(GEOID, STATEFP, NAME, year, pm25_mean) |>
  # make a complete year grid per county
  group_by(GEOID, STATEFP, NAME) |>
  tidyr::complete(year = 2005:2024) |>
  arrange(GEOID, year) |>
  # interpolate within each county: rule=2 lets us extrapolate ends (e.g., 2024)
  mutate(
    pm25_mean_interp = zoo::na.approx(pm25_mean, x = year, na.rm = FALSE, rule = 2),
    pm25_filled_2024 = if_else(year == 2024 & is.na(pm25_mean) & !is.na(pm25_mean_interp), TRUE, FALSE),
    pm25_mean = coalesce(pm25_mean, pm25_mean_interp)
  ) |>
  ungroup() |>
  select(GEOID, STATEFP, NAME, year, pm25_mean, pm25_filled_2024) |>
  arrange(GEOID, year)

# write out
readr::write_csv(pm25_all, "conus_county_pm25_2005_2024.csv")








