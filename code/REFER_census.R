
# install.packages(c("tidycensus", "tidyverse"))
library(tidycensus)
library(tidyverse)

# 0) Census API key (do once, then restart R or readRenviron)
# census_api_key("YOUR_KEY_HERE", install = TRUE)
# readRenviron("~/.Renviron")

# 1) Years (ACS 5-year "end" years)
years <- 2018:2023

# 2) CONUS states (48 + DC; exclude AK, HI, PR & territories)
conus_states <- setdiff(c(state.abb, "DC"), c("AK", "HI"))

# 3) Variables (estimates). We'll compute percents afterward.
# Notes:
# - Race alone counts (B02001); Hispanic (B03003)
# - Poverty (B17001)
# - Education 25+ (B15003)
# - Labor force & unemployment (B23025)
# - Insurance % from S2701 (may be unavailable in some years → we stub NA)
acs_vars <- c(
  total_pop      = "B01003_001",
  median_age     = "B01002_001",
  median_income  = "B19013_001",
  white_alone    = "B02001_002",
  black_alone    = "B02001_003",
  asian_alone    = "B02001_005",
  hispanic_any   = "B03003_003",
  pov_below      = "B17001_002",
  pov_universe   = "B17001_001",
  edu_total      = "B15003_001",
  lf_total       = "B23025_001",
  lf_in_lf       = "B23025_002",
  lf_unemployed  = "B23025_005",
  # Insurance subject table (percent insured, universe)
  insured_total  = "S2701_C01_001",
  insured_pct    = "S2701_C05_001"
)

# Helpers to build education group sums
expand_seq <- function(prefix, a, b) sprintf("%s_%03d", prefix, a:b)
edu_lt_hs_codes <- expand_seq("B15003", 2, 16)   # less than HS
edu_sc_aa_codes <- expand_seq("B15003", 18, 20)  # some college/AA
edu_bachp_codes <- expand_seq("B15003", 21, 25)  # bachelor's+

safe_pct <- function(num, den) ifelse(is.finite(num/den) & den > 0, 100 * num/den, NA_real_)

# 4) Fetch one (year, state) chunk
get_tract_chunk <- function(yr, st) {
  message("Fetching year=", yr, " state=", st, " ...")
  vars_needed <- unique(c(
    unname(acs_vars),
    edu_lt_hs_codes, edu_sc_aa_codes, edu_bachp_codes,
    "B15003_017" # HS grad
  ))
  
  acs <- get_acs(
    geography = "tract",
    state = st,
    variables = vars_needed,
    year = yr, survey = "acs5",
    output = "wide",
    cache_table = TRUE
  )
  
  # Helper to sum *_E columns
  sum_colsE <- function(df, codes) {
    cols <- paste0(codes, "E")
    rowSums(dplyr::select(df, dplyr::any_of(cols)), na.rm = TRUE)
  }
  
  # Check if S2701 is present for this (year,state)
  has_s2701 <- "S2701_C05_001E" %in% names(acs)
  
  out <- acs %>%
    dplyr::transmute(
      year = yr,
      geoid = GEOID,
      name  = NAME,
      
      total_pop     = .data[["B01003_001E"]],
      median_age    = .data[["B01002_001E"]],
      median_income = .data[["B19013_001E"]],
      
      pct_white = safe_pct(.data[["B02001_002E"]], .data[["B01003_001E"]]),
      pct_black = safe_pct(.data[["B02001_003E"]], .data[["B01003_001E"]]),
      pct_asian = safe_pct(.data[["B02001_005E"]], .data[["B01003_001E"]]),
      pct_hisp  = safe_pct(.data[["B03003_003E"]], .data[["B01003_001E"]]),
      
      pov_rate_pct = safe_pct(.data[["B17001_002E"]], .data[["B17001_001E"]]),
      
      # Education (25+)
      edu_totalE  = .data[["B15003_001E"]],
      edu_lt_hsE  = sum_colsE(cur_data_all(), edu_lt_hs_codes),
      edu_sc_aaE  = sum_colsE(cur_data_all(), edu_sc_aa_codes),
      edu_bachpE  = sum_colsE(cur_data_all(), edu_bachp_codes),
      pct_lt_hs   = safe_pct(edu_lt_hsE, edu_totalE),
      pct_hs      = safe_pct(.data[["B15003_017E"]], edu_totalE),
      pct_sc_aa   = safe_pct(edu_sc_aaE, edu_totalE),
      pct_bachplus= safe_pct(edu_bachpE, edu_totalE),
      
      labor_force    = .data[["B23025_002E"]],
      unemp_rate_pct = safe_pct(.data[["B23025_005E"]], .data[["B23025_002E"]]),
      
      # Insurance % (only if S2701 exists)
      pct_insured = if (has_s2701) .data[["S2701_C05_001E"]] else NA_real_
    ) %>%
    dplyr::select(-edu_totalE, -edu_lt_hsE, -edu_sc_aaE, -edu_bachpE)
  
  out
}


# 5) Run: all years x all CONUS states
all_tracts <- map_dfr(
  years,
  ~{
    yr <- .x
    map_dfr(conus_states, ~get_tract_chunk(yr, .x))
  }
)

# 6) Arrange & write CSV
all_tracts <- all_tracts %>%
  arrange(year, geoid)

out_path <- "conus_tract_sociodemographics_acs5_2018_2023.csv"
readr::write_csv(all_tracts, out_path)
message("Done! Wrote: ", out_path)

# Optional: sanity check counts (tracts per year)
all_tracts %>% count(year) %>% print(n = Inf)

# Filter out the 2023 data
tracts_2023 <- all_tracts %>% 
  dplyr::filter(year == 2023)

# Duplicate and change year to 2024
tracts_2024 <- tracts_2023 %>% 
  dplyr::mutate(year = 2024)

# Bind back together
all_tracts <- dplyr::bind_rows(all_tracts, tracts_2024) %>% 
  dplyr::arrange(year, geoid)

# Write out the updated file
out_path <- "conus_tract_sociodemographics_acs5_2018_2024.csv"
readr::write_csv(all_tracts, out_path)

message("Done! Wrote: ", out_path)



