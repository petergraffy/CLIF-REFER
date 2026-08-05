

library(tidyverse)
library(janitor)
library(lubridate)
library(sf)
library(tigris)
library(tmap)
library(viridisLite)
library(ggplot2)
library(dplyr)
library(grid)
library(rvest)
library(xml2)
library(flextable)
library(officer)
library(janitor)
library(glue)
library(metafor)
library(rlang)
library(forcats)
library(patchwork)
library(fs)
library(cowplot)

# ---- Paths ----
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(cmd_file)) normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = TRUE) else NA_character_
}
repo <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = FALSE) else normalizePath(getwd(), mustWork = TRUE)
dir_in  <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis"
dir_out <- dir_in

# ---- Ingest all site CSVs ----
files <- list.files(dir_in, pattern = "^arf_counts_by_county_year.*\\.csv$", full.names = TRUE)

read_one <- function(f){
  df <- readr::read_csv(f, show_col_types = FALSE) |> clean_names()
  # try to find columns by common aliases
  fips_col  <- intersect(c("county_fips","fips","geoid","geoid10","county"), names(df))[1]
  year_col  <- intersect(c("year","yr","yyyy","date"), names(df))[1]
  count_col <- intersect(c("arf_count","count","n","cases","arf_cases","arf_counts"), names(df))[1]
  
  if (is.na(fips_col) || is.na(year_col) || is.na(count_col)) {
    stop(paste("Could not find expected columns in", basename(f)))
  }
  
  df |>
    transmute(
      source_file = basename(f),
      fips  = str_pad(as.character(.data[[fips_col]]), width = 5, side = "left", pad = "0"),
      year  = if (inherits(df[[year_col]], "Date") || inherits(df[[year_col]], "POSIXt")) year(.data[[year_col]]) else as.integer(.data[[year_col]]),
      arf_count = as.numeric(.data[[count_col]])
    ) |>
    filter(!is.na(fips), !is.na(year), !is.na(arf_count))
}

raw_all <- purrr::map_dfr(files, read_one)

# ---- De-duplicate within a file/year/county if necessary (sum) ----
by_year <- raw_all |>
  dplyr::group_by(fips, year) |>
  dplyr::summarize(arf_count = sum(arf_count, na.rm = TRUE), .groups = "drop")

# Write per-year combined table
readr::write_csv(by_year, file.path(dir_out, "combined_arf_counts_by_county_year.csv"))

# ---- Aggregate across the entire study period ----
by_total <- by_year |>
  dplyr::group_by(fips) |>
  dplyr::summarize(arf_total = sum(arf_count, na.rm = TRUE), .groups = "drop")

# Write total counts table
readr::write_csv(by_total, file.path(dir_out, "arf_counts_by_county_total.csv"))

# ---- Build county map (CONUS + DC) and join ----
options(tigris_use_cache = TRUE)
drop_states <- c("02","15","72","60","66","69","78") # AK, HI, PR, territories
us_counties <- counties(cb = TRUE, year = 2020, class = "sf") |>
  filter(!STATEFP %in% drop_states)  # keep CONUS + DC

map_df <- us_counties |>
  rename(fips = GEOID) |>
  left_join(by_total, by = "fips") |>
  mutate(arf_total = replace_na(arf_total, 0)) |>
  st_transform(5070)  # NAD83 / Conus Albers

# --- Crop to CONUS bbox (lon/lat), then to Albers ---
conus_bb <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66.5, ymax = 49), crs = 4326))
map_df_conus <- map_df |>
  st_transform(4326) |>
  st_crop(conus_bb) |>
  st_transform(5070)

# us48 counties (no AK/HI/PR), in lon/lat for coord_sf limits
us48_ll <- counties(cb = TRUE, year = 2020, class = "sf") %>%
  filter(!STATEFP %in% drop_states) %>%
  rename(fips = GEOID)

# join totals
map_ll <- us48_ll %>% left_join(by_total, by = "fips")

# bins
min_count_to_show <- 10
pos_breaks <- c(10,25,50,100,250,500,1000,2500,5000,10000, Inf)
pos_labels <- c("10-25","26-50","51-100","101-250",
                "251-500","501-1k","1k-2.5k","2.5k-5k","5k-10k",">10k")

map_ll <- map_ll %>%
  mutate(
    arf_zero = !is.na(arf_total) & arf_total == 0,
    arf_suppressed = !is.na(arf_total) & arf_total > 0 & arf_total < min_count_to_show,
    arf_display = if_else(arf_suppressed, NA_real_, arf_total),
    arf_bin  = cut(arf_display, breaks = pos_breaks, labels = pos_labels,
                   include.lowest = TRUE, right = TRUE),
    arf_cat  = case_when(
      arf_zero       ~ "0 (No ARF)",
      arf_suppressed ~ "<10",
      !is.na(arf_bin) ~ as.character(arf_bin),
      TRUE           ~ NA_character_
    ),
    arf_cat = factor(arf_cat, levels = c("0 (No ARF)", "<10", pos_labels))
  )

stopifnot(!any(map_ll$arf_suppressed & !is.na(map_ll$arf_display)))

count_pal <- viridis(length(pos_labels))
pal <- c("#e6e6e6", "#24002f", count_pal)

p_arf <- ggplot(map_ll) +
  geom_sf(aes(fill = arf_cat), color = "white", linewidth = 0.05) +
  scale_fill_manual(
    name = "Number of Cases",
    values = pal,
    drop = FALSE,
    na.translate = FALSE,
    guide = guide_legend(reverse = TRUE)
  ) +
  coord_sf(
    xlim = c(-125, -66.5), ylim = c(24, 49), expand = FALSE
  ) +
  labs(
    title = "Total Acute Respiratory Failure Cases by U.S. County"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 14, hjust = 0, margin = margin(t = 4, b = 8)),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.margin = margin(8, 8, 8, 8)
  )
p_arf

ggsave(
  filename = file.path(dir_out, "arf_counts_total_by_county_map_gg.png"),
  plot = p_arf, width = 14, height = 9, dpi = 300, bg = "white"
)

sum(by_total$arf_total)


dir_tab1 <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/tab1"
out_dir  <- dir_tab1

# --- helpers --------------------------------------------------------------

site_from <- function(path) {
  m <- stringr::str_match(basename(path), "refer_([^_]+)_")[,2]
  ifelse(is.na(m), tools::file_path_sans_ext(basename(path)), m)
}

# --- replace your read_site() with this version ---
read_site <- function(path){
  # site label from filename
  m <- stringr::str_match(basename(path), "refer_([^_]+)_")[,2]
  site <- ifelse(is.na(m), tools::file_path_sans_ext(basename(path)), m)
  
  h   <- xml2::read_html(path)
  
  # keep leading spaces/NBSP -> trim = FALSE
  tab <- rvest::html_element(h, "table") |> rvest::html_table(header = TRUE, trim = FALSE)
  tab <- janitor::clean_names(tab)
  if (!"value" %in% names(tab)) names(tab)[2] <- "value"
  
  raw_label <- tab[[1]]
  value     <- tab[[2]]
  
  # Replace NBSP with normal space, but DO NOT trim yet
  raw_label <- gsub("\u00A0", " ", raw_label, fixed = TRUE)
  
  # Drop header/footnote rows by content
  kill <- grepl("^\\s*Characteristic\\s*$", raw_label) |
    grepl("^\\s*N\\s*=\\s*", raw_label) |
    grepl("Mean\\s*[±\\u00B1]|\\(Median;\\s*Q1,\\s*Q3\\)", raw_label)
  raw_label <- raw_label[!kill]
  value     <- value[!kill]
  
  # Subrow detection BEFORE squishing (2+ leading spaces)
  is_sub <- grepl("^\\s{2,}", raw_label)
  
  # Clean text without using a placeholder
  clean_text <- function(x){
    x <- sub("^[[:space:]]+", "", x)
    x <- gsub("\u00A0", " ", x, fixed = TRUE)
    x <- stringr::str_replace_all(x, "NO\\s*₂|NO₂", "NO2")  # normalize NO2
    x <- stringr::str_replace_all(x, "[\u2013\u2014]", "-") # en/em dash -> hyphen
    x <- stringr::str_squish(x)
    x
  }
  label_clean <- clean_text(raw_label)
  
  # Carry forward the last parent to build "parent — child"
  parent <- character(length(label_clean))
  last_parent <- NA_character_
  for (i in seq_along(label_clean)) {
    if (!is_sub[i]) last_parent <- label_clean[i]
    parent[i] <- last_parent
  }
  label <- ifelse(is_sub, paste0(parent, " — ", label_clean), label_clean)
  
  tibble::tibble(
    site  = site,
    label = label,
    value = stringr::str_squish(value)
  )
}


# --- read all HTMLs (unchanged) ---
files_html <- list.files(dir_tab1, pattern = "\\.html$", full.names = TRUE)

tbls <- purrr::map_dfr(files_html, read_site)

# Recode 0/1 → No/Yes; Vasopressor Unknown → No
tbls <- tbls %>%
  mutate(label = dplyr::case_when(
    grepl("^In-hospital death — 0$", label) ~ "In-hospital death — No",
    grepl("^In-hospital death — 1$", label) ~ "In-hospital death — Yes",
    grepl("^30-day death — 0$", label)      ~ "30-day death — No",
    grepl("^30-day death — 1$", label)      ~ "30-day death — Yes",
    grepl("^AKI — 0$", label)               ~ "AKI — No",
    grepl("^AKI — 1$", label)               ~ "AKI — Yes",
    grepl("^Vasopressor use — 1$", label)   ~ "Vasopressor use — Yes",
    grepl("^Vasopressor use — Unknown$", label) ~ "Vasopressor use — No",
    TRUE ~ label
  ))

# --- collapse duplicates (first non-empty) and pivot wide ---
tbls_clean <- tbls %>%
  group_by(site, label) %>%
  summarise(value = dplyr::first(value[value != "" & !is.na(value)]), .groups = "drop")

wide <- tbls_clean %>%
  mutate(site = factor(site, levels = sort(unique(site)))) %>%
  tidyr::pivot_wider(names_from = site, values_from = value)

# --- pooled Overall (same logic as before) ---
site_N <- purrr::map_dfr(files_html, function(p){
  doc    <- xml2::read_html(p)
  th_txt <- rvest::html_elements(doc, "th") |> rvest::html_text2()
  n_txt  <- th_txt[stringr::str_detect(th_txt, "^N\\s*=\\s*")]
  
  N <- if (length(n_txt) > 0) {
    num <- stringr::str_extract(n_txt[1], "\\d[\\d,]+")
    num <- stringr::str_remove_all(num, ",")
    as.numeric(num)
  } else {
    NA_real_
  }
  
  tibble::tibble(site = site_from(p), N = N)
})
overall_N <- sum(site_N$N, na.rm = TRUE)

parse_n <- function(x){
  if (is.na(x) || x == "") return(NA_real_)
  m <- stringr::str_match(x, "([0-9][0-9,]*)")[,2]
  ifelse(is.na(m), NA_real_, as.numeric(gsub(",", "", m)))
}
is_cont_row <- function(vs) any(stringr::str_detect(vs, "±|\\u00B1"), na.rm = TRUE)
site_cols <- setdiff(names(wide), "label")
site_N_vec <- setNames(site_N$N, site_N$site)

overall <- apply(wide[site_cols], 1, function(vals){
  if (is_cont_row(vals)) {
    msd <- t(vapply(vals, function(x){
      if (is.na(x) || x == "") return(c(NA_real_, NA_real_))
      m <- stringr::str_match(x, "([0-9.]+)\\s*[±\\u00B1]\\s*([0-9.]+)")
      if (all(is.na(m))) return(c(NA_real_, NA_real_))
      as.numeric(m[2:3])
    }, numeric(2)))
    df <- data.frame(site = names(vals), mean = msd[,1], sd = msd[,2]) %>%
      dplyr::left_join(site_N, by = "site") %>% tidyr::drop_na()
    if (nrow(df) == 0 || sum(df$N) <= 1) return(NA_character_)
    wmean <- sum(df$N * df$mean) / sum(df$N)
    num <- sum((df$N - 1) * (df$sd^2)) + sum(df$N * (df$mean - wmean)^2)
    wsd <- sqrt(pmax(num / (sum(df$N) - 1), 0))
    sprintf("%.1f ± %.1f", wmean, wsd)
  } else {
    tot <- sum(vapply(vals, parse_n, NA_real_), na.rm = TRUE)
    sprintf("%s (%.1f%%)", scales::comma(tot), 100 * tot / overall_N)
  }
})

wide_out <- wide %>% mutate(Overall = overall) %>% relocate(Overall, .after = dplyr::last_col())
readr::write_csv(wide_out, file.path(out_dir, "Table1_sites_plus_overall_from_html_fixed.csv"))





dir_in  <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis"
dir_out <- dir_in

`%||%` <- function(a,b) if (!is.null(a)) a else b

site_from <- function(path) {
  m <- stringr::str_match(basename(path), "refer_([^_]+)_")[,2]
  ifelse(is.na(m), tools::file_path_sans_ext(basename(path)), m)
}

# safer reader
read_one <- function(f) {
  x <- try(readr::read_csv(f, show_col_types = FALSE), silent = TRUE)
  if (inherits(x, "try-error") || is.null(x)) return(tibble())   # skip unreadable files
  df <- try(as_tibble(janitor::clean_names(x)), silent = TRUE)
  if (inherits(df, "try-error") || !is.data.frame(df)) return(tibble())
  
  n <- nrow(df)
  df <- mutate(df, site = site_from(f))
  
  get_chr <- function(cands) {
    nm <- intersect(cands, names(df))
    if (length(nm) == 0) return(rep(NA_character_, n))
    as.character(df[[nm[1]]])
  }
  pick_num <- function(cands) {
    nm <- intersect(cands, names(df))
    if (length(nm) == 0) return(rep(NA_real_, n))
    suppressWarnings(as.numeric(df[[nm[1]]]))
  }
  
  tibble(
    site      = df$site,
    outcome   = get_chr(c("outcome","cause","event")),
    term      = get_chr(c("term","label","model")),
    is_expo   = if ("exposure" %in% names(df)) as.character(df$exposure) else rep(NA_character_, n),
    SHR       = pick_num(c("shr","hr","estimate","exp_estimate")),
    conf_low  = pick_num(c("conf_low","ci_low","ci_lwr","lower","lcl","lo")),
    conf_high = pick_num(c("conf_high","ci_high","ci_upr","upper","ucl","hi")),
    p_value   = pick_num(c("p_value","p","pval","pvalue"))
  )
}

# file pattern: match “…cumul…” OR “…cumulative…”, then .csv
files <- list.files(
  dir_in,
  pattern = "(?:_shr_.*cumul|cumulative).*\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)
stopifnot(length(files) > 0)

raw <- map_dfr(files, read_one) |>
  filter(!is.na(SHR), !is.na(conf_low), !is.na(conf_high)) |>
  mutate(
    outcome = as.character(outcome),
    term    = as.character(term)
  )

# -------- keep only exposure rows and only NO2 / PM2.5
sel <- raw %>%
  # normalize is_expo -> logical
  mutate(
    is_expo_chr = str_to_lower(as.character(is_expo)),
    is_expo = case_when(
      is.na(is_expo_chr)                                  ~ NA,          # NA logical
      is_expo_chr %in% c("true","t","1","yes","y")        ~ TRUE,
      is_expo_chr %in% c("false","f","0","no","n")        ~ FALSE,
      TRUE                                                ~ NA
    )
  ) %>%
  select(-is_expo_chr) %>%
  filter(is.na(is_expo) | is_expo) %>%
  # detect pollutant + parse scaling
  mutate(
    term = as.character(term),
    term_l = str_to_lower(term),
    pollutant = case_when(
      str_detect(term_l, "pm\\s*2\\.?\\s*5|pm2\\.?5|pm25|pm₂\\.?₅|pm\\s*2,\\s*5") ~ "PM2.5",
      str_detect(term_l, "no\\s*2|no2|no₂") ~ "NO2",
      TRUE ~ NA_character_
    ),
    per_amount = suppressWarnings(as.numeric(str_match(term_l, "per\\s*([0-9]+\\.?[0-9]*)")[,2]))
  ) %>%
  filter(!is.na(pollutant))

# -------- harmonize scales: NO2 -> per 10 ppb, PM2.5 -> per 10 µg/m³
# If a row is already on that scale, it's unchanged. If per 5, we'll convert.
target_per <- tibble(pollutant = c("NO2","PM2.5"),
                     target_amt = c(10, 10))

sel2 <- sel %>%
  left_join(target_per, by = "pollutant") %>%
  mutate(
    # default: if we couldn't parse an amount, assume already on target to avoid dropping
    scale_from = ifelse(is.na(per_amount) | per_amount <= 0, target_amt, per_amount),
    scale_factor = target_amt / scale_from,
    # rescale SHR and CI: raise to power 'scale_factor' on the ratio scale
    SHR_rescaled       = SHR^scale_factor,
    conf_low_rescaled  = conf_low^scale_factor,
    conf_high_rescaled = conf_high^scale_factor,
    exposure_label = case_when(
      pollutant == "NO2"   ~ "NO2 (per 10 ppb)",
      pollutant == "PM2.5" ~ "PM2.5 (per 5 µg/m³)"
    )
  )

# -------- clean outcomes to friendly labels
sel2 <- sel2 %>%
  mutate(
    outcome_clean = case_when(
      str_detect(str_to_lower(outcome), "extub") | str_detect(outcome, "cause\\s*=\\s*1") ~ "Successful extubation",
      str_detect(str_to_lower(outcome), "prf|persist") | str_detect(outcome, "cause\\s*=\\s*3") ~ "Persistent respiratory failure",
      str_detect(str_to_lower(outcome), "death") | str_detect(outcome, "cause\\s*=\\s*2") ~ "Death",
      TRUE ~ outcome
    )
  )

# -------- meta-analysis per outcome x pollutant (REML)
meta_df <- sel2 %>%
  transmute(site, outcome = outcome_clean, exposure = exposure_label,
            est = SHR_rescaled, lo = conf_low_rescaled, hi = conf_high_rescaled) %>%
  filter(!is.na(est), !is.na(lo), !is.na(hi)) %>%
  mutate(yi = log(est),
         sei = (log(hi) - log(lo))/(2*1.96)) %>%
  filter(is.finite(yi), is.finite(sei), sei > 0)

pooled <- meta_df %>%
  group_by(outcome, exposure) %>%
  group_modify(~{
    if (nrow(.x) < 2) return(tibble(k = nrow(.x), pooled = NA_real_, lcl = NA_real_, ucl = NA_real_))
    fit <- rma(yi = .x$yi, sei = .x$sei, method = "REML")
    tibble(k = nrow(.x), pooled = exp(fit$b[1]), lcl = exp(fit$ci.lb), ucl = exp(fit$ci.ub))
  }) %>% ungroup()

# -------- final table: sites + overall
fmt <- function(e,l,u) sprintf("%.2f (%.2f, %.2f)", e,l,u)

tab_sites <- meta_df %>%
  transmute(Outcome = outcome, Exposure = exposure, Site = site,
            `SHR (95% CI)` = fmt(exp(yi), exp(yi - 1.96*sei), exp(yi + 1.96*sei)))

tab_overall <- pooled %>%
  transmute(Outcome = outcome, Exposure = exposure, Site = "Overall (REML)",
            `SHR (95% CI)` = fmt(pooled, lcl, ucl))

tab_out <- bind_rows(tab_sites, tab_overall) %>%
  arrange(Outcome, Exposure, Site)

readr::write_csv(tab_out, file.path(dir_out, "CIF_SHR_PM25_NO2_by_site_and_overall_rescaled.csv"))

# -------- forest plot (by outcome; columns = pollutant)

plot_df <- tab_sites %>%
  rename(outcome = Outcome, exposure = Exposure, site = Site) %>%
  left_join(meta_df %>% select(site, outcome, exposure, yi, sei), by = c("site","outcome","exposure")) %>%
  transmute(outcome, exposure, site,
            est = exp(yi), lo = exp(yi - 1.96*sei), hi = exp(yi + 1.96*sei), pooled = FALSE) %>%
  bind_rows(
    pooled %>% transmute(outcome, exposure, site = "Overall (REML)",
                         est = pooled, lo = lcl, hi = ucl, pooled = TRUE)
  )

site_levels <- c(sort(unique(plot_df$site[plot_df$site != "Overall (REML)"])), "Overall (REML)")
plot_df <- plot_df %>%
  mutate(site = factor(site, levels = site_levels),
         exposure = factor(exposure, levels = c("PM2.5 (per 5 µg/m³)", "NO2 (per 10 ppb)")),
         outcome  = factor(outcome, levels = c("Death","Successful extubation","Persistent respiratory failure")))

p <- ggplot(plot_df, aes(x = est, y = site, xmin = lo, xmax = hi,
                         shape = pooled, color = pooled)) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_pointrange(size = 0.32, na.rm = TRUE) +
  scale_x_log10(breaks = c(0.7, 1, 1.3, 1.6, 2.0),
                labels = scales::label_number(accuracy = 0.01)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18),
                     labels = c("Site","Overall (REML)"), name = NULL) +
  scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "#3366CC"),
                     labels = c("Site","Overall (REML)"), name = NULL) +
  labs(x = "Subdistribution Hazard Ratio (log scale)", y = NULL) +
  facet_grid(outcome ~ exposure, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(dir_out, "CIF_SHR_PM25_NO2_forest_rescaled.png"),
       p, width = 12, height = 8, dpi = 300)

p
# --- Recompute pooled stats with I2 (per outcome x exposure) ---
pooled_stats <- meta_df %>%
  group_by(outcome, exposure) %>%
  group_modify(~{
    if (nrow(.x) < 2) return(tibble(k = nrow(.x), pooled = NA_real_, lcl = NA_real_, ucl = NA_real_, tau2 = NA_real_, I2 = NA_real_))
    fit <- rma(yi = .x$yi, sei = .x$sei, method = "REML")
    tibble(
      k      = nrow(.x),
      pooled = exp(fit$b[1]),
      lcl    = exp(fit$ci.lb),
      ucl    = exp(fit$ci.ub),
      tau2   = fit$tau2,
      I2     = fit$I2
    )
  }) %>% ungroup()

# Save a clean pooled summary table (nice for supplement)
fmt <- function(e,l,u) sprintf("%.2f (%.2f, %.2f)", e,l,u)
pooled_table <- pooled_stats %>%
  transmute(Outcome = outcome, Exposure = exposure, Sites = k, `Pooled SHR (95% CI)` = fmt(pooled, lcl, ucl), `I^2 (%)` = round(I2, 1))
write_csv(pooled_table, file.path(dir_out, "CIF_SHR_PM25_NO2_overall_REML_I2.csv"))

# --- Build plotting data (sites + pooled) ---
sites_df <- meta_df %>%
  transmute(outcome, exposure, site, est = exp(yi), lo = exp(yi - 1.96*sei), hi = exp(yi + 1.96*sei), type = "Site")

pool_df  <- pooled_stats %>%
  transmute(outcome, exposure, site = "Overall (REML)", est = pooled, lo = lcl, hi = ucl, type = "Overall (REML)")

plot_df <- bind_rows(pool_df, sites_df)

# X limits based on data (log scale)
xr <- range(c(plot_df$lo, plot_df$hi), na.rm = TRUE)
lims <- c(max(0.6, xr[1]*0.95), xr[2]*1.05)
# nice wide, fixed axis for all facets
site_levels <- c("Overall (REML)", sort(unique(plot_df$site[plot_df$site != "Overall (REML)"])))
plot_df <- plot_df %>%
  mutate(
    site     = factor(site, levels = site_levels),
    exposure = factor(exposure, levels = c("PM2.5 (per 5 µg/m³)", "NO2 (per 10 ppb)")),
    outcome  = factor(outcome, levels = c("Death","Successful extubation","Persistent respiratory failure")),
    type     = factor(type, levels = c("Site","Overall (REML)"))
  )

# fixed wide axis for all facets
# 1) Clean & order for plotting (per outcome × exposure facet)
plot_df_clean <- plot_df %>%
  # rename here (handle case variants just in case)
  mutate(site = if_else(str_to_lower(site) == "hopkins", "JHU", site)) %>%
  filter(!is.na(site), !is.na(est), !is.na(lo), !is.na(hi)) %>%
  group_by(outcome, exposure) %>%
  arrange(desc(site == "Overall (REML)"), site, .by_group = TRUE) %>%
  mutate(site_f = factor(site, levels = rev(unique(site)))) %>%
  ungroup()

# 2) Facet-aware guide data for vertical lines
x_min <- 0.5
x_max <- 15.0
x_breaks <- c(0.5, 0.75, 1.00, 1.50, 2.00, 5.00, 15.00)

plot_df_clean <- plot_df_clean %>%
  mutate(
    exposure_facet = case_when(
      str_detect(exposure, "^NO2")   ~ "NO[2]~(per~10~ppb)",
      str_detect(exposure, "^PM2\\.5") ~ "PM[2.5]~(per~5~mu*g/m^3)",
      TRUE ~ exposure                 # fallback
    )
  )

pf <- plot_df_clean %>%
  mutate(
    lo_clip  = pmax(lo,  x_min),
    hi_clip  = pmin(hi,  x_max),
    est_clip = pmin(pmax(est, x_min), x_max),
    left_out  = lo < x_min,
    right_out = hi > x_max
  )

# facet keys for guide lines
panel_keys   <- dplyr::distinct(pf, outcome, exposure)
panel_guides <- tidyr::crossing(panel_keys, x = setdiff(x_breaks, 1))

# pf <- pf %>%
#   mutate(site = if_else(site == "hopkins", "JHU", site))
# pf <- pf %>%
#   mutate(site = if_else(site == "Overall (REML)", "Overall", site))

# ---- 3) Plot with clipped segments + edge arrows ----
p_pub <- ggplot(pf, aes(x = est_clip, y = site_f)) +
  facet_grid(
    outcome ~ exposure_facet,
    scales = "free_y",
    labeller = labeller(exposure_facet = label_parsed)
  ) +
  
  # reference lines
  geom_segment(data = dplyr::mutate(panel_keys, x = 1),
               aes(x = x, xend = x, y = -Inf, yend = Inf),
               linetype = 2, linewidth = 0.6, color = "grey35", inherit.aes = FALSE) +
  geom_segment(data = panel_guides,
               aes(x = x, xend = x, y = -Inf, yend = Inf),
               color = "grey92", linewidth = 0.4, inherit.aes = FALSE) +
  
  # --- SITE intervals (clipped) ---
  geom_segment(data = subset(pf, type == "Site"),
               aes(x = lo_clip, xend = hi_clip, yend = site_f),
               linewidth = 0.5, color = "grey30") +
  geom_point(data = subset(pf, type == "Site"),
             shape = 16, size = 1.9, color = "grey30") +
  # left/right arrows for truncated SITE CIs
  geom_segment(data = subset(pf, type == "Site" & left_out),
               aes(x = x_min, xend = x_min * 1.02, yend = site_f),
               arrow = arrow(length = unit(3, "pt"), ends = "first"),
               color = "grey30") +
  geom_segment(data = subset(pf, type == "Site" & right_out),
               aes(x = x_max, xend = x_max / 1.02, yend = site_f),
               arrow = arrow(length = unit(3, "pt"), ends = "last"),
               color = "grey30") +
  
  # --- POOLED intervals (clipped) ---
  geom_segment(data = subset(pf, type == "Overall (REML)"),
               aes(x = lo_clip, xend = hi_clip, yend = site_f),
               linewidth = 1.1, color = "#2C6BF6") +
  geom_point(data = subset(pf, type == "Overall (REML)"),
             shape = 23, size = 3.2, stroke = 0.6, color = "#2C6BF6", fill = "#2C6BF6") +
  # left/right arrows for truncated POOLED CIs
  geom_segment(data = subset(pf, type == "Overall (REML)" & left_out),
               aes(x = x_min, xend = x_min * 1.02, yend = site_f),
               arrow = arrow(length = unit(3, "pt"), ends = "first"),
               color = "#2C6BF6") +
  geom_segment(data = subset(pf, type == "Overall (REML)" & right_out),
               aes(x = x_max, xend = x_max / 1.02, yend = site_f),
               arrow = arrow(length = unit(3, "pt"), ends = "last"),
               color = "#2C6BF6") +
  
  scale_x_log10(limits = c(x_min, x_max), breaks = x_breaks,
                labels = scales::label_number(accuracy = 0.01)) +
  labs(x = "Subdistribution Hazard Ratio", y = NULL) +
  coord_cartesian(clip = "off") +  # ensure arrowheads at edges are visible
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 11),
    plot.margin = margin(6, 18, 6, 6)  # a little extra right space for arrows
  )

p_pub <- p_pub +
  theme(
    panel.border   = element_rect(colour = "#9aa0a6", fill = NA, linewidth = 0.7),
    panel.spacing  = unit(1.0, "lines"),           # a little breathing room between boxes
    strip.background = element_rect(fill = "#f7f7f7", colour = "#9aa0a6"),
    strip.text       = element_text(face = "bold")
  )

p_pub

ggsave(file.path(dir_out, "CIF_SHR_PM25_NO2_forest_pub.jpg"), p_pub, width = 10.5, height = 8.2, dpi = 600)
ggsave(file.path(dir_out, "CIF_SHR_PM25_NO2_forest_pub.pdf"), p_pub, width = 10.5, height = 8.2)

# ---- Build pooled stats if not present ----
if (!exists("pooled_stats")) {
  stopifnot(exists("meta_df"))  # meta_df must have yi, sei, outcome, exposure
  pooled_stats <- meta_df %>%
    group_by(outcome, exposure) %>%
    group_modify(~{
      if (nrow(.x) < 2) return(tibble(k = nrow(.x), pooled = NA_real_, lcl = NA_real_, ucl = NA_real_, I2 = NA_real_))
      fit <- rma(yi = .x$yi, sei = .x$sei, method = "REML")
      tibble(k = nrow(.x), pooled = exp(fit$b[1]), lcl = exp(fit$ci.lb), ucl = exp(fit$ci.ub), I2 = fit$I2)
    }) %>% ungroup()
}

# ---- Keep only Overall per outcome × pollutant and format labels ----
# Build the compact "overall only" data frame
overall_df <- pooled_stats %>%
  transmute(
    outcome,
    exposure,
    pollutant = ifelse(grepl("^PM", exposure), "PM2.5", "NO2"),
    est = pooled, lo = lcl, hi = ucl
  ) %>%
  filter(!is.na(est), !is.na(lo), !is.na(hi)) %>%
  mutate(
    outcome = factor(outcome,
                     levels = c("Death", "Successful extubation", "Persistent respiratory failure")),
    pollutant = factor(pollutant, levels = c("PM2.5", "NO2")),
    label_ci = sprintf("%.2f (%.2f–%.2f)", est, lo, hi)
  )

# 1) Make clean labels WITHOUT k / I²
# Subscripted labels
lab_PM25 <- "PM\u2082.\u2085"   # PM₂.₅
lab_NO2  <- "NO\u2082"         # NO₂

# (optional) keep a labeled version if you ever need to print the pollutant name in text
overall_lab <- overall_df %>%
  mutate(pollutant_lab = dplyr::recode(pollutant,
                                       "PM2.5" = lab_PM25,
                                       "NO2"   = lab_NO2))

# 2) Axis and guides
x_min <- 0.25
x_max <- max(overall_lab$hi, na.rm = TRUE) * 1.2
x_breaks_all <- c(0.25, 0.5, 1.0, 1.5, 2, 3, 5, 10)
x_breaks <- x_breaks_all[x_breaks_all >= x_min & x_breaks_all <= x_max]
guide_df <- data.frame(x = setdiff(x_breaks, 1))

# dodge for PM2.5 vs NO2
pd <- position_dodge(width = 0.55)

# 3) Single-panel JAMA-style plot
p_overall <- ggplot(overall_lab,
                    aes(x = est, y = outcome, xmin = lo, xmax = hi,
                        shape = pollutant, color = pollutant)) +
  # faint vertical guides + null line
  geom_vline(data = guide_df, aes(xintercept = x),
             color = "grey92", linewidth = 0.4, inherit.aes = FALSE) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.6, color = "grey35") +
  
  # CIs + points
  geom_errorbarh(height = 0, size = 0.9, lineend = "butt", position = pd) +
  geom_point(size = 3.1, position = pd) +
  
  # labels in the right gutter (never cover CIs)
  geom_text(
    data = overall_lab,
    aes(x = Inf, y = outcome, label = label_ci, group = pollutant),
    position = pd, hjust = -0.04, size = 3.6, color = "black"
  ) +
  
  # grayscale legend
  scale_shape_manual(values = c("PM2.5" = 16, "NO2" = 17),
                     labels = c("PM2.5" = lab_PM25, "NO2" = lab_NO2),
                     name = NULL) +
  scale_color_manual(values = c("PM2.5" = "black", "NO2" = "grey30"),
                     guide = "none") +
  
  # wider right side; avoid crowded tick labels
  scale_x_log10(limits = c(x_min, x_max), breaks = x_breaks,
                labels = scales::label_number(accuracy = 0.01), expand = c(0, 0)) +
  
  labs(x = "Subdistribution Hazard Ratio", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(margin = margin(t = 2)),
    axis.title.x = element_text(size = 11),
    # put legend at the top; horizontal; tight box so it won't touch labels
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box.margin = margin(b = 2),
    # give lots of room on the right for the gutter labels
    plot.margin  = margin(5, 210, 8, 5)
  ) +
  # if tick labels are still tight, let ggplot automatically drop overlaps
  guides(x = guide_axis(check.overlap = TRUE))

shape_labs <- list(expression(PM[2.5]), expression(NO[2]))

p_overall <- p_overall +
  scale_shape_manual(
    values = c("PM2.5" = 16, "NO2" = 17),
    labels = shape_labs,
    name = NULL
  ) +
  # keep color grayscale, no color legend
  scale_color_manual(values = c("PM2.5" = "black", "NO2" = "grey30"), guide = "none")

p_overall <- p_overall +
  theme(
    # put legend above the panel, tight spacing
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.justification = "left",
    legend.box.spacing = unit(2, "pt"),      # distance between legend and panel
    legend.box.margin  = margin(0, 0, 4, 0), # outside of legend box
    legend.margin      = margin(2, 6, 2, 6), # inside padding of the legend box
    
    # draw a box around the legend
    legend.box.background = element_rect(colour = "#6e7781", fill = "white", linewidth = 0.6),
    legend.background     = element_rect(fill = "white", colour = NA),
    legend.key            = element_rect(fill = "white", colour = NA),
    
    # slightly reduce top plot margin so the legend sits closer
    plot.margin = margin(2, 210, 8, 5)
  )

p_overall

ggsave(file.path(out_dir, "overall_SHR_PM25_NO2_JAMA_singlepanel.png"),
       p_overall, width = 8.2, height = 4.8, dpi = 600)
ggsave(file.path(out_dir, "overall_SHR_PM25_NO2_JAMA_singlepanel.pdf"),
       p_overall, width = 8.2, height = 4.8)


#### Model aggregation

# Point to your folder of CSVs
dir_models <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/models"

# --- 1) Ingest all CSVs and parse metadata from filenames ---------------------
# Expected filename pattern examples (order can vary slightly):
# refer_Emory_20251008_model_tidy_death30d_cumulative_no2_only_20251008_143656.csv
# refer_Michigan_20251009_model_tidy_inhosp_cumulative_pm25_only_20251009_161603.csv
# refer_OHSU_20251014_model_tidy_iculos_cumulative_no2_only_20251014_123455.csv
# refer_RUSH_20251007_model_tidy_venthrs_cumulative_no2_only_20251007_081351.csv

files <- list.files(dir_models, pattern = "\\.csv$", full.names = TRUE)

# Helper to parse site, outcome, pollutant from the filename
parse_meta <- function(path) {
  fn <- basename(path)
  
  site <- str_match(fn, "^refer_([^_]+)_")[,2]
  
  outcome_token <- tolower(str_match(fn, "model_tidy_([^_]+)")[,2])
  
  outcome <- case_when(
    outcome_token %in% c("death30d", "death_30d") ~ "death30d",
    outcome_token %in% c("inhosp", "in_hosp", "inhosp_death") ~ "inhosp_death",
    outcome_token %in% c("iculos", "icu_los") ~ "icu_los",
    outcome_token %in% c("venthrs", "vent_hours", "venthr") ~ "vent_hours",
    TRUE ~ outcome_token
  )
  
  pollutant <- case_when(
    str_detect(fn, regex("no2",  ignore_case = TRUE)) ~ "NO2",
    str_detect(fn, regex("pm25", ignore_case = TRUE)) ~ "PM25",
    TRUE ~ NA_character_
  )
  
  tibble(site, outcome, pollutant, file = path)
}

# Rebuild meta_df
files <- list.files(dir_models, pattern = "\\.csv$", full.names = TRUE)

meta_df <- map_dfr(files, parse_meta) %>%
  filter(!is.na(pollutant), !is.na(outcome))

# quick check
meta_df %>% count(pollutant, outcome)

# --- 2) Read each CSV and extract EXPOSURE row only ---------------------------
# We keep only the exposure term row (no2_* or pm25_*), ignoring factor terms
read_and_tag <- function(file, site, outcome, pollutant) {
  df <- suppressMessages(read_csv(file, show_col_types = FALSE))
  
  # Required columns in your screenshot:
  # term, estimate, std.error, statistic, p.value, conf.low, conf.high
  # We’ll be permissive about column names just in case:
  nm <- names(df)
  stopifnot(all(c("term","estimate","conf.low","conf.high","p.value") %in% nm))
  
  # pick exposure row
  exp_row <- df %>%
    filter(
      (pollutant == "NO2"  & str_detect(term, regex("no2",  ignore_case = TRUE))) |
        (pollutant == "PM25" & str_detect(term, regex("pm\\s*25|pm25", ignore_case = TRUE)))
    ) %>%
    # avoid factor() rows etc.
    filter(!str_detect(term, "^factor\\(")) %>%
    slice_head(n = 1)
  
  if (nrow(exp_row) == 0) return(NULL)
  
  exp_row %>%
    transmute(
      site, outcome, pollutant, term,
      estimate = as.numeric(estimate),
      conf.low = as.numeric(conf.low),
      conf.high= as.numeric(conf.high),
      p.value  = as.numeric(p.value)
    )
}

per_site <- pmap_dfr(
  meta_df,
  ~ read_and_tag(..4, ..1, ..2, ..3)  # file, site, outcome, pollutant
)

# Safety check
if (nrow(per_site) == 0) stop("No exposure rows found. Check filename patterns/term names.")

# --- 3) Compute log-effect and SE (works whether estimate is OR/HR or beta) ---
# If the CSV "estimate" is already exponentiated (as in your screenshot),
# using CI to back-calc SE on the log scale is safest.
per_site <- per_site %>%
  mutate(
    yi  = log(estimate),                                 # log effect
    sei = (log(conf.high) - log(conf.low)) / (2*1.96)    # SE from CI
  )

# --- 4) Meta-analysis across sites for each pollutant × outcome ---------------
pool_one <- function(df_group) {
  # Guard against bad rows
  df_group <- df_group %>% filter(is.finite(yi), is.finite(sei), sei > 0)
  
  if (nrow(df_group) < 2) {
    return(tibble(
      k           = nrow(df_group),
      method      = "insufficient_k",
      pooled_est  = NA_real_,
      pooled_lcl  = NA_real_,
      pooled_ucl  = NA_real_,
      tau2        = NA_real_,
      Q           = NA_real_,
      Q_p         = NA_real_
    ))
  }
  
  # Fixed-effect
  fe <- tryCatch(rma.uni(yi = yi, sei = sei, data = df_group, method = "FE"),
                 error = function(e) NULL)
  # Random-effects (DL)
  re <- tryCatch(rma.uni(yi = yi, sei = sei, data = df_group, method = "DL"),
                 error = function(e) NULL)
  
  bind_rows(
    if (!is.null(fe)) tibble(
      k          = fe$k,
      method     = "fixed",
      pooled_est = exp(as.numeric(fe$b)),
      pooled_lcl = exp(as.numeric(fe$ci.lb)),
      pooled_ucl = exp(as.numeric(fe$ci.ub)),
      tau2       = 0,
      Q          = as.numeric(fe$QE),
      Q_p        = as.numeric(fe$QEp)
    ) else NULL,
    if (!is.null(re)) tibble(
      k          = re$k,
      method     = "random_DL",
      pooled_est = exp(as.numeric(re$b)),
      pooled_lcl = exp(as.numeric(re$ci.lb)),
      pooled_ucl = exp(as.numeric(re$ci.ub)),
      tau2       = as.numeric(re$tau2),
      Q          = as.numeric(re$QE),
      Q_p        = as.numeric(re$QEp)
    ) else NULL
  )
}

pooled <- per_site %>%
  group_by(pollutant, outcome) %>%
  group_modify(~ pool_one(.x)) %>%
  ungroup()

# --- 5) Write outputs ----------------------------------------------------------
out_dir <- file.path(dir_models, "aggregated_outputs")
dir.create(out_dir, showWarnings = FALSE)

# Per-site (keeps original OR/HR + CI and the log-scale pieces)
write_csv(per_site, file.path(out_dir, "per_site_effects_exposure_only.csv"))

# Pooled results (fixed + random)
write_csv(pooled,   file.path(out_dir, "pooled_meta_by_pollutant_outcome.csv"))

# --- 6) Forest plots (one PDF per pollutant × outcome) ------------------------
# Optional quick forests using metafor
plot_dir <- file.path(out_dir, "forest_plots")
dir.create(plot_dir, showWarnings = FALSE)

make_forest <- function(df_group, pollutant, outcome) {
  df_group <- df_group %>% arrange(site)
  if (nrow(df_group) < 2) return(invisible(NULL))
  
  # Random-effects model (DL)
  re <- tryCatch(rma.uni(yi = yi, sei = sei, data = df_group, method = "DL"),
                 error = function(e) NULL)
  if (is.null(re)) return(invisible(NULL))
  
  pdf(file.path(plot_dir, sprintf("forest_%s_%s.pdf", pollutant, outcome)),
      width = 7, height = 6)
  par(mar = c(4, 4, 2, 2))
  metafor::forest(
    x = re,
    slab = df_group$site,
    xlab = "Effect (ratio scale)",
    atransf = exp
  )
  title(main = sprintf("%s – %s (random DL)", pollutant, outcome))
  dev.off()
}

per_site %>%
  group_by(pollutant, outcome) %>%
  group_walk(~ make_forest(.x, .y$pollutant, .y$outcome))

# --- 7) Nice, human-readable summary table in the console ---------------------
pooled %>%
  arrange(pollutant, outcome, method) %>%
  mutate(ci = sprintf("%.3f (%.3f–%.3f)", pooled_est, pooled_lcl, pooled_ucl)) %>%
  select(pollutant, outcome, method, k, ci, tau2, Q, Q_p) %>%
  print(n = Inf)





# --- Same ordering & offsets as before ---
order_levels <- c("Ventilation hours", "ICU length of stay", "In-hospital death", "30-day death")

dfr <- pooled %>%
  dplyr::filter(method == "random_DL") %>%
  dplyr::mutate(
    pollutant = dplyr::recode(pollutant, NO2 = "NO2", PM25 = "PM25"),
    outcome = dplyr::recode(
      outcome,
      death30d      = "30-day death",
      inhosp_death  = "In-hospital death",
      icu_los       = "ICU length of stay",
      vent_hours    = "Ventilation hours"
    ),
    outcome = factor(outcome, levels = order_levels),
    lab = sprintf("%.2f (%.2f–%.2f)", pooled_est, pooled_lcl, pooled_ucl)
  )

y_base <- tibble::tibble(outcome = factor(order_levels, levels = order_levels),
                         y0 = seq_along(order_levels))

dfr <- dfr %>%
  dplyr::left_join(y_base, by = "outcome") %>%
  dplyr::mutate(
    # keep your slight vertical separation
    y = y0 + if_else(pollutant == "PM25", -0.12, +0.12)
  )

# --- Guides like Panel A ---
x_min <- 0.85
x_max <- 1.25
x_breaks <- seq(0.90, 1.20, by = 0.05)              # light guide lines
guide_df <- tibble::tibble(x = x_breaks)

shape_labs <- list(expression(PM[2.5]), expression(NO[2]))

pB_jama <- ggplot(dfr, aes(x = pooled_est, y = y,
                           xmin = pooled_lcl, xmax = pooled_ucl,
                           shape = pollutant, color = pollutant)) +
  # faint vertical guides + null line (match p_overall)
  geom_vline(data = guide_df, aes(xintercept = x),
             color = "grey92", linewidth = 0.4, inherit.aes = FALSE) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.6, color = "grey35") +
  
  # horizontal CIs + points at manual y
  geom_errorbarh(height = 0, size = 0.9, lineend = "butt") +
  geom_point(size = 3.1) +
  
  # right-gutter labels (never cover CIs)
  geom_text(
    data = dfr,
    aes(x = Inf, y = y, label = lab, group = pollutant),
    hjust = -0.04, size = 3.6, color = "black", inherit.aes = FALSE
  ) +
  
  # y labels at the center baseline between the two markers
  scale_y_continuous(
    breaks = y_base$y0, labels = levels(y_base$outcome),
    limits = c(0.5, max(y_base$y0) + 0.5)
  ) +
  
  # grayscale legend (boxed), no color legend
  scale_shape_manual(values = c("PM25" = 16, "NO2" = 17),
                     labels = shape_labs, name = NULL) +
  scale_color_manual(values = c("PM25" = "black", "NO2" = "grey30"), guide = "none") +
  
  # x scale & limits like p_overall (linear here; use log10 if your ests span wide)
  scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks,
                     labels = scales::label_number(accuracy = 0.01),
                     expand = c(0, 0)) +
  
  labs(x = "Pooled Effects Ratio", y = NULL) +
  coord_cartesian(clip = "off") +  # allow right-gutter labels
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(margin = margin(t = 2)),
    axis.title.x = element_text(size = 11),
    
    # legend style: top, horizontal, boxed — matching p_overall
    legend.position       = "top",
    legend.direction      = "horizontal",
    legend.justification  = "left",
    legend.box.spacing    = unit(2, "pt"),
    legend.box.margin     = margin(0, 0, 4, 0),
    legend.margin         = margin(2, 6, 2, 6),
    legend.box.background = element_rect(colour = "#6e7781", fill = "white", linewidth = 0.6),
    legend.background     = element_rect(fill = "white", colour = NA),
    legend.key            = element_rect(fill = "white", colour = NA),
    
    # roomy right gutter for labels (same vibe as p_overall)
    plot.margin = margin(2, 210, 8, 5),
    
    # no internal grid (classic look)
    panel.grid = element_blank()
  )

# --- Save (same dims as your p_overall export if you like) ---
out_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/output"
ggsave(file.path(out_dir, "panelB_pooled_meta_JAMAmatched.png"),
       pB_jama, width = 8.2, height = 4.8, dpi = 600)
ggsave(file.path(out_dir, "panelB_pooled_meta_JAMAmatched.pdf"),
       pB_jama, width = 8.2, height = 4.8)


###### ONE BIG CIF

# --- paths ---
cif_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/cif"
out_dir <- file.path(cif_dir, "fig_out"); dir.create(out_dir, showWarnings = FALSE)

# --- helpers ---
detect_pollutant <- function(path_chr) {
  nm <- tolower(path_chr)
  if (str_detect(nm, "no2")) "NO2"
  else if (str_detect(nm, "pm25") || str_detect(nm, "pm2_5") || str_detect(nm, "pm2.5")) "PM2.5"
  else NA_character_
}

cause_map <- c(`1`="Successful extubation", `2`="Death", `3`="Persistent respiratory failure")
expo_pal  <- c("Low"="grey70","Mid"="grey45","High"="black")

clean_names <- function(x) {
  x %>% tolower() %>% str_replace_all("[^a-z0-9]+","_") %>% str_replace("^_|_$","")
}

read_one_cif <- function(fp) {
  df_raw <- readr::read_csv(fp, show_col_types = FALSE)
  names(df_raw) <- clean_names(names(df_raw))
  
  df <- df_raw %>%
    rename(
      site_id        = any_of("site_id"),
      site_name      = any_of("site_name"),
      exposure_var   = any_of("exposure_var"),
      exposure_lab   = any_of("exposure_label"),
      unit           = any_of("unit"),
      exposure_group = any_of(c("exposure_group","exposure_grp","exposure_grc")),
      legend_label   = any_of("legend_label"),
      cause          = any_of("cause"),
      day            = any_of("day"),
      cif            = any_of("cif"),
      risk_set       = any_of(c("risk_set","riskset")),
      d1             = any_of("d1"),
      d2             = any_of("d2"),
      d3             = any_of("d3"),
      d0             = any_of("d0")
    )
  
  # enforce single site_name per file
  sn <- df %>% pull(site_name) %>% unique() %>% na.omit()
  sn <- if (length(sn) >= 1L) sn[1] else NA_character_
  df  <- df %>% mutate(site_name = sn)
  
  # derive exposure bin from label prefix (prefer legend_label)
  lab_src <- coalesce(df$legend_label, df$exposure_group) %>% tolower()
  grp     <- stringr::str_match(lab_src, "^(low|mid|high)")[,2]
  grp_fact <- dplyr::recode(grp, low="Low", mid="Mid", high="High", .default = NA_character_)
  
  df %>%
    mutate(
      exposure_bin = factor(grp_fact, levels = c("Low","Mid","High")),
      cause_label  = dplyr::recode(as.character(cause), !!!cause_map),
      pollutant    = detect_pollutant(fp),
      source_file  = fs::path_file(fp)
    ) %>%
    filter(!is.na(cause_label), !is.na(exposure_bin))
}

# --- discover CIF files (allow optional " copy N" suffix) ---
cif_regex <- "site_cif_plotdf__.*_cumulative(?: copy \\d+)?\\.csv$"
plot_files <- fs::dir_ls(cif_dir, regexp = cif_regex, type = "file", recurse = TRUE)

# quick check of what we found
message("Found ", length(plot_files), " CIF files")
if (length(plot_files)) print(fs::path_file(plot_files))

if (length(plot_files) == 0) stop("No CIF files found. Check cif_dir and filename patterns.")

# --- read all CIFs, one by one (site_name enforced per file) ---
cif_all <- purrr::map_dfr(plot_files, read_one_cif)

# sanity: which sites/pollutants did we load?
cif_all %>% count(site_name, pollutant) %>% arrange(site_name) %>% print(n = 100)

# --- save master combined CSV ---
readr::write_csv(cif_all, file.path(out_dir, "combined_cif_all_sites.csv"))

# --- plotting function ---
plot_cif_all_sites <- function(df, pol,
                               which_cause = c("Death","Successful extubation","Persistent respiratory failure"),
                               y_max = NULL) {
  dd <- df %>% dplyr::filter(pollutant == pol, cause_label %in% which_cause)
  if (!nrow(dd)) stop(glue("No rows for pollutant = {pol}."))
  
  if (is.null(y_max)) {
    y_max <- dd %>% dplyr::summarize(mx = max(cif, na.rm = TRUE)) %>% pull(mx)
    y_max <- ifelse(is.finite(y_max), ceiling(y_max*20)/20, 0.2)
  }
  
  ggplot(dd, aes(day, cif, color = exposure_bin)) +
    geom_step(linewidth = 0.7) +
    scale_color_manual(values = expo_pal, name = NULL, breaks = c("Low","Mid","High")) +
    facet_grid(cause_label ~ site_name, scales = "fixed") +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title = glue("{pol} cumulative exposure: CIF by site"),
         x = "Days since ICU admission", y = "Cumulative incidence") +
    theme_classic(base_size = 11) +
    theme(
      panel.border    = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
      strip.background= element_rect(fill = "grey95", colour = "grey80"),
      legend.position = "top",
      legend.direction= "horizontal",
      legend.box.margin = margin(b = 4),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4))
    )
}

# --- make & save figures ---
p_no2  <- plot_cif_all_sites(cif_all, "NO2")
p_pm25 <- plot_cif_all_sites(cif_all, "PM2.5")

ggsave(file.path(out_dir, "CIF_all_sites_NO2.png"),  p_no2,  width = 14, height = 8.5, dpi = 300)
ggsave(file.path(out_dir, "CIF_all_sites_PM25.png"), p_pm25, width = 14, height = 8.5, dpi = 300)

# risk-set–weighted mean CIF at each day, by pollutant × exposure_bin × outcome
pooled_cif <- cif_all %>%
  group_by(pollutant, cause_label, exposure_bin, day) %>%
  dplyr::summarize(
    n_sites    = n_distinct(site_name),
    total_at_risk = sum(risk_set, na.rm = TRUE),
    cif        = weighted.mean(cif, w = pmax(risk_set, 0), na.rm = TRUE),
    .groups = "drop"
  )

# quick sanity: how many sites contributed per curve?
pooled_cif %>%
  distinct(pollutant, cause_label, exposure_bin, n_sites) %>%
  arrange(cause_label, pollutant, exposure_bin) %>%
  print(n = 50)

# ---------- PLOTTING (combined pollutants) ----------
pal_bins <- c("Low"="grey70","Mid"="grey45","High"="black")
lt_poll  <- c("NO2"="solid","PM2.5"="22")  # dashed for PM2.5, solid for NO2

plot_outcome_combined <- function(df, outcome, y_max = NULL) {
  dd <- df %>% filter(cause_label == outcome)
  
  if (is.null(y_max)) {
    y_max <- dd %>% dplyr::summarize(mx = max(cif, na.rm = TRUE)) %>% pull(mx)
    y_max <- ifelse(is.finite(y_max), ceiling(y_max*20)/20, 0.2)
  }
  
  ggplot(dd, aes(x = day, y = cif,
                 color = exposure_bin, linetype = pollutant, group = interaction(pollutant, exposure_bin))) +
    geom_step(linewidth = 1) +
    scale_color_manual(values = pal_bins, name = NULL, breaks = c("Low","Mid","High")) +
    scale_linetype_manual(values = lt_poll, name = NULL) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(
      title = outcome,
      x = "Days since ICU admission",
      y = "Cumulative incidence",
      caption = "Curves pooled across all sites; weighted by daily risk set"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box.margin = margin(b = 4),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4))
    )
}

p_death   <- plot_outcome_combined(pooled_cif, "Death")
p_extub   <- plot_outcome_combined(pooled_cif, "Successful extubation")
p_prf     <- plot_outcome_combined(pooled_cif, "Persistent respiratory failure")

# ---------- SAVE ----------
out_dir <- file.path(cif_dir, "fig_out"); dir.create(out_dir, showWarnings = FALSE)
ggsave(file.path(out_dir, "CIF_pooled_all_sites__Death_PM25_NO2.png"),   p_death, width = 7.5, height = 5, dpi = 300)
ggsave(file.path(out_dir, "CIF_pooled_all_sites__Extubation_PM25_NO2.png"), p_extub, width = 7.5, height = 5, dpi = 300)
ggsave(file.path(out_dir, "CIF_pooled_all_sites__PersistentRF_PM25_NO2.png"), p_prf, width = 7.5, height = 5, dpi = 300)


pooled_cif <- cif_all %>%
  dplyr::group_by(pollutant, cause_label, exposure_bin, day) %>%
  dplyr::summarize(
    n_sites        = dplyr::n_distinct(site_name),
    total_at_risk  = sum(risk_set, na.rm = TRUE),
    cif            = stats::weighted.mean(cif, w = pmax(risk_set, 0), na.rm = TRUE),
    .groups = "drop"
  )

smooth_pooled_cif <- function(df, span = 0.25, day_max = 30) {
  df %>%
    dplyr::filter(day <= day_max) %>%
    dplyr::arrange(pollutant, cause_label, exposure_bin, day) %>%
    dplyr::group_by(pollutant, cause_label, exposure_bin) %>%
    dplyr::group_modify(~{
      d <- .x
      # LOESS fit (skip if < 5 points)
      if (nrow(d) >= 5) {
        fit <- stats::loess(cif ~ day, data = d, span = span, degree = 2, surface = "direct")
        yhat <- stats::predict(fit, newdata = dplyr::tibble(day = d$day))
      } else {
        yhat <- d$cif
      }
      # clamp to [0,1] and enforce non-decreasing
      yhat <- pmin(pmax(yhat, 0), 1)
      yhat <- cummax(yhat)
      d$ci_smooth <- yhat
      d
    }) %>%
    dplyr::ungroup()
}

pooled30_sm <- smooth_pooled_cif(pooled_cif, span = 0.25, day_max = 30)

# color map: Low=blue, Mid=green, High=red (as requested)
# ----- nicer palette: High=red, Mid=green, Low=blue -----
pal_bins_pub <- c("High"="#D62728", "Mid"="#2CA02C", "Low"="#1F77B4")  # readable CMYK-ish

# ----- publication-ready plotting helper -----
plot_outcome_by_pollutant <- function(df_sm, df_raw, pol, outcome, day_max = 30, y_max = NULL) {
  dd_raw <- df_raw %>% dplyr::filter(pollutant == pol, cause_label == outcome, day <= day_max)
  dd     <- df_sm  %>% dplyr::filter(pollutant == pol, cause_label == outcome, day <= day_max)
  
  if (is.null(y_max)) {
    y_max <- dd %>% dplyr::summarize(mx = max(ci_smooth, na.rm = TRUE)) %>% dplyr::pull(mx)
    y_max <- ifelse(is.finite(y_max), ceiling(y_max*20)/20, 0.2)
  }
  
  ggplot() +
    # raw steps (faint, thin) for transparency
    geom_step(data = dd_raw, aes(day, cif, color = exposure_bin, group = exposure_bin),
              linewidth = 0.6, alpha = 0.25) +
    # smoothed monotone line (emphasized)
    geom_line(data = dd, aes(day, ci_smooth, color = exposure_bin, group = exposure_bin),
              linewidth = 1.4) +
    scale_color_manual(values = pal_bins_pub, breaks = c("High","Mid","Low"), name = NULL) +
    scale_x_continuous(limits = c(0, day_max), breaks = seq(0, day_max, by = 5), expand = c(0,0)) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title = glue::glue("{pol} — {outcome}"),
         x = "Days since ICU admission",
         y = "Cumulative incidence",
         caption = "Pooled across sites; LOESS-smoothed (monotone)") +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 14, hjust = 0),
      panel.border       = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
      legend.position    = "top",
      legend.direction   = "horizontal",
      legend.box.background = element_rect(colour = "#6e7781", fill = "white", linewidth = 0.6),
      legend.key         = element_rect(fill = "white", colour = NA),
      axis.title.x       = element_text(margin = margin(t = 6)),
      axis.title.y       = element_text(margin = margin(r = 6)),
      plot.margin        = margin(4, 8, 6, 6)
    )
}

# make and save (30-day window, smoothed)
outs <- c("Death", "Successful extubation", "Persistent respiratory failure")
for (pol in c("NO2","PM2.5")) {
  for (oc in outs) {
    p <- plot_outcome_by_pollutant(pooled30_sm, pooled_cif, pol = pol, outcome = oc, day_max = 30)
    stub <- glue::glue("{out_dir}/CIF_pooled30_smooth__{pol}__{gsub(' ', '', oc)}")
    ggsave(paste0(stub, ".png"), p, width = 6.8, height = 4.6, dpi = 600)
    ggsave(paste0(stub, ".pdf"), p, width = 6.8, height = 4.6, device = cairo_pdf)
  }
}


# ---------------------------------------------------------
# 0) INPUTS
# ---------------------------------------------------------
exposome_dir <- file.path(repo, "exposome")

# helper to pick the mean column robustly
pick_mean_col <- function(df, hints) {
  nm <- names(df)
  cand <- setdiff(nm, c("GEOID","geoid","YEAR","year"))
  # prioritize names containing the hint (e.g., "no2" or "pm25")
  cand <- c(grep(hints, cand, ignore.case = TRUE, value = TRUE), cand)
  # choose the first numeric column from candidates
  for (v in unique(cand)) if (is.numeric(df[[v]])) return(v)
  stop("Could not find a numeric mean column in ", deparse(substitute(df)))
}

# ---------------------------------------------------------
# 1) Read national exposure distributions & compute tertiles
# ---------------------------------------------------------
no2_us  <- readr::read_csv(file.path(exposome_dir, "no2_county_year.csv"), show_col_types = FALSE)
pm25_us <- readr::read_csv(file.path(exposome_dir, "pm25_county_year.csv"), show_col_types = FALSE)

no2_col  <- pick_mean_col(no2_us,  "no2")
pm25_col <- pick_mean_col(pm25_us, "pm25")

no2_vals  <- no2_us[[no2_col]]
pm25_vals <- pm25_us[[pm25_col]]

no2_cut  <- quantile(no2_vals,  probs = c(1/3, 2/3), na.rm = TRUE)
pm25_cut <- quantile(pm25_vals, probs = c(1/3, 2/3), na.rm = TRUE)

message("US tertiles — NO2:  ", paste(round(no2_cut, 3), collapse = " | "))
message("US tertiles — PM2.5:", paste(round(pm25_cut,3), collapse = " | "))

# ---------------------------------------------------------
# 2) Map each site's exposure group (lower/upper) to US Low/Mid/High
#     - use midpoint of [lower, upper]
# ---------------------------------------------------------


# ---------- (re)read ALL exposure-bins files, one-by-one ----------
bin_regex <- "site_exposure_bins__.*_cumulative(?: copy \\d+)?\\.csv$"
bin_files <- fs::dir_ls(cif_dir, regexp = bin_regex, type = "file", recurse = TRUE)

read_one_bins <- function(fp) {
  df_raw <- readr::read_csv(fp, show_col_types = FALSE)
  names(df_raw) <- names(df_raw) |>
    tolower() |> gsub("[^a-z0-9]+","_", x = _) |> gsub("^_|_$","", x = _)
  
  df <- df_raw %>%
    dplyr::rename(
      site_id       = dplyr::any_of("site_id"),
      site_name     = dplyr::any_of("site_name"),
      exposure_var  = dplyr::any_of("exposure_var"),
      exposure_lab  = dplyr::any_of("exposure_label"),
      unit          = dplyr::any_of("unit"),
      exposure_grp  = dplyr::any_of("exposure_group"),
      lower         = dplyr::any_of("lower"),
      upper         = dplyr::any_of("upper"),
      n_in_group    = dplyr::any_of("n_in_group")
    )
  
  # enforce single site_name from within the file
  sn <- df$site_name |> unique() |> na.omit()
  sn <- if (length(sn)>=1L) sn[1] else NA_character_
  df <- df %>% dplyr::mutate(site_name = sn)
  
  # add pollutant from filename
  pol <- if (grepl("no2", tolower(fp))) "NO2" else "PM2.5"
  
  df %>%
    dplyr::mutate(
      pollutant = pol,
      lower_num = suppressWarnings(as.numeric(lower)),
      upper_num = suppressWarnings(as.numeric(upper)),
      mid       = (lower_num + upper_num)/2,
      group_clean = stringr::str_match(tolower(exposure_grp), "^(low|mid|high)")[,2],
      group_clean = dplyr::recode(group_clean, low="Low", mid="Mid", high="High", .default = NA_character_),
      source_file = fs::path_file(fp)
    )
}

bins_all <- purrr::map_dfr(bin_files, read_one_bins)


bins_us <- bins_all %>%
  mutate(
    lower_num = as.numeric(lower),
    upper_num = as.numeric(upper),
    mid       = (lower_num + upper_num)/2,
    us_bin = case_when(
      pollutant == "NO2"  ~ as.character(cut(mid,  c(-Inf, no2_cut[1],  no2_cut[2],  Inf),
                                             labels = c("Low","Mid","High"), right = TRUE)),
      pollutant == "PM2.5"~ as.character(cut(mid,  c(-Inf, pm25_cut[1], pm25_cut[2], Inf),
                                             labels = c("Low","Mid","High"), right = TRUE)),
      TRUE ~ NA_character_
    ),
    us_bin = factor(us_bin, levels = c("Low","Mid","High"))
  ) %>%
  select(site_id, site_name, pollutant, exposure_var, "exposure_group" = exposure_grp, us_bin, mid, lower_num, upper_num) %>%
  distinct()

# ---------------------------------------------------------
# 3) Attach US-based bins to the CIF rows
#     - prefer US-based bin; if missing, fall back to original bin
# ---------------------------------------------------------
# ---- attach US-based bins to CIF rows ----
cif_us <- cif_all %>%
  dplyr::left_join(
    bins_us,
    by = c("site_id","site_name","pollutant","exposure_var","exposure_group")
  ) %>%
  dplyr::mutate(
    exposure_bin_us = dplyr::coalesce(as.character(us_bin), as.character(exposure_bin)),
    exposure_bin_us = factor(exposure_bin_us, levels = c("Low","Mid","High"))
  ) %>%
  dplyr::filter(!is.na(exposure_bin_us), !is.na(cause_label))

# risk-set weighted pool
pooled_us <- cif_us %>%
  dplyr::group_by(pollutant, cause_label, exposure_bin_us, day) %>%
  dplyr::summarize(
    n_sites       = dplyr::n_distinct(site_name),
    total_at_risk = sum(risk_set, na.rm = TRUE),
    cif           = stats::weighted.mean(cif, w = pmax(risk_set, 0), na.rm = TRUE),
    .groups = "drop"
  )

smooth_pooled_us <- function(df, span = 0.25, day_max = 30) {
  df %>%
    filter(day <= day_max) %>%
    arrange(pollutant, cause_label, exposure_bin_us, day) %>%
    group_by(pollutant, cause_label, exposure_bin_us) %>%
    group_modify(~{
      d <- .x
      if (nrow(d) >= 5) {
        fit  <- stats::loess(cif ~ day, data = d, span = span, degree = 2, surface = "direct")
        yhat <- stats::predict(fit, newdata = tibble(day = d$day))
      } else yhat <- d$cif
      yhat <- pmin(pmax(yhat, 0), 1)
      yhat <- cummax(yhat)
      d$ci_smooth <- yhat
      d
    }) %>%
    ungroup()
}

# (reuse your smoothing + plotting functions, but point them at exposure_bin_us)
pooled30_us_sm <- smooth_pooled_us(pooled_us, span = 0.25, day_max = 30)

pal_bins_pub <- c("High"="#D62728", "Mid"="#2CA02C", "Low"="#1F77B4")

plot_us_bins <- function(df_sm, df_raw, pol, outcome, day_max = 30, y_max = NULL) {
  dd_raw <- df_raw %>% dplyr::filter(pollutant == pol, cause_label == outcome, day <= day_max)
  dd     <- df_sm  %>% dplyr::filter(pollutant == pol, cause_label == outcome, day <= day_max)
  
  thr_txt <- if (pol == "NO2")
    glue::glue("US tertiles (ppb): ≤{round(no2_cut[1],2)}, {round(no2_cut[1],2)}–{round(no2_cut[2],2)}, ≥{round(no2_cut[2],2)}")
  else
    glue::glue("US tertiles (µg/m³): ≤{round(pm25_cut[1],2)}, {round(pm25_cut[1],2)}–{round(pm25_cut[2],2)}, ≥{round(pm25_cut[2],2)}")
  
  ggplot() +
    geom_step(data = dd_raw, aes(day, cif, color = exposure_bin_us, group = exposure_bin_us),
              linewidth = 0.6, alpha = 0.25) +
    geom_line(data = dd, aes(day, ci_smooth, color = exposure_bin_us, group = exposure_bin_us),
              linewidth = 1.4) +
    scale_color_manual(values = pal_bins_pub, breaks = c("High","Mid","Low"), name = NULL) +
    scale_x_continuous(limits = c(0, day_max), breaks = seq(0, day_max, by = 5), expand = c(0,0)) +
    labs(title = glue::glue("{pol} — {outcome}"), x = "Days since ICU admission", y = "Cumulative incidence",
         caption = glue::glue("Pooled across sites (risk-set weighted). Bins by US distribution — {thr_txt}")) +
    theme_classic(base_size = 12) +
    theme(panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
          legend.position = "top", legend.box.background = element_rect(colour = "#6e7781", fill = "white", linewidth = 0.6))
}

# save figures
outs <- c("Death", "Successful extubation", "Persistent respiratory failure")
for (pol in c("NO2","PM2.5")) {
  for (oc in outs) {
    p <- plot_us_bins(pooled30_us_sm, pooled_us, pol = pol, outcome = oc, day_max = 30)
    stub <- glue::glue("{out_dir}/CIF_USbins30_smooth__{pol}__{gsub(' ', '', oc)}")
    ggsave(paste0(stub, ".png"), p, width = 6.8, height = 4.6, dpi = 600)
    ggsave(paste0(stub, ".pdf"), p, width = 6.8, height = 4.6, device = cairo_pdf)
  }
}

# -------- legend labels with US thresholds --------
# legend labels with US thresholds (unchanged)
# ---------- legend labels with thresholds ----------
lab_strings <- function(pol) {
  if (pol == "NO2") {
    c(High = paste0("High (≥", round(no2_cut[2], 2), " ppb)"),
      Mid  = paste0("Mid (", round(no2_cut[1], 2), "–", round(no2_cut[2], 2), " ppb)"),
      Low  = paste0("Low (≤", round(no2_cut[1], 2), " ppb)"))
  } else {
    c(High = paste0("High (≥", round(pm25_cut[2], 2), " µg/m³)"),
      Mid  = paste0("Mid (", round(pm25_cut[1], 2), "–", round(pm25_cut[2], 2), " µg/m³)"),
      Low  = paste0("Low (≤", round(pm25_cut[1], 2), " µg/m³)"))
  }
}

pal_bins <- c("High"="#D62728","Mid"="#2CA02C","Low"="#1F77B4")
outs     <- c("Death","Persistent respiratory failure","Successful extubation")

# outcome-wise y-limits (same across columns)
ymax_by_outcome <- pooled30_us_sm |>
  dplyr::group_by(cause_label) |>
  dplyr::summarise(ymax = ceiling(max(ci_smooth, na.rm = TRUE)*20)/20, .groups="drop") |>
  tibble::deframe()

# base panel builder (no legend, no axis titles by default)
panel_base <- function(pol, outcome, show_x = FALSE, show_y = FALSE, day_max = 30) {
  labs_vec <- lab_strings(pol)
  dd <- pooled30_us_sm |>
    dplyr::filter(pollutant == pol, cause_label == outcome, day <= day_max)
  
  ggplot(dd, aes(day, ci_smooth, color = exposure_bin_us, group = exposure_bin_us)) +
    geom_line(linewidth = 1.3) +
    scale_color_manual(values = pal_bins,
                       breaks = c("High","Mid","Low"),
                       labels = labs_vec[c("High","Mid","Low")],
                       name = NULL) +
    scale_x_continuous(
      limits = c(0, day_max),
      breaks = seq(0, day_max, 5),
      expand = expansion(mult = c(0, 0.02))   # 2% right padding so "30" shows
    ) +
    coord_cartesian(ylim = c(0, ymax_by_outcome[[outcome]])) +
    labs(
      x = if (show_x) "Days since ICU admission" else NULL,
      y = if (show_y) "Cumulative incidence"     else NULL
    ) +
    theme_classic(base_size = 11.5) +
    theme(
      panel.border  = element_rect(colour="grey35", fill=NA, linewidth=0.6),
      axis.title.x  = element_text(margin = margin(t = 4)),
      axis.title.y  = element_text(margin = margin(r = 4)),
      plot.margin   = margin(2, 10, 2, 2),     # +right margin
      legend.position = "none"
    )
}

# build 3x2 grid of panels (only left column has y, only bottom row has x)
# left column = NO2, right column = PM2.5
panels_left  <- mapply(\(oc, i) panel_base("NO2",  oc, show_x = (i==3), show_y = TRUE),
                       outs, seq_along(outs), SIMPLIFY = FALSE)
panels_right <- mapply(\(oc, i) panel_base("PM2.5", oc, show_x = (i==3), show_y = FALSE),
                       outs, seq_along(outs), SIMPLIFY = FALSE)

# column headers (tight, no extra rows)
# ----- compact column headers (unchanged) -----
col_header <- function(txt) {
  ggplot() + theme_void() +
    labs(title = txt) +
    theme(plot.title = element_text(face="bold", size=12, hjust=0.5),
          plot.margin = margin(0,0,2,0))
}

hdr_no2  <- ggplot() + theme_void() +
  labs(title = expression(NO[2])) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.margin = margin(0,0,2,0))

hdr_pm25 <- ggplot() + theme_void() +
  labs(title = expression(PM[2.5])) +  # <- dot is subscript, not baseline
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.margin = margin(0,0,2,0))


# ----- make one legend per pollutant (single row), then put them side-by-side -----
# Legend that always shows all three levels (no sampling from data)
legend_for <- function(pol) {
  labs_vec <- lab_strings(pol)
  legend_df <- tibble::tibble(
    exposure_bin_us = factor(c("High","Mid","Low"), levels = c("High","Mid","Low")),
    day = c(1,2,3), ci_smooth = c(1,2,3)
  )
  p_legend <- ggplot(legend_df, aes(day, ci_smooth, color = exposure_bin_us)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = pal_bins,
                       breaks = c("High","Mid","Low"),
                       labels = labs_vec[c("High","Mid","Low")],
                       name = NULL) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 10),
          legend.key.width = unit(16, "pt"))
  cowplot::get_legend(p_legend + guides(color = guide_legend(nrow = 1)))
}

# -------------------- rebuild column stacks (no headers inside) --------------------
# stacks of panels (no headers, no legends in panels)
left_col_stack  <- cowplot::plot_grid(panels_left[[1]], panels_left[[2]], panels_left[[3]],
                                      ncol = 1, rel_heights = c(1,1,1), align = "v")
right_col_stack <- cowplot::plot_grid(panels_right[[1]], panels_right[[2]], panels_right[[3]],
                                      ncol = 1, rel_heights = c(1,1,1), align = "v")

# --- tiny vertical row-label strip, tight to y-axes ---
# 1) tiny vertical label (same as before)
vrow_lab <- function(txt) {
  cowplot::ggdraw() +
    cowplot::draw_label(txt, angle = 90, fontface = "bold",
                        hjust = 0.5, vjust = 0.5, size = 11) +
    theme_void()
}

# 2) build each row: [label | (NO2 panel , PM2.5 panel)]
row1 <- cowplot::plot_grid(
  vrow_lab("Death"),
  cowplot::plot_grid(panels_left[[1]],  panels_right[[1]],  ncol = 2, rel_widths = c(1,1), align = "h"),
  ncol = 2, rel_widths = c(0.055, 1), align = "h"
)
row2 <- cowplot::plot_grid(
  vrow_lab("Persistent respiratory failure"),
  cowplot::plot_grid(panels_left[[2]],  panels_right[[2]],  ncol = 2, rel_widths = c(1,1), align = "h"),
  ncol = 2, rel_widths = c(0.055, 1), align = "h"
)
row3 <- cowplot::plot_grid(
  vrow_lab("Successful extubation"),
  cowplot::plot_grid(panels_left[[3]],  panels_right[[3]],  ncol = 2, rel_widths = c(1,1), align = "h"),
  ncol = 2, rel_widths = c(0.055, 1), align = "h"
)

# stack the three outcome rows (labels are now centered per row)
panel_grid_labeled <- cowplot::plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1,1,1), align = "v")

# 3) column headers (leave a spacer in the label column)
hdr_no2  <- ggplot() + theme_void() + labs(title = expression(NO[2]))    +
  theme(plot.title = element_text(face="bold", size=12, hjust=0.5))
hdr_pm25 <- ggplot() + theme_void() + labs(title = expression(PM[2.5])) +
  theme(plot.title = element_text(face="bold", size=12, hjust=0.5))

header_row <- cowplot::plot_grid(
  ggplot() + theme_void(),  # spacer where label strip would be
  hdr_no2, hdr_pm25,
  ncol = 3, rel_widths = c(0.055, 1, 1)
)

# 4) per-column legends (again with a left spacer)
leg_no2  <- legend_for("NO2")
leg_pm25 <- legend_for("PM2.5")
legend_row <- cowplot::plot_grid(
  ggplot() + theme_void(),          # spacer for label strip
  cowplot::ggdraw(leg_no2),
  cowplot::ggdraw(leg_pm25),
  ncol = 3, rel_widths = c(0.055, 1, 1)
)

# 5) title + final assembly
top_title <- ggplot() + theme_void() +
  labs(title   = "Cumulative incidence functions by US exposure bins",
       subtitle= "Pooled across sites (risk-set weighted); curves shown for 0–30 days") +
  theme(plot.title   = element_text(face="bold", size=14, hjust=0),
        plot.subtitle= element_text(size=11, hjust=0, margin = margin(t=1, b=4)))

final_fig <- cowplot::plot_grid(
  top_title,
  header_row,
  panel_grid_labeled,
  legend_row,
  ncol = 1,
  rel_heights = c(0.10, 0.06, 1, 0.12),
  align = "v"
)

ggsave(file.path(out_dir, "CIF_USbins30_smooth__six_panel_CENTERED-ROW-LABELS.png"),
       final_fig, width = 10.6, height = 8.9, dpi = 600)
ggsave(file.path(out_dir, "CIF_USbins30_smooth__six_panel_CENTERED-ROW-LABELS.pdf"),
       final_fig, width = 10.6, height = 8.9, device = cairo_pdf)




# --- vertical row label (same as before) ---
vrow_lab <- function(txt) {
  cowplot::ggdraw() +
    cowplot::draw_label(txt, angle = 90, fontface = "bold",
                        hjust = 0.5, vjust = 0.5, size = 11) +
    theme_void()
}

# --- build only the two outcome rows ---
row1 <- cowplot::plot_grid(
  vrow_lab("Death"),
  cowplot::plot_grid(panels_left[[1]], panels_right[[1]],
                     ncol = 2, rel_widths = c(1,1), align = "h"),
  ncol = 2, rel_widths = c(0.055, 1), align = "h"
)

row2 <- cowplot::plot_grid(
  vrow_lab("Successful extubation"),
  cowplot::plot_grid(panels_left[[3]], panels_right[[3]],  # skip the middle one
                     ncol = 2, rel_widths = c(1,1), align = "h"),
  ncol = 2, rel_widths = c(0.055, 1), align = "h"
)

# --- stack the two rows ---
panel_grid_labeled <- cowplot::plot_grid(
  row1, row2,
  ncol = 1, rel_heights = c(1,1),
  align = "v"
)

# --- headers over columns (leave left spacer for label strip) ---
hdr_no2  <- ggplot() + theme_void() +
  labs(title = expression(NO[2])) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
hdr_pm25 <- ggplot() + theme_void() +
  labs(title = expression(PM[2.5])) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

header_row <- cowplot::plot_grid(
  ggplot() + theme_void(),  # spacer for row labels
  hdr_no2, hdr_pm25,
  ncol = 3, rel_widths = c(0.055, 1, 1)
)

# --- legends per column ---
leg_no2  <- legend_for("NO2")
leg_pm25 <- legend_for("PM2.5")

legend_row <- cowplot::plot_grid(
  ggplot() + theme_void(),          # spacer for label strip
  cowplot::ggdraw(leg_no2),
  cowplot::ggdraw(leg_pm25),
  ncol = 3, rel_widths = c(0.055, 1, 1)
)

# --- title + subtitle ---
top_title <- ggplot() + theme_void() +
  labs(
    title   = "Cumulative incidence functions by US exposure bins",
    subtitle= "Pooled across sites (risk-set weighted); curves shown for 0–30 days"
  ) +
  theme(
    plot.title   = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle= element_text(size = 11, hjust = 0, margin = margin(t = 1, b = 4))
  )

# --- final 2x2 layout ---
final_fig <- cowplot::plot_grid(
  top_title,
  header_row,
  panel_grid_labeled,
  legend_row,
  ncol = 1,
  rel_heights = c(0.10, 0.06, 1, 0.14),
  align = "v"
)

ggsave(file.path(out_dir, "CIF_USbins30_smooth__FOUR-PANEL.png"),
       final_fig, width = 11.6, height = 7.6, dpi = 600)
ggsave(file.path(out_dir, "CIF_USbins30_smooth__FOUR-PANEL.pdf"),
       final_fig, width = 11.6, height = 7.6, device = cairo_pdf)



# --- combined CIF plot per pollutant: Death (CIF) vs Extubation (1 − CIF) ---
# expects columns: pollutant, cause_label, exposure_bin_us, day, ci_smooth
# ---- publication styling helpers ----


# if you already have: brks <- get_breaks(no2_vals, pm25_vals, scheme = "tertiles" or similar)
get_or_make_cuts <- function(no2_vals, pm25_vals, brks = NULL) {
  if (!is.null(brks)) return(list(NO2 = brks$NO2, PM25 = brks$PM25))
  list(
    NO2  = quantile(no2_vals,  probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    PM25 = quantile(pm25_vals, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  )
}

load_us_vals <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE) #%>% filter(year %in% 2018:2024)
  # grab the first numeric column that isn't GEOID or year
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, c("year"))
  # GEOID is often character; if it’s numeric in your file, exclude it too:
  num_cols <- setdiff(num_cols, c("GEOID","geoid"))
  df[[num_cols[1]]]
}

no2_vals  <- load_us_vals(file.path(exposome_dir, "no2_county_year.csv"))
pm25_vals <- load_us_vals(file.path(exposome_dir, "pm25_county_year.csv"))



# tertiles (change probs if using another scheme)
cuts <- list(
  NO2  = quantile(no2_vals,  probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  PM25 = quantile(pm25_vals, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
)

legend_labels_from_cuts <- function(cuts, unit) {
  c(
    "Low"  = sprintf("Low (\u2264%.2f %s)", cuts[2], unit),
    "Mid"  = sprintf("Mid (%.2f–%.2f %s)", cuts[2], cuts[3], unit),
    "High" = sprintf("High (\u2265%.2f %s)", cuts[3], unit)
  )
}

labs_NO2  <- legend_labels_from_cuts(cuts$NO2,  "ppb")
labs_PM25 <- legend_labels_from_cuts(cuts$PM25, "\u00B5g/m\u00B3")


pal_bins_pub <- c("High"="#D62728", "Mid"="#2CA02C", "Low"="#1F77B4")  # keep your spec
ln_death <- 1.4   # line width
ln_ext   <- 1.4

prep_combined <- function(df, pol, day_max = 30) {
  dd <- df %>%
    dplyr::filter(pollutant == pol,
                  cause_label %in% c("Death","Successful extubation"),
                  day <= day_max)
  
  death <- dd %>%
    dplyr::filter(cause_label == "Death") %>%
    dplyr::mutate(outcome = "Death (CIF)")
  
  ext1m <- dd %>%
    dplyr::filter(cause_label == "Successful extubation") %>%
    dplyr::mutate(ci_smooth = 1 - ci_smooth,
                  outcome   = "Extubation (1 − CIF)")
  
  dplyr::bind_rows(death, ext1m)
}

make_combined_plot_pub <- function(df_pol, pol, labs_color) {
  y_top <- ceiling(max(df_pol$ci_smooth, na.rm = TRUE) * 20) / 20
  
  ggplot(df_pol, aes(
    x = day, y = ci_smooth,
    color = exposure_bin_us,
    linetype = outcome,
    group = interaction(exposure_bin_us, outcome)
  )) +
    geom_line(linewidth = 1.4) +
    scale_color_manual(
      values = c("High"="#D62728","Mid"="#2CA02C","Low"="#1F77B4"),
      breaks = c("Low","Mid","High"),                 # enforce order
      labels = labs_color[c("Low","Mid","High")],     # add cutpoints
      name   = NULL
    ) +
    scale_linetype_manual(
      values = c("Death (CIF)" = "solid", "Extubation (1 − CIF)" = "11"),
      breaks = c("Death (CIF)", "Extubation (1 − CIF)"),
      labels = c("Death (CIF)", "Extubation (1 − CIF)"),
      name   = NULL
    ) +
    scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, y_top)) +
    labs(
      subtitle = if (pol == "NO2") expression(NO[2]) else expression(PM[2.5]),
      x = "Days since ICU admission",
      y = expression(paste("Cumulative Incidence"))
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      panel.border   = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      axis.title.x   = element_text(margin = margin(t = 8)),
      axis.title.y   = element_text(margin = margin(r = 10)),
      plot.subtitle  = element_text(face = "bold", hjust = 0.5, margin = margin(b = 2)),
      legend.position = "bottom",
      legend.box      = "vertical",           # <-- stacked like PM2.5
      legend.text     = element_text(size = 10.5),
      legend.key.width = unit(28, "pt"),   # was 16 pt
      legend.key.height= unit(8,  "pt"),
      legend.spacing.y= unit(2, "pt"),
      plot.margin     = margin(4, 6, 4, 6)
    ) +
    guides(
      color    = guide_legend(order = 1, ncol = 3, byrow = TRUE,
                              override.aes = list(linetype = "solid", linewidth = 1.4)),
      linetype = guide_legend(order = 2, ncol = 3,                      # keep your layout
                              override.aes = list(color = "black",
                                                  linewidth = 1.6,
                                                  linetype = c("solid", "11")))
    )
}

# ---- build polished two-panel figure (collect legends) ----
df_no2   <- prep_combined(pooled30_us_sm, "NO2")
df_pm25  <- prep_combined(pooled30_us_sm, "PM2.5")

p_no2_pub  <- make_combined_plot_pub(df_no2,  "NO2",  labs_color = labs_NO2)
p_pm25_pub <- make_combined_plot_pub(df_pm25, "PM2.5", labs_color = labs_PM25)


final_two_panel <- p_no2_pub | p_pm25_pub +
  plot_annotation(
    title    = "Competing risks on a common scale",
    subtitle = "Solid: CIF(Death)   ·   Dashed: 1 − CIF(Extubation)   ·   Color: US exposure bin (cutpoints shown)",
    theme = theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle= element_text(size = 11, margin = margin(t = 2, b = 6))
    )
  ) &
  theme(legend.position = "bottom")    # same stacked block under each panel

ggsave(file.path(out_dir, "CIF_competing__NO2_PM25_overlaid__pub_MATCHED-LEGENDS_studyprd.png"),
       final_two_panel, width = 11.8, height = 5.4, dpi = 600)
ggsave(file.path(out_dir, "CIF_competing__NO2_PM25_overlaid__pub_MATCHED-LEGENDS_studyprd.pdf"),
       final_two_panel, width = 11.8, height = 5.4, device = cairo_pdf)

#### low/high cuts

# --- MEDIAN-BASED CUTS (add after you load no2_us / pm25_us and pick_mean_col) ---
no2_vals  <- no2_us[[no2_col]]
pm25_vals <- pm25_us[[pm25_col]]

no2_med  <- median(no2_vals,  na.rm = TRUE)
pm25_med <- median(pm25_vals, na.rm = TRUE)

message("US medians — NO2: ", round(no2_med, 3), " ppb;  PM2.5: ", round(pm25_med, 3), " µg/m³")

# --- Reuse/read the same exposure-bins files you already collected in `bins_all` ---
# (If `bins_all` exists from your tertile code, skip re-reading and just run this.)

bins_us_med <- bins_all %>%
  mutate(
    lower_num = as.numeric(lower),
    upper_num = as.numeric(upper),
    mid       = (lower_num + upper_num)/2,
    us_bin_med = case_when(
      pollutant == "NO2"   ~ if_else(mid <= no2_med, "Low", "High"),
      pollutant == "PM2.5" ~ if_else(mid <= pm25_med, "Low", "High"),
      TRUE ~ NA_character_
    ),
    us_bin_med = factor(us_bin_med, levels = c("Low","High"))
  ) %>%
  select(
    site_id, site_name, pollutant, exposure_var,
    exposure_group = exposure_grp, mid, us_bin_med
  ) %>%
  distinct()

# Attach Low/High (median) bins to CIF rows
cif_median <- cif_all %>%
  left_join(
    bins_us_med,
    by = c("site_id","site_name","pollutant","exposure_var","exposure_group")
  ) %>%
  mutate(
    # If a row has no median bin (rare), fall back by computing from `mid` if available
    exposure_bin_us_med = coalesce(
      as.character(us_bin_med),
      case_when(
        pollutant == "NO2"   & !is.na(mid) ~ if_else(mid <= no2_med, "Low", "High"),
        pollutant == "PM2.5" & !is.na(mid) ~ if_else(mid <= pm25_med, "Low", "High"),
        TRUE ~ NA_character_
      )
    ),
    exposure_bin_us_med = factor(exposure_bin_us_med, levels = c("Low","High"))
  ) %>%
  filter(!is.na(exposure_bin_us_med), !is.na(cause_label))

# Risk-set–weighted pooling by median bins
pooled_med <- cif_median %>%
  group_by(pollutant, cause_label, exposure_bin_us_med, day) %>%
  dplyr::summarize(
    n_sites        = n_distinct(site_name),
    total_at_risk  = sum(risk_set, na.rm = TRUE),
    cif            = weighted.mean(cif, w = pmax(risk_set, 0), na.rm = TRUE),
    .groups = "drop"
  )

smooth_pooled_med <- function(df, span = 0.25, day_max = 30) {
  df %>%
    filter(day <= day_max) %>%
    arrange(pollutant, cause_label, exposure_bin_us_med, day) %>%
    group_by(pollutant, cause_label, exposure_bin_us_med) %>%
    group_modify(~{
      d <- .x
      if (nrow(d) >= 5) {
        fit  <- stats::loess(cif ~ day, data = d, span = span, degree = 2, surface = "direct")
        yhat <- stats::predict(fit, newdata = tibble(day = d$day))
      } else yhat <- d$cif
      yhat <- pmin(pmax(yhat, 0), 1)
      yhat <- cummax(yhat)
      d$ci_smooth <- yhat
      d
    }) %>%
    ungroup()
}
pooled30_med_sm <- smooth_pooled_med(pooled_med, span = 0.25, day_max = 30)

# Two-color palette for median split
pal_bins_med <- c("Low"="#1F77B4", "High"="#D62728")

lab_strings_median <- function(pol) {
  if (pol == "NO2") {
    c("Low"  = paste0("≤ ", round(no2_med, 2), " ppb"),
      "High" = paste0("> ",  round(no2_med, 2), " ppb"))
  } else {
    c("Low"  = paste0("≤ ", round(pm25_med, 2), " µg/m³"),
      "High" = paste0("> ",  round(pm25_med, 2), " µg/m³"))
  }
}

plot_median_bins <- function(df_sm, df_raw, pol, outcome, day_max = 30, y_max = NULL) {
  dd_raw <- df_raw %>% filter(pollutant == pol, cause_label == outcome, day <= day_max)
  dd     <- df_sm  %>% filter(pollutant == pol, cause_label == outcome, day <= day_max)
  labs_vec <- lab_strings_median(pol)
  
  if (is.null(y_max)) {
    y_max <- dd %>% dplyr::summarize(mx = max(ci_smooth, na.rm = TRUE)) %>% pull(mx)
    y_max <- ifelse(is.finite(y_max), ceiling(y_max*20)/20, 0.2)
  }
  
  ggplot() +
    geom_step(data = dd_raw,
              aes(day, cif, color = exposure_bin_us_med, group = exposure_bin_us_med),
              linewidth = 0.6, alpha = 0.25) +
    geom_line(data = dd,
              aes(day, ci_smooth, color = exposure_bin_us_med, group = exposure_bin_us_med),
              linewidth = 1.4) +
    scale_color_manual(values = pal_bins_med,
                       breaks = c("Low","High"),
                       labels = labs_vec[c("Low","High")],
                       name = NULL) +
    scale_x_continuous(limits = c(0, day_max), breaks = seq(0, day_max, by = 5), expand = c(0,0)) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title = paste0(pol, " — ", outcome),
         x = "Days since ICU admission",
         y = "Cumulative incidence",
         caption = "Pooled across sites (risk-set weighted). Bins: at national median") +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
      legend.position = "top",
      legend.box.background = element_rect(colour = "#6e7781", fill = "white", linewidth = 0.6)
    )
}

# --- SAVE median Low/High CIF figures ---
outs <- c("Death", "Successful extubation", "Persistent respiratory failure")
out_dir <- file.path(cif_dir, "fig_out"); dir.create(out_dir, showWarnings = FALSE)

for (pol in c("NO2","PM2.5")) {
  for (oc in outs) {
    p <- plot_median_bins(pooled30_med_sm, pooled_med, pol = pol, outcome = oc, day_max = 30)
    stub <- glue::glue("{out_dir}/CIF_USMEDIAN30_smooth__{pol}__{gsub(' ', '', oc)}")
    ggsave(paste0(stub, ".png"), p, width = 6.8, height = 4.6, dpi = 600)
    ggsave(paste0(stub, ".pdf"), p, width = 6.8, height = 4.6, device = cairo_pdf)
  }
}

# Prepare combined Death vs 1 − Extubation for median bins
prep_combined_med <- function(df, pol, day_max = 30) {
  dd <- df %>%
    filter(pollutant == pol,
           cause_label %in% c("Death","Successful extubation"),
           day <= day_max)
  
  death <- dd %>%
    filter(cause_label == "Death") %>%
    mutate(outcome = "Death (CIF)")
  
  ext1m <- dd %>%
    filter(cause_label == "Successful extubation") %>%
    mutate(ci_smooth = 1 - ci_smooth,
           outcome   = "Extubation (1 − CIF)")
  
  bind_rows(death, ext1m)
}

make_combined_plot_med <- function(df_pol, pol) {
  labs_color <- lab_strings_median(pol)
  y_top <- ceiling(max(df_pol$ci_smooth, na.rm = TRUE) * 20) / 20
  
  ggplot(df_pol, aes(
    x = day, y = ci_smooth,
    color = exposure_bin_us_med,
    linetype = outcome,
    group = interaction(exposure_bin_us_med, outcome)
  )) +
    geom_line(linewidth = 1.4) +
    scale_color_manual(
      values = pal_bins_med,
      breaks = c("Low","High"),
      labels = labs_color[c("Low","High")],
      name   = NULL
    ) +
    scale_linetype_manual(
      values = c("Death (CIF)" = "solid", "Extubation (1 − CIF)" = "11"),
      breaks = c("Death (CIF)", "Extubation (1 − CIF)"),
      labels = c("Death (CIF)", "Extubation (1 − CIF)"),
      name   = NULL
    ) +
    scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, y_top)) +
    labs(
      subtitle = if (pol == "NO2") expression(NO[2]) else expression(PM[2.5]),
      x = "Days since ICU admission",
      y = "Cumulative Incidence"
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      panel.border   = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      plot.subtitle  = element_text(face = "bold", hjust = 0.5, margin = margin(b = 2)),
      legend.position = "bottom",
      legend.box      = "vertical",
      legend.key.width = unit(28, "pt"),
      legend.key.height= unit(8,  "pt"),
      legend.spacing.y= unit(2, "pt"),
      plot.margin     = margin(4, 6, 4, 6)
    ) +
    guides(
      color    = guide_legend(order = 1, ncol = 2, byrow = TRUE,
                              override.aes = list(linetype = "solid", linewidth = 1.4)),
      linetype = guide_legend(order = 2, ncol = 2,
                              override.aes = list(color = "black",
                                                  linewidth = 1.6,
                                                  linetype = c("solid", "11")))
    )
}

df_no2_med  <- prep_combined_med(pooled30_med_sm, "NO2")
df_pm25_med <- prep_combined_med(pooled30_med_sm, "PM2.5")

p_no2_med  <- make_combined_plot_med(df_no2_med,  "NO2")
p_pm25_med <- make_combined_plot_med(df_pm25_med, "PM2.5")

final_two_panel_med <- p_no2_med | p_pm25_med +
  patchwork::plot_annotation(
    title    = "Competing risks on a common scale (median split)",
    subtitle = "Solid: CIF(Death)   ·   Dashed: 1 − CIF(Extubation)   ·   Color: exposure ≤ vs > national median",
    theme = theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle= element_text(size = 11, margin = margin(t = 2, b = 6))
    )
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "CIF_competing__NO2_PM25_overlaid__MEDIAN.png"),
       final_two_panel_med, width = 11.8, height = 5.4, dpi = 600)
ggsave(file.path(out_dir, "CIF_competing__NO2_PM25_overlaid__MEDIAN.pdf"),
       final_two_panel_med, width = 11.8, height = 5.4, device = cairo_pdf)


#### quartiles

# --- QUARTILE CUTS (Q1, Q2/median, Q3) ---
no2_col  <- pick_mean_col(no2_us,  "no2")
pm25_col <- pick_mean_col(pm25_us, "pm25")

no2_vals  <- no2_us[[no2_col]]
pm25_vals <- pm25_us[[pm25_col]]

no2_qu  <- quantile(no2_vals,  probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
pm25_qu <- quantile(pm25_vals, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

message("US quartiles — NO2 (ppb): Q1≤", round(no2_qu[1],3), " | Q2=", round(no2_qu[2],3),
        " | Q3=", round(no2_qu[3],3))
message("US quartiles — PM2.5 (µg/m³): Q1≤", round(pm25_qu[1],3), " | Q2=", round(pm25_qu[2],3),
        " | Q3=", round(pm25_qu[3],3))

# Expect bins_all to contain per-site [lower, upper] with site/pollutant metadata (as in your code)
bins_us_q4 <- bins_all %>%
  mutate(
    lower_num = as.numeric(lower),
    upper_num = as.numeric(upper),
    mid       = (lower_num + upper_num)/2,
    us_bin_q4 = case_when(
      pollutant == "NO2"   ~ as.character(cut(mid,  c(-Inf, no2_qu[1], no2_qu[2], no2_qu[3], Inf),
                                              labels = c("Q1","Q2","Q3","Q4"), right = TRUE)),
      pollutant == "PM2.5" ~ as.character(cut(mid,  c(-Inf, pm25_qu[1], pm25_qu[2], pm25_qu[3], Inf),
                                              labels = c("Q1","Q2","Q3","Q4"), right = TRUE)),
      TRUE ~ NA_character_
    ),
    us_bin_q4 = factor(us_bin_q4, levels = c("Q1","Q2","Q3","Q4"))
  ) %>%
  select(
    site_id, site_name, pollutant, exposure_var,
    exposure_group = exposure_grp, mid, us_bin_q4
  ) %>%
  distinct()

# Attach Q1–Q4 bins to CIF rows
cif_q4 <- cif_all %>%
  left_join(
    bins_us_q4,
    by = c("site_id","site_name","pollutant","exposure_var","exposure_group")
  ) %>%
  mutate(
    exposure_bin_us_q4 = us_bin_q4
  ) %>%
  filter(!is.na(exposure_bin_us_q4), !is.na(cause_label))

pooled_q4 <- cif_q4 %>%
  group_by(pollutant, cause_label, exposure_bin_us_q4, day) %>%
  dplyr::summarize(
    n_sites        = n_distinct(site_name),
    total_at_risk  = sum(risk_set, na.rm = TRUE),
    cif            = weighted.mean(cif, w = pmax(risk_set, 0), na.rm = TRUE),
    .groups = "drop"
  )

smooth_pooled_q4 <- function(df, span = 0.25, day_max = 30) {
  df %>%
    filter(day <= day_max) %>%
    arrange(pollutant, cause_label, exposure_bin_us_q4, day) %>%
    group_by(pollutant, cause_label, exposure_bin_us_q4) %>%
    group_modify(~{
      d <- .x
      if (nrow(d) >= 5) {
        fit  <- stats::loess(cif ~ day, data = d, span = span, degree = 2, surface = "direct")
        yhat <- stats::predict(fit, newdata = tibble(day = d$day))
      } else yhat <- d$cif
      yhat <- pmin(pmax(yhat, 0), 1)
      yhat <- cummax(yhat)   # enforce non-decreasing CIF
      d$ci_smooth <- yhat
      d
    }) %>%
    ungroup()
}
pooled30_q4_sm <- smooth_pooled_q4(pooled_q4, span = 0.25, day_max = 30)

legend_labels_quartiles <- function(pol) {
  if (pol == "NO2") {
    c(
      "Q1" = paste0("Q1 (≤", round(no2_qu[1],2), " ppb)"),
      "Q2" = paste0("Q2 (",  round(no2_qu[1],2), "–", round(no2_qu[2],2), " ppb)"),
      "Q3" = paste0("Q3 (",  round(no2_qu[2],2), "–", round(no2_qu[3],2), " ppb)"),
      "Q4" = paste0("Q4 (≥", round(no2_qu[3],2), " ppb)")
    )
  } else {
    c(
      "Q1" = paste0("Q1 (≤", round(pm25_qu[1],2), " µg/m³)"),
      "Q2" = paste0("Q2 (",  round(pm25_qu[1],2), "–", round(pm25_qu[2],2), " µg/m³)"),
      "Q3" = paste0("Q3 (",  round(pm25_qu[2],2), "–", round(pm25_qu[3],2), " µg/m³)"),
      "Q4" = paste0("Q4 (≥", round(pm25_qu[3],2), " µg/m³)")
    )
  }
}

# Four-level palette (colorblind-friendly-ish; adjust if you prefer)
pal_q4 <- c("Q1"="#1F77B4", "Q2"="#2CA02C", "Q3"="#FF7F0E", "Q4"="#D62728")

plot_quartiles <- function(df_sm, df_raw, pol, outcome, day_max = 30, y_max = NULL) {
  dd_raw <- df_raw %>% filter(pollutant == pol, cause_label == outcome, day <= day_max)
  dd     <- df_sm  %>% filter(pollutant == pol, cause_label == outcome, day <= day_max)
  labs_vec <- legend_labels_quartiles(pol)
  
  if (is.null(y_max)) {
    y_max <- dd %>% dplyr::summarize(mx = max(ci_smooth, na.rm = TRUE)) %>% pull(mx)
    y_max <- ifelse(is.finite(y_max), ceiling(y_max*20)/20, 0.2)
  }
  
  ggplot() +
    geom_step(data = dd_raw,
              aes(day, cif, color = exposure_bin_us_q4, group = exposure_bin_us_q4),
              linewidth = 0.6, alpha = 0.22) +
    geom_line(data = dd,
              aes(day, ci_smooth, color = exposure_bin_us_q4, group = exposure_bin_us_q4),
              linewidth = 1.35) +
    scale_color_manual(values = pal_q4,
                       breaks = c("Q1","Q2","Q3","Q4"),
                       labels = labs_vec[c("Q1","Q2","Q3","Q4")],
                       name = NULL) +
    scale_x_continuous(limits = c(0, day_max), breaks = seq(0, day_max, by = 5), expand = c(0,0)) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title = paste0(pol, " — ", outcome),
         x = "Days since ICU admission",
         y = "Cumulative incidence",
         caption = "Pooled across sites (risk-set weighted). Bins: national quartiles") +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box.background = element_rect(colour = "#6e7781", fill = "white", linewidth = 0.6)
    )
}

outs    <- c("Death", "Successful extubation", "Persistent respiratory failure")
out_dir <- file.path(cif_dir, "fig_out"); dir.create(out_dir, showWarnings = FALSE)

for (pol in c("NO2","PM2.5")) {
  for (oc in outs) {
    p <- plot_quartiles(pooled30_q4_sm, pooled_q4, pol = pol, outcome = oc, day_max = 30)
    stub <- glue::glue("{out_dir}/CIF_USQ4_30_smooth__{pol}__{gsub(' ', '', oc)}")
    ggsave(paste0(stub, ".png"), p, width = 6.8, height = 4.6, dpi = 600)
    ggsave(paste0(stub, ".pdf"), p, width = 6.8, height = 4.6, device = cairo_pdf)
  }
}

prep_combined_q4 <- function(df, pol, day_max = 30) {
  dd <- df %>%
    filter(pollutant == pol,
           cause_label %in% c("Death","Successful extubation"),
           day <= day_max)
  
  death <- dd %>%
    filter(cause_label == "Death") %>%
    mutate(outcome = "Death (CIF)")
  
  ext1m <- dd %>%
    filter(cause_label == "Successful extubation") %>%
    mutate(ci_smooth = 1 - ci_smooth,
           outcome   = "Extubation (1 − CIF)")
  
  bind_rows(death, ext1m)
}

make_competing_q4 <- function(df_pol, pol) {
  labs_color <- legend_labels_quartiles(pol)
  y_top <- ceiling(max(df_pol$ci_smooth, na.rm = TRUE) * 20) / 20
  
  ggplot(df_pol, aes(
    x = day, y = ci_smooth,
    color = exposure_bin_us_q4,
    linetype = outcome,
    group = interaction(exposure_bin_us_q4, outcome)
  )) +
    geom_line(linewidth = 1.35) +
    scale_color_manual(
      values = pal_q4,
      breaks = c("Q1","Q2","Q3","Q4"),
      labels = labs_color[c("Q1","Q2","Q3","Q4")],
      name   = NULL
    ) +
    scale_linetype_manual(
      values = c("Death (CIF)" = "solid", "Extubation (1 − CIF)" = "11"),
      breaks = c("Death (CIF)", "Extubation (1 − CIF)"),
      labels = c("Death (CIF)", "Extubation (1 − CIF)"),
      name   = NULL
    ) +
    scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, y_top)) +
    labs(
      subtitle = if (pol == "NO2") expression(NO[2]) else expression(PM[2.5]),
      x = "Days since ICU admission",
      y = "Cumulative Incidence"
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      panel.border   = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      plot.subtitle  = element_text(face = "bold", hjust = 0.5, margin = margin(b = 2)),
      legend.position = "bottom",
      legend.box      = "vertical",
      legend.key.width = unit(26, "pt"),
      legend.key.height= unit(8,  "pt"),
      legend.spacing.y= unit(2, "pt"),
      plot.margin     = margin(4, 6, 4, 6)
    ) +
    guides(
      color = guide_legend(
        order = 1,
        ncol = 2,           # << force two columns
        byrow = TRUE,
        override.aes = list(
          linetype = "solid",
          linewidth = 1.35
        )
      ),
      linetype = guide_legend(
        order = 2,
        ncol = 2,
        override.aes = list(
          color = "black",
          linewidth = 1.6,
          linetype = c("solid", "11")
        )
      )
    )
}

df_no2_q4  <- prep_combined_q4(pooled30_q4_sm, "NO2")
df_pm25_q4 <- prep_combined_q4(pooled30_q4_sm, "PM2.5")

p_no2_q4  <- make_competing_q4(df_no2_q4,  "NO2")
p_pm25_q4 <- make_competing_q4(df_pm25_q4, "PM2.5")

final_two_panel_q4 <- p_no2_q4 | p_pm25_q4 +
  patchwork::plot_annotation(
    title    = "Competing risks on a common scale (quartiles)",
    subtitle = "Solid: CIF(Death) · Dashed: 1 − CIF(Extubation) · Color: national quartiles (Q1–Q4)",
    theme = theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle= element_text(size = 11, margin = margin(t = 2, b = 6))
    )
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "CIF_competing__NO2_PM25_overlaid__QUARTILES.png"),
       final_two_panel_q4, width = 11.8, height = 5.4, dpi = 600)
ggsave(file.path(out_dir, "CIF_competing__NO2_PM25_overlaid__QUARTILES.pdf"),
       final_two_panel_q4, width = 11.8, height = 5.4, device = cairo_pdf)




#### flow chart inclusions

# ---- Packages ----
# install.packages(c("tidyverse","DiagrammeR","DiagrammeRsvg","rsvg","scales"))  # if needed
library(tidyverse)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(scales)

# ---- 1) Data: unified flow (Option B, multi-site) ----
flow_nodes <- tribble(
  ~id,          ~label,                                       ~n,       ~type,
  "candidates", "ICU 2018–2024 (candidates)",                 665074L,  "include",
  "age_in",     "≥18 years",                                  665042L,  "include",
  "demo_in",    "Demographics present",                       663876L,  "include",
  "icu_ge24",   "ICU stay ≥24h",                              532905L,  "include",
  "geo_present","Geography present (county, census tract)",         522090L,  "include",
  "abg_in",     "ABG or continuous SpO₂ in 24h",              522004L,  "include",
  "arf_yes",    "Meets ARF criteria in 24h",                  166152L,  "include",
  # Exclusion branches
  "icu_lt24",   "ICU stay < 24h",                             130971L,  "exclude",
  "geo_missing","Missing geography",                           10815L,  "exclude",
  "no_abg",     "No ABG or continuous SpO₂",                      86L,  "exclude",
  "arf_no",     "No ARF criteria in ±24h",                    355852L,  "exclude"
)

flow_edges <- tribble(
  ~from,        ~to,
  "candidates", "age_in",
  "age_in",     "demo_in",
  "demo_in",    "icu_ge24",
  "demo_in",    "icu_lt24",
  "icu_ge24",   "geo_present",
  "icu_ge24",   "geo_missing",
  "geo_present","abg_in",
  "geo_present","no_abg",
  "abg_in",     "arf_yes",
  "abg_in",     "arf_no"
)

# ---- 2) Build node labels with counts ----
node_lab <- flow_nodes %>%
  mutate(
    label_full = paste0(label, "\n", "n = ", comma(n))
  )

# ---- 3) Compose Graphviz DOT ----
# Styling
fill_include <- "#F5F5F5"  # light gray
fill_exclude <- "#FDECEA"  # light red
border_excl  <- "#D93025"  # red border for exclusions

# Node lines
node_lines <- node_lab %>%
  mutate(
    fill = if_else(type == "include", fill_include, fill_exclude),
    color = if_else(type == "include", "black", border_excl),
    line  = sprintf('%s [label="%s", shape=box, style="rounded,filled", color="%s", fillcolor="%s"];',
                    id, label_full, color, fill)
  ) %>%
  pull(line) %>%
  paste(collapse = "\n  ")

# Edge lines
edge_lines <- flow_edges %>%
  transmute(line = sprintf("%s -> %s;", from, to)) %>%
  pull(line) %>%
  paste(collapse = "\n  ")

# Ranks to align siblings horizontally (optional but prettier)
rank_same_blocks <- paste(
  "{rank=same; icu_ge24; icu_lt24}",
  "{rank=same; geo_present; geo_missing}",
  "{rank=same; abg_in; no_abg}",
  "{rank=same; arf_yes; arf_no}",
  sep = "\n  "
)

dot <- sprintf('digraph cohort_flow {
  graph [rankdir=TB, nodesep="0.35", ranksep="0.5", fontsize=10];
  node  [fontname="Helvetica", fontsize=10, margin="0.05,0.05"];
  edge  [color="black", penwidth=1.2, arrowsize=0.7, arrowhead=normal];

  %s

  %s

  %s
}', node_lines, edge_lines, rank_same_blocks)

# ---- 4) Render in RStudio/Viewer ----
g <- DiagrammeR::grViz(dot)
g

# ---- 5) Export high-resolution PNG ----
svg_txt <- DiagrammeRsvg::export_svg(g)
rsvg::rsvg_png(charToRaw(svg_txt),
               file = "unified_cohort_flow.png",
               width = 2400, height = 3000)  # ~8x10 in at 300 dpi

# The PNG "unified_cohort_flow.png" is now in your working directory.
# If you also want an SVG:
rsvg::rsvg_svg(charToRaw(svg_txt), file = "unified_cohort_flow.svg")


# ---- Packages ----

options(tigris_use_cache = TRUE)
theme_set(theme_minimal(base_size = 11))

# ---- Get CONUS counties (unchanged) ----
# Replace your earlier call with a pre-2023 vintage (e.g., 2019 or 2020)
get_conus_counties <- function(year = 2019) {
  drop_states <- c("02","15","72","60","66","69","78") # AK, HI, PR, territories
  tigris::counties(cb = TRUE, year = year, class = "sf") |>
    dplyr::filter(!STATEFP %in% drop_states) |>
    dplyr::mutate(GEOID = as.character(GEOID)) |>
    sf::st_transform(5070)
}

us_counties <- get_conus_counties(year = 2019)

# ---- Summarize to time-average per county (with safe GEOID handling) ----
summarize_county_mean <- function(df, years = 2005:2024,
                                  value_candidates = c("value","no2_mean","pm25_mean")) {
  stopifnot(all(c("GEOID","year") %in% names(df)))
  # ensure 5-char, zero-padded character GEOID
  df <- df |>
    mutate(
      GEOID = str_pad(as.character(GEOID), width = 5, side = "left", pad = "0")
    )
  
  val_col <- value_candidates[value_candidates %in% names(df)][1]
  if (is.na(val_col)) stop("No value column found in: ", paste(names(df), collapse=", "))
  
  df |>
    filter(year %in% years) |>
    group_by(GEOID) |>
    summarise(
      mean_value = mean(.data[[val_col]], na.rm = TRUE),
      n_years    = sum(!is.na(.data[[val_col]]) & year %in% years),
      .groups = "drop"
    )
}

# ---- Mapper with taller legend and dense breaks ----
make_pollutant_map <- function(summary_df, title_expr, legend_expr,
                               filename = NULL,
                               cmap = inferno(256), na_fill = "grey85",
                               legend_breaks = NULL) {
  
  # join (GEOID already character in both)
  plot_df <- us_counties |>
    left_join(summary_df, by = "GEOID")
  
  # choose breaks if not supplied
  if (is.null(legend_breaks)) {
    rng <- range(plot_df$mean_value, na.rm = TRUE)
    legend_breaks <- pretty(rng, n = 9)  # more ticks
  }
  
  p <- ggplot(plot_df) +
    geom_sf(aes(fill = mean_value), color = NA) +
    scale_fill_gradientn(
      colors = cmap,
      na.value = na_fill,
      breaks = legend_breaks,
      labels = label_number(accuracy = 0.1),    # ~1 decimal
      name   = legend_expr
    ) +
    labs(
      title    = title_expr,
      subtitle = "County-level average, 2018–2024"     # adjust if you map a different span
    ) +
    guides(
      fill = guide_colorbar(
        barheight = unit(7.5, "cm"),   # longer legend bar
        ticks.colour = "grey30",
        frame.colour = "grey30"
      )
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      plot.title   = element_text(face = "bold")
    )
  
  if (!is.null(filename)) {
    ggsave(filename, p, width = 12, height = 7, dpi = 600)
  }
  p
}

# ---- Build the summaries (your existing data frames) ----
no2_cnty_mean  <- summarize_county_mean(no2_us,  years = 2018:2024, value_candidates = c("no2_mean","value"))
pm25_cnty_mean <- summarize_county_mean(pm25_us, years = 2018:2024, value_candidates = c("pm25_mean","value"))

# ---- Draw PM2.5 map with correct math labels & longer legend ----
p_pm25 <- make_pollutant_map(
  pm25_cnty_mean,
  title_expr  = expression("Average PM"[2.5]*" by County, 2018–2024"),
  legend_expr = expression(PM[2.5]~"(" * mu * "g/" * m^3 * ")"),
  filename    = "map_pm25_conus_2018_2024.png"
)

# ---- Draw NO2 map with math labels & longer legend ----
p_no2 <- make_pollutant_map(
  no2_cnty_mean,
  title_expr  = expression("Average NO"[2]*" by County, 2018–2024"),
  legend_expr = expression(NO[2]~"(ppb)"),
  filename    = "map_no2_conus_2018_2024.png"
)

# Print to viewer
p_pm25
p_no2

# ---- Optional: quickly check for any lingering NA joins (should be 0 except places with no data) ----
# Which counties are still NA after join?
 anti_df <- us_counties %>%
   left_join(pm25_cnty_mean, by = "GEOID") %>%
   st_drop_geometry() %>%
   filter(is.na(mean_value)) %>% select(STATEFP, COUNTYFP, GEOID, NAME)
 nrow(anti_df); head(anti_df)

## hard cutpoints for no2 and pm2.5

 # ---------------------------------------------------------
 # 0) Paths and basic output dir
 # ---------------------------------------------------------
 cif_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/cif"
 out_dir <- file.path(cif_dir, "fig_out")
 dir.create(out_dir, showWarnings = FALSE)
 
 # ---------------------------------------------------------
 # 1) Read per-site exposure-bin CSVs (bins_all)
 # ---------------------------------------------------------
 bin_regex <- "site_exposure_bins__.*_cumulative(?: copy \\d+)?\\.csv$"
 bin_files <- fs::dir_ls(cif_dir, regexp = bin_regex, type = "file", recurse = TRUE)
 
 read_one_bins <- function(fp) {
   df_raw <- readr::read_csv(fp, show_col_types = FALSE)
   names(df_raw) <- names(df_raw) |>
     tolower() |>
     gsub("[^a-z0-9]+","_", x = _) |>
     gsub("^_|_$","", x = _)
   
   df <- df_raw %>%
     dplyr::rename(
       site_id       = dplyr::any_of("site_id"),
       site_name     = dplyr::any_of("site_name"),
       exposure_var  = dplyr::any_of("exposure_var"),
       exposure_lab  = dplyr::any_of("exposure_label"),
       unit          = dplyr::any_of("unit"),
       exposure_grp  = dplyr::any_of("exposure_group"),
       lower         = dplyr::any_of("lower"),
       upper         = dplyr::any_of("upper"),
       n_in_group    = dplyr::any_of("n_in_group")
     )
   
   # enforce single site_name per file
   sn <- df$site_name |> unique() |> na.omit()
   sn <- if (length(sn) >= 1L) sn[1] else NA_character_
   df <- df %>% dplyr::mutate(site_name = sn)
   
   pol <- if (grepl("no2", tolower(fp))) "NO2" else "PM2.5"
   
   df %>%
     dplyr::mutate(
       pollutant  = pol,
       lower_num  = suppressWarnings(as.numeric(lower)),
       upper_num  = suppressWarnings(as.numeric(upper)),
       mid        = (lower_num + upper_num)/2,
       source_file = fs::path_file(fp)
     )
 }
 
 bins_all <- purrr::map_dfr(bin_files, read_one_bins)
 
 # ---------------------------------------------------------
 # 2) Define HARD CUTPOINTS for NO2 & PM2.5
 #    (you can tweak these if you prefer different thresholds)
 # ---------------------------------------------------------
 
 # NO2 (ppb): Low <10, Mid 10–20, High >20
 no2_breaks  <- c(-Inf, 2.5, 7, Inf)
 
 # PM2.5 (µg/m³): Low <7, Mid 7–12, High >12
 pm25_breaks <- c(-Inf, 6, 9, Inf)
 
 # ---------------------------------------------------------
 # 3) Map each site's mid-bin value to Low / Mid / High using these cutpoints
 # ---------------------------------------------------------
 
 bins_us_fixed <- bins_all %>%
   dplyr::mutate(
     exposure_bin_fixed = dplyr::case_when(
       pollutant == "NO2" ~ as.character(
         cut(mid,
             breaks = no2_breaks,
             labels = c("Low","Mid","High"),
             right  = TRUE)
       ),
       pollutant == "PM2.5" ~ as.character(
         cut(mid,
             breaks = pm25_breaks,
             labels = c("Low","Mid","High"),
             right  = TRUE)
       ),
       TRUE ~ NA_character_
     ),
     exposure_bin_fixed = factor(exposure_bin_fixed,
                                 levels = c("Low","Mid","High"))
   ) %>%
   dplyr::select(
     site_id, site_name, pollutant, exposure_var,
     exposure_group = exposure_grp, mid, exposure_bin_fixed
   ) %>%
   dplyr::distinct()
 
 # ---------------------------------------------------------
 # 4) Attach fixed bins to CIF rows and pool across sites
 # ---------------------------------------------------------
 # cif_all is assumed to already exist from your earlier code.
 
 cif_fixed <- cif_all %>%
   dplyr::left_join(
     bins_us_fixed,
     by = c("site_id","site_name","pollutant","exposure_var","exposure_group")
   ) %>%
   dplyr::mutate(
     exposure_bin_fixed = factor(exposure_bin_fixed,
                                 levels = c("Low","Mid","High"))
   ) %>%
   dplyr::filter(!is.na(exposure_bin_fixed), !is.na(cause_label))
 
 pooled_fixed <- cif_fixed %>%
   dplyr::group_by(pollutant, cause_label, exposure_bin_fixed, day) %>%
   dplyr::summarise(
     n_sites        = dplyr::n_distinct(site_name),
     total_at_risk  = sum(risk_set, na.rm = TRUE),
     cif            = stats::weighted.mean(cif, w = pmax(risk_set, 0), na.rm = TRUE),
     .groups = "drop"
   )
 
 # ---------------------------------------------------------
 # 5) Smooth CIFs (LOESS, monotone) up to 30 days
 # ---------------------------------------------------------
 smooth_pooled_fixed <- function(df, span = 0.25, day_max = 30) {
   df %>%
     dplyr::filter(day <= day_max) %>%
     dplyr::arrange(pollutant, cause_label, exposure_bin_fixed, day) %>%
     dplyr::group_by(pollutant, cause_label, exposure_bin_fixed) %>%
     dplyr::group_modify(~{
       d <- .x
       if (nrow(d) >= 5) {
         fit  <- stats::loess(cif ~ day, data = d, span = span,
                              degree = 2, surface = "direct")
         yhat <- stats::predict(fit, newdata = tibble(day = d$day))
       } else {
         yhat <- d$cif
       }
       # clamp to [0,1] and enforce non-decreasing
       yhat <- pmin(pmax(yhat, 0), 1)
       yhat <- cummax(yhat)
       d$ci_smooth <- yhat
       d
     }) %>%
     dplyr::ungroup()
 }
 
 pooled30_fixed_sm <- smooth_pooled_fixed(pooled_fixed, span = 0.25, day_max = 30)
 
 # ---------------------------------------------------------
 # 6) Legend labels from cutpoints
 # ---------------------------------------------------------
 legend_labels_fixed <- function(pol) {
   if (pol == "NO2") {
     c(
       "Low"  = sprintf("Low (<%g ppb)",  no2_breaks[2]),
       "Mid"  = sprintf("Mid (%g–%g ppb)", no2_breaks[2], no2_breaks[3]),
       "High" = sprintf("High (≥%g ppb)", no2_breaks[3])
     )
   } else {
     c(
       "Low"  = sprintf("Low (<%g \u00B5g/m\u00B3)",  pm25_breaks[2]),
       "Mid"  = sprintf("Mid (%g–%g \u00B5g/m\u00B3)", pm25_breaks[2], pm25_breaks[3]),
       "High" = sprintf("High (≥%g \u00B5g/m\u00B3)", pm25_breaks[3])
     )
   }
 }
 
 # color palette: Low=blue, Mid=green, High=red
 pal_fixed <- c("Low" = "#1F77B4", "Mid" = "#2CA02C", "High" = "#D62728")
 
 # ---------------------------------------------------------
 # 7) Build combined dataset: Death (CIF) & Extubation (1 − CIF)
 # ---------------------------------------------------------
 prep_combined_fixed <- function(df, pol, day_max = 30) {
   dd <- df %>%
     dplyr::filter(
       pollutant   == pol,
       cause_label %in% c("Death", "Successful extubation"),
       day <= day_max
     )
   
   death <- dd %>%
     dplyr::filter(cause_label == "Death") %>%
     dplyr::mutate(outcome = "Death (CIF)")
   
   ext1m <- dd %>%
     dplyr::filter(cause_label == "Successful extubation") %>%
     dplyr::mutate(
       ci_smooth = 1 - ci_smooth,
       outcome   = "Extubation (1 \u2212 CIF)"
     )
   
   dplyr::bind_rows(death, ext1m)
 }
 
 # ---------------------------------------------------------
 # 8) Plot function: competing risks on common scale, fixed cutpoints
 # ---------------------------------------------------------
 make_competing_fixed <- function(df_pol, pol) {
   labs_color <- legend_labels_fixed(pol)
   y_top <- ceiling(max(df_pol$ci_smooth, na.rm = TRUE) * 20) / 20
   
   ggplot(df_pol,
          aes(x = day, y = ci_smooth,
              color    = exposure_bin_fixed,
              linetype = outcome,
              group    = interaction(exposure_bin_fixed, outcome))) +
     geom_line(linewidth = 1.35) +
     scale_color_manual(
       values = pal_fixed,
       breaks = c("Low","Mid","High"),
       labels = labs_color[c("Low","Mid","High")],
       name   = NULL
     ) +
     scale_linetype_manual(
       values = c("Death (CIF)" = "solid",
                  "Extubation (1 \u2212 CIF)" = "11"),
       breaks = c("Death (CIF)", "Extubation (1 \u2212 CIF)"),
       labels = c("Death (CIF)", "Extubation (1 \u2212 CIF)"),
       name   = NULL
     ) +
     scale_x_continuous(
       limits = c(0, 30),
       breaks = seq(0, 30, 5),
       expand = expansion(mult = c(0, 0.02))
     ) +
     coord_cartesian(ylim = c(0, y_top)) +
     labs(
       subtitle = if (pol == "NO2") expression(NO[2]) else expression(PM[2.5]),
       x = "Days since ICU admission",
       y = "Cumulative Incidence"
     ) +
     theme_classic(base_size = 12.5) +
     theme(
       panel.border        = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
       panel.grid.major.y  = element_line(colour = "grey90", linewidth = 0.35),
       plot.subtitle       = element_text(face = "bold", hjust = 0.5, margin = margin(b = 2)),
       legend.position     = "bottom",
       legend.box          = "vertical",
       legend.key.width    = unit(26, "pt"),
       legend.key.height   = unit(8,  "pt"),
       legend.spacing.y    = unit(2,  "pt"),
       plot.margin         = margin(4, 6, 4, 6)
     ) +
     guides(
       color = guide_legend(
         order = 1,
         ncol  = 2,
         byrow = TRUE,
         override.aes = list(
           linetype = "solid",
           linewidth = 1.35
         )
       ),
       linetype = guide_legend(
         order = 2,
         ncol  = 2,
         override.aes = list(
           color    = "black",
           linewidth= 1.6,
           linetype = c("solid", "11")
         )
       )
     )
 }
 
 # ---------------------------------------------------------
 # 9) Build and save the final two-panel figure
 # ---------------------------------------------------------
 df_no2_fixed   <- prep_combined_fixed(pooled30_fixed_sm, "NO2")
 df_pm25_fixed  <- prep_combined_fixed(pooled30_fixed_sm, "PM2.5")
 
 p_no2_fixed  <- make_competing_fixed(df_no2_fixed,  "NO2")
 p_pm25_fixed <- make_competing_fixed(df_pm25_fixed, "PM2.5")
 
 final_two_panel_fixed <- p_no2_fixed | p_pm25_fixed +
   patchwork::plot_annotation(
     title    = "Competing risks on a common scale (fixed NO\u2082 and PM\u2082.\u2085 cutpoints)",
     subtitle = "Solid: CIF(Death) \u00b7 Dashed: 1 \u2212 CIF(Extubation) \u00b7 Color: national exposure categories",
     theme = theme(
       plot.title   = element_text(face = "bold", size = 14, hjust = 0),
       plot.subtitle= element_text(size = 11, margin = margin(t = 2, b = 6))
     )
   ) &
   theme(legend.position = "bottom")
 
 ggsave(file.path(out_dir, "CIF_competing__NO2_PM25_overlaid__FIXEDCUTS.png"),
        final_two_panel_fixed, width = 11.8, height = 5.4, dpi = 600)
 ggsave(file.path(out_dir, "CIF_competing__NO2_PM25_overlaid__FIXEDCUTS.pdf"),
        final_two_panel_fixed, width = 11.8, height = 5.4, device = cairo_pdf)
 


##### subtype pooling
 
 library(tidyverse)
 
 # ---- inverse-variance pooling on link scale ----
 pool_link <- function(pred, lo, hi, link = c("logit", "log")) {
   link <- match.arg(link)
   
   if (link == "logit") {
     est <- qlogis(pred)
     se  <- (qlogis(hi) - qlogis(lo)) / (2 * 1.96)
     inv <- plogis
   } else {
     est <- log(pred)
     se  <- (log(hi) - log(lo)) / (2 * 1.96)
     inv <- exp
   }
   
   w <- 1 / se^2
   mu <- sum(w * est, na.rm = TRUE) / sum(w, na.rm = TRUE)
   se_mu <- sqrt(1 / sum(w, na.rm = TRUE))
   
   tibble(
     pred = inv(mu),
     lo   = inv(mu - 1.96 * se_mu),
     hi   = inv(mu + 1.96 * se_mu)
   )
 }
 
 pool_curves <- function(files, group_var, xvar, link) {
   
   read_and_tag <- function(f) {
     read_csv(f, show_col_types = FALSE) %>%
       mutate(site = tools::file_path_sans_ext(basename(f)))
   }
   
   dat <- map_dfr(files, read_and_tag)
   
   dat %>%
     group_by(.data[[group_var]], .data[[xvar]]) %>%
     summarise(
       pool = list(pool_link(pred, lo, hi, link)),
       .groups = "drop"
     ) %>%
     unnest(pool)
 }
 
 
 # -------------------------
 # Cleaning helpers
 # -------------------------
 clean_sex <- function(x) {
   x_chr <- as.character(x)
   x_chr <- str_trim(x_chr)
   x_chr <- str_to_lower(x_chr)
   
   dplyr::case_when(
     x_chr %in% c("m", "male") ~ "Male",
     x_chr %in% c("f", "female") ~ "Female",
     x_chr %in% c("unknown", "unk", "u", "na", "n/a", "", "missing") ~ "Unknown",
     is.na(x_chr) ~ "Unknown",
     TRUE ~ str_to_title(x_chr)
   )
 }
 
 clean_race <- function(x) {
   x_chr <- as.character(x)
   x_chr <- str_trim(x_chr)
   x_chr <- str_replace_all(x_chr, "\\s+", " ")
   
   # standardize common patterns before title-casing
   x_chr <- str_replace_all(x_chr, regex("^non[- ]?hispanic", ignore_case = TRUE), "Non-Hispanic")
   x_chr <- str_replace_all(x_chr, regex("^hispanic", ignore_case = TRUE), "Hispanic")
   
   str_to_title(x_chr)
 }
 
 # -------------------------
 # Read + pool helpers
 # -------------------------
 read_site_preds <- function(file) {
   readr::read_csv(file, show_col_types = FALSE) %>%
     mutate(site_file = basename(file))
 }
 
 # This assumes your pool_curves() already exists and returns columns:
 #   no2_10, pred, lo, hi, and the group variable
 # If not, adjust below accordingly.
 
prep_ribbon_df <- function(df, group_var, group_clean_fn = NULL, levels = NULL) {
   group_var <- rlang::ensym(group_var)
   group_nm  <- rlang::as_name(group_var)
   
   out <- df %>%
     mutate(
       no2_10 = as.numeric(.data$no2_10),
       pred   = as.numeric(.data$pred),
       lo     = as.numeric(.data$lo),
       hi     = as.numeric(.data$hi)
     ) %>%
     filter(is.finite(no2_10), is.finite(pred), is.finite(lo), is.finite(hi)) %>%
     mutate(
       # enforce monotone CI ordering
       lo = pmin(lo, hi),
       hi = pmax(lo, hi),
       # keep pred within bounds (helps if rounding/transform artifacts exist)
       pred = pmax(lo, pmin(pred, hi))
     )
   
   if (!is.null(group_clean_fn)) {
     out[[group_nm]] <- group_clean_fn(out[[group_nm]])
   }
   
   # CRITICAL: ensure ONE row per group × x (ribbon cannot behave well otherwise)
   out <- out %>%
     group_by(.data[[group_nm]], no2_10) %>%
     summarise(
       pred = mean(pred, na.rm = TRUE),
       lo   = mean(lo,   na.rm = TRUE),
       hi   = mean(hi,   na.rm = TRUE),
       .groups = "drop"
     ) %>%
     arrange(.data[[group_nm]], no2_10)
   
   if (!is.null(levels)) {
     out[[group_nm]] <- factor(out[[group_nm]], levels = levels)
   }
   
   out
 }
 
 plot_ribbon_pub <- function(df, group_var, title_expr, ylab, is_count = FALSE) {
   group_var <- rlang::ensym(group_var)
   
   ggplot(df, aes(x = no2_10, y = pred, color = !!group_var, fill = !!group_var, group = !!group_var)) +
     geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.07, color = NA) +
     geom_line(linewidth = 1.1) +
     labs(
       title = title_expr,
       x     = expression(NO[2] * " (per 10 ppb)"),
       y     = ylab,
       color = NULL,
       fill  = NULL
     ) +
     theme_classic(base_size = 13) +
     theme(
       plot.title = element_text(face = "bold"),
       axis.title = element_text(face = "bold"),
       legend.position = "right"
     )
 }
 
 # -------------------------
 # FILES
 # -------------------------
 subtype_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype"
 
 sex_30d_files <- list.files(subtype_dir, pattern = "site_preds_30d_by_sex.*\\.csv$", full.names = TRUE)
 race_vent_files <- list.files(subtype_dir, pattern = "site_preds_vent_by_race.*\\.csv$", full.names = TRUE)
 
 # -------------------------
 # POOL + PLOT: 30d mortality by sex
 # -------------------------
 pool_30d_sex_raw <- pool_curves(
   files     = sex_30d_files,
   group_var = "sex_category",
   xvar      = "no2_10",
   link      = "logit"
 )
 
 pool_30d_sex <- prep_ribbon_df(
   pool_30d_sex_raw,
   group_var = sex_category,
   group_clean_fn = clean_sex,
   levels = c("Female", "Male", "Unknown")
 )
 
 p_sex_30d <- plot_ribbon_pub(
   pool_30d_sex,
   group_var = sex_category,
   title_expr = expression("30-day mortality vs " * NO[2] * " by sex"),
   ylab = "Predicted probability"
 )
 
 print(p_sex_30d)
 
 # -------------------------
 # POOL + PLOT: Ventilation hours by race/ethnicity (log link)
 # -------------------------
 pool_vent_race_raw <- pool_curves(
   files     = race_vent_files,
   group_var = "race_ethnicity_simple",
   xvar      = "no2_10",
   link      = "log"
 )
 
 pool_vent_race <- prep_ribbon_df(
   pool_vent_race_raw,
   group_var = race_ethnicity_simple,
   group_clean_fn = clean_race
 )
 
 p_race_vent <- plot_ribbon_pub(
   pool_vent_race,
   group_var = race_ethnicity_simple,
   title_expr = expression("Ventilation hours vs " * NO[2] * " by race/ethnicity"),
   ylab = "Predicted mean ventilation hours"
 )
 
 print(p_race_vent)
 
 
 
 
 
 
 library(dplyr)
 library(stringr)
 library(readr)
 library(purrr)
 library(ggplot2)
 library(rlang)
 
 # -------------------------
 # Cleaning helpers
 # -------------------------
 clean_sex <- function(x) {
   x_chr <- as.character(x)
   x_chr <- stringr::str_trim(x_chr)
   x_chr <- stringr::str_to_lower(x_chr)
   
   dplyr::case_when(
     x_chr %in% c("m", "male") ~ "Male",
     x_chr %in% c("f", "female") ~ "Female",
     x_chr %in% c("unknown", "unk", "u", "na", "n/a", "", "missing") ~ "Unknown",
     is.na(x) ~ "Unknown",
     TRUE ~ stringr::str_to_title(x_chr)
   )
 }
 
 clean_race <- function(x) {
   x_chr <- as.character(x)
   x_chr <- stringr::str_trim(x_chr)
   x_chr <- stringr::str_replace_all(x_chr, "\\s+", " ")
   x_chr <- stringr::str_replace_all(
     x_chr,
     stringr::regex("^non[- ]?hispanic", ignore_case = TRUE),
     "Non-Hispanic"
   )
   x_chr <- stringr::str_replace_all(
     x_chr,
     stringr::regex("^hispanic", ignore_case = TRUE),
     "Hispanic"
   )
   stringr::str_to_title(x_chr)
 }
 
 # -------------------------
 # Pool/plot helpers
 # -------------------------
 prep_ribbon_df <- function(df, group_var, group_clean_fn = NULL, levels = NULL,
                            drop_levels = NULL) {
   group_var <- rlang::ensym(group_var)
   group_nm  <- rlang::as_name(group_var)
   
   out <- df %>%
     dplyr::mutate(
       no2_10 = as.numeric(.data$no2_10),
       pred   = as.numeric(.data$pred),
       lo     = as.numeric(.data$lo),
       hi     = as.numeric(.data$hi)
     ) %>%
     dplyr::filter(is.finite(no2_10), is.finite(pred), is.finite(lo), is.finite(hi)) %>%
     dplyr::mutate(
       lo = pmin(lo, hi),
       hi = pmax(lo, hi),
       pred = pmax(lo, pmin(pred, hi))
     )
   
   if (!is.null(group_clean_fn)) {
     out[[group_nm]] <- group_clean_fn(out[[group_nm]])
   }
   
   # DROP unwanted categories (e.g., Unknown for sex)
   if (!is.null(drop_levels)) {
     out <- out %>% dplyr::filter(!(.data[[group_nm]] %in% drop_levels))
   }
   
   # critical: one row per group × x for clean ribbons
   out <- out %>%
     dplyr::group_by(.data[[group_nm]], no2_10) %>%
     dplyr::summarise(
       pred = mean(pred, na.rm = TRUE),
       lo   = mean(lo,   na.rm = TRUE),
       hi   = mean(hi,   na.rm = TRUE),
       .groups = "drop"
     ) %>%
     dplyr::arrange(.data[[group_nm]], no2_10)
   
   if (!is.null(levels)) {
     out[[group_nm]] <- factor(out[[group_nm]], levels = levels)
   }
   
   out
 }
 
 library(RColorBrewer)
 # Brewer palettes
 pal_sex  <- brewer.pal(3, "Set2")   # we'll only use first 2 after dropping Unknown
 pal_race <- brewer.pal(8, "Set1")   # up to 8 race categories
 
 plot_ribbon_pub <- function(df, group_var, title_expr, ylab, palette = NULL) {
   group_var <- rlang::ensym(group_var)
   group_nm  <- rlang::as_name(group_var)
   
   # If palette is named or longer than needed, align it to factor levels
   if (!is.null(palette)) {
     levs <- levels(df[[group_nm]])
     if (is.null(levs)) levs <- sort(unique(df[[group_nm]]))
     if (is.null(names(palette))) {
       palette_use <- palette[seq_len(min(length(palette), length(levs)))]
       names(palette_use) <- levs[seq_len(length(palette_use))]
     } else {
       palette_use <- palette[levs]
     }
   }
   
   p <- ggplot2::ggplot(
     df,
     ggplot2::aes(
       x = no2_10,
       y = pred,
       color = !!group_var,
       fill  = !!group_var,
       group = !!group_var
     )
   ) +
     ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.07, color = NA) +
     ggplot2::geom_line(linewidth = 1.15) +
     ggplot2::labs(
       title = title_expr,
       x     = expression(NO[2] * " (per 10 ppb)"),
       y     = ylab,
       color = NULL,
       fill  = NULL
     ) +
     ggplot2::theme_classic(base_size = 16) +
     ggplot2::theme(
       plot.title   = ggplot2::element_text(face = "bold"),
       axis.title   = ggplot2::element_text(face = "bold"),
       legend.position = "right"
     )
   
   if (!is.null(palette)) {
     p <- p +
       ggplot2::scale_color_manual(values = palette_use, drop = TRUE) +
       ggplot2::scale_fill_manual(values  = palette_use, drop = TRUE)
   }
   
   p
 }
 
 # -------------------------
 # Configuration
 # -------------------------
 subtype_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype"
 
 outcomes <- tibble::tribble(
   ~outcome_key, ~file_key, ~ylab,                              ~link,   ~title_stub,
   "30d",        "30d",     "Predicted probability",            "logit", "30-day mortality vs ",
   "inhosp",     "inhosp",  "Predicted probability",            "logit", "In-hospital mortality vs ",
   "vent",       "vent",    "Predicted mean ventilation hours", "log",   "Ventilation hours vs "
 )
 
 group_specs <- tibble::tribble(
   ~group_key, ~group_var,               ~clean_fn,   ~levels,                     ~drop_levels,
   "sex",      "sex_category",           clean_sex,   list(c("Female", "Male")),   list(c("Unknown")),
   "race",     "race_ethnicity_simple",  clean_race,  list(NULL),                  list(NULL)
 )
 
 # -------------------------
 # Run all: sex/race × 30d/inhosp/vent
 # -------------------------
 plots <- list()
 
 for (g in seq_len(nrow(group_specs))) {
   g_key    <- group_specs$group_key[g]
   g_var    <- group_specs$group_var[g]
   g_clean  <- group_specs$clean_fn[[g]]
   g_levels <- group_specs$levels[[g]][[1]]
   g_drop   <- group_specs$drop_levels[[g]][[1]]
   
   for (o in seq_len(nrow(outcomes))) {
     o_key  <- outcomes$outcome_key[o]
     f_key  <- outcomes$file_key[o]
     ylab   <- outcomes$ylab[o]
     link   <- outcomes$link[o]
     t_stub <- outcomes$title_stub[o]
     
     files <- list.files(
       subtype_dir,
       pattern = paste0("^site_preds_", f_key, "_by_", g_key, ".*\\.csv$"),
       full.names = TRUE
     )
     
     if (length(files) == 0) {
       message("No files found for ", g_key, " x ", o_key,
               " (pattern: site_preds_", f_key, "_by_", g_key, "...)")
       next
     }
     
     pooled_raw <- pool_curves(
       files     = files,
       group_var = g_var,
       xvar      = "no2_10",
       link      = link
     )
     
     pooled <- prep_ribbon_df(
       pooled_raw,
       group_var = !!rlang::sym(g_var),
       group_clean_fn = g_clean,
       levels = g_levels,
       drop_levels = g_drop
     )
     
     title_expr <- if (g_key == "sex") {
       as.expression(bquote(.(t_stub) * NO[2] * " by sex"))
     } else {
       as.expression(bquote(.(t_stub) * NO[2] * " by race/ethnicity"))
     }
     
     p <- plot_ribbon_pub(
       pooled,
       group_var = !!rlang::sym(g_var),
       title_expr = title_expr,
       ylab = ylab,
       palette = if (g_key == "sex") pal_sex else pal_race
     )
     
     plot_name <- paste0(g_key, "_", o_key)
     plots[[plot_name]] <- p
     
     print(p)
   }
 }
 
 # -------------------------
 # Save figures (high-res)
 # -------------------------
 output_dir <- file.path(subtype_dir, "output")
 dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
 
 png_w <- 7.5   # inches
 png_h <- 5.0
 dpi   <- 600
 
 for (nm in names(plots)) {
   p <- plots[[nm]]
   
   ggsave(
     filename = file.path(output_dir, paste0(nm, ".png")),
     plot     = p,
     width    = png_w,
     height   = png_h,
     units    = "in",
     dpi      = dpi,
     bg       = "white"
   )
   
   ggsave(
     filename = file.path(output_dir, paste0(nm, ".pdf")),
     plot     = p,
     width    = png_w,
     height   = png_h,
     units    = "in",
     device   = cairo_pdf
   )
 }
 
 
 
####### subtype pooled 
 
 # ---- 1) Locate STRICT x_subtype model tidy CSVs (and validate) ----
 all_csvs <- list.files(
   subtype_dir,
   pattern = "\\.csv$",
   full.names = TRUE
 )
 
 # Strict: filename must include x_subtype
 xsub_files <- all_csvs[str_detect(basename(all_csvs), "x_subtype")]
 
 # Keep only the three outcomes of interest
 xsub_files <- xsub_files[
   str_detect(basename(xsub_files), "death30d|inhosp_death|vent_hours")
 ]
 
 stopifnot(length(xsub_files) > 0)
 
 # Validate file contents so we don't accidentally include wrong csvs
 is_valid_xsub <- function(f) {
   df <- tryCatch(
     suppressMessages(readr::read_csv(f, show_col_types = FALSE)),
     error = function(e) NULL
   )
   if (is.null(df)) return(FALSE)
   
   # Required columns
   needed <- c("term", "estimate", "std.error")
   if (!all(needed %in% names(df))) return(FALSE)
   
   # Must contain main NO2 effect
   if (!("no2_10" %in% df$term)) return(FALSE)
   
   # Must contain at least one NO2 x subtype interaction
   if (!any(str_detect(df$term, "^no2_10:arf_subtype"))) return(FALSE)
   
   TRUE
 }
 
 subtype_files <- xsub_files[vapply(xsub_files, is_valid_xsub, logical(1))]
 
 stopifnot(length(subtype_files) > 0)
 
 message("Using ", length(subtype_files), " validated x_subtype files.")
 
 
 # ---- Helpers ----
 clean_subtype <- function(x) {
   x %>%
     as.character() %>%
     str_replace_all("_", " ") %>%
     str_trim() %>%
     str_to_title()
 }
 
 # Parse outcome from filename
 # Expected patterns (examples):
 # refer_Emory_..._model_tidy_death30d_x_subtype_....csv
 # refer_Emory_..._model_tidy_inhosp_death_x_subtype_....csv
 # refer_Emory_..._model_tidy_vent_hours_x_subtype_....csv
 infer_outcome <- function(fn) {
   fn <- basename(fn)
   if (str_detect(fn, "death30d")) return("30-day mortality")
   if (str_detect(fn, "inhosp_death")) return("In-hospital mortality")
   if (str_detect(fn, "vent_hours")) return("Ventilation hours")
   return(NA_character_)
 }
 
 # Link type per outcome
 # (death outcomes: logit => OR; vent_hours: log => IRR)
 infer_link <- function(outcome) {
   if (outcome %in% c("30-day mortality","In-hospital mortality")) return("logit")
   if (outcome %in% c("Ventilation hours")) return("log")
   return(NA_character_)
 }
 
 # Try to infer site from filename (assumes refer_<SITE>_...)
 infer_site <- function(fn) {
   fn <- basename(fn)
   # e.g., "refer_Emory_20251008_..." -> "Emory"
   m <- str_match(fn, "^refer_([^_]+)_")
   if (!is.na(m[,2])) return(m[,2])
   # fallback: first chunk
   return(str_split(fn, "_", simplify = TRUE)[1])
 }
 
 # Extract a single term row from tidy df; returns 1-row tibble or empty tibble
 get_term <- function(df, term_name) {
   df %>% filter(.data$term == term_name) %>% slice(1)
 }
 
 # Compute subtype-specific NO2 effect (beta + se) for one site/outcome
 # Assumes:
 # - main effect term: "no2_10"
 # - interaction terms: "no2_10:arf_subtypeHypercapnic", "no2_10:arf_subtypeMixed", etc.
 # - reference subtype is the baseline (not appearing as arf_subtype... term)
 compute_site_subtype_effects <- function(tidy_df, site_id, outcome) {
   
   # Pull main NO2
   no2_row <- get_term(tidy_df, "no2_10")
   if (nrow(no2_row) == 0) return(tibble())  # cannot compute anything
   
   b0  <- no2_row$estimate[[1]]
   se0 <- no2_row$std.error[[1]]
   
   # Find which subtype levels exist in this model output
   # Main-effect subtype terms look like: "arf_subtypeHypercapnic", "arf_subtypeMixed"
   subtype_terms <- tidy_df %>%
     filter(str_detect(term, "^arf_subtype")) %>%
     pull(term)
   
   # Infer non-reference subtype names from those terms
   nonref_subtypes <- subtype_terms %>%
     str_replace("^arf_subtype", "") %>%
     unique()
   
   # Also consider interaction terms might exist even if main-effect subtype term isn't present
   int_terms <- tidy_df %>%
     filter(str_detect(term, "^no2_10:arf_subtype")) %>%
     pull(term)
   
   nonref_from_int <- int_terms %>%
     str_replace("^no2_10:arf_subtype", "") %>%
     unique()
   
   nonref_all <- unique(c(nonref_subtypes, nonref_from_int))
   nonref_all <- nonref_all[!is.na(nonref_all) & nzchar(nonref_all)]
   
   # Reference subtype label (assumed)
   # If your reference is NOT Hypoxemic, change this label for display only.
   ref_label <- "Hypoxemic"
   
   # Start with reference subtype effect = main NO2
   out <- tibble(
     site_id   = site_id,
     outcome   = outcome,
     subtype   = ref_label,
     beta      = b0,
     se        = se0
   )
   
   # Add each non-reference subtype: beta = b0 + b_int ; se = sqrt(se0^2 + se_int^2)
   if (length(nonref_all) > 0) {
     add <- map_dfr(nonref_all, function(st) {
       int_name <- paste0("no2_10:arf_subtype", st)
       int_row  <- get_term(tidy_df, int_name)
       
       # If interaction term missing, skip (cannot compute subtype-specific effect)
       if (nrow(int_row) == 0) return(tibble())
       
       b1  <- int_row$estimate[[1]]
       se1 <- int_row$std.error[[1]]
       
       tibble(
         site_id = site_id,
         outcome = outcome,
         subtype = clean_subtype(st),
         beta    = b0 + b1,
         se      = sqrt(se0^2 + se1^2)  # covariance unavailable in tidy csv
       )
     })
     
     out <- bind_rows(out, add)
   }
   
   # Clean subtype display for ref as well
   out <- out %>%
     mutate(subtype = if_else(subtype == "Hypoxemic", "Hypoxemic", subtype))
   
   out
 }
 
 # Random-effects meta-analysis pooling (DerSimonian-Laird by default in metafor::rma.uni)
 pool_by_outcome_subtype <- function(df) {
   # df has one row per site per subtype per outcome with beta & se on link scale
   df %>%
     group_by(outcome, subtype) %>%
     group_modify(~{
       dat <- .x %>% filter(is.finite(beta), is.finite(se), se > 0)
       
       if (nrow(dat) == 0) {
         return(tibble(k = 0, beta_pool = NA_real_, se_pool = NA_real_,
                       ci_lo = NA_real_, ci_hi = NA_real_, tau2 = NA_real_))
       }
       
       # If only 1 site, just pass through
       if (nrow(dat) == 1) {
         b  <- dat$beta[[1]]
         se <- dat$se[[1]]
         return(tibble(
           k = 1,
           beta_pool = b,
           se_pool   = se,
           ci_lo     = b - 1.96*se,
           ci_hi     = b + 1.96*se,
           tau2      = 0
         ))
       }
       
       fit <- metafor::rma.uni(yi = dat$beta, sei = dat$se, method = "DL")
       
       tibble(
         k         = nrow(dat),
         beta_pool = as.numeric(fit$b),
         se_pool   = as.numeric(fit$se),
         ci_lo     = as.numeric(fit$ci.lb),
         ci_hi     = as.numeric(fit$ci.ub),
         tau2      = as.numeric(fit$tau2)
       )
     }) %>%
     ungroup()
 }
 
 # Forest-style plot (pooled only)
 plot_pooled_forest <- function(pooled_df, outcome, link) {
   
   # Exponentiate for OR/IRR and build CIs
   lab <- if (link == "logit") "Pooled OR per 10 ppb NO\u2082"
   else if (link == "log") "Pooled IRR per 10 ppb NO\u2082"
   else "Pooled effect per 10 ppb NO\u2082"
   
   dfp <- pooled_df %>%
     filter(outcome == !!outcome) %>%
     mutate(
       eff    = exp(beta_pool),
       eff_lo = exp(ci_lo),
       eff_hi = exp(ci_hi),
       subtype = fct_relevel(subtype, "Hypoxemic")
     ) %>%
     arrange(subtype)
   
   ggplot(dfp, aes(y = subtype, x = eff)) +
     geom_vline(xintercept = 1, linetype = 3, linewidth = 0.5) +
     geom_errorbarh(aes(xmin = eff_lo, xmax = eff_hi), height = 0.18, linewidth = 0.7) +
     geom_point(size = 2.6) +
     scale_x_log10() +
     labs(
       title = paste0(outcome, ": pooled NO\u2082 effect by ARF subtype"),
       x = lab,
       y = NULL
     ) +
     theme_classic(base_size = 13) +
     theme(
       plot.title = element_text(face = "bold"),
       axis.title.x = element_text(face = "bold")
     )
 }
 
 
 
 # ---- 2) Read + compute site-level subtype-specific NO2 effects ----
 site_subtype_effects <- map_dfr(subtype_files, function(f) {
   outcome <- infer_outcome(f)
   site_id <- infer_site(f)
   
   df <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
   
   # Defensive: ensure expected columns exist
   needed <- c("term","estimate","std.error")
   if (!all(needed %in% names(df))) return(tibble())
   
   compute_site_subtype_effects(df, site_id = site_id, outcome = outcome)
 })
 
 stopifnot(nrow(site_subtype_effects) > 0)
 
 # ---- 3) Pool across sites within each outcome x subtype ----
 pooled_subtype <- pool_by_outcome_subtype(site_subtype_effects)
 
 # ---- 4) Make and save plots for each outcome ----
 outcomes <- c("30-day mortality","In-hospital mortality","Ventilation hours")
 
 for (oc in outcomes) {
   link <- infer_link(oc)
   
   p <- plot_pooled_forest(pooled_subtype, outcome = oc, link = link)
   
   fn <- paste0(
     "pooled_no2_effect_by_subtype_",
     case_when(
       oc == "30-day mortality" ~ "death30d",
       oc == "In-hospital mortality" ~ "inhosp_death",
       oc == "Ventilation hours" ~ "vent_hours",
       TRUE ~ "outcome"
     ),
     ".png"
   )
   
   ggsave(
     filename = file.path(out_dir, fn),
     plot     = p,
     width    = 7.5,
     height   = 4.2,
     dpi      = 600
   )
 }
 
 # ---- 5) Optional: also save a tidy pooled table ----
 pooled_table <- pooled_subtype %>%
   mutate(
     effect = exp(beta_pool),
     lo     = exp(ci_lo),
     hi     = exp(ci_hi),
     measure = case_when(
       outcome %in% c("30-day mortality","In-hospital mortality") ~ "OR",
       outcome == "Ventilation hours" ~ "IRR",
       TRUE ~ "Effect"
     )
   ) %>%
   select(outcome, subtype, k, measure, effect, lo, hi, tau2)
 
 readr::write_csv(pooled_table, file.path(out_dir, "pooled_no2_by_subtype_summary.csv"))
 
 message("Done. Outputs written to: ", out_dir)
 
 
 
 
 
 
 
 # ---- Libraries ----
 library(dplyr)
 library(stringr)
 library(readr)
 library(tidyr)
 library(purrr)
 library(ggplot2)
 
 # ---- Paths ----
 subtype_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype"
 out_dir <- file.path(subtype_dir, "output")
 dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
 
 # ---- Helper: validate file is an x_subtype tidy model with NO2 main + interaction ----
 is_valid_xsub <- function(f) {
   df <- tryCatch(
     suppressMessages(readr::read_csv(f, show_col_types = FALSE)),
     error = function(e) NULL
   )
   if (is.null(df)) return(FALSE)
   needed <- c("term", "estimate", "std.error")
   if (!all(needed %in% names(df))) return(FALSE)
   if (!("no2_10" %in% df$term)) return(FALSE)
   if (!any(str_detect(df$term, "^no2_10:arf_subtype"))) return(FALSE)
   TRUE
 }
 
 # ---- Helper: extract a readable site id from filename ----
 site_from_file <- function(f) {
   # e.g., refer_Emory_20251008_...csv -> Emory
   b <- basename(f)
   m <- str_match(b, "^refer_([^_]+)_")
   if (!is.na(m[,2])) return(m[,2])
   # fallback: first chunk
   str_split(b, "_", simplify = TRUE)[1]
 }
 
 # ---- Helper: fixed-effect meta-analysis on log scale ----
 meta_fixed <- function(beta, se) {
   w <- 1 / (se^2)
   mu <- sum(w * beta) / sum(w)
   se_mu <- sqrt(1 / sum(w))
   tibble(beta = mu, se = se_mu, lo = mu - 1.96 * se_mu, hi = mu + 1.96 * se_mu)
 }
 
 # ---- (Optional) random-effects DL (kept here if you want it) ----
 meta_random_DL <- function(beta, se) {
   w <- 1 / (se^2)
   mu_fe <- sum(w * beta) / sum(w)
   Q <- sum(w * (beta - mu_fe)^2)
   df <- length(beta) - 1
   C <- sum(w) - (sum(w^2) / sum(w))
   tau2 <- max(0, (Q - df) / C)
   w_re <- 1 / (se^2 + tau2)
   mu <- sum(w_re * beta) / sum(w_re)
   se_mu <- sqrt(1 / sum(w_re))
   tibble(beta = mu, se = se_mu, lo = mu - 1.96 * se_mu, hi = mu + 1.96 * se_mu, tau2 = tau2)
 }
 
 # ---- Core: compute subtype-specific NO2 betas per site ----
 compute_site_subtype_effects <- function(df, file) {
   
   site_id <- site_from_file(file)
   
   # Pull NO2 main effect (OR scale in CSV)
   no2_row <- df %>% filter(term == "no2_10") %>% slice(1)
   if (nrow(no2_row) == 0) return(NULL)
   
   beta_no2 <- log(no2_row$estimate)
   se_no2   <- (log(no2_row$conf.high) - log(no2_row$conf.low)) / (2 * 1.96)
   
   # Interaction terms
   inter <- df %>%
     filter(str_detect(term, "^no2_10:arf_subtype")) %>%
     mutate(
       subtype = str_replace(term, "^no2_10:arf_subtype", ""),
       beta_int = log(estimate),
       se_int   = (log(conf.high) - log(conf.low)) / (2 * 1.96)
     )
   
   # Reference subtype = model reference (Hypoxemic)
   out <- tibble(
     site_id    = site_id,
     arf_subtype = "Hypoxemic",
     beta       = beta_no2,
     se         = se_no2
   )
   
   # Add non-reference subtypes
   if (nrow(inter) > 0) {
     out <- bind_rows(
       out,
       inter %>%
         transmute(
           site_id = site_id,
           arf_subtype = str_to_title(subtype),
           beta = beta_no2 + beta_int,
           se   = sqrt(se_no2^2 + se_int^2)
         )
     )
   }
   
   out
 }
 
 # ---- Read & filter files ----
 all_csvs <- list.files(subtype_dir, pattern = "\\.csv$", full.names = TRUE)
 xsub_files <- all_csvs[str_detect(basename(all_csvs), "x_subtype")]
 xsub_files <- xsub_files[str_detect(basename(xsub_files), "death30d|inhosp_death|vent_hours")]
 xsub_files <- xsub_files[vapply(xsub_files, is_valid_xsub, logical(1))]
 
 stopifnot(length(xsub_files) > 0)
 message("Using ", length(xsub_files), " validated x_subtype files.")
 
 # ---- Attach outcome labels ----
 file_map <- tibble(
   file = xsub_files,
   outcome = case_when(
     str_detect(basename(file), "death30d")     ~ "30-day mortality",
     str_detect(basename(file), "inhosp_death") ~ "In-hospital mortality",
     str_detect(basename(file), "vent_hours")   ~ "Ventilation hours",
     TRUE ~ NA_character_
   )
 ) %>% filter(!is.na(outcome))
 
 # ---- Build per-site subtype-specific effects ----
 site_effects <- file_map %>%
   mutate(df = map(file, ~ suppressMessages(read_csv(.x, show_col_types = FALSE)))) %>%
   mutate(effects = map2(df, file, compute_site_subtype_effects)) %>%
   select(outcome, file, effects) %>%
   unnest(effects) %>%
   mutate(
     arf_subtype = str_replace_all(arf_subtype, "_", " "),
     arf_subtype = str_trim(arf_subtype)
   )
 
 stopifnot(nrow(site_effects) > 0)
 
 # ---- Pool across sites (log scale), then exponentiate once ----
 pooled <- site_effects %>%
   group_by(outcome, arf_subtype) %>%
   summarise(
     k = n(),
     meta = list(meta_fixed(beta, se)),  # swap to meta_random_DL(beta,se) if desired
     .groups = "drop"
   ) %>%
   unnest(meta) %>%
   mutate(
     effect = exp(beta),
     lo = exp(lo),
     hi = exp(hi),
     measure = case_when(
       outcome %in% c("30-day mortality", "In-hospital mortality") ~ "OR per 10 ppb NO\u2082",
       outcome == "Ventilation hours" ~ "IRR per 10 ppb NO\u2082"
     )
   )
 
 # ---- Plot: forest by outcome ----
 plot_forest <- function(df_out) {
   df_out <- df_out %>%
     arrange(effect) %>%
     mutate(arf_subtype = factor(arf_subtype, levels = arf_subtype))
   
   ggplot(df_out, aes(y = arf_subtype, x = effect)) +
     geom_vline(xintercept = 1, linetype = 3, linewidth = 0.6) +
     geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.8) +
     geom_point(size = 3.2) +
     scale_x_log10() +
     labs(
       x = unique(df_out$measure),
       y = NULL,
       title = paste0(unique(df_out$outcome), ": pooled NO\u2082 effect by ARF subtype"),
     ) +
     theme_classic(base_size = 14) +
     theme(
       plot.title = element_text(face = "bold", size = 18),
       axis.title.x = element_text(face = "bold"),
       axis.text.y = element_text(size = 13)
     )
 }
 
 # Save high-res
 for (o in unique(pooled$outcome)) {
   p <- plot_forest(pooled %>% filter(outcome == o))
   fn <- file.path(out_dir, paste0("pooled_no2_effect_by_subtype_", gsub("[^A-Za-z0-9]+","_", tolower(o)), ".png"))
   ggsave(fn, p, width = 11, height = 6.5, dpi = 600)
   message("Saved: ", fn)
 }
 
 # Optional: write pooled table
 write_csv(pooled, file.path(out_dir, "pooled_no2_effect_by_subtype_all_outcomes.csv"))
 
 
 
###### facet wrapped 
 
 
 subtype_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype"
 out_dir <- file.path(subtype_dir, "output")
 dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
 
 # --- helpers ---
 site_from_file <- function(f) {
   # if your site_id is in the file itself you can ignore this.
   # otherwise, derive from filename prefix like "refer_Emory_..."
   bn <- basename(f)
   str_match(bn, "^refer_([^_]+)_")[,2] %||% "UnknownSite"
 }
 
 calc_se_from_ci <- function(lo, hi) (log(hi) - log(lo)) / (2 * 1.96)
 
 # Correct subtype-specific NO2 effects from a tidy model table (OR/IRR scale)
 compute_site_subtype_effects <- function(df, file) {
   
   site_id <- site_from_file(file)
   
   # main NO2 effect (ratio scale)
   no2_row <- df %>% filter(term == "no2_10") %>% slice(1)
   if (nrow(no2_row) == 0) return(NULL)
   
   beta_no2 <- log(no2_row$estimate)
   se_no2   <- calc_se_from_ci(no2_row$conf.low, no2_row$conf.high)
   
   inter <- df %>%
     filter(str_detect(term, "^no2_10:arf_subtype")) %>%
     mutate(
       arf_subtype = str_replace(term, "^no2_10:arf_subtype", ""),
       arf_subtype = str_to_title(arf_subtype),
       beta_int    = log(estimate),
       se_int      = calc_se_from_ci(conf.low, conf.high)
     )
   
   out <- tibble(
     site_id = site_id,
     arf_subtype = "Hypoxemic",
     beta = beta_no2,
     se   = se_no2
   )
   
   if (nrow(inter) > 0) {
     out <- bind_rows(
       out,
       inter %>%
         transmute(
           site_id,
           arf_subtype,
           beta = beta_no2 + beta_int,
           se   = sqrt(se_no2^2 + se_int^2)   # assumes 0 covariance
         )
     )
   }
   
   out
 }
 
 # Fixed-effect pooling on log scale
 pool_log_effects <- function(dat) {
   dat <- dat %>% filter(is.finite(beta), is.finite(se), se > 0)
   
   if (nrow(dat) == 0) return(tibble(est = NA_real_, lo = NA_real_, hi = NA_real_, k = 0L))
   
   w <- 1 / (dat$se^2)
   beta_hat <- sum(w * dat$beta) / sum(w)
   se_hat   <- sqrt(1 / sum(w))
   
   tibble(
     est = exp(beta_hat),
     lo  = exp(beta_hat - 1.96 * se_hat),
     hi  = exp(beta_hat + 1.96 * se_hat),
     k   = nrow(dat)
   )
 }
 
 # --- read + compute pooled effects for one outcome file pattern ---
 pool_one_outcome <- function(pattern, outcome_label) {
   
   files <- list.files(subtype_dir, pattern = pattern, full.names = TRUE)
   files <- files[str_detect(basename(files), "x_subtype")]  # safety
   
   if (length(files) == 0) return(NULL)
   
   per_site <- map_dfr(files, \(f) {
     df <- suppressMessages(read_csv(f, show_col_types = FALSE))
     compute_site_subtype_effects(df, f)
   })
   
   pooled <- per_site %>%
     mutate(
       arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed"))
     ) %>%
     group_by(arf_subtype) %>%
     group_modify(~ pool_log_effects(.x)) %>%
     ungroup() %>%
     mutate(outcome = outcome_label)
   
   pooled
 }
 
 # --- run pooling for all 3 outcomes (ONLY x_subtype files) ---
 pooled_all <- bind_rows(
   pool_one_outcome("model_tidy_death30d.*x_subtype.*\\.csv$",   "30-day mortality (OR)"),
   pool_one_outcome("model_tidy_inhosp_death.*x_subtype.*\\.csv$", "In-hospital mortality (OR)"),
   pool_one_outcome("model_tidy_vent_hours.*x_subtype.*\\.csv$", "Ventilation hours (IRR)")
 ) %>%
   filter(!is.na(est)) %>%
   mutate(
     outcome = factor(outcome, levels = c("30-day mortality (OR)",
                                          "In-hospital mortality (OR)",
                                          "Ventilation hours (IRR)"))
   )
 
 # --- publication-ready faceted forest plot ---
 p <- ggplot(pooled_all, aes(x = est, y = 1)) +
   geom_vline(xintercept = 1, linetype = 3, linewidth = 0.5) +
   geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.8) +
   geom_point(size = 2.6) +
   scale_x_log10() +
   facet_grid(outcome ~ arf_subtype, scales = "free_x") +
   labs(
     title = expression("Pooled NO"[2]*" effect per 10 ppb by ARF subtype"),
     x = expression("Effect estimate per 10 ppb NO"[2]*" (log scale)"),
     y = NULL
   ) +
   theme_classic(base_size = 13) +
   theme(
     strip.background = element_blank(),
     strip.text = element_text(face = "bold"),
     axis.text.y = element_blank(),
     axis.ticks.y = element_blank(),
     plot.title = element_text(face = "bold", size = 16),
     panel.spacing = unit(1.0, "lines")
   )
 
 ggsave(
   filename = file.path(out_dir, "pooled_no2_effect_by_subtype_all_outcomes_facet.png"),
   plot = p,
   width = 11.5, height = 7.0, dpi = 600
 )
 
 ggsave(
   filename = file.path(out_dir, "pooled_no2_effect_by_subtype_all_outcomes_facet.pdf"),
   plot = p,
   width = 11.5, height = 7.0
 )
 
 p



 
 # Ensure clean ordering
 pooled_plot <- pooled_all %>%
   mutate(
     arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed")),
     outcome = factor(outcome, levels = c(
       "30-day mortality (OR)",
       "In-hospital mortality (OR)",
       "Ventilation hours (IRR)"
     ))
   )
 
 # Global x limits (common axis across all panels)
 x_min <- min(pooled_plot$lo, na.rm = TRUE)
 x_max <- max(pooled_plot$hi, na.rm = TRUE)
 
 p2 <- ggplot(pooled_plot, aes(x = est, y = outcome, color = outcome)) +
   geom_vline(xintercept = 1, linetype = 3, linewidth = 0.5, color = "grey30") +
   geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.9) +
   geom_point(size = 2.6) +
   facet_wrap(~ arf_subtype, nrow = 1) +
   scale_x_log10(limits = c(x_min, x_max)) +
   labs(
     title = expression("Pooled NO"[2]*" effect per 10 ppb by ARF subtype"),
     x = expression("Effect estimate per 10 ppb NO"[2]*" (log scale)"),
     y = NULL,
     color = "Outcome"
   ) +
   theme_classic(base_size = 13) +
   theme(
     plot.title = element_text(face = "bold", size = 15),
     strip.background = element_blank(),
     strip.text = element_text(face = "bold"),
     legend.position = "bottom",
     legend.title = element_text(face = "bold"),
     panel.spacing = unit(1.0, "lines")
   )
 
 # Save high-res
 out_dir <- file.path("/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype", "output")
 dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
 
 ggsave(
   filename = file.path(out_dir, "pooled_no2_effect_by_subtype_all_outcomes_color_commonx.png"),
   plot = p2,
   width = 10.5, height = 4.2, dpi = 600
 )
 
 ggsave(
   filename = file.path(out_dir, "pooled_no2_effect_by_subtype_all_outcomes_color_commonx.pdf"),
   plot = p2,
   width = 10.5, height = 4.2
 )
 
 p2
 
 
 # pooled_all must contain: outcome, arf_subtype, est, lo, hi
 # outcome should already be labeled like:
 # "30-day mortality (OR)", "In-hospital mortality (OR)", "Ventilation hours (IRR)"
 plot_dat <- pooled_all %>%
   mutate(
     arf_subtype = as.character(arf_subtype),
     arf_subtype = case_when(
       str_detect(str_to_lower(arf_subtype), "hypox") ~ "Hypoxemic",
       str_detect(str_to_lower(arf_subtype), "hypercap") ~ "Hypercapnic",
       str_detect(str_to_lower(arf_subtype), "mixed") ~ "Mixed",
       TRUE ~ str_to_title(arf_subtype)
     ),
     arf_subtype = factor(arf_subtype, levels = c("Hypoxemic", "Hypercapnic", "Mixed")),
     outcome = factor(outcome, levels = c(
       "30-day mortality (OR)",
       "In-hospital mortality (OR)",
       "Ventilation hours (IRR)"
     ))
   )
 
 # Common x-limits across everything
 x_min <- min(plot_dat$lo, na.rm = TRUE)
 x_max <- max(plot_dat$hi, na.rm = TRUE)
 
 pd <- position_dodge(width = 0.55)
 
 p <- ggplot(plot_dat, aes(x = est, y = outcome, color = arf_subtype)) +
   geom_vline(xintercept = 1, linetype = 3, linewidth = 0.5, color = "grey35") +
   geom_errorbarh(aes(xmin = lo, xmax = hi), position = pd, height = 0, linewidth = 0.95) +
   geom_point(position = pd, size = 2.8) +
   scale_x_log10(limits = c(x_min, x_max)) +
   labs(
     title = expression("Pooled NO"[2]*" effect per 10 ppb by outcome and ARF subtype"),
     x = expression("Effect estimate per 10 ppb NO"[2]*" (log scale)"),
     y = NULL,
     color = "ARF subtype"
   ) +
   theme_classic(base_size = 13) +
   theme(
     plot.title = element_text(face = "bold", size = 15),
     legend.position = "bottom",
     legend.title = element_text(face = "bold"),
     axis.title.x = element_text(face = "bold"),
     axis.text.y = element_text(face = "bold"),
     panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
     panel.grid.minor = element_blank()
   )
 
 # Save high-res
 out_dir <- file.path(
   "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype",
   "output"
 )
 dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
 
 ggsave(
   filename = file.path(out_dir, "pooled_no2_effect_all_outcomes_grouped_by_outcome_color_subtype.png"),
   plot = p,
   width = 9.2, height = 4.6, dpi = 600
 )
 
 ggsave(
   filename = file.path(out_dir, "pooled_no2_effect_all_outcomes_grouped_by_outcome_color_subtype.pdf"),
   plot = p,
   width = 9.2, height = 4.6
 )
 
 p

 
 library(cowplot)
 out_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype/output"
 
 race_30d_png    <- file.path(out_dir, "race_30d.png")
 race_inhosp_png <- file.path(out_dir, "race_inhosp.png")
 race_vent_png   <- file.path(out_dir, "race_vent.png")
 
 sex_30d_png     <- file.path(out_dir, "sex_30d.png")
 sex_inhosp_png  <- file.path(out_dir, "sex_inhosp.png")
 sex_vent_png    <- file.path(out_dir, "sex_vent.png")
 
 # ---- READ PNGs AS GROBS ----
 g_race_30d    <- ggdraw() + draw_image(race_30d_png,    scale = 1)
 g_race_inhosp <- ggdraw() + draw_image(race_inhosp_png, scale = 1)
 g_race_vent   <- ggdraw() + draw_image(race_vent_png,   scale = 1)
 
 g_sex_30d     <- ggdraw() + draw_image(sex_30d_png,     scale = 1)
 g_sex_inhosp  <- ggdraw() + draw_image(sex_inhosp_png,  scale = 1)
 g_sex_vent    <- ggdraw() + draw_image(sex_vent_png,    scale = 1)
 
 # ---- BUILD GRID (Race column | Sex column) ----
 panel <- plot_grid(
   g_race_30d,    g_sex_30d,
   g_race_inhosp, g_sex_inhosp,
   g_race_vent,   g_sex_vent,
   ncol = 2,
   align = "hv",
   axis = "tblr",
   rel_widths = c(1, 1),
   rel_heights = c(1, 1, 1),
   labels = c("A", "B", "C", "D", "E", "F"),
   label_fontface = "bold",
   label_size = 16,
   label_x = 0.02,
   label_y = 0.98,
   hjust = 0,
   vjust = 1
 )
 
 # Optional: add column headers
 panel2 <- ggdraw(panel) +
   draw_label("Race/ethnicity", x = 0.25, y = 0.995, hjust = 0.5, vjust = 1,
              fontface = "bold", size = 14) +
   draw_label("Sex",            x = 0.75, y = 0.995, hjust = 0.5, vjust = 1,
              fontface = "bold", size = 14)
 
 # ---- SAVE HIGH-RES ----
 ggsave(
   filename = file.path(out_dir, "Figure_pooled_curves_race_vs_sex_A-F.png"),
   plot = panel2,
   width = 12.5, height = 14, dpi = 600, bg = "white"
 )
 
 ggsave(
   filename = file.path(out_dir, "Figure_pooled_curves_race_vs_sex_A-F.pdf"),
   plot = panel2,
   width = 12.5, height = 14, useDingbats = FALSE
 )
 
 panel2

 ###### covid sensitivity analysis
 
 
 subtype_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype"
 
 # Output folder
 out_dir <- file.path(subtype_dir, "output", "covid_sensitivity")
 dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
 
 # -------------------------
 # 1) Find ONLY the covid-adjusted tidy model files
 # -------------------------
 covid_files <- list.files(
   subtype_dir,
   pattern = "_adj_plus_covid_.*\\.csv$",
   full.names = TRUE
 )
 
 if (length(covid_files) == 0) stop("No *_adj_plus_covid_*.csv files found in subtype_dir.")
 
 # Keep only the outcomes you care about (include ICU LOS here)
 # NOTE: adjust patterns to match your filenames if needed
 outcome_map <- tibble::tribble(
   ~pattern,        ~outcome_key,   ~outcome_label,                 ~measure_label,
   "death30d",      "death30d",     "30-day mortality",             "OR",
   "inhosp_death",  "inhosp_death", "In-hospital mortality",        "OR",
   "vent_hours",    "vent_hours",   "Ventilation hours",            "IRR",
   "icu_los",       "icu_los",      "ICU length of stay",           "IRR"
 )
 
 tag_outcome <- function(path) {
   fn <- basename(path)
   hit <- outcome_map %>% filter(str_detect(fn, fixed(pattern)))
   if (nrow(hit) == 0) return(NULL)
   hit[1,] %>% mutate(file = path)
 }
 
 files_tagged <- covid_files %>%
   map(tag_outcome) %>%
   compact() %>%
   bind_rows()
 
 if (nrow(files_tagged) == 0) {
   stop("Found *_adj_plus_covid_* files, but none matched outcome patterns in outcome_map.")
 }
 
 # Attempt to infer site name from filename (edit if your naming differs)
 infer_site <- function(fn) {
   # e.g., refer_Emory_20251008_... => "Emory"
   x <- str_match(fn, "^refer_([^_]+)_")[,2]
   ifelse(is.na(x), "UnknownSite", x)
 }
 
 files_tagged <- files_tagged %>%
   mutate(
     filename = basename(file),
     site = infer_site(filename)
   )
 
 # -------------------------
 # 2) Read each CSV and extract NO2 term row (no2_10)
 # -------------------------
 read_tidy <- function(path) {
   # robust read for tab or comma
   dat <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
   if (!("term" %in% names(dat))) {
     dat <- suppressWarnings(readr::read_tsv(path, show_col_types = FALSE))
   }
   dat
 }
 
 extract_no2 <- function(path) {
   dat <- read_tidy(path)
   
   # standardize column names (in case of capitalization differences)
   names(dat) <- tolower(names(dat))
   
   needed <- c("term", "estimate", "conf.low", "conf.high")
   if (!all(needed %in% names(dat))) {
     stop(glue("Missing required columns in {basename(path)}. Need: {paste(needed, collapse=', ')}"))
   }
   
   row <- dat %>%
     filter(term == "no2_10") %>%
     slice(1)
   
   if (nrow(row) == 0) return(NULL)
   
   row %>%
     transmute(
       estimate   = as.numeric(estimate),
       conf_low   = as.numeric(conf.low),
       conf_high  = as.numeric(conf.high)
     )
 }
 
 no2_site_effects <- files_tagged %>%
   mutate(no2 = map(file, extract_no2)) %>%
   filter(!map_lgl(no2, is.null)) %>%
   unnest(no2)
 
 if (nrow(no2_site_effects) == 0) {
   stop("No files contained a 'no2_10' term row. Confirm term naming in the tidy outputs.")
 }
 
 # -------------------------
 # 3) Meta-analysis prep:
 # pool on log scale (log(OR) or log(IRR))
 # derive SE from CI if needed
 # -------------------------
 no2_site_effects <- no2_site_effects %>%
   mutate(
     yi = log(estimate),
     sei = (log(conf_high) - log(conf_low)) / (2 * 1.96)
   ) %>%
   filter(is.finite(yi), is.finite(sei), sei > 0)
 
 # -------------------------
 # 4) Random-effects pooling per outcome
 # -------------------------
 pool_one_outcome <- function(df) {
   # REML random-effects meta-analysis
   fit <- metafor::rma(yi = yi, sei = sei, method = "REML", data = df)
   
   tibble(
     k = nrow(df),
     pooled_log = as.numeric(fit$b[1]),
     pooled_se  = as.numeric(fit$se[1]),
     pooled_ci_lo_log = as.numeric(fit$ci.lb),
     pooled_ci_hi_log = as.numeric(fit$ci.ub),
     tau2 = as.numeric(fit$tau2),
     i2   = as.numeric(fit$I2),
     pval = as.numeric(fit$pval)
   )
 }
 
 pooled <- no2_site_effects %>%
   group_by(outcome_key, outcome_label, measure_label) %>%
   group_modify(~ pool_one_outcome(.x)) %>%
   ungroup() %>%
   mutate(
     pooled_est = exp(pooled_log),
     pooled_lo  = exp(pooled_ci_lo_log),
     pooled_hi  = exp(pooled_ci_hi_log)
   )
 
 # Save pooled table
 readr::write_csv(pooled, file.path(out_dir, "pooled_no2_adj_plus_covid_summary.csv"))
 
 # Also create a long table that includes site-specific + pooled rows (handy for plotting)
 plot_dat <- no2_site_effects %>%
   transmute(
     outcome_key, outcome_label, measure_label,
     site,
     est = estimate, lo = conf_low, hi = conf_high,
     type = "Site"
   ) %>%
   bind_rows(
     pooled %>%
       transmute(
         outcome_key, outcome_label, measure_label,
         site = "Pooled (REML)",
         est = pooled_est, lo = pooled_lo, hi = pooled_hi,
         type = "Pooled"
       )
   ) %>%
   mutate(
     site = factor(site, levels = c(sort(unique(no2_site_effects$site)), "Pooled (REML)")),
     outcome_label = factor(outcome_label, levels = outcome_map$outcome_label)
   )
 
 readr::write_csv(plot_dat, file.path(out_dir, "no2_adj_plus_covid_site_and_pooled.csv"))
 
 # -------------------------
 # 5) Publication-ready forest plot (facet by outcome)
 # -------------------------
 p_forest <- ggplot(plot_dat, aes(x = est, y = site)) +
   geom_vline(xintercept = 1, linetype = 3, linewidth = 0.4) +
   geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.6) +
   geom_point(aes(size = type), shape = 16) +
   scale_x_log10() +
   scale_size_manual(values = c("Site" = 2.2, "Pooled" = 3.2), guide = "none") +
   facet_grid(outcome_label ~ ., scales = "free_y", space = "free_y") +
   labs(
     title = expression("NO"[2] * " effect per 10 ppb (models adjusted for COVID period)"),
     x = expression("Effect estimate per 10 ppb NO"[2] * " (log scale)"),
     y = NULL
   ) +
   theme_classic(base_size = 13) +
   theme(
     plot.title = element_text(face = "bold"),
     strip.text.y = element_text(face = "bold"),
     axis.title.x = element_text(face = "bold")
   )
 
 ggsave(
   filename = file.path(out_dir, "forest_no2_adj_plus_covid_by_outcome.png"),
   plot = p_forest,
   width = 8.5, height = 9.0, units = "in", dpi = 600, bg = "white"
 )
 ggsave(
   filename = file.path(out_dir, "forest_no2_adj_plus_covid_by_outcome.pdf"),
   plot = p_forest,
   width = 8.5, height = 9.0, units = "in", device = cairo_pdf
 )
 
 # -------------------------
 # 6) A concise “results” table for manuscript text
 # -------------------------
 pooled_for_text <- pooled %>%
   transmute(
     outcome = outcome_label,
     measure = measure_label,
     k,
     est = pooled_est,
     lo = pooled_lo,
     hi = pooled_hi,
     p = pval,
     I2 = i2
   ) %>%
   mutate(
     effect_ci = sprintf("%.2f (%.2f–%.2f)", est, lo, hi),
     p_fmt = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
     I2_fmt = sprintf("%.1f%%", I2)
   ) %>%
   select(outcome, measure, k, effect_ci, p_fmt, I2_fmt)
 
 readr::write_csv(pooled_for_text, file.path(out_dir, "pooled_no2_adj_plus_covid_for_text.csv"))
 
 print(pooled_for_text)
 

##### pooled ARF subtypes

 # =========================
 # Pooled NO2 curves by ARF subtype
 # using site_preds_*_by_{sex|race} files that already contain arf_subtype
 # =========================
 

 # ---- REQUIRED: pool_curves() must already exist in your session ----
 # (you used it earlier; this script assumes it is defined)
 
 # -------------------------
 # Cleaning helpers
 # -------------------------
 clean_subtype <- function(x) {
   x_chr <- as.character(x)
   x_chr <- str_trim(x_chr)
   x_chr <- str_replace_all(x_chr, "\\s+", " ")
   x_chr <- str_to_lower(x_chr)
   
   dplyr::case_when(
     x_chr %in% c("hypoxemic", "hypoxaemic") ~ "Hypoxemic",
     x_chr %in% c("hypercapnic")            ~ "Hypercapnic",
     x_chr %in% c("mixed")                  ~ "Mixed",
     TRUE                                   ~ str_to_title(x_chr)
   )
 }
 
 # -------------------------
 # Pool/plot helpers
 # -------------------------
 prep_ribbon_df <- function(df, group_var, group_clean_fn = NULL, levels = NULL) {
   group_var <- rlang::ensym(group_var)
   group_nm  <- rlang::as_name(group_var)
   
   out <- df %>%
     mutate(
       no2_10 = as.numeric(.data$no2_10),
       pred   = as.numeric(.data$pred),
       lo     = as.numeric(.data$lo),
       hi     = as.numeric(.data$hi)
     ) %>%
     filter(is.finite(no2_10), is.finite(pred), is.finite(lo), is.finite(hi)) %>%
     mutate(
       lo = pmin(lo, hi),
       hi = pmax(lo, hi),
       pred = pmax(lo, pmin(pred, hi))
     )
   
   if (!is.null(group_clean_fn)) {
     out[[group_nm]] <- group_clean_fn(out[[group_nm]])
   }
   
   # one row per group × x for clean ribbons
   out <- out %>%
     group_by(.data[[group_nm]], no2_10) %>%
     summarise(
       pred = mean(pred, na.rm = TRUE),
       lo   = mean(lo,   na.rm = TRUE),
       hi   = mean(hi,   na.rm = TRUE),
       .groups = "drop"
     ) %>%
     arrange(.data[[group_nm]], no2_10)
   
   if (!is.null(levels)) {
     out[[group_nm]] <- factor(out[[group_nm]], levels = levels)
   }
   
  out
}

prep_rug_df <- function(files, group_var, group_clean_fn = NULL, levels = NULL) {
  group_var <- rlang::ensym(group_var)
  group_nm  <- rlang::as_name(group_var)
  
  out <- purrr::map_dfr(files, function(f) {
    readr::read_csv(f, show_col_types = FALSE) %>%
      dplyr::mutate(site_file = basename(f))
  }) %>%
    dplyr::mutate(no2_10 = as.numeric(.data$no2_10)) %>%
    dplyr::filter(is.finite(no2_10), !is.na(.data[[group_nm]]))
  
  if (!is.null(group_clean_fn)) {
    out[[group_nm]] <- group_clean_fn(out[[group_nm]])
  }
  
  if (!is.null(levels)) {
    out[[group_nm]] <- factor(out[[group_nm]], levels = levels)
  }
  
  out %>%
    dplyr::distinct(site_file, no2_10, .data[[group_nm]])
}

plot_ribbon_pub <- function(df, group_var, title_expr, ylab, palette = NULL, rug_df = NULL) {
  group_var <- rlang::ensym(group_var)
  group_nm  <- rlang::as_name(group_var)
  
  p <- ggplot(
     df,
     aes(
       x = no2_10,
       y = pred,
       color = !!group_var,
       fill  = !!group_var,
       group = !!group_var
     )
   ) +
     geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.08, color = NA) +
     geom_line(linewidth = 1.15) +
     labs(
       title = title_expr,
       x     = expression(NO[2] * " (per 10 ppb)"),
       y     = ylab,
       color = NULL,
       fill  = NULL
     ) +
     theme_classic(base_size = 13) +
     theme(
       plot.title       = element_text(face = "bold"),
       axis.title       = element_text(face = "bold"),
       legend.position  = "right",
       legend.key.width = unit(0.9, "lines")
     )
  
  if (!is.null(rug_df) && nrow(rug_df) > 0) {
    y_rng <- range(c(df$pred, df$lo, df$hi), na.rm = TRUE, finite = TRUE)
    y_span <- diff(y_rng)
    if (!is.finite(y_span) || y_span <= 0) y_span <- max(abs(y_rng), 1, na.rm = TRUE)
    rug_levels <- levels(df[[group_nm]])
    if (is.null(rug_levels)) rug_levels <- unique(as.character(df[[group_nm]]))
    
    rug_marks <- rug_df %>%
      dplyr::mutate(
        rug_row = match(as.character(.data[[group_nm]]), rug_levels),
        rug_y = y_rng[1] + y_span * (0.018 + 0.022 * (rug_row - 1)),
        rug_yend = rug_y + y_span * 0.018
      ) %>%
      dplyr::filter(!is.na(rug_row))
    
    p <- p +
      geom_segment(
        data = rug_marks,
        aes(x = no2_10, xend = no2_10, y = rug_y, yend = rug_yend, color = !!group_var),
        inherit.aes = FALSE,
        alpha = 0.45,
        linewidth = 0.35
      )
  }
   
   if (!is.null(palette)) {
     p <- p +
       scale_color_manual(values = palette) +
       scale_fill_manual(values  = palette)
   }
   
   p
 }
 
 # -------------------------
 # Configuration
 # -------------------------
 subtype_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype"
 
 # Choose which set of site_preds files to use as the source:
 # "sex" uses:  site_preds_{outcome}_by_sex*.csv
 # "race" uses: site_preds_{outcome}_by_race*.csv
 source_group <- "sex"  # <-- change to "race" if desired
 
 outcomes <- tibble::tribble(
   ~outcome_key, ~file_key, ~ylab,                              ~link,   ~title_stub,
   "30d",        "30d",     "Predicted probability",            "logit", "30-day mortality vs ",
   "inhosp",     "inhosp",  "Predicted probability",            "logit", "In-hospital mortality vs ",
   "vent",       "vent",    "Predicted mean\nventilation hours", "log",   "Ventilation hours vs "
 )
 
 # ARF subtype ordering + palette
 subtype_levels <- c("Hypoxemic", "Hypercapnic", "Mixed")
 pal_subtype <- brewer.pal(3, "Dark2")  # clean, colorblind-friendly-ish
 
 # -------------------------
 # Run: pooled curves by ARF subtype for each outcome
 # -------------------------
 plots <- list()
 
 for (o in seq_len(nrow(outcomes))) {
   o_key  <- outcomes$outcome_key[o]
   f_key  <- outcomes$file_key[o]
   ylab   <- outcomes$ylab[o]
   link   <- outcomes$link[o]
   t_stub <- outcomes$title_stub[o]
   
   files <- list.files(
     subtype_dir,
     pattern = paste0("^site_preds_", f_key, "_by_", source_group, ".*\\.csv$"),
     full.names = TRUE
   )
   
   if (length(files) == 0) {
     message("No files found for outcome=", o_key, " using by_", source_group,
             " (pattern: site_preds_", f_key, "_by_", source_group, "...)")
     next
   }
   
   pooled_raw <- pool_curves(
     files     = files,
     group_var = "arf_subtype",
     xvar      = "no2_10",
     link      = link
   )
   
  pooled <- prep_ribbon_df(
    pooled_raw,
    group_var = arf_subtype,
    group_clean_fn = clean_subtype,
    levels = subtype_levels
  )
  
  rug_df <- prep_rug_df(
    files,
    group_var = arf_subtype,
    group_clean_fn = clean_subtype,
    levels = subtype_levels
  )
  
  title_expr <- as.expression(bquote(.(t_stub) * NO[2] * ""))
  
  p <- plot_ribbon_pub(
    pooled,
    group_var  = arf_subtype,
    title_expr = title_expr,
    ylab       = ylab,
    palette    = pal_subtype,
    rug_df     = rug_df
  )
   
   plot_name <- paste0("subtype_", o_key, "_from_", source_group)
   plots[[plot_name]] <- p
   print(p)
 }
 
 # -------------------------
 # Save figures (high-res)
 # -------------------------
 output_dir <- file.path(subtype_dir, "output")
 dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
 
png_w <- 11.0
png_h <- 5.0
dpi   <- 600
 
 for (nm in names(plots)) {
   p <- plots[[nm]]
   
   ggsave(
     filename = file.path(output_dir, paste0(nm, ".png")),
     plot     = p,
     width    = png_w,
     height   = png_h,
     units    = "in",
     dpi      = dpi,
     bg       = "white"
   )
   
   ggsave(
     filename = file.path(output_dir, paste0(nm, ".pdf")),
     plot     = p,
     width    = png_w,
     height   = png_h,
     units    = "in",
     device   = cairo_pdf
   )
 }
 
 message("Done. Saved to: ", output_dir)
 
 
 # If you used source_group <- "sex" in the script:
 pA <- plots[["subtype_30d_from_sex"]]
 pB <- plots[["subtype_inhosp_from_sex"]]
 pC <- plots[["subtype_vent_from_sex"]]
 
 # If you used source_group <- "race" instead, swap to:
 # pA <- plots[["subtype_30d_from_race"]]
 # pB <- plots[["subtype_inhosp_from_race"]]
 # pC <- plots[["subtype_vent_from_race"]]
 
 # Build a single-row, 3-panel figure with A–C labels
 fig_ABC <- cowplot::plot_grid(
   pA, pB, pC,
   nrow = 3,
   labels = c("A", "B", "C"),
   label_size = 16,
   label_fontface = "bold",
   align = "hv",
   axis = "tb"
 )
 
 # Save high-res outputs
 output_dir <- file.path(subtype_dir, "output")
 dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
 
ggsave(
  filename = file.path(output_dir, "NO2_ARFsubtype_pooled_ABC.png"),
  plot     = fig_ABC,
  width    = 10.5, height = 9, units = "in",
  dpi      = 600,
  bg       = "white"
)

ggsave(
  filename = file.path(output_dir, "NO2_ARFsubtype_pooled_ABC.pdf"),
  plot     = fig_ABC,
  width    = 10.5, height = 9, units = "in",
  device   = cairo_pdf
)
 
 fig_ABC

 
 library(tidyverse)
 library(sf)
 library(tigris)
 library(scales)
 library(viridisLite)
 
 options(tigris_use_cache = TRUE)
 
 # ----------------------------
 # 1) Counties (CONUS + DC)
 # ----------------------------
 get_conus_counties <- function(year = 2020, remove_holes = TRUE) {
   drop_states <- c("02","15","72","60","66","69","78") # AK, HI, PR, territories
   
   x <- tigris::counties(cb = TRUE, year = year, class = "sf") |>
     dplyr::filter(!STATEFP %in% drop_states) |>
     dplyr::mutate(GEOID = as.character(GEOID)) |>
     sf::st_transform(4326) |>
     sf::st_make_valid()
   
   if (remove_holes) {
     # Rebuild polygons to drop interior rings (water bodies)
     x <- x |>
       dplyr::group_by(GEOID) |>
       dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
   }
   
   x
 }
 
 us_counties_ll <- get_conus_counties(year = 2020, remove_holes = TRUE)
 
 conus_xlim <- c(-125, -66.5)
 conus_ylim <- c(24, 49)
 
 # ----------------------------
 # 2) Summarize county mean
 # ----------------------------
 summarize_county_mean <- function(df, years = 2018:2024,
                                   value_candidates = c("value","no2_mean","pm25_mean")) {
   stopifnot(all(c("GEOID","year") %in% names(df)))
   
   df <- df |>
     mutate(
       GEOID = str_pad(as.character(GEOID), width = 5, side = "left", pad = "0"),
       year  = as.integer(year)
     )
   
   val_col <- value_candidates[value_candidates %in% names(df)][1]
   if (is.na(val_col)) stop("No value column found in: ", paste(names(df), collapse=", "))
   
   df |>
     filter(year %in% years) |>
     group_by(GEOID) |>
     summarise(
       mean_value = mean(.data[[val_col]], na.rm = TRUE),
       n_years    = sum(!is.na(.data[[val_col]])),
       .groups = "drop"
     )
 }
 
 # ----------------------------
 # 3) Bins + labels
 # ----------------------------
 make_bins <- function(x, probs = c(0.10, 0.25, 0.40, 0.55, 0.70, 0.82, 0.90, 0.95, 0.975, 0.99)) {
   qs <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
   unique(c(-Inf, qs, Inf))
 }
 
 fmt_range_labels <- function(brks, digits = 2) {
   b <- brks
   n <- length(b) - 1
   labs <- character(n)
   for (i in seq_len(n)) {
     lo <- b[i]; hi <- b[i+1]
     if (is.infinite(lo) && is.finite(hi)) {
       labs[i] <- paste0("≤", round(hi, digits))
     } else if (is.finite(lo) && is.infinite(hi)) {
       labs[i] <- paste0("≥", round(lo, digits))
     } else {
       labs[i] <- paste0(round(lo, digits), "–", round(hi, digits))
     }
   }
   labs
 }
 
 # ----------------------------
 # 4) ARF-style binned map with per-map palette
 # ----------------------------
 legend_expr_pm25 <- expression(bold(PM[2.5]~"(" * mu * "g/" * m^3 * ")"))
legend_expr_no2  <- expression(bold(NO[2]~"(ppb)"))
 
 make_pollutant_map_arf_style <- function(summary_df,
                                          title,
                                          legend_title,
                                          filename = NULL,
                                          palette_option = "magma",   # distinct per pollutant
                                          na_fill = "#bdbdbd",
                                          border_col = "white",
                                          border_lwd = 0.05,
                                          bins = NULL,
                                          digits = 2) {
   
   plot_df <- us_counties_ll |>
     left_join(summary_df, by = "GEOID")
   
   vals <- plot_df$mean_value
   if (is.null(bins)) bins <- make_bins(vals)
   labels <- fmt_range_labels(bins, digits = digits)
   
   plot_df <- plot_df |>
     mutate(
       bin = cut(mean_value, breaks = bins, include.lowest = TRUE, right = TRUE, labels = labels),
       bin = factor(bin, levels = labels)
     )
   
   pal <- viridisLite::viridis(length(labels), option = palette_option)
   names(pal) <- labels
   
   gg <- ggplot(plot_df) +
     geom_sf(aes(fill = bin), color = border_col, linewidth = border_lwd) +
     scale_fill_manual(
       name = legend_title,
       values = pal,
       drop = FALSE,
       na.value = na_fill,
       na.translate = FALSE,
       guide = guide_legend(reverse = TRUE)  # high at top
     ) +
     coord_sf(xlim = conus_xlim, ylim = conus_ylim, expand = FALSE) +
     labs(
       title = title,
       subtitle = "County-level average, 2018–2024"
     ) +
     theme_void(base_size = 12) +
     theme(
       plot.title = element_text(size = 20, face = "bold", hjust = 0),
       plot.subtitle = element_text(size = 14, hjust = 0, margin = margin(t = 4, b = 8)),
       legend.position = "right",
       legend.title = element_text(size = 18, face = "bold"),
       legend.text  = element_text(size = 12),
       plot.margin = margin(8, 8, 8, 8)
     )
   
   if (!is.null(filename)) {
     ggsave(filename, gg, width = 14, height = 9, dpi = 300)
   }
   gg
 }
 
 # ----------------------------
 # 5) Build summaries (your existing dfs)
 # ----------------------------
 no2_cnty_mean  <- summarize_county_mean(no2_us,  years = 2018:2024, value_candidates = c("no2_mean","value"))
 pm25_cnty_mean <- summarize_county_mean(pm25_us, years = 2018:2024, value_candidates = c("pm25_mean","value"))
 
 # ----------------------------
 # 6) Draw maps with distinct palettes
 # ----------------------------
 p_pm25 <- make_pollutant_map_arf_style(
   pm25_cnty_mean,
   title        = "Average fine particular matter (2.5) concentration estimate by U.S. County",
   legend_title = expression(PM[2.5]~"("*mu*"g/"*m^3*")"),
   palette_option = "magma",  # PM2.5 palette
   filename     = "map_pm25_conus_2018_2024_ARFstyle.png"
 )
 
 p_no2 <- make_pollutant_map_arf_style(
   no2_cnty_mean,
   title        = "Average nitrogen dioxide concentration estimate by U.S. County",
   legend_title = expression(NO[2]~"(ppb)"),
   palette_option = "cividis", # NO2 palette (distinct from PM2.5)
   filename     = "map_no2_conus_2018_2024_ARFstyle.png"
 )
 
 p_pm25
 p_no2
 
 ggsave(
   filename = file.path(dir_out, "map_pm25_conus_2018_2024.png"),
   plot     = p_pm25,
   width    = 14,
   height   = 9,
   dpi      = 600
 )
 
 ggsave(
   filename = file.path(dir_out, "map_no2_conus_2018_2024.png"),
   plot     = p_no2,
   width    = 14,
   height   = 9,
   dpi      = 600
 )


 library(patchwork)
 
 # Rename for clarity if needed

 
 p_3panel_stacked <- (p_arf / p_pm25 / p_no2) +
   plot_annotation(tag_levels = "A") &
   theme(
     plot.tag = element_text(face = "bold", size = 18),
     plot.tag.position = c(0.01, 0.99)  # top-left of each panel
   )
 
 # Save: same WIDTH as ARF, triple HEIGHT
 ggsave(
   filename = file.path(dir_out, "FIG_3panel_ARF_PM25_NO2_STACKED_ABC.png"),
   plot     = p_3panel_stacked,
   width    = 14,
   height   = 27,   # 3 × 9 inches
   dpi      = 300
 )
 
 ggsave(
   filename = file.path(dir_out, "FIG_3panel_ARF_PM25_NO2_STACKED_ABC.pdf"),
   plot     = p_3panel_stacked,
   width    = 14,
   height   = 27,
   device   = cairo_pdf
 )
 
 p_3panel_stacked
 
 library(cowplot)
 
 p_3panel_stacked <- plot_grid(
   p_arf, p_pm25, p_no2,
   ncol = 1,
   labels = c("A", "B", "C"),
   label_fontface = "bold",
   label_size = 38,
   label_x = 0.02,
   label_y = 0.98,
   hjust = 0,
   vjust = 1,
   align = "v",
   axis = "l"
 )
 
 p_3panel_stacked
 
 ggsave(
   filename = file.path(dir_out, "FIG_3panel_ARF_PM25_NO2_STACKED_ABC.png"),
   plot     = p_3panel_stacked,
   width    = 14,
   height   = 27,
   dpi      = 600
 )
 
 
 
 
 # tighter margins for ALL three plots (keeps ARF-like look, just less whitespace)
 tight_theme <- theme(
   plot.margin   = margin(2, 6, 2, 6),  # top, right, bottom, left (shrink these)
   plot.subtitle = element_text(margin = margin(t = 2, b = 4))
 )
 
 p_3panel_stacked_tight <-
   (p_arf / p_pm25 / p_no2) +
   plot_annotation(tag_levels = "A") &
   tight_theme &
   theme(
     plot.tag = element_text(face = "bold", size = 18),
     plot.tag.position = c(0.01, 0.99)
   )
 
 # Reduce the spacing between patchwork panels
 p_3panel_stacked_tight <- p_3panel_stacked_tight +
   plot_layout(ncol = 1, heights = c(1, 1, 1)) &
   theme(panel.spacing = unit(0.15, "lines"))  # smaller gap
 
 p_3panel_stacked_tight
 
 p_3panel_stacked_nogap <-
   (p_arf + theme(plot.margin = margin(0, 6, 0, 6))) /
   (p_pm25 + theme(plot.margin = margin(0, 6, 0, 6))) /
   (p_no2 + theme(plot.margin = margin(0, 6, 0, 6))) +
   plot_annotation(tag_levels = "A") &
   theme(plot.tag = element_text(face="bold", size=18),
         plot.tag.position = c(0.01, 0.99))
 
 p_3panel_stacked_nogap
 
 
 ggsave(
   filename = file.path(dir_out, "FIG_3panel_ARF_PM25_NO2_STACKED_ABC_tight.png"),
   plot     = p_3panel_stacked_nogap,
   width    = 14,
   height   = 27,
   dpi      = 600
 )
 
 
 library(tidyverse)
 library(scales)
 
 # --- Enter your table + effect-type mapping ---
 df <- tribble(
   ~Outcome,                          ~PM25,                     ~NO2,                      ~EffectType,
   "30-day mortality",                "1.06 (1.04–1.09)",        "1.11 (1.05–1.17)",         "OR",
   "In-hospital mortality",           "1.06 (1.04–1.08)",        "1.09 (1.04–1.14)",         "OR",
   "ICU length of stay",              "0.99 (0.98–1.01)",        "0.99 (0.96–1.01)",         "IRR",
   "Ventilation hours",               "1.02 (0.97–1.07)",        "1.01 (1.00–1.03)",         "IRR",
   "Persistent respiratory failure", "1.67 (0.71–3.97)",        "1.32 (0.87–1.99)",         "SHR",
   "Successful extubation",          "0.66 (0.40–1.09)",        "0.84 (0.69–1.03)",         "SHR",
   "Death",                          "2.31 (1.22–4.38)",        "1.50 (1.16–1.94)",         "SHR"
 )
 
 parse_est <- function(x){
   x <- str_replace_all(x, "–", "-")  # en-dash -> hyphen
   est <- as.numeric(str_extract(x, "^[0-9.]+"))
   lo  <- as.numeric(str_extract(x, "(?<=\\()[0-9.]+"))
   hi  <- as.numeric(str_extract(x, "(?<=-)[0-9.]+(?=\\))"))
   tibble(est = est, lo = lo, hi = hi)
 }
 
 # your df tribble as you created it above ...
 
 plot_df <- df %>%
   pivot_longer(cols = c(PM25, NO2), names_to = "Pollutant", values_to = "Estimate") %>%
   mutate(
     Pollutant = recode(Pollutant,
                        PM25 = "PM[2.5]",   # store as parseable text
                        NO2  = "NO[2]")     # store as parseable text
   ) %>%
   bind_cols(map_dfr(.$Estimate, parse_est)) %>%
   mutate(
     Outcome = factor(Outcome, levels = rev(unique(df$Outcome))),
     EffectType = factor(EffectType, levels = c("OR", "IRR", "SHR"))
   )
 
 shape_map <- c(OR = 16, IRR = 17, SHR = 15)
 color_map <- c(OR = "#1b9e77", IRR = "#d95f02", SHR = "#7570b3")
 
 xmin <- min(plot_df$lo, na.rm = TRUE)
 xmax <- max(plot_df$hi, na.rm = TRUE)
 
 p <- ggplot(plot_df, aes(x = est, y = Outcome, color = EffectType, shape = EffectType)) +
   geom_vline(xintercept = 1, linewidth = 0.4, linetype = "dashed", alpha = 0.7) +
   geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.6) +
   geom_point(size = 2.6, stroke = 0.4) +
   facet_wrap(~ Pollutant, ncol = 2, labeller = label_parsed) +
   scale_x_log10(
     limits = c(xmin * 0.95, xmax * 1.05),
     breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4),
     labels = label_number(accuracy = 0.01)
   ) +
   scale_shape_manual(values = shape_map, name = "Effect measure") +
   scale_color_manual(values = color_map, name = "Effect measure") +
   labs(
     title = "Associations of PM₂.₅ and NO₂ with Acute Respiratory Failure Outcomes",
     subtitle = "Point estimates with 95% confidence intervals (log scale; dashed line indicates null)",
     x = "Effect estimate (ratio scale)",
     y = NULL
   ) +
   theme_minimal(base_size = 12) +
   theme(
     panel.grid.major.y = element_blank(),
     panel.grid.minor = element_blank(),
     strip.text = element_text(face = "bold"),
     legend.position = "bottom",
     legend.title = element_text(face = "bold"),
     plot.title = element_text(face = "bold")
   )
 
 p
 ggsave("forest_pm25_no2_outcomes.png", p, width = 14, height = 9, dpi = 300)
 
 
 
 # ============================================================
 # Combine already-saved subtype/output panels into ONE 2x2 figure
 # Assumes you already have:
 #   sex_inhosp.(png/pdf), sex_vent.(png/pdf),
 #   race_inhosp.(png/pdf), race_vent.(png/pdf)
 # in: subtype_dir/output
 # ============================================================
 
 library(cowplot)
 library(ggplot2)
 
 subtype_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-ARFVI/sites/analysis/subtype"
 out_dir     <- file.path(subtype_dir, "output")
 
 # ---- file paths (match your naming convention) ----
 f_sex_inhosp  <- file.path(out_dir, "sex_inhosp.png")
 f_sex_vent    <- file.path(out_dir, "sex_vent.png")
 f_race_inhosp <- file.path(out_dir, "race_inhosp.png")
 f_race_vent   <- file.path(out_dir, "race_vent.png")
 
 # ---- basic existence checks ----
 stopifnot(file.exists(f_sex_inhosp))
 stopifnot(file.exists(f_sex_vent))
 stopifnot(file.exists(f_race_inhosp))
 stopifnot(file.exists(f_race_vent))
 
 # ---- read in as grobs ----
 g_sex_inhosp  <- cowplot::ggdraw() + cowplot::draw_image(f_sex_inhosp)
 g_sex_vent    <- cowplot::ggdraw() + cowplot::draw_image(f_sex_vent)
 g_race_inhosp <- cowplot::ggdraw() + cowplot::draw_image(f_race_inhosp)
 g_race_vent   <- cowplot::ggdraw() + cowplot::draw_image(f_race_vent)
 
 # Optional: add panel labels (A–D). Remove label_* args if you don't want.
 p_4panel <- cowplot::plot_grid(
   g_sex_inhosp, g_sex_vent,
   g_race_inhosp, g_race_vent,
   ncol = 2,
   align = "hv",
   axis = "tblr",
   labels = c("A", "B", "C", "D"),
   label_size = 18,
   label_fontface = "bold"
 )
 
 # ---- display ----
 print(p_4panel)
 
 # ---- save ----
 ggsave(
   filename = file.path(out_dir, "sex_race_inhosp_vent_4panel.png"),
   plot     = p_4panel,
   width    = 14.5,
   height   = 10.5,
   units    = "in",
   dpi      = 600,
   bg       = "white"
 )
 
 ggsave(
   filename = file.path(out_dir, "sex_race_inhosp_vent_4panel.pdf"),
   plot     = p_4panel,
   width    = 14.5,
   height   = 10.5,
   units    = "in",
   device   = cairo_pdf
 )
 
 
 
 
 
 
 
 
 
