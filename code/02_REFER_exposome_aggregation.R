

suppressPackageStartupMessages({
  library(terra)
  library(ncdf4)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(glue)
  library(scales)
  library(lubridate)
  library(sf)
  library(tigris)
  library(exactextractr)
})


# ---- Point this to your folder of .nc files ----
nc_dir <- "~/Desktop/Peter/Postdoc/CLIF"   # <-- change me

# ---- Naming helpers ([RESULT]_[SITE]_[SYSTEM_TIME]) ----
sanitize_tag <- function(x) {
  x <- if (is.null(x)) "SITE" else as.character(x)
  x <- iconv(x, to = "ASCII//TRANSLIT"); x <- gsub("[^A-Za-z0-9]+","_",x); x <- gsub("^_+|_+$","",x)
  if (!nzchar(x)) "SITE" else x
}
SITE_NAME   <- sanitize_tag(if (exists("config")) config$site_name else Sys.getenv("SITE_NAME","SITE"))
SYSTEM_TIME <- format(Sys.time(), "%Y%m%dT%H%M%S")
make_name   <- function(result, ext=NULL) { nm <- paste(result, SITE_NAME, SYSTEM_TIME, sep="_"); if (is.null(ext)) nm else paste0(nm,".",ext) }
out_dir <- file.path("output", paste0("no2_preview_", SITE_NAME, "_", SYSTEM_TIME))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Find files & years ----
pat <- "^annual_mean_tropomi_lur_conus_surface_no2_(\\d{4})\\.v\\d+\\.\\d+\\.nc$"
nc_files <- list.files(nc_dir, pattern = pat, full.names = TRUE)
stopifnot(length(nc_files) > 0)

years <- as.integer(str_match(basename(nc_files), pat)[,2])
ord <- order(years)
nc_files <- nc_files[ord]; years <- years[ord]

# ---- Peek NetCDF variable/units (first file) ----
peek_nc <- function(path){
  nc <- nc_open(path)
  on.exit(nc_close(nc), add = TRUE)
  vars <- names(nc$var)
  dims <- names(nc$dim)
  cat("Variables:", paste(vars, collapse=", "), "\n")
  cat("Dims:", paste(dims, collapse=", "), "\n")
  for (v in vars) {
    u <- nc$var[[v]]$units %||% ""
    ln <- nc$var[[v]]$longname %||% ""
    cat(sprintf("  - %s  units: %s  longname: %s\n", v, u, ln))
  }
}
peek_nc(nc_files[1])

# ---- Read one file with terra & print metadata ----
r1 <- rast(nc_files[1])  # terra picks the main data variable
names(r1) <- names(r1) %||% "no2"

cat("\nCRS:", crs(r1), "\n")
print(ext(r1)); print(res(r1)); print(dim(r1))
print(global(r1, fun = c("min","max","mean","sd"), na.rm = TRUE))

# ---- Quicklook map (downsampled) ----
r1_q <- aggregate(r1, fact = 4)  # coarser for speed (adjust as you like)
png_path <- file.path(out_dir, make_name(glue("no2_map_{years[1]}"), "png"))
png(png_path, width=1400, height=900, res=150)
plot(r1_q, main = glue("Annual mean NO2 ({years[1]})"))
dev.off(); message("Saved: ", png_path)

# ---- Stack all years & basic summaries ----
s <- rast(nc_files)           # lazy stack
names(s) <- paste0("y", years)

stats <- lapply(seq_along(years), function(i){
  g <- global(s[[i]], fun = c("min","max","mean","sd"), na.rm = TRUE)
  tibble(year = years[i],
         min = g[1,1], max = g[1,2], mean = g[1,3], sd = g[1,4])
}) %>% bind_rows()

stats_path <- file.path(out_dir, make_name("no2_yearly_stats","csv"))
write_csv(stats, stats_path); message("Saved: ", stats_path)

# ---- Tiny trend preview (fast, non-spatial): annual means line ----
p_line <- ggplot(stats, aes(year, mean)) + geom_line() + geom_point() +
  labs(title="CONUS annual mean NO2 (file-wide mean)", x=NULL, y="NO2 (mean of grid cells)") +
  theme_minimal(base_size = 12)
line_png <- file.path(out_dir, make_name("no2_yearly_means","png"))
ggsave(line_png, p_line, width=7, height=4, dpi=150); message("Saved: ", line_png)

# ---- Optional: pixel-wise trend (slope) map (can be heavy; try on coarse aggregate) ----
# Start coarse to be safe:
s_coarse <- aggregate(s, fact = 4)  # increase fact if needed
yrs_num  <- years
trend_fun <- function(v) {
  if (all(is.na(v))) return(NA_real_)
  # simple linear slope per pixel
  coef(lm(v ~ yrs_num))[2]
}
slope <- app(s_coarse, trend_fun)
slope_png <- file.path(out_dir, make_name("no2_trend_slope_coarse","png"))
png(slope_png, width=1400, height=900, res=150)
plot(slope, main="NO2 trend slope (value/year), coarse grid")
dev.off(); message("Saved: ", slope_png)

options(tigris_use_cache = TRUE)

# ---------------------------
# 1) Get CONUS county shapes
# ---------------------------
# Exclude non-CONUS by STATEFP (FIPS): AK=02, HI=15, AS=60, GU=66, MP=69, PR=72, VI=78
non_conus_statefp <- c("02","15","60","66","69","72","78")

# Counties (cartographic boundary = generalized)
counties_sf <- tigris::counties(cb = TRUE, year = 2020, class = "sf") %>%
  filter(!STATEFP %in% non_conus_statefp) %>%
  select(GEOID, NAME, STATEFP, geometry)

# (Optional) Add STUSPS if available, but don’t rely on it
states_min <- tigris::states(cb = TRUE, year = 2020, class = "sf") %>%
  st_drop_geometry() %>%
  select(STATEFP, STUSPS) %>%
  distinct()

if ("STATEFP" %in% names(states_min) && "STUSPS" %in% names(states_min)) {
  counties_sf <- counties_sf %>%
    left_join(states_min, by = "STATEFP")
}

# Reproject to match your raster stack `s`
counties_sf <- st_transform(counties_sf, crs = crs(s))

# terra vectors for optional terra::extract fallback
counties_v <- vect(counties_sf)

# ------------------------------------------
# 2) County means for each NO2 raster layer
# ------------------------------------------
use_exact <- requireNamespace("exactextractr", quietly = TRUE)

if (use_exact) {
  # exactextractr does area-weighted means and handles partial overlaps
  message("Using exactextractr (area-weighted means) ...")
  # Iterate over layers to keep memory use moderate
  res_list <- lapply(seq_len(nlyr(s)), function(i) {
    vals <- exactextractr::exact_extract(s[[i]], counties_sf, "mean")  # coverage-fraction weighted
    tibble(GEOID = counties_sf$GEOID, year = years[i], no2_mean = vals)
  })
  no2_county_year <- bind_rows(res_list)
} else {
  message("exactextractr not installed; using terra::extract (set exact=TRUE if available).")
  # terra::extract can compute means; exact=TRUE gives partial cell weighting (terra >= 1.7)
  ex <- terra::extract(s, counties_v, fun = mean, na.rm = TRUE, exact = TRUE)
  # ex has columns ID + one per layer; bring GEOID and pivot longer
  no2_county_year <- ex %>%
    as_tibble() %>%
    mutate(GEOID = counties_sf$GEOID[ID]) %>%
    select(-ID) %>%
    pivot_longer(cols = starts_with("y"), names_to = "band", values_to = "no2_mean") %>%
    mutate(year = as.integer(str_remove(band, "^y"))) %>%
    select(GEOID, year, no2_mean)
}

# Optional sanity: how many counties per year have values?
no2_county_year %>% count(year) %>% print(n=Inf)



### =================================================== PM 2.5







# ---------- 0) Naming helpers (reuse if you already defined them) ----------
if (!exists("make_name")) {
  sanitize_tag <- function(x) {
    x <- if (is.null(x)) "SITE" else as.character(x)
    x <- iconv(x, to = "ASCII//TRANSLIT")
    x <- gsub("[^A-Za-z0-9]+","_",x)
    x <- gsub("^_+|_+$","",x)
    if (!nzchar(x)) "SITE" else x
  }
  SITE_NAME   <- sanitize_tag(if (exists("config")) config$site_name else Sys.getenv("SITE_NAME","SITE"))
  SYSTEM_TIME <- format(Sys.time(), "%Y%m%dT%H%M%S")
  make_name   <- function(result, ext=NULL) {
    base <- paste(result, SITE_NAME, SYSTEM_TIME, sep="_")
    if (is.null(ext)) base else paste0(base, ".", ext)
  }
}
out_dir_pm25 <- file.path("output", paste0("pm25_preview_", SITE_NAME, "_", SYSTEM_TIME))
if (!dir.exists(out_dir_pm25)) dir.create(out_dir_pm25, recursive = TRUE)

# ---------- 1) Find files & parse years ----------
pm_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF"  # <-- set your folder
pat_pm <- "^V6GL\\d+\\.\\d+\\.CNNPM25\\.NA\\.(\\d{6})-(\\d{6})\\.nc$"
pm_files_all <- list.files(pm_dir, pattern = pat_pm, full.names = TRUE)
stopifnot(length(pm_files_all) > 0)

m <- str_match(basename(pm_files_all), pat_pm)
start_ym <- m[,2]  # e.g., "201801"
years_pm <- as.integer(substr(start_ym, 1, 4))
ord <- order(years_pm)
pm_files <- pm_files_all[ord]
years_pm <- years_pm[ord]

# ---------- 2) Inspect first file (optional peek) ----------
pm_peek <- function(path){
  nc <- nc_open(path); on.exit(nc_close(nc), add=TRUE)
  cat("Vars:", paste(names(nc$var), collapse=", "), "\n")
  for (v in names(nc$var)) {
    u <- nc$var[[v]]$units %||% ""
    ln <- nc$var[[v]]$longname %||% ""
    cat(sprintf("  - %s  units: %s  longname: %s\n", v, u, ln))
  }
}
pm_peek(pm_files[1])

# ---------- 3) Build annual-mean rasters (memory-friendly) ----------
annual_list <- vector("list", length(pm_files))
for (i in seq_along(pm_files)) {
  r <- rast(pm_files[i])             # may have 12 monthly layers
  # If multiple variables exist, terra picks one; to force, use rast(pm_files[i], subds="PM25")
  ann <- mean(r, na.rm = TRUE)       # annual mean across layers in the file
  names(ann) <- paste0("y", years_pm[i])
  annual_list[[i]] <- ann
}
pm_stack <- rast(annual_list)         # one layer per year
names(pm_stack) <- paste0("y", years_pm)

# ---------- 4) Quick stats + a preview map ----------
pm_stats <- lapply(seq_along(years_pm), function(i) {
  g <- global(pm_stack[[i]], fun=c("min","max","mean","sd"), na.rm=TRUE)
  tibble(year = years_pm[i], min=g[1,1], max=g[1,2], mean=g[1,3], sd=g[1,4])
}) %>% bind_rows()

stats_path <- file.path(out_dir_pm25, make_name("pm25_yearly_stats","csv"))
write_csv(pm_stats, stats_path); message("Saved: ", stats_path)

pm1 <- aggregate(pm_stack[[1]], fact = 4)  # coarse preview
png_path <- file.path(out_dir_pm25, make_name(glue("pm25_map_{years_pm[1]}"), "png"))
png(png_path, width=1400, height=900, res=150)
plot(pm1, main = glue("Annual mean PM2.5 ({years_pm[1]})"))
dev.off(); message("Saved: ", png_path)

# ---------- 5) County polygons (CONUS) ----------
options(tigris_use_cache = TRUE)
non_conus_statefp <- c("02","15","60","66","69","72","78")  # AK, HI, AS, GU, MP, PR, VI
counties_sf <- tigris::counties(cb=TRUE, year=2020, class="sf") %>%
  filter(!STATEFP %in% non_conus_statefp) %>%
  select(GEOID, NAME, STATEFP, geometry) %>%
  st_transform(crs = crs(pm_stack))
counties_v <- vect(counties_sf)

# ---------- 6) County means per year ----------
use_exact <- requireNamespace("exactextractr", quietly = TRUE)

if (use_exact) {
  message("Using exactextractr (area-weighted means) ...")
  res_list <- lapply(seq_len(nlyr(pm_stack)), function(i) {
    vals <- exactextractr::exact_extract(pm_stack[[i]], counties_sf, "mean")
    tibble(GEOID = counties_sf$GEOID, year = years_pm[i], pm25_mean = vals)
  })
  pm25_county_year <- bind_rows(res_list)
} else {
  message("Using terra::extract (exact=TRUE if supported) ...")
  ex <- terra::extract(pm_stack, counties_v, fun=mean, na.rm=TRUE, exact=TRUE)
  pm25_county_year <- ex %>%
    as_tibble() %>%
    mutate(GEOID = counties_sf$GEOID[ID]) %>%
    select(-ID) %>%
    pivot_longer(cols = starts_with("y"), names_to="band", values_to="pm25_mean") %>%
    mutate(year = as.integer(str_remove(band, "^y"))) %>%
    select(GEOID, year, pm25_mean)
}

pm_csv <- file.path(out_dir_pm25, make_name("pm25_county_year","csv"))
write_csv(pm25_county_year, pm_csv); message("Saved county PM2.5: ", pm_csv)



# -----------------------------------
# 1) Pixel-wise slope of PM2.5
# -----------------------------------

# Trend resolution (use a coarser preview too)
SLOPE_AGG_FACT <- 4   # increase to 8/16 if memory is tight for the preview

# -----------------------------------
# 1) Ensure we have years & PM stack
# -----------------------------------
stopifnot(exists("pm_stack") && inherits(pm_stack, "SpatRaster"))

if (!exists("years_pm")) {
  years_pm <- as.integer(gsub("^y", "", names(pm_stack)))
}

# -----------------------------------
# 2) Build CONUS mask in the raster CRS
# -----------------------------------
options(tigris_use_cache = TRUE)
non_conus_statefp <- c("02","15","60","66","69","72","78")  # AK, HI, AS, GU, MP, PR, VI

counties_sf <- tigris::counties(cb = TRUE, year = 2020, class = "sf") |>
  dplyr::filter(!STATEFP %in% non_conus_statefp) |>
  dplyr::select(GEOID, geometry) |>
  st_transform(crs = crs(pm_stack))              # match raster CRS

# Dissolve to a single CONUS polygon (fast union)
counties_v <- vect(counties_sf)
conus_v    <- terra::aggregate(counties_v)       # dissolve counties into one multipart polygon

# -----------------------------------
# 3) Crop & mask the PM stack to CONUS
# -----------------------------------
pm_conus <- terra::crop(pm_stack, terra::ext(conus_v), snap = "out")
pm_conus <- terra::mask(pm_conus, conus_v)

# -----------------------------------
# 4) Pixel-wise trend slope (units / year)
#     slope = cov(x,y)/var(x), x = years, y = values
# -----------------------------------
yrs <- years_pm
slope_fun <- function(v) {
  idx <- !is.na(v)
  if (sum(idx) < 2) return(NA_real_)
  x <- yrs[idx]; y <- v[idx]
  x <- x - mean(x); y <- y - mean(y)
  sum(x * y) / sum(x * x)
}

# (a) Full-resolution slope (may take a bit—terra streams by chunks)
pm25_slope_conus <- terra::app(pm_conus, slope_fun)

# Save full-res GeoTIFF
slope_tif <- file.path(out_dir_pm25, make_name("pm25_trend_slope_conus", "tif"))
terra::writeRaster(pm25_slope_conus, slope_tif, overwrite = TRUE)
message("Saved full-res slope GeoTIFF: ", slope_tif)

# -----------------------------------
# 5) Quicklook PNG (coarse) for easy viewing
# -----------------------------------
pm_conus_coarse <- terra::aggregate(pm_conus, fact = SLOPE_AGG_FACT)
slope_coarse    <- terra::app(pm_conus_coarse, slope_fun)

# Choose a symmetric range around 0 for a diverging palette
rng <- quantile(values(slope_coarse), probs = c(0.02, 0.98), na.rm = TRUE)
mx  <- max(abs(rng), na.rm = TRUE); brks <- c(-mx, mx)

png_path <- file.path(out_dir_pm25, make_name("pm25_trend_slope_conus_quicklook", "png"))
png(png_path, width = 1400, height = 900, res = 150)
plot(slope_coarse,
     main = glue("PM2.5 trend slope (units/year), CONUS\nYears: {min(yrs)}–{max(yrs)}  (coarse preview)"),
     col = hcl.colors(100, "Blue-Red 3"),
     type = "continuous",
     range = brks)
dev.off()
message("Saved slope quicklook PNG: ", png_path)
