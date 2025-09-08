# ================================================================================================
# ICU REspiratory Failure Environmental Risk (REFER) Index | PI: Peter Graffy (graffy@uchicago.edu)
# Exposome Aggregation and Linkage to Cohort
# Years: 2018–2023 (PM2.5) and 2019-2024 (NO2)
# Outputs:
#   pm/no2_county_year       : (2) aggregated air pollutants by county geoID
#   pm/no2_map.              : (2) single year maps of air pollutant
#   pm/no2_ trend_slope      : (2) pixel-wise trends slopes of air pollutant
#   pm/no2_yearly_stats      : (2) min/max/mean/sd summary stats of air pollutant
# =================================================================================================

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
  library(daymetr)
  library(fs)
})

# ---- Point this to your folder of .nc files ----
nc_dir <- "~/Desktop/Peter/Postdoc/CLIF"   ################## <-- change me

# ---- Naming helpers ([RESULT]_[SITE]_[SYSTEM_TIME]) ----
sanitize_tag <- function(x) {
  x <- if (is.null(x)) "SITE" else as.character(x)
  x <- iconv(x, to = "ASCII//TRANSLIT"); x <- gsub("[^A-Za-z0-9]+","_",x); x <- gsub("^_+|_+$","",x)
  if (!nzchar(x)) "SITE" else x
}
SITE_NAME   <- sanitize_tag(if (exists("config")) config$site_name else Sys.getenv("SITE_NAME","SITE"))
SYSTEM_TIME <- format(Sys.Date(), "%Y%m%d")
make_name   <- function(result, ext=NULL) { nm <- paste(result, SITE_NAME, SYSTEM_TIME, sep="_"); if (is.null(ext)) nm else paste0(nm,".",ext) }
out_dir <- file.path("output", paste0("no2_preview_", SITE_NAME, "_", SYSTEM_TIME))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ========================================
# NO2
# ========================================

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

# infer years if not already defined
if (!exists("years")) {
  # if you still have the NO2 yearly stack `s`, pull from its layer names like "y2019"
  if (exists("s") && inherits(s, "SpatRaster")) {
    years <- as.integer(gsub("^y", "", names(s)))
  } else {
    years <- NA_integer_
  }
}

# symmetric range around 0 using central 2%–98% quantiles (like the others)
q  <- quantile(values(slope), probs = c(0.02, 0.98), na.rm = TRUE)
mx <- max(abs(q), na.rm = TRUE)
rng <- c(-mx, mx)

# title (uses NO₂; switch to "NO2" if your graphics device dislikes Unicode)
ttl <- if (all(!is.na(years))) {
  glue("NO2 trend slope (units/year), CONUS\nYears: {min(years)}–{max(years)} (coarse preview)")
} else {
  "NO2 trend slope (units/year), CONUS (coarse preview)"
}

# choose output path (falls back to a temp file if naming helpers aren't defined)
outfile <- if (exists("make_name") && exists("out_dir")) {
  file.path(out_dir, make_name("no2_trend_slope_conus_quicklook", "png"))
} else {
  "no2_trend_slope_conus_quicklook.png"
}

png(outfile, width = 1400, height = 900, res = 150)
plot(slope,
     main = ttl,
     col  = hcl.colors(100, "Blue-Red 3"),
     type = "continuous",
     range = rng)
dev.off()
message("Saved: ", outfile)


# ---------------------------
# 1) Get CONUS county shapes
# ---------------------------

options(tigris_use_cache = TRUE)
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

no2_csv <- file.path(out_dir, make_name("no2_county_year","csv"))
write_csv(no2_county_year, no2_csv); message("Saved county NO2: ", no2_csv)

# ========================================
# PM2.5
# ========================================

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
  SYSTEM_TIME <- format(Sys.Date(), "%Y%m%d")
  make_name   <- function(result, ext=NULL) {
    base <- paste(result, SITE_NAME, SYSTEM_TIME, sep="_")
    if (is.null(ext)) base else paste0(base, ".", ext)
  }
}
out_dir_pm25 <- file.path("output", paste0("pm25_preview_", SITE_NAME, "_", SYSTEM_TIME))
if (!dir.exists(out_dir_pm25)) dir.create(out_dir_pm25, recursive = TRUE)

# ---------- 1) Find files & parse years ----------
pm_dir <- "~/Desktop/Peter/Postdoc/CLIF"   ################## <-- make sure it's pointed to the same directory that had NO2

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

# =======================================================
# Daymet (tmax, tmin, vp, and prcp) from 2018-2024
# =======================================================

####################################################
#### IF YOU WANT TO DOWNLOAD THE DATA YOURSELF #####
####################################################

# # --- Naming helpers (reuse yours if already defined) ---
# sanitize_tag <- function(x) {
#   x <- if (is.null(x)) "SITE" else as.character(x)
#   x <- iconv(x, to = "ASCII//TRANSLIT"); x <- gsub("[^A-Za-z0-9]+","_",x); x <- gsub("^_+|_+$","",x)
#   if (!nzchar(x)) "SITE" else x
# }
# SITE_NAME   <- if (exists("SITE_NAME")) SITE_NAME else sanitize_tag(if (exists("config")) config$site_name else Sys.getenv("SITE_NAME","SITE"))
# SYSTEM_TIME <- if (exists("SYSTEM_TIME")) SYSTEM_TIME else format(Sys.time(), "%Y%m%dT%H%M%S")
# make_name   <- function(result, ext=NULL){ nm <- paste(result, SITE_NAME, SYSTEM_TIME, sep="_"); if (is.null(ext)) nm else paste0(nm,".",ext) }
# 
# # --- Where to put things ---
# base_dir <- file.path("env_data", paste0("daymet_conus_", SITE_NAME, "_", SYSTEM_TIME))
# dir_create(base_dir)
# 
# # --- What to download ---
# years   <- 2018:2024
# params  <- c("tmax","tmin","vp","prcp")  # Daymet variable keys
# # CONUS bbox (top-left lat, lon; bottom-right lat, lon) for NCSS:
# bbox_conus <- c(49.5, -125.0, 24.0, -66.0)
# 
# # --- Helper: safe NCSS download for 1 param-year, with rename ---
# dl_ncss_one <- function(param, yr){
#   message(glue("NCSS: {param} {yr} ..."))
#   nc_path <- tryCatch(
#     daymetr::download_daymet_ncss(
#       location  = bbox_conus,
#       start     = yr, end = yr,
#       param     = param,
#       frequency = "annual",      # annual keeps file sizes reasonable
#       mosaic    = "na",          # North America mosaic
#       path      = base_dir,
#       silent    = FALSE,
#       force     = TRUE           # allow newest full year
#     ),
#     error = function(e) e
#   )
#   if (inherits(nc_path, "error")) return(nc_path)  # signal failure
#   
#   # Rename to your convention
#   out_name <- make_name(glue("daymet_{param}_annual_CONUS_{yr}"), "nc")
#   out_path <- file.path(base_dir, out_name)
#   file_move(nc_path, out_path)
#   message("  -> saved: ", out_path)
#   out_path
# }
# 
# # --- Helper: tile fallback (daily tiles) if NCSS fails for size/other reasons ---
# # NOTE: This downloads daily tile NetCDFs (bigger volume). Consider running later
# # a reducer to produce annual summaries from the tiles if you hit this path.
# dl_tiles_one <- function(param, yr){
#   message(glue("Tiles (fallback): {param} {yr} ..."))
#   daymetr::download_daymet_tiles(
#     location = bbox_conus,  # selects all intersecting tiles
#     start    = yr, end = yr,
#     param    = param,
#     path     = file.path(base_dir, "tiles"),
#     silent   = FALSE,
#     force    = TRUE
#   )
#   # Files keep Daymet’s original names; we keep them in base_dir/tiles/
#   message(glue("  -> downloaded daily tiles to {file.path(base_dir,'tiles')}"))
#   invisible(TRUE)
# }
# 
# # --- Drive the downloads ---
# for (p in params) {
#   for (y in years) {
#     res <- dl_ncss_one(p, y)
#     if (inherits(res, "error")) {
#       message("  NCSS failed: ", conditionMessage(res), " — switching to tiles …")
#       dl_tiles_one(p, y)
#     }
#   }
# }
# 
# message("\nDone. Annual NC files (when NCSS succeeded) live in:\n  ", base_dir,
#         "\nTile fallbacks (daily) live in:\n  ", file.path(base_dir, "tiles"))

###############################################
#### IF YOU WANT TO USE PREDOWNLOADED DATA ####
###############################################

# ---------- 0) Where the NetCDFs live ----------
dm_dir <- "~/Desktop/Peter/Postdoc/CLIF"   # <-- set this folder

# ---------- 1) Naming helpers ----------
sanitize_tag <- function(x) {
  x <- if (is.null(x)) "SITE" else as.character(x)
  x <- iconv(x, to = "ASCII//TRANSLIT"); x <- gsub("[^A-Za-z0-9]+","_",x); x <- gsub("^_+|_+$","",x)
  if (!nzchar(x)) "SITE" else x
}
SITE_NAME   <- if (exists("SITE_NAME")) SITE_NAME else sanitize_tag(if (exists("config")) config$site_name else Sys.getenv("SITE_NAME","SITE"))
SYSTEM_TIME <- if (exists("SYSTEM_TIME")) SYSTEM_TIME else format(Sys.time(), "%Y%m%dT%H%M%S")
make_name   <- function(result, ext=NULL){ nm <- paste(result, SITE_NAME, SYSTEM_TIME, sep="_"); if (is.null(ext)) nm else paste0(nm,".",ext) }

out_dir_daymet <- file.path("output", paste0("daymet_conus_", SITE_NAME, "_", SYSTEM_TIME))
if (!dir.exists(out_dir_daymet)) dir.create(out_dir_daymet, recursive = TRUE)

# ---------- 2) Build CONUS mask (once) ----------
options(tigris_use_cache = TRUE)
non_conus_statefp <- c("02","15","60","66","69","72","78")  # AK, HI, AS, GU, MP, PR, VI

# We'll project to raster CRS per-variable later; keep in WGS84 for now
counties_sf_wgs <- tigris::counties(cb = TRUE, year = 2020, class = "sf") |>
  dplyr::filter(!STATEFP %in% non_conus_statefp) |>
  dplyr::select(GEOID, geometry)

# ---------- 3) File discovery helpers ----------
find_files <- function(var) {
  # annavg for tmax/tmin/vp; annttl for prcp
  if (var == "prcp") {
    pat <- "^daymet_v4_prcp_annttl_na_(\\d{4})\\.nc$"
  } else {
    pat <- paste0("^daymet_v4_", var, "_annavg_na_(\\d{4})\\.nc$")
  }
  fs <- list.files(dm_dir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
  stopifnot(length(fs) > 0)
  yrs <- as.integer(str_match(basename(fs), pat)[,2])
  o <- order(yrs)
  list(files = fs[o], years = yrs[o])
}

# ---------- 4) County aggregation helper ----------
get_county_means <- function(rstack, years_vec, counties_sf) {
  use_exact <- requireNamespace("exactextractr", quietly = TRUE)
  if (use_exact) {
    # area-weighted mean per county per year
    res_list <- lapply(seq_len(nlyr(rstack)), function(i) {
      vals <- exactextractr::exact_extract(rstack[[i]], counties_sf, "mean")
      tibble(GEOID = counties_sf$GEOID, year = years_vec[i], value = vals)
    })
    bind_rows(res_list)
  } else {
    # terra::extract exact=TRUE if available in your terra version
    counties_v <- vect(counties_sf)
    ex <- terra::extract(rstack, counties_v, fun = mean, na.rm = TRUE, exact = TRUE)
    as_tibble(ex) |>
      mutate(GEOID = counties_sf$GEOID[ID]) |>
      select(-ID) |>
      mutate(.band = names(rstack)) |>
      tidyr::pivot_longer(cols = starts_with("y"), names_to = "band", values_to = "value") |>
      mutate(year = as.integer(str_remove(band, "^y"))) |>
      select(GEOID, year, value)
  }
}

# ---------- 5) Trend slope helper ----------
trend_slope <- function(years_num) {
  force(years_num)
  function(v){
    idx <- !is.na(v)
    if (sum(idx) < 2) return(NA_real_)
    x <- years_num[idx]; y <- v[idx]
    x <- x - mean(x); y <- y - mean(y)
    sum(x*y) / sum(x*x)  # units per year
  }
}

# ---------- 6) One variable pipeline ----------
run_var <- function(var, sample_year = NULL, quick_fact = 4) {
  # Discover files/years
  ff <- find_files(var)
  files <- ff$files
  years <- ff$years
  if (is.null(sample_year)) sample_year <- max(years)  # map the latest year we have
  if (!sample_year %in% years) sample_year <- years[1L]
  
  # Build stack (one layer per year)
  s <- rast(files)
  names(s) <- paste0("y", years)
  
  # Reproject counties to raster CRS and dissolve to CONUS polygon
  counties_sf <- st_transform(counties_sf_wgs, crs = crs(s))
  conus_v <- terra::aggregate(vect(counties_sf))
  
  # Crop/mask to CONUS
  s_conus <- crop(s, ext(conus_v), snap = "out")
  s_conus <- mask(s_conus, conus_v)
  
  # ---- Yearly stats ----
  stats <- lapply(seq_along(years), function(i) {
    g <- global(s_conus[[i]], fun = c("min","max","mean","sd"), na.rm = TRUE)
    tibble(year = years[i], min = g[1,1], max = g[1,2], mean = g[1,3], sd = g[1,4])
  }) |> bind_rows()
  write_csv(stats, file.path(out_dir_daymet, make_name(glue("daymet_{var}_yearly_stats"), "csv")))
  
  # ---- Single-year map (quicklook) ----
  lyr <- which(years == sample_year)[1]
  r_year <- aggregate(s_conus[[lyr]], fact = quick_fact)
  png_map <- file.path(out_dir_daymet, make_name(glue("daymet_{var}_map_{sample_year}"), "png"))
  png(png_map, width = 1400, height = 900, res = 150)
  plot(r_year, main = glue("Daymet {var} {sample_year} (annual) — CONUS"))
  dev.off()
  
  # ---- Trend slope (full-res GeoTIFF + coarse PNG) ----
  slope_fun <- trend_slope(years)
  slope_full <- app(s_conus, slope_fun)
  tif_path <- file.path(out_dir_daymet, make_name(glue("daymet_{var}_trend_slope_conus"), "tif"))
  writeRaster(slope_full, tif_path, overwrite = TRUE)
  
  slope_coarse <- app(aggregate(s_conus, fact = quick_fact), slope_fun)
  q <- quantile(values(slope_coarse), probs = c(0.02, 0.98), na.rm = TRUE)
  mx <- max(abs(q), na.rm = TRUE)
  png_slope <- file.path(out_dir_daymet, make_name(glue("daymet_{var}_trend_slope_conus_quicklook"), "png"))
  png(png_slope, width = 1400, height = 900, res = 150)
  plot(slope_coarse,
       main = glue("Daymet {var} trend slope (units/year), CONUS\n{min(years)}–{max(years)} (coarse preview)"),
       col = hcl.colors(100, "Blue-Red 3"),
       type = "continuous", range = c(-mx, mx))
  dev.off()
  
  # ---- County means per year ----
  county_year <- get_county_means(s_conus, years, counties_sf) |>
    rename(!!paste0(var, if (var == "prcp") "_mean" else "_mean") := value)
  # Note: prcp files are annual totals per pixel; county aggregation here is the area-weighted mean of pixel totals
  
  out_csv <- file.path(out_dir_daymet, make_name(glue("daymet_{var}_county_year"), "csv"))
  write_csv(county_year, out_csv)
  
  list(
    years = years,
    stats = stats,
    county_year = county_year,
    map_png = png_map,
    slope_tif = tif_path,
    slope_png = png_slope
  )
}

# ---------- 7) Run for all variables & also write a combined county×year table ----------
vars <- c("tmax","tmin","vp","prcp")
results <- list()
for (v in vars) {
  cat("\n=== Processing", v, "===\n")
  results[[v]] <- run_var(v, sample_year = NULL, quick_fact = 4)
}

# Combined county×year (wide, one row per county-year with all vars)
combine_cy <- function(lst) {
  # lst: named list of per-var county_year tibbles with columns GEOID, year, <var>_mean
  ky <- Reduce(function(x, y) full_join(x, y, by = c("GEOID","year")),
               lapply(names(lst), function(v) results[[v]]$county_year))
  arrange(ky, GEOID, year)
}
daymet_allvars <- combine_cy(results)

write_csv(daymet_allvars, file.path(out_dir_daymet, make_name("daymet_county_year_allvars", "csv")))

cat("\nDone. Outputs in:\n", normalizePath(out_dir_daymet), "\n", sep = "")


######################################################
########### SOCIAL VULNERABILITY INDEX ###############
######################################################

# ---- Naming helpers (reuse yours if present) ----
sanitize_tag <- function(x){x <- if (is.null(x)) "SITE" else as.character(x); x <- iconv(x,"ASCII//TRANSLIT"); x <- gsub("[^A-Za-z0-9]+","_",x); x <- gsub("^_+|_+$","",x); if (!nzchar(x)) "SITE" else x}
SITE_NAME   <- if (exists("SITE_NAME")) SITE_NAME else sanitize_tag(if (exists("config")) config$site_name else Sys.getenv("SITE_NAME","SITE"))
SYSTEM_TIME <- if (exists("SYSTEM_TIME")) SYSTEM_TIME else format(Sys.time(), "%Y%m%dT%H%M%S")
make_name   <- function(result, ext=NULL){ nm <- paste(result, SITE_NAME, SYSTEM_TIME, sep="_"); if (is.null(ext)) nm else paste0(nm,".",ext) }
out_dir_svi <- file.path("output", paste0("svi_conus_", SITE_NAME, "_", SYSTEM_TIME))
if (!dir.exists(out_dir_svi)) dir.create(out_dir_svi, recursive = TRUE)

# ---------- Config ----------
YEARS_TARGET <- 2018:2024
nearest_release <- function(y) dplyr::case_when(
  y %in% c(2018, 2019) ~ 2018L,
  y %in% c(2020, 2021) ~ 2020L,
  TRUE                 ~ 2022L            # 2022 for 2022–2024
)
non_conus <- c("02","15","60","66","69","72","78")  # AK, HI, AS, GU, MP, PR, VI
map_year  <- 2022L

# ---------- Read local CSVs & standardize ----------
read_local_svi <- function(year){
  # be tolerant of case differences
  pat <- paste0("^SVI[_-]?", year, "[_-]US[_-]county\\.csv$")
  fp  <- list.files(".", pattern = pat, ignore.case = TRUE, full.names = TRUE)
  stopifnot(length(fp) == 1)
  suppressMessages(read_csv(fp, show_col_types = FALSE))
}

to_std <- function(df, yr){
  # normalize names
  nml <- tolower(names(df))
  names(df) <- nml
  fips_col <- if ("fips" %in% nml) "fips" else if ("geoid" %in% nml) "geoid" else NULL
  stopifnot(!is.null(fips_col))
  
  # some files use rpl_themes / rpl_theme1..4; keep those
  pick <- function(nm) { if (nm %in% names(df)) df[[nm]] else NA_real_ }
  out <- tibble(
    GEOID       = str_pad(as.character(df[[fips_col]]), 5, "left", "0"),
    statefp     = substr(GEOID, 1, 2),
    svi_overall = pick("rpl_themes"),
    svi_theme1  = pick("rpl_theme1"),
    svi_theme2  = pick("rpl_theme2"),
    svi_theme3  = pick("rpl_theme3"),
    svi_theme4  = pick("rpl_theme4"),
    year_release = yr
  )
  out
}

svi_release <- bind_rows(
  to_std(read_local_svi(2018), 2018),
  to_std(read_local_svi(2020), 2020),
  to_std(read_local_svi(2022), 2022)
) %>% filter(!statefp %in% non_conus)

# ---------- Build county×year table (2018–2024 using nearest release) ----------
svi_cy <- tidyr::expand_grid(GEOID = unique(svi_release$GEOID), year = YEARS_TARGET) %>%
  mutate(year_release = nearest_release(year)) %>%
  left_join(svi_release %>% select(GEOID, year_release, starts_with("svi_")),
            by = c("GEOID","year_release")) %>%
  arrange(GEOID, year)

# ---------- Yearly stats (overall SVI across counties) ----------
svi_yearly_stats <- svi_cy %>%
  group_by(year) %>%
  summarize(
    min  = min(svi_overall, na.rm = TRUE),
    max  = max(svi_overall, na.rm = TRUE),
    mean = mean(svi_overall, na.rm = TRUE),
    sd   = sd(svi_overall,   na.rm = TRUE),
    .groups = "drop"
  )
write_csv(svi_yearly_stats, file.path(out_dir_svi, make_name("svi_yearly_stats", "csv")))

# ---------- Map (single year) ----------
options(tigris_use_cache = TRUE)
counties_sf <- tigris::counties(cb = TRUE, year = 2020, class = "sf") %>%
  mutate(GEOID = str_pad(GEOID, 5, "left", "0"),
         statefp = substr(GEOID,1,2)) %>%
  filter(!statefp %in% non_conus) %>%
  select(GEOID, geometry)

svi_map_df <- svi_cy %>% filter(year == map_year) %>% select(GEOID, svi_overall)
map_sf <- counties_sf %>% left_join(svi_map_df, by = "GEOID")

p_map <- ggplot(map_sf) +
  geom_sf(aes(fill = svi_overall), linewidth = 0) +
  scale_fill_viridis_c(option = "magma", limits = c(0,1), oob = squish,
                       na.value = "grey90", name = glue("SVI overall ({map_year})")) +
  labs(title = glue("CDC/ATSDR SVI overall — counties, {map_year}"),
       caption = "SVI is a percentile (0–1); higher = more vulnerable") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major = element_line(color = "grey92"))
ggsave(file.path(out_dir_svi, make_name(glue("svi_overall_map_{map_year}"), "png")),
       p_map, width = 11, height = 7, dpi = 150)

# ---------- (Optional) Trend slope map from releases (2018,2020,2022) ----------
rel_years <- c(2018L, 2020L, 2022L)
svi_rel <- svi_release %>% select(GEOID, year_release, svi_overall) %>%
  filter(year_release %in% rel_years)

slope_df <- svi_rel %>%
  group_by(GEOID) %>%
  summarize(
    slope = {
      x <- year_release; y <- svi_overall
      if (sum(!is.na(y)) < 2) NA_real_ else { x <- x - mean(x); y <- y - mean(y); sum(x*y)/sum(x*x) }
    },
    .groups = "drop"
  )

map_slope_sf <- counties_sf %>% left_join(slope_df, by = "GEOID")
q  <- quantile(map_slope_sf$slope, c(0.02, 0.98), na.rm = TRUE); mx <- max(abs(q), na.rm = TRUE)
p_slope <- ggplot(map_slope_sf) +
  geom_sf(aes(fill = slope), linewidth = 0) +
  scale_fill_distiller(palette = "RdBu", limits = c(-mx, mx), oob = squish,
                       name = "SVI slope\n(Δ per year)") +
  labs(title = "CDC/ATSDR SVI county trend slope (2018–2022)",
       caption = "Caution: SVI percentiles are within-year; slopes reflect change in relative standing.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major = element_line(color = "grey92"))
ggsave(file.path(out_dir_svi, make_name("svi_overall_trend_slope_conus", "png")),
       p_slope, width = 11, height = 7, dpi = 150)

# ---------- Save county×year table ----------
write_csv(svi_cy, file.path(out_dir_svi, make_name("svi_county_year", "csv")))



































