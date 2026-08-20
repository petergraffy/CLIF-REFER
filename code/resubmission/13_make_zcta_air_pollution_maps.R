#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(ragg)
  library(scales)
  library(sf)
  library(viridis)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

find_repo <- function() {
  env_repo <- Sys.getenv("REFER_REPO", unset = "")
  if (nzchar(env_repo) && dir.exists(env_repo)) {
    return(normalizePath(env_repo, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg)) {
    candidate <- normalizePath(dirname(file_arg[[1]]), mustWork = TRUE)
    repeat {
      if (file.exists(file.path(candidate, "renv.lock")) &&
          dir.exists(file.path(candidate, "exposome", "zcta"))) {
        return(candidate)
      }
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  normalizePath(getwd(), mustWork = TRUE)
}

repo <- find_repo()
out_dir <- normalizePath(
  Sys.getenv(
    "REFER_MAP_OUTPUT_DIR",
    unset = file.path(repo, "figures", "air_pollution_maps_zcta")
  ),
  mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zcta_dir <- file.path(repo, "exposome", "zcta")
arf_county_counts_path <- Sys.getenv(
  "REFER_ARF_COUNTY_COUNTS_PATH",
  unset = file.path(repo, "sites", "analysis", "arf_counts_by_county_total.csv")
)
shape_path <- Sys.getenv(
  "REFER_ZCTA_SHAPEFILE",
  unset = file.path(repo, "data", "geography", "zcta_2020", "cb_2020_us_zcta520_500k.shp")
)
county_shape_path <- Sys.getenv(
  "REFER_COUNTY_SHAPEFILE",
  unset = file.path(repo, "data", "geography", "county_2022", "cb_2022_us_county_500k.shp")
)
zip_county_path <- Sys.getenv(
  "REFER_ZIP_COUNTY_CROSSWALK",
  unset = "/Users/saborpete/Desktop/Peter/Postdoc/Data/zip_county_latest_by_year.rds"
)
if (!file.exists(shape_path)) {
  stop(
    "ZCTA shapefile not found at ", shape_path, ". ",
    "Download cb_2020_us_zcta520_500k.zip from the Census cartographic boundary files ",
    "and unzip it into data/geography/zcta_2020/."
  )
}
if (!file.exists(county_shape_path)) {
  stop(
    "County shapefile not found at ", county_shape_path, ". ",
    "Expected the Census cb_2022_us_county_500k shapefile in data/geography/county_2022/."
  )
}
if (!file.exists(arf_county_counts_path)) {
  stop("ARF county count file not found at ", arf_county_counts_path, ".")
}

message("Reading ZCTA boundaries...")
zcta_ll <- sf::st_read(shape_path, quiet = TRUE) %>%
  rename(zip = ZCTA5CE20) %>%
  mutate(zip = as.character(zip)) %>%
  sf::st_crop(xmin = -126, xmax = -65, ymin = 22, ymax = 51)

zcta <- zcta_ll %>%
  sf::st_transform(5070)

conus_bbox <- sf::st_bbox(zcta)
conus_pad_x <- as.numeric(conus_bbox[["xmax"]] - conus_bbox[["xmin"]]) * 0.015
conus_pad_y <- as.numeric(conus_bbox[["ymax"]] - conus_bbox[["ymin"]]) * 0.035

message("Reading county boundaries...")
county <- sf::st_read(county_shape_path, quiet = TRUE) %>%
  transmute(county_fips = GEOID, geometry) %>%
  sf::st_crop(xmin = -126, xmax = -65, ymin = 22, ymax = 51) %>%
  sf::st_transform(5070)

no2_path <- file.path(zcta_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")
pm25_path <- file.path(zcta_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet")

message("Summarizing NO2, 2015-2024...")
no2_map <- arrow::read_parquet(no2_path) %>%
  transmute(
    zip = as.character(zip),
    year = as.integer(year),
    value = as.numeric(no2)
  ) %>%
  filter(year >= 2015, year <= 2024, is.finite(value)) %>%
  group_by(zip) %>%
  summarise(
    mean_exposure = mean(value, na.rm = TRUE),
    years_observed = n_distinct(year),
    .groups = "drop"
  ) %>%
  mutate(
    pollutant = "Nitrogen dioxide",
    abbreviation = "NO2",
    units = "ppb",
    requested_year_start = 2015L,
    requested_year_end = 2024L,
    mapped_year_start = 2015L,
    mapped_year_end = 2024L
  )

message("Summarizing PM2.5, 2015-2023 available ZCTA data...")
pm25_map <- arrow::read_parquet(pm25_path) %>%
  transmute(
    zip = as.character(zip),
    year = as.integer(year),
    month = as.integer(month),
    value = as.numeric(pm25_ug_m3)
  ) %>%
  filter(year >= 2015, year <= 2023, is.finite(value)) %>%
  group_by(zip) %>%
  summarise(
    mean_exposure = mean(value, na.rm = TRUE),
    years_observed = n_distinct(year),
    months_observed = n_distinct(paste(year, month)),
    .groups = "drop"
  ) %>%
  mutate(
    pollutant = "Fine particulate matter",
    abbreviation = "PM2.5",
    units = "ug/m3",
    requested_year_start = 2015L,
    requested_year_end = 2024L,
    mapped_year_start = 2015L,
    mapped_year_end = 2023L
  )

readr::write_csv(no2_map, file.path(out_dir, "zcta_mean_no2_2015_2024.csv"))
readr::write_csv(pm25_map, file.path(out_dir, "zcta_mean_pm25_2015_2023_available.csv"))

message("Reading pooled ARF patient counts by county...")
county_counts <- readr::read_csv(arf_county_counts_path, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    county_fips = sprintf("%05s", as.character(fips)),
    arf_patients = as.integer(arf_total)
  ) %>%
  filter(!is.na(county_fips), !is.na(arf_patients)) %>%
  arrange(desc(arf_patients))

county_count_qc <- tibble::tibble(
  source_file = arf_county_counts_path,
  counties_in_source = nrow(county_counts),
  counties_with_any_arf_patients = sum(county_counts$arf_patients > 0, na.rm = TRUE),
  total_arf_patients = sum(county_counts$arf_patients, na.rm = TRUE),
  max_county_arf_patients = max(county_counts$arf_patients, na.rm = TRUE)
)

readr::write_csv(county_counts, file.path(out_dir, "county_arf_patient_counts_total.csv"))
readr::write_csv(county_count_qc, file.path(out_dir, "county_arf_patient_counts_qc.csv"))

theme_map <- function(base_size = 12) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(8, 10, 8, 10)
    )
}

make_map <- function(summary_df, legend_title, filename_stem, option) {
  limits <- as.numeric(quantile(summary_df$mean_exposure, probs = c(0.01, 0.99), na.rm = TRUE))
  breaks <- pretty(limits, n = 5)
  map_df <- zcta %>%
    left_join(summary_df, by = "zip") %>%
    mutate(fill_value = pmin(pmax(mean_exposure, limits[[1]]), limits[[2]]))

  p <- ggplot(map_df) +
    geom_sf(aes(fill = fill_value), color = NA, linewidth = 0) +
    scale_fill_viridis_c(
      option = option,
      direction = 1,
      limits = limits,
      breaks = breaks,
      labels = label_number(accuracy = 0.1),
      oob = squish,
      na.value = "grey92",
      name = legend_title,
      guide = guide_colorbar(
        barwidth = unit(4.8, "in"),
        barheight = unit(0.18, "in"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    coord_sf(
      xlim = c(conus_bbox[["xmin"]] - conus_pad_x, conus_bbox[["xmax"]] + conus_pad_x),
      ylim = c(conus_bbox[["ymin"]] - conus_pad_y, conus_bbox[["ymax"]] + conus_pad_y),
      expand = FALSE
    ) +
    theme_map()

  png_path <- file.path(out_dir, paste0(filename_stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(filename_stem, ".pdf"))
  ragg::agg_png(png_path, width = 10.5, height = 5.6, units = "in", res = 450)
  print(p)
  dev.off()
  ggsave(pdf_path, p, width = 10.5, height = 5.6, units = "in", device = cairo_pdf)
  invisible(list(plot = p, png = png_path, pdf = pdf_path))
}

make_county_count_map <- function(counts_df) {
  map_df <- county %>%
    left_join(counts_df, by = "county_fips") %>%
    mutate(
      arf_patients = if_else(is.na(arf_patients), NA_integer_, as.integer(arf_patients)),
      arf_count_bin = case_when(
        is.na(arf_patients) | arf_patients <= 0L ~ "NA",
        arf_patients < 30L ~ "<30",
        arf_patients < 50L ~ "30-49",
        arf_patients < 100L ~ "50-99",
        arf_patients < 250L ~ "100-249",
        arf_patients < 500L ~ "250-499",
        arf_patients < 1000L ~ "500-999",
        arf_patients < 2500L ~ "1,000-2,499",
        arf_patients < 5000L ~ "2,500-4,999",
        arf_patients < 10000L ~ "5,000-9,999",
        arf_patients < 25000L ~ "10,000-24,999",
        TRUE ~ ">=25,000"
      )
    )

  bin_levels <- c(
    "NA", "<30", "30-49", "50-99", "100-249", "250-499", "500-999",
    "1,000-2,499", "2,500-4,999", "5,000-9,999", "10,000-24,999", ">=25,000"
  )
  observed_levels <- bin_levels[bin_levels %in% map_df$arf_count_bin]
  map_df <- map_df %>%
    mutate(arf_count_bin = factor(arf_count_bin, levels = observed_levels))

  arf_cols <- c(
    "NA" = "grey82",
    stats::setNames(
      viridis::viridis(length(bin_levels) - 1, option = "cividis", direction = 1),
      bin_levels[-1]
    )
  )[observed_levels]

  p <- ggplot(map_df) +
    geom_sf(aes(fill = arf_count_bin), color = NA, linewidth = 0) +
    scale_fill_manual(
      values = arf_cols,
      drop = FALSE,
      name = "ARF patients (n)",
      guide = guide_legend(
        nrow = 2,
        byrow = TRUE,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    coord_sf(
      xlim = c(conus_bbox[["xmin"]] - conus_pad_x, conus_bbox[["xmax"]] + conus_pad_x),
      ylim = c(conus_bbox[["ymin"]] - conus_pad_y, conus_bbox[["ymax"]] + conus_pad_y),
      expand = FALSE
    ) +
    theme_map()

  png_path <- file.path(out_dir, "county_map_arf_patient_counts_total.png")
  pdf_path <- file.path(out_dir, "county_map_arf_patient_counts_total.pdf")
  ragg::agg_png(png_path, width = 10.5, height = 5.6, units = "in", res = 450)
  print(p)
  dev.off()
  ggsave(pdf_path, p, width = 10.5, height = 5.6, units = "in", device = cairo_pdf)
  invisible(list(plot = p, png = png_path, pdf = pdf_path))
}

no2_out <- make_map(
  no2_map,
  expression("Nitrogen dioxide (NO"[2]*", ppb)"),
  "zcta_map_no2_mean_2015_2024",
  "magma"
)

pm25_out <- make_map(
  pm25_map,
  expression("Fine particulate matter (PM"[2.5]*", "*mu*"g/m"^3*")"),
  "zcta_map_pm25_mean_2015_2023_available",
  "inferno"
)

county_count_out <- make_county_count_map(county_counts)

ab_air_pollution_plot <- (no2_out$plot / pm25_out$plot) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.015, 0.98)
  )

ab_air_pollution_png <- file.path(out_dir, "figure_ab_air_pollution_maps.png")
ab_air_pollution_pdf <- file.path(out_dir, "figure_ab_air_pollution_maps.pdf")
ragg::agg_png(ab_air_pollution_png, width = 10.5, height = 11.2, units = "in", res = 450)
print(ab_air_pollution_plot)
dev.off()
ggsave(ab_air_pollution_pdf, ab_air_pollution_plot, width = 10.5, height = 11.2, units = "in", device = cairo_pdf)

abc_plot <- (county_count_out$plot / no2_out$plot / pm25_out$plot) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.015, 0.98)
  )

abc_png <- file.path(out_dir, "figure_abc_air_pollution_and_arf_count_maps.png")
abc_pdf <- file.path(out_dir, "figure_abc_air_pollution_and_arf_count_maps.pdf")
ragg::agg_png(abc_png, width = 10.5, height = 16.8, units = "in", res = 450)
print(abc_plot)
dev.off()
ggsave(abc_pdf, abc_plot, width = 10.5, height = 16.8, units = "in", device = cairo_pdf)

metadata <- tibble::tibble(
  map = c("NO2", "PM2.5", "ARF patient counts", "AB air pollution composite", "ABC composite"),
  png = c(no2_out$png, pm25_out$png, county_count_out$png, ab_air_pollution_png, abc_png),
  pdf = c(no2_out$pdf, pm25_out$pdf, county_count_out$pdf, ab_air_pollution_pdf, abc_pdf),
  source_data = c(no2_path, pm25_path, arf_county_counts_path, paste(c(no2_path, pm25_path), collapse = ";"), paste(c(arf_county_counts_path, no2_path, pm25_path), collapse = ";")),
  boundary_data = c(shape_path, shape_path, county_shape_path, shape_path, paste(c(shape_path, county_shape_path), collapse = ";")),
  mapped_years = c("2015-2024", "2015-2023", "all available pooled ARF county counts", "mixed", "mixed"),
  note = c(
    "NO2 annual ZCTA data available through 2025.",
    "Bundled ZCTA PM2.5 data available through 2023; 2024 was not mapped.",
    "County counts read from sites/analysis/arf_counts_by_county_total.csv.",
    "Composite contains NO2 and PM2.5 maps.",
    "Composite contains ARF patient county counts, NO2, and PM2.5."
  )
)
readr::write_csv(metadata, file.path(out_dir, "zcta_air_pollution_map_manifest.csv"))

message("Wrote maps to: ", out_dir)
