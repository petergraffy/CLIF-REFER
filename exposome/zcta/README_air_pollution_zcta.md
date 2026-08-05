# Air Pollution ZCTA Parquet Release Assets

These files contain final filled air-pollution exposure surfaces aggregated to 2020 ZCTA5 polygons.

Assets:

- `air_pollution_zcta_pm25_monthly_2005_2023.parquet`: monthly PM2.5, `pm25_ug_m3`.
- `air_pollution_zcta_o3_monthly_2005_2023.parquet`: monthly ozone, `o3_ppb`.
- `air_pollution_zcta_no2_annual_2005_2025.parquet`: annual NO2, `no2`.
- `air_pollution_zcta_parquet_manifest.csv`: sizes, row counts, completeness, and SHA-256 checksums.
- `fill_summary.csv`: nearest-fill audit by source file.
- `fill_audit_summary.csv`: additional fill audit summary when available.

Columns include `zip`, `pollutant`, `year`, `month`, `temporal_resolution`, the pollutant value column, `value_source`, and `fill_distance_km`.
PM2.5 and ozone currently cover 2005-2023 in the local source files; NO2 covers 2005-2025.
