# Monthly NO2 ZCTA Parquet Release Assets

These files contain monthly TROPOMI/LUR CONUS surface NO2 estimates aggregated to 2020 ZCTA5 polygons.

Assets:

- `no2_zcta_monthly_YYYY.parquet`: yearly monthly ZCTA NO2 values, one row per ZCTA-month.
- `no2_zcta_monthly_parquet_manifest.csv`: row counts, month coverage, completeness, file sizes, and SHA-256 checksums.
- `no2_zcta_monthly_qc_all_months.csv`: monthly completeness and distribution QC from the source aggregation.
- `README.md`: variables and source notes.

Columns:

- `zip`
- `pollutant`
- `year`
- `month`
- `temporal_resolution`
- `no2_ppbv`
- `value_source`
- `fill_distance_m`
- `fill_cell`

The `no2_ppbv` column is monthly surface NO2 in parts per billion by volume. `value_source`, `fill_distance_m`, and `fill_cell` provide audit metadata for nearest-cell fills; the completed 2019-2025 local outputs have zero missing ZCTA-month NO2 values.

Citation:

Nawaz, M. Omar. Monthly and Annual US TROPOMI Surface NO2 Estimates (~1km x 1km), Version v2. Zenodo, 2026. https://doi.org/10.5281/zenodo.18919769.

Zenodo record: https://zenodo.org/records/18919769.
Related publication DOI listed by Zenodo: https://doi.org/10.1021/acsestair.4c00153.
Also cite this repository/release for the ZCTA aggregation workflow.

The reproducible aggregation script is `code/56_aggregate_monthly_no2_to_zcta.R` and this packaging script is `code/57_package_monthly_no2_release_assets.R`.
