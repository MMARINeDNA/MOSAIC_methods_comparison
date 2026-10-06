These are the processed data files after the bioinformatics pipeline and the associated metadata files for the samples used in this study.


**Revision round 2 — sample locations.** Positions are validated and repaired from the ship GPS tracks (`ship_tracks/`, not tracked in git: `SR2311_combined_gps_coordinates.csv`, `SR2408_combined_gps_coordinates.csv`; UTC time, latitude, longitude) by `scripts/1.5_metadata_geocorrections.R`. `reference/CalCOFIStationOrder.csv` is the official CalCOFI station position table (https://calcofi.org/downloads/maps/CalCOFIStationOrder.csv), used for CTD position fill and checks.
