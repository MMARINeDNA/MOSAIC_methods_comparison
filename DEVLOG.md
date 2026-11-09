# DEVLOG — MOSAIC1 (Paper 1)

## 2026-10-05 — Revision round 2: corrected sample locations through the whole pipeline
- Did:
  - New `scripts/1.5_metadata_geocorrections.R` (runs between 01 and 02). Reads the original phyloseq object,
    `data/MOSAIC_taxnamed_sample_data_corrected-1.xlsx` and the hand-made Inshore/Offshore sheet; applies
    per-method position rules; runs QC; writes `outputs/geodata_corrected_MOSAIC_metadata.csv`
    (+ `outputs/geodata_QC_report.csv`). Settings (paths, method names, thresholds, optional study-area box
    and optional CalCOFI station file) are at the top of the script. `02` now reads that CSV.
  - `04` no longer uses hard-coded drop lists: it drops MOSAIC-A/B filters with implied speed > 14.8 kt and any
    filter missing a required covariate (incl. filtration end time). The F021 end-time correction is now
    written back to the data (previously only applied to a temporary copy).
  - Retained 163 of 165 filters (101 A, 20 B, 6 C, 36 CTD); the committed pre-round-2 output had 138.
    Excluded: F076 (no volume), F062 (no filtration end time).
  - Re-ran every model (08a-d, 09, v3/v4 comparison models) and the dependent figures/tables; max Gelman-Rubin
    point estimate across the main models <= 1.019. Primary MOSAIC-A vs CTD odds ratio 4.6 (95% CrI 1.9-11.8);
    transect-only model per SD: volume 1.31 (0.86-2.06), distance 1.08 (0.69-1.69), filtration time 0.90
    (0.75-1.07). Total distance sampled (overlap counted once): 861 km (Fall 2023), 2,991 km (Summer 2024).
  - Figure 6 layout changes (species legend as dots, larger names, panel letters A-D, no subtitle).
  - Supporting checks: scripts 12-17 (non-overlapping distance, MOSAIC-C tracks before/after, CTD vs CalCOFI
    station grid, neighbouring-transect stations, CTD evidence table, F025 map).
- Decided (with the data owner):
  - Workbook positions come from the ship log at the filtration start/end times, so they are used only for
    MOSAIC-A/B filters with a time at both ends and for MOSAIC-D (CTD, filtered on a CalCOFI station).
    MOSAIC-C (carboy) and filters with a missing time keep their original coordinates.
  - If the workbook's own note says "kept original coordinates" but its position is > 1 km from the original,
    the original is kept (rule selects F101, F102, F087, F031, F057).
  - A transect replicate whose end point equals that of a shorter filter with the same start time
    (F025-1 inherited F026's end) is repaired from the filter's other replicates.
- Open / next:
  - CTD casts still flagged in `outputs/geodata_QC_report.csv` (workbook position not on the labelled station):
    C47, C35/C36, C68/C69, C10/C11, C29, C76/C77, C15/C16 — decide positions.
  - F062: needs its filtration end time (and end position).
  - Manuscript text and tables still carry pre-round-2 numbers; a tracked-changes version was built locally
    (not committed). The text's 107 MOSAIC-A filters = 103 real + 4 field-blank filters labelled MOSAIC-A (Fall23, Fall45, Fall62, Fall79); the tracked-changes version says 103.
  - Model result files (`outputs/v4`, `outputs/v5`, `outputs/v6`, ~1.7 GB) are not tracked.
  - `.gitignore` contains unresolved merge-conflict markers (pre-existing).

## 2026-10-05 (later) — Locations from the ship GPS track instead of the workbook
- Did:
  - Found that `data/MOSAIC_taxnamed_sample_data_corrected-1.xlsx` looked up ship positions with a fixed UTC-7 offset, so every Fall 2023
    filter from 5 Nov onward (after daylight saving ended) was one hour off; Fall34 showed 1 km instead of ~10 km. The workbook is no longer used.
  - Rewrote `scripts/1.5_metadata_geocorrections.R`: MOSAIC-A/B starts and ends are checked against the ship track (`ship_tracks/`, not in git)
    at the logged clock time (local time -> UTC via America/Los_Angeles); positions > 2 km from the track, or missing, are replaced by the
    track position. F021's end time (12:22) is retried 12 h later. MOSAIC-C keeps its field coordinates (carboy filled along the transect).
    CTD casts are placed at the mean ship position of the latest stationary period (< 0.5 kt, >= 8 min) before filtering began; the CalCOFI
    station label is only a check (and the fallback). Controls (field blanks) stay in the output table with no collection method;
    `02` separates them. `outputs/controls_MOSAIC_metadata.csv` removed.
  - CTD labels contradicted by the ship stops: C35/C36 (labelled 86.7 40, ship stopped at 86.7 35), C47 (labelled 86.7 80, stopped at 86.7 90);
    F010 (no label) stopped at 93.3 50; F057's field position (86.7 110) differs from its ship stop (83.3 110). See
    `outputs/ctd_label_evidence_ship_stops_neighbours.csv` (script 21).
  - Re-ran 02, 04, 08b builder, all models (08a-d, 09, v3, v4) and every figure/table script. Analysis set unchanged (163 filters: 101 A, 20 B,
    6 C, 36 CTD). New supporting scripts 18-21.
- Results after the re-run (read from the output files): Table 2 odds ratios MOSAIC-A vs CTD 4.6 (1.9-12.2) primary, 4.7 (1.8-13.3) balanced,
  4.6 (1.9-12.0) all surface methods, 6.7 (2.4-20.3) Fall pilot. Transect-only (v5 stage 2) per SD: volume 1.32 (0.87-2.05), distance 0.99
  (0.67-1.47), filtration time 0.90 (0.75-1.08). Non-overlapping distance sampled: 1,190 km (Fall 2023), 2,987 km (Summer 2024).
  Convergence (3 chains x 5000 saved draws, computed from the saved chains): max Gelman-Rubin over all monitored parameters 1.016 (primary), 1.059 (balanced), 1.026 (all surface), 1.032 (Fall pilot), 1.005 (transect-only v5), 1.007 (v3), 1.041 (v4); the only values > 1.05 are two latent site x family occupancy states in the balanced model (phi[1,5], phi[1,6]); every non-latent parameter <= 1.007; min ESS 532-1549 (mu_alpha).
- Decided: ship stops (not station labels) place CTDs; a field position is kept for A/B only if within 2 km of the track.
- Open / next:
  - `manuscript_revision2/MOSAIC manuscript Revision2_tracked_changes.docx` was rebuilt from these results (88 tracked edits, 11 comments; reject-all reproduces the submitted Revision 1 text exactly). `viz/final_figures` PDFs were refreshed from the new outputs. Check the Table 2 'Max R-hat' wording (the Methods text says < 1.05 for all parameters; the balanced model has two latent values at 1.059).
  - F062 has no filtration end time (dropped); F076 has no volume (dropped).
  - Summer 2024 data-sheet times were compared with the track only through script 1.5; scripts 13, 17 and 19 still read the old workbook for comparison only.
  - `.gitignore` still contains unresolved merge-conflict markers (pre-existing).

## 2026-11-09 — release-sync: run_all.R and clean-copy test
- Did: added `scripts/run_all.R`; rewrote `scripts/README.md`, updated top-level and `processed_data/` READMEs; `.gitignore` now ignores `outputs/`, `viz/`, `logs/`, `model_evolution_tracking/`. Small fixes: `01` handles taxonomizr >= 0.11 (`domain` replaces `superkingdom`) and takes `ACCESSION_TAXA_SQL`; `03` same env var; `11` warns instead of stopping when pandoc is too old for the HTML export.
- Tested the full pipeline in a clean copy (R 4.6.1, NIMBLE 1.4.3): all 25 scripts ran. Regenerated data, Table 1 and distance table match the rev2 reference; taxa-table family labels for 5 taxa differ (NCBI taxonomy version), which removes Alopiidae detections (2 -> 0) and moves the balanced-subset OR to 4.4. Started from `processed_data/processed_phyloseq/`, Table 2 reproduces exactly (4.6, 4.7, 4.6, 6.7).
- Open / next: decide how to pin the taxonomy so `01`-`04` reproduce the shipped taxa table; confirm `10_rare_taxa_panels_final` is the Fig 6 script.
- Update: per the taxonomy mismatch above, `01` now reads the pinned `processed_data/bioinformatic_outputs/taxonomy_lookup_NCBI.csv` (tax_table of MOSAIC.ps from the rev2 run) and the live lookup is commented out. Re-ran 01, 02, 04 in a clean copy: the three filtered CSVs and the `.ps.RDS` taxa table are identical to `processed_phyloseq/`. Models were not re-run after this change (inputs are now identical to the files that reproduced Table 2).
