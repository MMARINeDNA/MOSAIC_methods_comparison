# Analysis scripts

Everything in the paper that is computed (Tables 1–2, Figs 3–7, the supplementary figures, the visual-vs-eDNA comparison and the three MOSAIC-vs-CTD sensitivity analyses) can be re-created from the processed data in this repository. Figs 1 and 2 (station map, apparatus photo) have no script here.

```bash
Rscript scripts/run_all.R                 # whole pipeline, from the repository root
Rscript scripts/run_all.R --skip-models   # everything except the five NIMBLE models
Rscript scripts/run_all.R --from=08e      # restart at the first script whose name starts with "08e"
Rscript scripts/run_all.R --only=08a,08c  # just these scripts (prefix match)
```

`run_all.R` runs each script in its own `Rscript` session from the repository root, writes `logs/<script>.log`, creates the output folders, and stops at the first error. Every script can also be run by itself (`Rscript scripts/05_mosaic_figure3_familylevel.R`) once the scripts before it have run. Outputs go to `outputs/` (intermediate data, figures, model results), `viz/`, `model_evolution_tracking/` and `logs/`; all except `viz/` are git-ignored because they are regenerated. The figures from the paper are in `viz/`; the scripts also write their regenerated figures there.

## Inputs

| Input | Used by |
|---|---|
| `processed_data/bioinformatic_outputs/` (three taxon tables, `taxonomy_lookup_NCBI.csv`, `SR2311_SR2408_eDNA_metadata.csv` and `geodata_corrected_MOSAIC_metadata.csv`) | `01`, `02` |
| `visual_data/` (CalCOFI sightings and segments workbooks, species codes, `visuals_for_MOSAIC.csv`) | `03` |
| `visual_data/visuals_coalesced.csv` | `03`: the **one manual input**, the visual sightings coalesced by hand into events. It is read, not generated. |

`processed_data/processed_phyloseq/` is the reference copy of what `01`–`04` produce (see `processed_data/README.md`). Correcting the sample positions against the ship tracks is not part of this pipeline; it was done once, outside the scripts in this repository (the ship tracks are in `ship_tracks/`), and its result is `geodata_corrected_MOSAIC_metadata.csv`.

Taxonomy in `01` is pinned to the lookup used for the paper (`processed_data/bioinformatic_outputs/taxonomy_lookup_NCBI.csv`); the live NCBI lookup is kept in `01` as commented-out code. `03` still looks up the visual-survey species codes with `taxonomizr`, which needs the NCBI `accessionTaxa.sql` database (large; downloaded to the working directory on first use). To reuse a copy: `export ACCESSION_TAXA_SQL=/path/to/accessionTaxa.sql`.

## Run order

| Script | Produces |
|---|---|
| `01_csvs_to_phyloseq.R` | `outputs/MOSAIC.ps.RDS`: merges the three taxon tables with the metadata and adds the pinned taxonomy. |
| `02_prep_detection_prob_metadata.R` | Reads the corrected metadata; removes blanks and contaminants; `outputs/MOSAIC_taxnamed*` (phyloseq and CSVs), mammal-only phyloseq objects and tables. |
| `03_prep_visual_data.R` | Visual survey tables joined to eDNA detections: `outputs/visveDNA.csv`, `visvMOSAIC_eDNA.csv`, `visv4methods_eDNA.csv`, `Cetacean_MOSAIC_vs_visuals.csv`. **Visual vs eDNA comparison.** |
| `04_MOSAIC_metadata_cleaning_for_modeling.R` | `outputs/MOSAIC_taxnamed_{otu_table,taxa_table,sample_data}_metadata_filtered.csv` and `MOSAIC_taxnamed_metadata_filtered.ps.RDS`: the analysis data set (163 filters, 491 libraries, 365 taxa). All later scripts read these. |
| `05_mosaic_figure3_familylevel.R` | **Fig 3** (`viz/mosaic_figure3_familylevel_metadata_filtered.pdf`) |
| `06_suppfig_species-accumulation.R` | Supplementary species-accumulation curves (`viz/species_accumulation_curves_*.pdf`) |
| `07_suppfig_nested-venn.R` | Supplementary nested Venn of visual vs eDNA species (`viz/Figure5_nested_venn_marine_mammals.pdf`; the file name is a leftover, it is a supplementary figure). Uses base graphics only. |
| `08a_model_primary_MOSAIC-A_vs_CTD.R` | **Primary model** (MOSAIC-A vs CTD, all data): `outputs/v6/*_AvsD_v6_results.RData` and diagnostics in `outputs/v6/figures/stage1_v6/` |
| `08b_build_balanced-subset-pairs.R` | `outputs/v6/matched_pairs_v6.csv`: 1:1 CTD–A pairs within cruise × inshore/offshore |
| `08b_model_A_vs_CTD_balanced-subset.R` | **Sensitivity 1**, equal-number matched subset: `*_v6_matched_results.RData`, `figures/stage1_v6_matched/` |
| `08c_model_all-MOSAIC_vs_CTD.R` | **Sensitivity 2**, all MOSAIC methods (A/B/C) vs CTD: `*_allMOSAIC_v6_results.RData`, `figures/stage1_v6_allMOSAIC/` |
| `08d_model_A_vs_CTD_Fall-2023-pilot.R` | **Sensitivity 3**, Fall 2023 single-cruise comparison: `outputs/v5/*_stage1_AvsD_results.RData`, `outputs/v5/figures/stage1_AvsD/` |
| `08e_figure4_capture-probability.R` | **Fig 4** (`outputs/v6/figures/fig4_v6primary.pdf`), from the `08a` posterior |
| `08f_table2_MOSAIC_vs_CTD.R` | **Table 2**: `model_evolution_tracking/Table2_MOSAIC_vs_CTD.csv` and `.docx`, from the `08a`–`08d` posteriors |
| `09_model_transect-only.R` | Surface-only (A/B/C) model with volume and distance covariates: `outputs/v5/*_v5_results.RData`, `figures/v5_surface/` |
| `09a_figure5_transect-only-coefficients.R` | Coefficient panels for Fig 5 (`outputs/v5/figures/panels/`) |
| `08.6_A_B_taxa_detection_map.R`, `08.6_A_transect_source_ordination.R`, `08.6_transect_PCoA_dbRDA.R`, `08.6_transect_species_corroboration.R` | Panels and supplementary figures for the within-transect analysis (`outputs/v5/figures/transect_source/`, `panels/`, plus `outputs/v5/transect_source_*.csv`; the dbRDA statistics are printed in its log) |
| `assemble_manuscript_figures.R` | **Fig 5** (`outputs/v5/figures/v5_surface/fig5_model.pdf`) and **Fig 7** (`outputs/v5/figures/transect_source/fig7_q3_transect_source.pdf`) from the panels above |
| `09.5_A-D_cetacean_detection_chisq.R` | MOSAIC-A vs CTD cetacean detection rate and chi-squared test (printed in `logs/09.5_*.log`) |
| `10_rare_taxa_panels_final.R` | **Fig 6** rare-taxa maps (`viz/rare_taxa_panels_SR2408.pdf`, `viz/rare_taxa_panels_RS2311.pdf`, individual panels as PNG) |
| `11_make_summary_Table_1.R` | **Table 1** (`outputs/summary_table1.docx`; the `.html` versions need a recent pandoc) |
| `12_total_distance_sampled.R` | `outputs/total_distance_sampled.csv`: distance sampled per cruise, overlapping water counted once |

Dependencies that are not obvious from the numbering: `08e`, `08f` need the `08a`–`08d` RData; `09a` needs `09`; `assemble_manuscript_figures.R` needs `09a` and the four `08.6_*` scripts; `07` needs `03`; `10` reads the unfiltered `outputs/MOSAIC_taxnamed_*.csv` written by `02`.

## Runtimes

Measured on a 10-core Apple-silicon Mac (R 4.6.1, NIMBLE 1.4.3), in a clean copy of the repository. All scripts except the models take under a minute each (`01` about 15 s, `10` about 30 s). Each NIMBLE model (3 chains) took 4–10 minutes when the five were run at the same time: `08a` 8.4 min, `08b_model` 4.6, `08c` 9.7, `08d` 4.2, `09_model` 8.4. Run one after another that is roughly 35 minutes; the NIMBLE scripts are single-threaded, so they can be run in parallel with `--only=`. The complete pipeline is therefore well under an hour. The first call to NIMBLE compiles C++ code, so a working compiler toolchain is required.

## R packages

CRAN: `tidyverse` (dplyr, tidyr, purrr, tibble, ggplot2, readr, stringr, lubridate), `readxl`, `geosphere`, `geodist`, `vegan`, `maps`, `ggrepel`, `scales`, `cowplot`, `patchwork`, `nimble` (needs a C++ compiler), `coda`, `bayesplot`, `flextable`, `officer`, `kableExtra`, `taxonomizr`, `clue` (optional; `08b_build` falls back to a greedy match without it).
Bioconductor / GitHub: `phyloseq`, `microViz`, `metagMisc`, `speedyseq`.
No Venn package is needed. Saving Table 1 as `.html` needs a recent pandoc; if it is too old the script warns and still writes the `.docx`.
