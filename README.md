# MOSAIC_methods_comparison [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23942563.svg)](https://doi.org/10.5281/zenodo.23942563)
protocols and data scripts associated with MOSAIC eDNA sampling methods comparison paper

## Reproducing the analyses

| Folder | Contents |
|---|---|
| `processed_data/` | Processed eDNA data and corrected sample metadata (the starting point; see its README) |
| `data/` | Sample metadata (`SR2311_SR2408_eDNA_metadata.csv`) |
| `visual_data/` | CalCOFI visual survey data used for the visual-vs-eDNA comparison |
| `scripts/` | Numbered analysis scripts and `run_all.R` (see `scripts/README.md` for run order, outputs and R packages) |
| `ship_tracks/`, `protocols/`, `supplement/` | Ship GPS tracks, laboratory protocols, supplementary files |
| `viz/` | Figures from the paper (the scripts also write regenerated figures here) |
| `bioinformatics_code/` | Sequence-processing pipeline code (reads to taxon tables) |
| `archive/` | Superseded and diagnostic scripts, including the metadata position-correction script; not needed to reproduce the paper |

From the repository root, `Rscript scripts/run_all.R` re-creates every table and figure except Figs 1 and 2 (station map and apparatus photograph, no script). Results go to `outputs/`, `viz/` and `logs/` (`outputs/` and `logs/` are git-ignored).
