# MOSAIC_methods_comparison [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23942562.svg)](https://doi.org/10.5281/zenodo.23942562)
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

From the repository root, `Rscript scripts/run_all.R` re-creates every table and figure except Figs 1 and 2 (station map and apparatus photograph, no script). Results go to `outputs/`, `viz/` and `logs/` (`outputs/` and `logs/` are git-ignored).

## License

Code: MIT (`LICENSE`). Data in `processed_data/`, `data/` and `visual_data/`: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/), which requires attribution.

## How to cite

Wisely E, Patin N, Alksne MN, Dinasquet J, Swalethorp R, Bowlin N, McKaughan G, Coden F, Satterthwaite EV, Baumann-Pickering S, Semmens BX. Continuous Large-Volume Active Transect Sampling Improves Detection of Rare Marine Vertebrate eDNA Across Large Spatial Scales. Methods in Ecology and Evolution (2026, in press). Please also cite this repository (DOI above); see `CITATION.cff`.
