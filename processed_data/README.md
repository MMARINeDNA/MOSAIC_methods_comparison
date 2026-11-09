# MOSAIC eDNA metabarcoding data: bioinformatic outputs and processed phyloseq data

Processed eDNA metabarcoding data and corrected sample metadata used in the analyses of a study of active transect-based eDNA sampling. Seawater was filtered on CalCOFI cruises RS2311 (Fall 2023) and SR2408 (Summer 2024) with four methods (MOSAIC-A, -B, -C, -D). Each library is one PCR replicate of a vertebrate 12S rRNA amplicon (MarVer1 assay), sequenced on an Illumina NextSeq 2000 (2x150 bp). Raw reads are archived separately in NCBI SRA under BioProject PRJNA1556353.

Analysis code: https://github.com/MMARINeDNA/MOSAIC_methods_comparison

License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). Please cite the paper and this repository; see `CITATION.cff`.

## Raw sequence data

The raw sequencing reads (FASTQ) are archived in NCBI SRA under BioProject PRJNA1556353.

## Sampling methods

All filters are 5 um mixed cellulose ester. The method is recorded in `Collection method` / `Collection.method`.

| Label | Description |
|---|---|
| MOSAIC-A | Surface underway water filtered continuously through one filter over a whole transect (the ship's path between two stations). |
| MOSAIC-B | The same underway filtration, one filter per hour of transect, so a transect is covered by consecutive filters. SR2408 only. |
| MOSAIC-C | Underway water collected along a whole transect into a 50 L carboy, mixed, and 5 L filtered. SR2408 only. |
| MOSAIC-D | Vertically integrated CTD/Niskin samples (labelled "CTD" in the paper's figures). One filter, F102, is a discrete single-depth sample; see `Sample type`. |

## Folder contents

The data are in two folders of this repository (archived with it on Zenodo, DOI 10.5281/zenodo.23942562), described below. `bioinformatic_outputs/` is the input to `scripts/01`-`04`, which rebuild `processed_phyloseq/` (written to `outputs/`); `processed_phyloseq/` is the reference copy used for the paper. See `scripts/README.md`.

### `processed_phyloseq/`: the data analysed in the paper

163 filters and 491 PCR-replicate libraries, after removal of human and human-associated farm-animal reads. The three CSVs together reconstitute the phyloseq object.

| File | Contents |
|---|---|
| `MOSAIC_taxnamed_otu_table_metadata_filtered.csv` | 365 taxa (rows) x 491 samples (columns). Integer counts of reads assigned to each taxon. Row names are the lowest assigned rank (for example "Belonidae family", "Harengula thrissina"). Counts are not abundance estimates. |
| `MOSAIC_taxnamed_taxa_table_metadata_filtered.csv` | 365 taxa x 6 ranks (phylum, class, order, family, genus, species). Row names match the OTU table. |
| `MOSAIC_taxnamed_sample_data_metadata_filtered.csv` | 491 samples x 41 fields, one row per PCR-replicate library. First column is the sample name (for example `C47-1` = filter C47, replicate 1) and matches the OTU table columns. |
| `MOSAIC_taxnamed_metadata_filtered.ps.RDS` | The three tables as one phyloseq object (R package `phyloseq`). |

In these CSVs the first column holds the row names: it is labelled `Taxon` in the OTU and taxonomy tables and `SampleName` in the sample data (and in `geodata_corrected_MOSAIC_metadata.csv`). Read them with the first column as row names (for example `row.names = 1` in R) to rebuild the phyloseq object.

Notes for the sample data:
- The unit of analysis is the filter (`FilterID`). Four PCRs were run per sample, and replicates with no visible gel band were not sequenced, so not every filter has four libraries.
- Filtration positions and times were checked and corrected against the research vessel's GPS track and the CalCOFI station positions.
- `Inshore.Offshore` classifies each sample as inshore (on the continental shelf) or offshore (off the shelf). `dist_km` is the transect distance over which the water was filtered (0 for MOSAIC-D CTD samples). `datetime_start`, `datetime_end` and `filtration_time_min` are derived from the filtration times.

### `bioinformatic_outputs/`: sequence-processing outputs and corrected metadata

Taxon tables as produced by the bioinformatic pipeline, before human and human-associated farm-animal reads were removed and before restriction to the paper's samples. They contain additional samples and non-target taxa (including human and livestock taxa such as *Homo sapiens*). Values are integer read counts.

| File | Contents |
|---|---|
| `taxon_table_lineage_wide.csv` | 228 taxa x 369 samples. Columns: `BestTaxon`, Kingdom, Phylum, Class, Order, Family, Genus, Species, then one column per sample. |
| `taxon_table_wide - taxon_table_wide.csv` | 339 taxa x 331 samples. Columns: `BestTaxon`, `Class`, then one column per sample. |
| `taxon_table_wide-SR2408-GGBC2.csv` | 354 taxa x 188 samples. Columns: `BestTaxon`, `Class`, then one column per sample. |
| `geodata_corrected_MOSAIC_metadata.csv` | Corrected sample metadata, 508 rows x 37 fields (plus the `SampleName` column), covering the samples in the study including some not analysed in the paper. Filtration positions and times were checked and corrected against the research vessel's GPS track and the CalCOFI station positions. |
| `SR2311_SR2408_eDNA_metadata.csv` | Sample metadata (field, filtration and extraction records) for the filters and libraries of cruises RS2311 and SR2408, as recorded before the position corrections. Read by `scripts/01`. |
| `taxonomy_lookup_NCBI.csv` | Taxonomy (phylum to species) for each of the 415 OTUs in the taxon tables, as looked up from NCBI for the paper. Row names are the OTU ids. `scripts/01` reads it so that the analysis does not depend on the current NCBI taxonomy. |

## Taxonomic assignment

Reads were trimmed and denoised into amplicon sequence variants (Cutadapt, DADA2) and assigned to taxa with BLAST+ against the NCBI eukaryote nucleotide database, using a custom accession exclusion list and a lowest-common-ancestor step (taxonkit). Taxon labels are given as assigned. A detection is evidence of presence for the assigned taxon, not a confirmed observation. Full methods are in the paper and its supplement. The scripts that perform these steps are in the `bioinformatics_code/` directory of the GitHub repository (https://github.com/MMARINeDNA/MOSAIC_methods_comparison).
