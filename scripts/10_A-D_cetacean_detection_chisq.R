library(tidyverse)

# Load tables
sample <- read_csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv")
taxa   <- read_csv("outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv")
otu    <- read_csv("outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv")

# Cetacean OTUs at family level or finer
cet_taxa <- taxa |>
  filter(order == "Artiodactyla", !str_detect(family, "order")) |>
  pull(`...1`)

# Per-replicate detection: any cetacean reads > 0
cet_detected <- otu |>
  filter(`...1` %in% cet_taxa) |>
  select(-`...1`) |>
  summarise(across(everything(), \(x) any(x > 0))) |>
  pivot_longer(everything(), names_to = "rep_id", values_to = "detected")

# Join metadata, filter to MOSAIC-A and MOSAIC-D, collapse to filter level (OR across PCR replicates)
filter_det <- cet_detected |>
  left_join(sample |> select(rep_id = `...1`, FilterID, Collection.method), by = "rep_id") |>
  filter(Collection.method %in% c("MOSAIC-A", "MOSAIC-D")) |>
  summarise(detected = any(detected), .by = c(FilterID, Collection.method))

# Detection rates by method
filter_det |>
  summarise(n_detected = sum(detected), n_filters = n(),
            pct = 100 * mean(detected), .by = Collection.method)

# Chi-squared: MOSAIC-A vs MOSAIC-D (CTD)
ct <- table(filter_det$Collection.method, filter_det$detected)
chisq.test(ct)
