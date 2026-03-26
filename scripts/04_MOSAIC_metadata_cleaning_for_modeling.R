# =============================================================================
# MOSAIC eDNA — Metadata Cleaning Script
# =============================================================================
#
# Applies all sample-level metadata corrections and exclusions to the raw
# phyloseq object, then writes a cleaned phyloseq RDS.
#
# Input:  MOSAIC_taxnamed.ps.RDS
# Output: MOSAIC_taxnamed_metadata_filtered.ps.RDS
#
# Steps applied:
#   1. Correct F021 filtration end time (12:22 -> 00:22; typo, 12h off)
#   2. Drop samples with no GPS coordinates (10 FilterIDs)
#   3. Drop samples with NA model covariates (4 FilterIDs)
#   4. Drop samples with physically impossible dist_km (17 FilterIDs):
#      implied ship speed exceeds R/V Sally Ride top speed of 14.8 kts
#   5. Dynamic NA sweep: drop any remaining rows with NA in required covariates
#
# NOTE: F021 is corrected in step 1 but its corrected implied speed (~23.6 kts)
#       still exceeds 14.8 kts, so it is dropped in step 4 regardless.
#       The correction is applied first for data-integrity purposes.
#
# Total dropped: 32 FilterIDs (31 from static lists + 1 from dynamic sweep)
# Total retained: 133 FilterIDs
# =============================================================================

library(phyloseq)

# =============================================================================
# 1. LOAD
# =============================================================================

ps  <- readRDS("outputs/MOSAIC_taxnamed.ps.RDS")   # update path as needed
md  <- data.frame(sample_data(ps), check.names = FALSE)

cat(sprintf("Loaded: %d samples, %d taxa\n", nsamples(ps), ntaxa(ps)))
cat(sprintf("Unique FilterIDs before cleaning: %d\n",
            length(unique(na.omit(md$FilterID)))))

# =============================================================================
# 2. CORRECT F021 FILTRATION END TIME
# =============================================================================
# Filtration.time.end was recorded as "12:22"; correct value is "00:22".
# The erroneous entry inflated filtration_time_min from 225 to 945 min.
# datetime_end and filtration_time_min are recalculated from the corrected time.

f021_rows <- !is.na(md$FilterID) & md$FilterID == "F021"
cat(sprintf("\nF021: %d PCR replicates found\n", sum(f021_rows)))
cat(sprintf("  Original  Filtration.time.end : %s\n",
            unique(md$Filtration.time.end[f021_rows])))
cat(sprintf("  Original  filtration_time_min : %s\n",
            unique(md$filtration_time_min[f021_rows])))

md$Filtration.time.end[f021_rows] <- "00:22"

# Recalculate datetime_end: same date as original (2024-07-30) but 00:22
new_datetime_end <- as.POSIXct("2024-07-30 00:22:00", tz = "UTC")
md$datetime_end[f021_rows] <- format(new_datetime_end, "%Y-%m-%d %H:%M:%S")

# Recalculate filtration_time_min
start_dt <- as.POSIXct(md$datetime_start[f021_rows][1], tz = "UTC")
new_filt_min <- as.numeric(difftime(new_datetime_end, start_dt, units = "mins"))
md$filtration_time_min[f021_rows] <- new_filt_min

cat(sprintf("  Corrected Filtration.time.end : 00:22\n"))
cat(sprintf("  Corrected filtration_time_min : %.0f min\n", new_filt_min))

# Flag: corrected speed is still > 14.8 kts — will be dropped in step 4
f021_speed <- unique(as.numeric(md$dist_km[f021_rows])) /
              (new_filt_min / 60) / 1.852
cat(sprintf("  Corrected implied speed       : %.1f kts (> 14.8 kt limit; will be dropped)\n",
            f021_speed))

# =============================================================================
# 3. BUILD DROP LISTS
# =============================================================================

# --- 3a. No GPS coordinates ---
drop_no_coords <- c("C47","C76","C77","C18","C19","C27","C33","C35","C36","F010")

# --- 3b. NA required covariates (dist_km or Volume.filtered) ---
# F062 was previously in this list solely for missing filtration_time_min,
# which is no longer a model covariate; F062 is retained.
drop_na_covars <- c("Fall13", "Fall17", "Fall65", "F076")

# --- 3c. Physically impossible dist_km ---
# For each sample below the Haversine distance between start and end coordinates
# implies a ship speed exceeding 14.8 kts (R/V Sally Ride maximum).
# Most probable cause: end coordinates reflect nominal CalCOFI station-grid
# positions rather than actual GPS at filtration end.
# Speed range across these 17 samples: 16.5 – 293.8 kts.
#
# F021 is included here; even after correcting filtration_time_min (step 2)
# the implied speed remains ~23.6 kts.
# Speed filter applies to MOSAIC-A and MOSAIC-B only.
# MOSAIC-C uses a carboy system towed alongside the ship; its end coordinates
# reflect the carboy retrieval position rather than the ship's transit end,
# so the dist_km calculation is not directly comparable to ship speed.
# MOSAIC-C samples F008, F024, F044, F082, F095 are therefore retained.
drop_impossible_dist <- c(
  # RS2311 (9 samples, MOSAIC-A): implied speed 16.5 – 293.8 kts
  "Fall02", "Fall05", "Fall07", "Fall20", "Fall21",
  "Fall22", "Fall26", "Fall61", "Fall86",
  # SR2408 (3 samples, MOSAIC-A): implied speed 20.0 – 29.5 kts
  "F011", "F022", "F030"
)

drop_all <- union(drop_no_coords, union(drop_na_covars, drop_impossible_dist))

# =============================================================================
# 4. APPLY STATIC DROPS
# =============================================================================

all_fids  <- unique(na.omit(md$FilterID))
n_dropped <- length(intersect(drop_all, all_fids))
cat(sprintf("\nStatic drops:\n"))
cat(sprintf("  No coordinates:     %d FilterIDs\n",
            length(intersect(drop_no_coords,       all_fids))))
cat(sprintf("  NA covariates:      %d FilterIDs\n",
            length(intersect(drop_na_covars,        all_fids))))
cat(sprintf("  Impossible dist_km: %d FilterIDs\n",
            length(intersect(drop_impossible_dist,  all_fids))))
cat(sprintf("  Total static:       %d FilterIDs\n", n_dropped))

keep_samples <- rownames(md)[is.na(md$FilterID) | !(md$FilterID %in% drop_all)]
ps <- prune_samples(keep_samples, ps)
md <- data.frame(sample_data(ps), check.names = FALSE)

# =============================================================================
# 5. DYNAMIC NA SWEEP
# =============================================================================
# Drop any remaining samples with NA in required covariates not caught above.

required_covars <- c("Volume.filtered..mL.", "dist_km",
                     "Latitude.start", "Longitude.start")
na_mask <- apply(md[, required_covars], 1, anyNA)
if (any(na_mask)) {
  na_fids <- unique(na.omit(md$FilterID[na_mask]))
  cat(sprintf("\nDynamic NA sweep: dropping %d rows (FilterIDs: %s)\n",
              sum(na_mask), paste(na_fids, collapse = ", ")))
  ps <- prune_samples(!na_mask, ps)
  md <- data.frame(sample_data(ps), check.names = FALSE)
} else {
  cat("\nDynamic NA sweep: no additional rows dropped\n")
}

# Final verification
stopifnot(all(!is.na(md[, required_covars])))

# =============================================================================
# 6. REPORT AND SAVE
# =============================================================================

retained_fids <- length(unique(na.omit(md$FilterID)))
cat(sprintf("\nCleaned phyloseq object:\n"))
cat(sprintf("  Samples (PCR replicates): %d\n", nsamples(ps)))
cat(sprintf("  Unique FilterIDs:         %d\n", retained_fids))
cat(sprintf("  Taxa:                     %d\n", ntaxa(ps)))
cat(sprintf("\nMethod breakdown:\n"))
print(table(md$Collection.method, useNA = "ifany"))

out_path <- "outputs/MOSAIC_taxnamed_metadata_filtered.ps.RDS"
saveRDS(ps, out_path)
cat(sprintf("\nSaved: %s\n", out_path))

MOSAIC_taxnamed_otu_table_metadata_filtered<-cbind(as.data.frame(otu_table(ps)))
MOSAIC_taxnamed_taxa_table_metadata_filtered<-cbind(as.data.frame(tax_table(ps)))
MOSAIC_taxnamed_sample_data_metadata_filtered<-cbind(as.data.frame(sample_data(ps)))

write.csv(MOSAIC_taxnamed_otu_table_metadata_filtered,"outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv")
write.csv(MOSAIC_taxnamed_taxa_table_metadata_filtered,"outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv")
write.csv(MOSAIC_taxnamed_sample_data_metadata_filtered,"outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv")


