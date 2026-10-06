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
#   2. Drop MOSAIC-A/B samples whose implied ship speed
#      (dist_km / filtration time) exceeds the R/V Sally Ride top speed of 14.8 kts
#   3. Drop any samples with NA in the required covariates
#
# Revision round 2: start/end locations (and therefore dist_km) now come from
# the corrected-locations workbook (applied in 02_prep_detection_prob_metadata.R).
# The hard-coded drop lists used in the original submission (10 "no GPS" FilterIDs,
# 4 "NA covariate" FilterIDs, 12 "impossible dist_km" FilterIDs; this header said 32 dropped /
# 133 retained, but the committed output had 138) existed only because of the uncorrected locations, so exclusions
# are now applied as rules to the corrected data. The old lists are kept below
# (legacy_static_drops) purely to report which of those samples are now retained.
# =============================================================================
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

# Corrected implied speed (informational): with corrected locations F021 is no
# longer a speed outlier, so it is retained
f021_speed <- unique(as.numeric(md$dist_km[f021_rows])) /
              (new_filt_min / 60) / 1.852
cat(sprintf("  Corrected implied speed       : %.1f kts\n", f021_speed))

# Write the corrected F021 values back into the phyloseq object (previously they
# were only applied to the md copy, which was fine while F021 was always dropped)
sample_data(ps) <- sample_data(md)

# =============================================================================
# 3. SPEED-BASED EXCLUSION (corrected locations)
# =============================================================================
# Implied ship speed = dist_km / hours filtered. A value above the R/V Sally Ride
# maximum (14.8 kts) is physically impossible and indicates bad end coordinates.
# Speed filter applies to MOSAIC-A and MOSAIC-B only: MOSAIC-C uses a carboy
# system (water collected, homogenized, then a quick 5 L filtration), so
# dist_km / time filtered is not comparable to ship speed; CTD (MOSAIC-D) has
# dist_km = 0 by definition.
max_speed_kts <- 14.8
filt_min      <- as.numeric(md$filtration_time_min)
filt_min[!is.na(filt_min) & filt_min <= 0] <- NA   # 0 min = missing, not instantaneous
speed_kts     <- as.numeric(md$dist_km) / (filt_min / 60) / 1.852
too_fast      <- md$Collection.method %in% c("MOSAIC-A", "MOSAIC-B") &
                 !is.na(speed_kts) & speed_kts > max_speed_kts
drop_impossible_speed <- unique(na.omit(md$FilterID[too_fast]))

cat(sprintf("\nSpeed rule (> %.1f kts, MOSAIC-A/B): %d FilterIDs dropped%s\n",
            max_speed_kts, length(drop_impossible_speed),
            if (length(drop_impossible_speed))
              paste0(" (", paste(drop_impossible_speed, collapse = ", "), ")") else ""))

# Original-submission static drops, kept only for the "now retained" report below
legacy_static_drops <- c(
  "C47","C76","C77","C18","C19","C27","C33","C35","C36","F010",   # no GPS
  "Fall13", "Fall17", "Fall65", "F076",                           # NA covariates
  "Fall02", "Fall05", "Fall07", "Fall20", "Fall21", "Fall22",     # impossible dist_km
  "Fall26", "Fall61", "Fall86", "F011", "F022", "F030",
  "F021")

# =============================================================================
# 4. APPLY SPEED EXCLUSION
# =============================================================================

keep_samples <- rownames(md)[is.na(md$FilterID) | !(md$FilterID %in% drop_impossible_speed)]
ps <- prune_samples(keep_samples, ps)
md <- data.frame(sample_data(ps), check.names = FALSE)

# =============================================================================
# 5. DYNAMIC NA SWEEP
# =============================================================================
# Drop any remaining samples with NA in required covariates not caught above.

# datetime_end / filtration_time_min are required by the models (temporal-overlap
# site structure and the filtration-time covariate); F062 has no filtration end
# time ("na") and previously dropped out via its missing end coordinates instead.
required_covars <- c("Volume.filtered..mL.", "dist_km",
                     "Latitude.start", "Longitude.start",
                     "datetime_end", "filtration_time_min")
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

legacy_now_retained <- intersect(legacy_static_drops, unique(na.omit(md$FilterID)))
cat(sprintf("\nOriginal-submission static drops now retained (%d of %d): %s\n",
            length(legacy_now_retained), length(legacy_static_drops),
            paste(legacy_now_retained, collapse = ", ")))
cat(sprintf("Original-submission static drops still excluded: %s\n",
            paste(setdiff(legacy_static_drops, legacy_now_retained), collapse = ", ")))

out_path <- "outputs/MOSAIC_taxnamed_metadata_filtered.ps.RDS"
saveRDS(ps, out_path)
cat(sprintf("\nSaved: %s\n", out_path))

MOSAIC_taxnamed_otu_table_metadata_filtered<-cbind(as.data.frame(otu_table(ps)))
MOSAIC_taxnamed_taxa_table_metadata_filtered<-cbind(as.data.frame(tax_table(ps)))
MOSAIC_taxnamed_sample_data_metadata_filtered<-cbind(as.data.frame(sample_data(ps)))

write.csv(MOSAIC_taxnamed_otu_table_metadata_filtered,"outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv")
write.csv(MOSAIC_taxnamed_taxa_table_metadata_filtered,"outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv")
write.csv(MOSAIC_taxnamed_sample_data_metadata_filtered,"outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv")


