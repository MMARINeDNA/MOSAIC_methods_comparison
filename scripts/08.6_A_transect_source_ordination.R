# =============================================================================
# Q3 — Where along the A transect does the eDNA come from?
# =============================================================================
# The concurrent B-chain tiles each parent-A transect into consecutive segments
# (seg1 = origin end). Using the WHOLE community (all 365 taxa) we ask where the
# integrated A sample's composition originates, then contrast that with the rare
# marine-mammal detections specifically.
#   Figure: whole-community (all-taxa) relative-abundance Bray-Curtis
#       dissimilarity between each B segment and its parent-A filter, vs the
#       segment's distance from the transect origin. Computed as a direct
#       pairwise vegdist on the relative-abundance vectors (NOT from the NMDS or
#       PCoA ordinations).
#   (The earlier NMDS panel and the per-family detection-position panel were
#    dropped from the figure; the family-position data is still written to CSV
#    and reported to the console below.)
# ONLY MOSAIC-A and MOSAIC-B are used here. MOSAIC-C (carboy) is excluded from
# this question entirely — it is not part of an A/B transect series.
# Exploratory: 7 multi-segment chains (SR2408, offshore-heavy). Read depth is
# controlled by relative abundance.
# Inputs : outputs/MOSAIC_taxnamed_{sample_data,taxa_table,otu_table}_metadata_filtered.csv
# Outputs: outputs/v5/figures/transect_source/B_series_vs_A_braycurtis_dissimilarity.pdf
#          outputs/v5/transect_source_A_vs_segments.csv
#          outputs/v5/transect_source_family_positions.csv
# =============================================================================

suppressMessages({library(vegan); library(ggplot2); library(patchwork)})
dir.create("outputs/v5/figures/transect_source", recursive = TRUE, showWarnings = FALSE)
set.seed(42)

md  <- read.csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv",
                row.names = 1, stringsAsFactors = FALSE)
tax <- read.csv("outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv",
                row.names = 1, stringsAsFactors = FALSE)
otu <- read.csv("outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv",
                row.names = 1, check.names = FALSE)
otu <- otu[, rownames(md), drop = FALSE]

otu_filt <- t(rowsum(t(otu), group = md$FilterID))            # taxa x filters
relab    <- sweep(otu_filt, 2, pmax(colSums(otu_filt), 1), "/")

focal_families <- c("Delphinidae", "Phocoenidae", "Balaenopteridae",
                    "Otariidae", "Phocidae", "Alopiidae")
fam_col <- setNames(c("#0072B2", "#CC79A7", "#009E73",
                      "#E69F00", "#56B4E9", "#D55E00"), focal_families)
famreads <- function(fid, fam) {
  ft <- rownames(tax)[!is.na(tax$family) & tax$family == fam]; ft <- ft[ft %in% rownames(otu_filt)]
  if (!length(ft)) return(0); sum(otu_filt[ft, fid])
}
hav <- function(a, b, c, d) { R <- 6371; p <- pi/180
  x <- sin((c-a)*p/2)^2 + cos(a*p)*cos(c*p)*sin((d-b)*p/2)^2; 2*R*asin(pmin(1, sqrt(x))) }

fm <- md[!duplicated(md$FilterID), ]
fm$t0 <- as.numeric(as.POSIXct(fm$datetime_start, tz = "UTC"))
fm$t1 <- as.numeric(as.POSIXct(fm$datetime_end,   tz = "UTC"))
fm <- fm[order(fm$Cruise, fm$t0), ]

# --- B-chains and their parent A (temporal overlap); A and B ONLY (never C) ---
isB <- fm$Collection.method == "MOSAIC-B"
isA <- fm$Collection.method == "MOSAIC-A"
Bidx <- which(isB)
# Group each MOSAIC-B segment with its parent MOSAIC-A (same cruise, greatest temporal
# overlap), keeping parent-A transects with >= 2 B segments. Grouping by parent-A rather
# than by runs of consecutive B rows is robust to a parent-A that shares a start time with
# its first B segment (which otherwise interleaves the A between two B's, splits the chain,
# and drops segments; e.g. parent-A F005 and F065's first segment F064).
pA <- vapply(Bidx, function(b) {
  cand <- which(isA & fm$Cruise == fm$Cruise[b] &
                pmin(fm$t1, fm$t1[b]) > pmax(fm$t0, fm$t0[b]))
  if (!length(cand)) return(NA_character_)
  fm$FilterID[cand[which.max(pmin(fm$t1[cand], fm$t1[b]) - pmax(fm$t0[cand], fm$t0[b]))]]
}, character(1))
grp  <- split(Bidx, pA)
grp  <- grp[vapply(grp, length, 1L) >= 2]
mult <- lapply(grp, function(ix) ix[order(fm$t0[ix])])
chains <- mult

seg_rows <- list(); fam_rows <- list()
for (ci in seq_along(mult)) {
  ix <- mult[[ci]]; cw0 <- min(fm$t0[ix]); cw1 <- max(fm$t1[ix])
  Ai <- which(fm$Collection.method == "MOSAIC-A" & fm$Cruise == fm$Cruise[ix[1]] &
              pmin(fm$t1, cw1) > pmax(fm$t0, cw0))
  if (!length(Ai)) next
  aid <- fm$FilterID[Ai[1]]
  total <- sum(sapply(ix, function(b) hav(fm$Latitude.start[b], fm$Longitude.start[b],
                                          fm$Latitude.end[b],   fm$Longitude.end[b])))
  cum <- 0
  for (j in seq_along(ix)) {
    b <- ix[j]
    seglen <- hav(fm$Latitude.start[b], fm$Longitude.start[b],
                  fm$Latitude.end[b],   fm$Longitude.end[b])
    seg_rows[[length(seg_rows) + 1]] <- data.frame(
      chain = ci, seg = j, filter = fm$FilterID[b], parentA = aid,
      dist_mid_km = cum + seglen / 2)
    for (fam in focal_families) if (famreads(fm$FilterID[b], fam) > 0)
      fam_rows[[length(fam_rows) + 1]] <- data.frame(
        chain = ci, family = fam, seg = j, frac = (cum + seglen / 2) / total,
        # was this B-segment detection also seen in the paired parent-A transect?
        corroborated = famreads(aid, fam) > 0)
    cum <- cum + seglen
  }
}
segdf <- do.call(rbind, seg_rows)
rr    <- do.call(rbind, fam_rows)
rr$family <- factor(rr$family, levels = focal_families)

# For each series x family where the paired A detected the taxon (corroborated),
# flag the FIRST (earliest / origin-most) B segment that corroborated it.
rr$first_corrob <- FALSE
for (g in split(seq_len(nrow(rr)), list(rr$chain, rr$family), drop = TRUE)) {
  if (rr$corroborated[g][1]) rr$first_corrob[g[which.min(rr$seg[g])]] <- TRUE
}

# --- (b) whole-community Bray-Curtis: parent A vs each B segment ---
segdf$bc_A <- NA_real_
for (r in seq_len(nrow(segdf))) {
  segdf$bc_A[r] <- as.numeric(vegdist(t(relab[, c(segdf$parentA[r], segdf$filter[r])]), "bray"))
}
write.csv(segdf, "outputs/v5/transect_source_A_vs_segments.csv", row.names = FALSE)
write.csv(rr,    "outputs/v5/transect_source_family_positions.csv", row.names = FALSE)

closest <- tapply(seq_len(nrow(segdf)), segdf$chain,
                  function(idx) segdf$seg[idx][which.min(segdf$bc_A[idx])])
rho <- suppressWarnings(cor(segdf$dist_mid_km, segdf$bc_A, method = "spearman"))
cat(sprintf("Whole community: A closest to origin seg in %d/%d chains; Spearman(dist,BC-to-A)=%.2f\n",
            sum(closest == 1), length(closest), rho))
# report all six focal families, and (for contrast) the non-Delphinidae subset
non_d <- rr[rr$family != "Delphinidae", ]
cat(sprintf("Focal-family detections in chains: n=%d, mean position=%.2f\n",
            nrow(rr), mean(rr$frac)))
cat(sprintf("  Delphinidae: n=%d mean=%.2f (near-ubiquitous) | non-Delphinidae: n=%d mean=%.2f, distal third=%d/%d\n",
            sum(rr$family == "Delphinidae"), mean(rr$frac[rr$family == "Delphinidae"]),
            nrow(non_d), mean(non_d$frac), sum(non_d$frac > 0.67), nrow(non_d)))

bc_plot <- ggplot(segdf, aes(dist_mid_km, bc_A, colour = factor(chain))) +
  geom_line(alpha = 0.6) + geom_point(size = 2.4) +
  labs(title = "B-series vs parent-A Bray-Curtis dissimilarity",
       subtitle = sprintf("whole-community relative-abundance (all taxa); A closest to origin segment in %d/%d chains, Spearman(distance, BC) = %.2f",
                          sum(closest == 1), length(closest), rho),
       x = "B segment distance from transect origin (km)",
       y = "Bray-Curtis dissimilarity\n(B segment vs parent A)", colour = "transect") + theme_bw()

ggsave("outputs/v5/figures/transect_source/B_series_vs_A_braycurtis_dissimilarity.pdf",
       bc_plot, width = 8, height = 5, device = grDevices::cairo_pdf)
cat("Saved: outputs/v5/figures/transect_source/B_series_vs_A_braycurtis_dissimilarity.pdf\n")

# clean title-less panel for manuscript assembly (Fig 6)
dir.create("outputs/v5/figures/panels", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/v5/figures/panels/panel_q3_bc.pdf",
       bc_plot + labs(title = NULL, subtitle = NULL, caption = NULL),
       width = 5, height = 4, device = grDevices::cairo_pdf)
saveRDS(bc_plot + labs(title = NULL, subtitle = NULL, caption = NULL),
        "outputs/v5/figures/panels/bc_obj.rds")   # for assembly
