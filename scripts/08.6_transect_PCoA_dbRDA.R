# =============================================================================
# Q3 ordination — PCoA of the SR2408 A+B community (Bray-Curtis)
# =============================================================================
# Figure: PCoA (metric MDS) on Bray-Curtis of the SR2408 A+B community — the
#     metric, interpretable-axes ordination; each chain's B segments connected
#     in order, parent A marked.
# Also computed and reported to the console (NOT plotted): partial distance-based
#     RDA on the B SEGMENTS only, community ~ distance_from_origin +
#     Condition(chain) — does composition change with position ALONG the transect
#     after removing between-transect (chain / parent-A) differences? "chain" is
#     the parent-A grouping (the condition). Kept as a stat, dropped from the figure.
# Read depth controlled by relative abundance. Exploratory: 20 B segments, 7
# chains (SR2408, offshore-heavy).
# Inputs : outputs/MOSAIC_taxnamed_{sample_data,otu_table}_metadata_filtered.csv
# Output : outputs/v5/figures/transect_source/transect_PCoA_ordination.pdf
# =============================================================================

suppressMessages({library(vegan); library(ggplot2); library(patchwork)})
dir.create("outputs/v5/figures/transect_source", recursive = TRUE, showWarnings = FALSE)
set.seed(42)

md  <- read.csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv",
                row.names = 1, stringsAsFactors = FALSE)
otu <- read.csv("outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv",
                row.names = 1, check.names = FALSE)
otu <- otu[, rownames(md), drop = FALSE]
otu_filt <- t(rowsum(t(otu), group = md$FilterID))
relab    <- sweep(otu_filt, 2, pmax(colSums(otu_filt), 1), "/")
hav <- function(a, b, c, d) { R <- 6371; p <- pi/180
  x <- sin((c-a)*p/2)^2 + cos(a*p)*cos(c*p)*sin((d-b)*p/2)^2; 2*R*asin(pmin(1, sqrt(x))) }

fm <- md[!duplicated(md$FilterID), ]
fm$t0 <- as.numeric(as.POSIXct(fm$datetime_start, tz = "UTC"))
fm$t1 <- as.numeric(as.POSIXct(fm$datetime_end,   tz = "UTC"))
fm <- fm[order(fm$Cruise, fm$t0), ]
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

# B-segment table: distance from origin + chain (= parent-A grouping)
segrows <- list(); pa_rows <- list()
for (ci in seq_along(mult)) {
  ix <- mult[[ci]]; cw0 <- min(fm$t0[ix]); cw1 <- max(fm$t1[ix])
  Ai <- which(fm$Collection.method == "MOSAIC-A" & fm$Cruise == fm$Cruise[ix[1]] &
              pmin(fm$t1, cw1) > pmax(fm$t0, cw0))
  if (length(Ai)) pa_rows[[length(pa_rows) + 1]] <-
    data.frame(chain = factor(ci), parentA = fm$FilterID[Ai[1]])
  cum <- 0
  for (j in seq_along(ix)) {
    b <- ix[j]
    seglen <- hav(fm$Latitude.start[b], fm$Longitude.start[b], fm$Latitude.end[b], fm$Longitude.end[b])
    segrows[[length(segrows) + 1]] <- data.frame(
      filter = fm$FilterID[b], chain = factor(ci), seg = j, dist_km = cum + seglen / 2)
    cum <- cum + seglen
  }
}
segdf  <- do.call(rbind, segrows)
pa_map <- do.call(rbind, pa_rows)      # chain -> its parent-A transect

# ---------------------------------------------------------------------------
# (a) PCoA of SR2408 A + B (Bray-Curtis)
# ---------------------------------------------------------------------------
ab   <- fm$FilterID[fm$Cruise == "SR2408" & fm$Collection.method %in% c("MOSAIC-A", "MOSAIC-B")]
comm <- t(relab[, ab, drop = FALSE]); comm <- comm[, colSums(comm) > 0, drop = FALSE]
D    <- vegdist(comm, "bray")
pco  <- cmdscale(D, k = 2, eig = TRUE)
ve   <- round(100 * pco$eig[1:2] / sum(pco$eig[pco$eig > 0]), 1)
sc   <- data.frame(PCo1 = pco$points[, 1], PCo2 = pco$points[, 2], filter = ab)
sc$chain <- segdf$chain[match(sc$filter, segdf$filter)]
sc$seg   <- segdf$seg[match(sc$filter, segdf$filter)]
sc$isB   <- sc$filter %in% segdf$filter
# assign each PARENT A its chain so it plots (as a triangle) with its own B series
sc$chain[match(pa_map$parentA, sc$filter)] <- pa_map$chain
sc$isParentA <- sc$filter %in% pa_map$parentA
sc$isA <- fm$Collection.method[match(sc$filter, fm$FilterID)] == "MOSAIC-A"
seg_pts <- sc[sc$isB, ]; seg_pts <- seg_pts[order(seg_pts$chain, seg_pts$seg), ]
pa_pts  <- sc[sc$isParentA, ]
# shape encodes MOSAIC method (A = parent transect vs B = segment) -> its own legend
seg_pts$role <- "B (segment)"; pa_pts$role <- "A (parent transect)"

pa <- ggplot(sc, aes(PCo1, PCo2)) +
  geom_point(data = subset(sc, is.na(chain) & !isA), colour = "grey80", size = 1.4) +
  geom_point(data = subset(sc, isA & !isParentA), shape = 17, colour = "grey55", size = 2) +
  geom_path(data = seg_pts, aes(group = chain, colour = chain),
            arrow = grid::arrow(length = grid::unit(0.1, "cm")), alpha = 0.85) +
  geom_point(data = seg_pts, aes(colour = chain, shape = role), size = 5.4) +
  geom_text(data = seg_pts, aes(label = seg), size = 4.0, fontface = "bold", colour = "white") +
  # parent A of each chain: large triangle in the chain colour, labelled A
  geom_point(data = pa_pts, aes(colour = chain, shape = role), size = 5.6) +
  geom_text(data = pa_pts, aes(label = "A", colour = chain),
            vjust = -1.1, size = 4.2, fontface = "bold", show.legend = FALSE) +
  scale_shape_manual(name = "MOSAIC method",
                     values = c(`A (parent transect)` = 17, `B (segment)` = 16)) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)),
         shape  = guide_legend(order = 2, override.aes = list(colour = "black"))) +
  labs(title = "PCoA (Bray-Curtis, SR2408 A+B)",
       subtitle = sprintf("triangle = parent A; numbered dots = its B segments; axes = %.1f%% / %.1f%% variation", ve[1], ve[2]),
       caption = "grey triangles = other (non-chain) A samples; grey dots = other B",
       x = sprintf("PCo1 (%.1f%%)", ve[1]), y = sprintf("PCo2 (%.1f%%)", ve[2]),
       colour = "transect") + theme_bw()

# ---------------------------------------------------------------------------
# (b) partial db-RDA: B segments ~ dist_km + Condition(chain)
# ---------------------------------------------------------------------------
bcomm <- t(relab[, segdf$filter, drop = FALSE]); bcomm <- bcomm[, colSums(bcomm) > 0, drop = FALSE]
dbr   <- dbrda(bcomm ~ dist_km + Condition(chain), data = segdf, distance = "bray")
ant   <- anova(dbr, permutations = 999)
prop_constr <- 100 * dbr$CCA$tot.chi / dbr$tot.chi
Fval <- ant$F[1]; pval <- ant$`Pr(>F)`[1]
cat(sprintf("db-RDA (B segments ~ distance + Condition(chain)): constrained = %.1f%% of variation | F = %.2f, p = %.3f\n",
            prop_constr, Fval, pval))

ggsave("outputs/v5/figures/transect_source/transect_PCoA_ordination.pdf",
       pa, width = 8, height = 6, device = grDevices::cairo_pdf)
cat("Saved: outputs/v5/figures/transect_source/transect_PCoA_ordination.pdf\n")

# clean title-less panel for manuscript assembly (Fig 6)
dir.create("outputs/v5/figures/panels", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/v5/figures/panels/panel_q3_pcoa.pdf",
       pa + labs(title = NULL, subtitle = NULL, caption = NULL),
       width = 5.5, height = 4.5, device = grDevices::cairo_pdf)
saveRDS(pa + labs(title = NULL, subtitle = NULL, caption = NULL),
        "outputs/v5/figures/panels/pcoa_obj.rds")   # for assembly
