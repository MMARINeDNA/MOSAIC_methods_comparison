# =============================================================================
# Map — target-taxa detections along the 7 A+B transects
# =============================================================================
# Reuses the repo's mapping idiom (scripts/09: maps::map_data("world") coastline,
# CalCOFI extent, transects as start->end geom_segments over an ocean-blue panel).
# Restricted to the 6 concurrent A/B CHAINS (each parent-A transect + its B
# segment series). Each focal FAMILY has its own colour (repo rare-taxa style).
# Every segment is LABELLED "A" (whole transect) or "B" (segment) at its midpoint;
# segments that detected NONE of the target families are drawn as dashed grey.
# A small per-family perpendicular offset keeps co-located detections legible.
# No C, no D.
# Inputs : outputs/MOSAIC_taxnamed_{sample_data,taxa_table,otu_table}_metadata_filtered.csv
# Output : outputs/v5/figures/transect_source/AB_chain_taxa_map.pdf
# =============================================================================

suppressMessages({library(ggplot2); library(maps); library(ggrepel)})
dir.create("outputs/v5/figures/transect_source", recursive = TRUE, showWarnings = FALSE)
A_OFF <-  0.06   # A transects drawn this far NORTH of true position
B_OFF <- -0.06   # B segments drawn this far SOUTH (so A is always on top of B)

md  <- read.csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv",
                row.names = 1, stringsAsFactors = FALSE)
tax <- read.csv("outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv",
                row.names = 1, stringsAsFactors = FALSE)
otu <- read.csv("outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv",
                row.names = 1, check.names = FALSE)
otu <- otu[, rownames(md), drop = FALSE]
otu_filt <- t(rowsum(t(otu), group = md$FilterID))

# Phocoenidae dropped: 0 detections across the AB-chain samples in this map.
focal_families <- c("Delphinidae", "Balaenopteridae",
                    "Otariidae", "Phocidae", "Alopiidae")
fam_col <- setNames(c("#0072B2", "#009E73",
                      "#E69F00", "#56B4E9", "#D55E00"), focal_families)
famdet <- function(fid, fam) {
  ft <- rownames(tax)[!is.na(tax$family) & tax$family == fam]; ft <- ft[ft %in% rownames(otu_filt)]
  if (!length(ft)) return(0L); as.integer(sum(otu_filt[ft, fid]) > 0)
}

# --- identify the A/B chains (7 in SR2408) (parent A + its B segment series) ---
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

rows <- list()
for (ci in seq_along(mult)) {
  ix <- mult[[ci]]; cw0 <- min(fm$t0[ix]); cw1 <- max(fm$t1[ix])
  Ai <- which(fm$Collection.method == "MOSAIC-A" & fm$Cruise == fm$Cruise[ix[1]] &
              pmin(fm$t1, cw1) > pmax(fm$t0, cw0))
  for (m in c(if (length(Ai)) Ai[1], ix)) rows[[length(rows) + 1]] <- data.frame(
    chain = ci, FilterID = fm$FilterID[m],
    AB = ifelse(fm$Collection.method[m] == "MOSAIC-A", "A", "B"),
    lon_start = fm$Longitude.start[m], lat_start = fm$Latitude.start[m],
    lon_end   = fm$Longitude.end[m],   lat_end   = fm$Latitude.end[m])
}
chain_geom <- do.call(rbind, rows)
chain_geom$lon_mid <- (chain_geom$lon_start + chain_geom$lon_end) / 2
chain_geom$lat_mid <- (chain_geom$lat_start + chain_geom$lat_end) / 2
chain_geom$n_target <- vapply(chain_geom$FilterID,
  function(f) sum(vapply(focal_families, function(x) famdet(f, x), integer(1))), integer(1))
stopifnot(!any(fm$Collection.method[match(chain_geom$FilterID, fm$FilterID)] %in% c("MOSAIC-C", "MOSAIC-D")))

no_det <- chain_geom[chain_geom$n_target == 0, ]     # dashed grey
cat(sprintf("%d chains | %d A + %d B segments | %d segments with no target taxa (dashed grey)\n",
            length(mult), sum(chain_geom$AB == "A"), sum(chain_geom$AB == "B"), nrow(no_det)))

# --- detections (long) + per-family perpendicular offset lanes ---
det <- do.call(rbind, lapply(seq_along(focal_families), function(fi) {
  f <- focal_families[fi]
  d <- chain_geom
  d$family <- f
  d <- d[vapply(d$FilterID, famdet, integer(1), fam = f) == 1, ]
  if (!nrow(d)) return(NULL)
  # A drawn NORTH (+lat), B drawn SOUTH (-lat); small per-family stacking within band
  sh <- ifelse(d$AB == "A", A_OFF, B_OFF) + (fi - (length(focal_families) + 1) / 2) * 0.008
  d$lat_start <- d$lat_start + sh; d$lat_end <- d$lat_end + sh
  d
}))
det$family <- factor(det$family, levels = focal_families)

# no-detection segments and A/B labels share the same A-north / B-south banding
no_det$lat_start <- no_det$lat_start + ifelse(no_det$AB == "A", A_OFF, B_OFF)
no_det$lat_end   <- no_det$lat_end   + ifelse(no_det$AB == "A", A_OFF, B_OFF)
chain_geom$lab_lat <- chain_geom$lat_mid + ifelse(chain_geom$AB == "A", A_OFF, B_OFF)

world_df <- map_data("world")
ca_lon <- c(-125, -116); ca_lat <- c(29, 36)

p <- ggplot() +
  geom_polygon(data = world_df, aes(long, lat, group = group),
               fill = "#c8d8b0", colour = "#7a9a5a", linewidth = 0.3) +
  # segments with no target taxa: dashed grey
  geom_segment(data = no_det,
               aes(x = lon_start, y = lat_start, xend = lon_end, yend = lat_end),
               colour = "grey55", linewidth = 0.5, linetype = "dashed") +
  # detections: coloured by family; A north of B
  geom_segment(data = det,
               aes(x = lon_start, y = lat_start, xend = lon_end, yend = lat_end, colour = family),
               linewidth = 1.1, lineend = "round", alpha = 0.95) +
  # every segment labelled A / B (ggrepel, no overlap); A labels north, B south
  ggrepel::geom_text_repel(data = chain_geom, aes(x = lon_mid, y = lab_lat, label = AB),
                           size = 4.2, fontface = "bold", colour = "grey15",
                           max.overlaps = Inf, seed = 42, min.segment.length = 0.15,
                           segment.size = 0.2, segment.colour = "grey60") +
  # place-name labels over the land (green), not the ocean
  annotate("text", x = -116.55, y = 32.55, label = "San Diego",
           hjust = 1, fontface = "italic", size = 4.4, colour = "black") +
  annotate("text", x = -120.30, y = 34.62, label = "Point Conception",
           hjust = 0, fontface = "italic", size = 4.4, colour = "black") +
  scale_colour_manual(name = "Family", values = fam_col, drop = FALSE) +
  coord_fixed(ratio = 1, xlim = ca_lon, ylim = ca_lat, expand = FALSE) +
  labs(title = sprintf("Target-taxa detections along the %d A+B transects (SR2408)", length(chains)),
       subtitle = "A transect drawn north of its B segments; dashed grey = no target taxa; colour = family",
       x = NULL, y = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2.5))) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_rect(fill = "#d6e8f5"),
        panel.grid = element_blank(), legend.position = "right")

ggsave("outputs/v5/figures/transect_source/AB_chain_taxa_map.pdf",
       p, width = 8.5, height = 8, device = grDevices::cairo_pdf, bg = "white")
cat("Saved: outputs/v5/figures/transect_source/AB_chain_taxa_map.pdf\n")

# clean title-less panel for manuscript assembly (Fig 6)
dir.create("outputs/v5/figures/panels", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/v5/figures/panels/panel_q3_abmap.pdf",
       p + labs(title = NULL, subtitle = NULL, caption = NULL),
       width = 6, height = 6.5, device = grDevices::cairo_pdf, bg = "white")
saveRDS(p + labs(title = NULL, subtitle = NULL, caption = NULL),
        "outputs/v5/figures/panels/abmap_obj.rds")   # for assembly
print(table(det$family, det$AB))
