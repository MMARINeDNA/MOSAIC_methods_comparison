# =============================================================================
# Q3 — origin-to-distal corroboration of the parent-A detection, by B segment
# =============================================================================
# For each B segment in a chain, mark where along the transect (0 = origin) it
# detected a focal-mammal taxon, and ring the FIRST B segment to corroborate a
# detection the paired parent-A transect also made (same identity, earliest B).
# Two identity resolutions, both faceted by family:
#   (1) TAXON level  -> matches the ordination's taxon identities (species
#       binomials AND genus/family placeholders each as their own row).
#       outputs/v5/figures/transect_source/taxon_corrob_byidentity.pdf
#   (2) SPECIES level -> real binomials only (placeholders excluded, per the
#       project's taxonomic conventions), with the ring also encoding HOW CLOSE
#       the paired A got (species / family-only / missed).
#       outputs/v5/figures/transect_source/species_corrob_byspecies.pdf
# A and B only (no C, no D). Exploratory: 6 chains, SR2408.
# =============================================================================

suppressMessages({library(ggplot2)})
dir.create("outputs/v5/figures/transect_source", recursive = TRUE, showWarnings = FALSE)

md  <- read.csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv", row.names = 1, stringsAsFactors = FALSE)
tax <- read.csv("outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv", row.names = 1, stringsAsFactors = FALSE)
otu <- read.csv("outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv", row.names = 1, check.names = FALSE)
otu <- otu[, rownames(md), drop = FALSE]
of  <- t(rowsum(t(otu), group = md$FilterID))          # taxon x filter (rownames = assigned identity)

focal_families <- c("Delphinidae", "Phocoenidae", "Balaenopteridae", "Otariidae", "Phocidae", "Alopiidae")
is_placeholder <- function(s) grepl("genus|family|\\bsp\\.|\\bspp\\b|environmental|uncultured|mgcode|cf\\.|order|class|phylum",
                                    s, ignore.case = TRUE)
tr        <- function(fid, taxon) of[taxon, fid]                       # reads for a taxon identity
fam_reads <- function(fid, fam) { ft <- rownames(tax)[!is.na(tax$family) & tax$family == fam]; ft <- ft[ft %in% rownames(of)]
  if (!length(ft)) return(0); sum(of[ft, fid]) }
gen_reads <- function(fid, ge) { if (is.na(ge)) return(0); ft <- rownames(tax)[!is.na(tax$genus) & tax$genus == ge]; ft <- ft[ft %in% rownames(of)]
  if (!length(ft)) return(0); sum(of[ft, fid]) }
hav <- function(a, b, c, d) { R <- 6371; p <- pi/180
  x <- sin((c-a)*p/2)^2 + cos(a*p)*cos(c*p)*sin((d-b)*p/2)^2; 2*R*asin(pmin(1, sqrt(x))) }

# --- taxonomic-resolution machinery (drives the ring/square in the taxon panel) --------
# rank order none < family < genus < species; a_resolution() = finest level a filter
# reached within the lineage of identity `t`, so we can compare A vs each B segment.
foc_ids <- rownames(tax)[tax$family %in% focal_families & rownames(tax) %in% rownames(of)]
sp_ids  <- foc_ids[!is_placeholder(foc_ids)]
lvl_of  <- function(t) if (grepl("family", t)) "family" else if (grepl("genus", t)) "genus" else "species"
rank_of <- c(none = 0, family = 1, genus = 2, species = 3)
a_resolution <- function(fid, t) {
  fa <- tax$family[match(t, rownames(tax))]; g <- tax$genus[match(t, rownames(tax))]; lt <- lvl_of(t)
  if (lt == "species") {
    if (tr(fid, t) > 0) return("species")
    if (gen_reads(fid, g) > 0) return("genus")
    if (fam_reads(fid, fa) > 0) return("family"); return("none")
  }
  if (lt == "genus") {
    spg <- sp_ids[tax$genus[match(sp_ids, rownames(tax))] == g]
    if (length(spg) && sum(of[spg, fid]) > 0) return("species")
    if (gen_reads(fid, g) > 0) return("genus")
    if (fam_reads(fid, fa) > 0) return("family"); return("none")
  }
  spf <- sp_ids[tax$family[match(sp_ids, rownames(tax))] == fa]
  if (length(spf) && sum(of[spf, fid]) > 0) return("species")
  gph <- foc_ids[grepl("genus", foc_ids) & tax$family[match(foc_ids, rownames(tax))] == fa]
  if (length(gph) && sum(of[gph, fid]) > 0) return("genus")
  if (fam_reads(fid, fa) > 0) return("family"); return("none")
}

fm <- md[!duplicated(md$FilterID), ]
fm$t0 <- as.numeric(as.POSIXct(fm$datetime_start, tz = "UTC")); fm$t1 <- as.numeric(as.POSIXct(fm$datetime_end, tz = "UTC"))
fm <- fm[order(fm$Cruise, fm$t0), ]
isB <- fm$Collection.method == "MOSAIC-B"; isA <- fm$Collection.method == "MOSAIC-A"; Bidx <- which(isB)
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
grp  <- split(Bidx, pA); grp <- grp[vapply(grp, length, 1L) >= 2]
mult <- lapply(grp, function(ix) ix[order(fm$t0[ix])])

# precompute per-chain: parent A, B filters, segment fractions (origin -> distal)
CH <- lapply(seq_along(mult), function(ci) {
  ix <- mult[[ci]]; cw0 <- min(fm$t0[ix]); cw1 <- max(fm$t1[ix])
  Ai <- which(fm$Collection.method == "MOSAIC-A" & fm$Cruise == fm$Cruise[ix[1]] & pmin(fm$t1, cw1) > pmax(fm$t0, cw0))
  if (!length(Ai)) return(NULL)
  seglen <- sapply(ix, function(b) hav(fm$Latitude.start[b], fm$Longitude.start[b], fm$Latitude.end[b], fm$Longitude.end[b]))
  list(ci = ci, aid = fm$FilterID[Ai[1]], bids = fm$FilterID[ix],
       frac = (cumsum(seglen) - seglen/2) / sum(seglen))
})
CH <- Filter(Negate(is.null), CH)

# --- build detection table for an arbitrary taxon identity set ---
build_det <- function(taxa) {
  rows <- list()
  for (h in CH) for (t in taxa) {
    Adet <- tr(h$aid, t) > 0
    bdet <- which(vapply(h$bids, tr, numeric(1), taxon = t) > 0)
    if (!length(bdet)) next
    a_res <- a_resolution(h$aid, t); lt <- rank_of[[lvl_of(t)]]
    matched <- a_res != "none"        # only count when A also detected the lineage (>= family)
    ge_A  <- matched && lt >= rank_of[[a_res]]   # first B meets-or-beats A's resolution
    gt_A  <- matched && lt >  rank_of[[a_res]]   # first B strictly out-resolves A
    jfirst <- bdet[which.min(h$frac[bdet])]
    for (j in bdet) rows[[length(rows) + 1]] <- data.frame(
      chain = h$ci, seg = j, taxon = t,
      family = tax$family[match(t, rownames(tax))], frac = h$frac[j],
      corrob = Adet, first_corrob = Adet && j == min(bdet),
      ring = ge_A && j == jfirst, hires = gt_A && j == jfirst, stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows); d$family <- factor(d$family, levels = focal_families); d
}

plot_corrob <- function(d, title, sub) {
  ggplot(d, aes(frac, taxon, colour = factor(chain))) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_point(data = subset(d, hires), shape = 0, colour = "black", size = 5.6, stroke = 1) +
    geom_point(data = subset(d, ring),  shape = 1, colour = "black", size = 4.0, stroke = 1) +
    facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = title, subtitle = sub,
         x = "position along transect (0 = origin, 1 = distal)", y = NULL, colour = "transect") +
    theme_bw() + theme(strip.text.y.left = element_text(angle = 0, face = "bold"),
                       strip.placement = "outside")
}

# ============================ (1) TAXON LEVEL ===============================
taxa_all <- rownames(tax)[tax$family %in% focal_families & rownames(tax) %in% rownames(of)]
dt <- build_det(taxa_all)
cat(sprintf("TAXON level: %d B-segment detections across %d taxon identities; %d first-corroborations\n",
            nrow(dt), length(unique(dt$taxon)), sum(dt$first_corrob)))
p_tax <- plot_corrob(dt,
  "B-segment detections along the transect — TAXON-level identities",
  "circle = first B segment matching or beating the integrated A's taxonomic resolution; square = first B that out-resolves A (finer ID than A reached for that lineage)")
ggsave("outputs/v5/figures/transect_source/taxon_corrob_byidentity.pdf", p_tax,
       width = 9, height = 7.5, device = grDevices::cairo_pdf)
cat("Saved: taxon_corrob_byidentity.pdf\n")

# clean title-less panel for manuscript assembly (Fig 5D)
dir.create("outputs/v5/figures/panels", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/v5/figures/panels/panel_fig5d_taxon_corrob.pdf",
       p_tax + labs(title = NULL, subtitle = NULL, caption = NULL),
       width = 6.5, height = 6, device = grDevices::cairo_pdf)
saveRDS(p_tax + labs(title = NULL, subtitle = NULL, caption = NULL),
        "outputs/v5/figures/panels/corrob_obj.rds")   # for assembly

# ============================ (2) SPECIES LEVEL =============================
taxa_sp <- taxa_all[!is_placeholder(taxa_all)]
ds <- build_det(taxa_sp)
# how close did A get for each species detection: species / family-only / missed
ds$A_reached <- mapply(function(ci, sp) {
  aid <- CH[[which(sapply(CH, `[[`, "ci") == ci)]]$aid
  fa  <- tax$family[match(sp, rownames(tax))]
  if (tr(aid, sp) > 0) "species" else if (fam_reads(aid, fa) > 0) "family only" else "missed"
}, ds$chain, ds$taxon)
ds$ring <- ifelse(ds$first_corrob, "species (first B)",
            ifelse(ds$A_reached == "species", "species (later B)",
             ifelse(ds$A_reached == "family only", "A: family only", NA)))
cat(sprintf("SPECIES level: %d detections; A reached -> %s\n", nrow(ds),
            paste(names(table(ds$A_reached)), table(ds$A_reached), sep = ":", collapse = "  ")))

ring_lv <- c("species (first B)", "species (later B)", "A: family only")
rings <- subset(ds, !is.na(ring)); rings$ring <- factor(rings$ring, levels = ring_lv)
p_sp <- ggplot(ds, aes(frac, taxon)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = factor(chain)), size = 2.8, alpha = 0.85) +
  geom_point(data = rings, aes(shape = ring, size = ring),
             colour = "black", fill = NA, stroke = 1) +
  scale_shape_manual(name = "paired A", values = c("species (first B)" = 1,
                     "species (later B)" = 1, "A: family only" = 5), drop = FALSE) +
  scale_size_manual(name = "paired A", values = c("species (first B)" = 5,
                    "species (later B)" = 3.6, "A: family only" = 4.2), drop = FALSE) +
  facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "B-segment detections along the transect — SPECIES level",
       subtitle = "ring = how close paired A got: bold circle = species (first B) | thin circle = species (later B) | diamond = family only | no ring = A missed it",
       x = "position along transect (0 = origin, 1 = distal)", y = NULL, colour = "transect") +
  theme_bw() + theme(strip.text.y.left = element_text(angle = 0, face = "bold"),
                     strip.placement = "outside")
ggsave("outputs/v5/figures/transect_source/species_corrob_byspecies.pdf", p_sp,
       width = 9.5, height = 6, device = grDevices::cairo_pdf)
cat("Saved: species_corrob_byspecies.pdf\n")

file.remove("outputs/v5/figures/transect_source/species_corrob_pooled.pdf")  # superseded by taxon-level
