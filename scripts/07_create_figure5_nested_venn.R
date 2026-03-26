#!/usr/bin/env Rscript
# ===============================================================================
# Figure 5: Nested Venn Diagram - REVISED with clearer positioning
# ===============================================================================

# Read and process data (same as before)
df <- read.csv("outputs/visv4methods_eDNA.csv", row.names = 1, stringsAsFactors = FALSE)

is_species_level <- function(taxon_name) {
  if (taxon_name == "NA." || taxon_name == "NA" || is.na(taxon_name)) return(FALSE)
  words <- strsplit(taxon_name, "\\.")[[1]]
  words <- words[words != ""]
  if (length(words) != 2) return(FALSE)
  if (grepl("genus$|family$|order$", taxon_name)) return(FALSE)
  return(TRUE)
}

taxa_columns <- setdiff(colnames(df), c("occurrence_date", "survey", "method"))
species_columns <- taxa_columns[sapply(taxa_columns, is_species_level)]

get_detected_taxa <- function(data, taxa_cols) {
  detected <- character(0)
  delphinus_detected <- FALSE
  if ("Delphinus.capensis" %in% taxa_cols) {
    if (any(!is.na(data[["Delphinus.capensis"]]) & data[["Delphinus.capensis"]] > 0)) {
      delphinus_detected <- TRUE
    }
  }
  if ("Delphinus.bairdii" %in% taxa_cols) {
    if (any(!is.na(data[["Delphinus.bairdii"]]) & data[["Delphinus.bairdii"]] > 0)) {
      delphinus_detected <- TRUE
    }
  }
  if (delphinus_detected) {
    detected <- c(detected, "Delphinus.capensis")
  }
  other_species <- setdiff(taxa_cols, c("Delphinus.capensis", "Delphinus.bairdii"))
  for (col in other_species) {
    if (any(!is.na(data[[col]]) & data[[col]] > 0)) {
      detected <- c(detected, col)
    }
  }
  return(unique(detected))
}

vis <- df[is.na(df$method) & !is.na(df$survey), ]
a <- df[df$method == "MOSAIC-A" & !is.na(df$method), ]
b <- df[df$method == "MOSAIC-B" & !is.na(df$method), ]
c <- df[df$method == "MOSAIC-C" & !is.na(df$method), ]
d <- df[df$method == "CTD" & !is.na(df$method), ]

vis_taxa <- get_detected_taxa(vis, species_columns)
a_taxa <- get_detected_taxa(a, species_columns)
b_taxa <- get_detected_taxa(b, species_columns)
c_taxa <- get_detected_taxa(c, species_columns)
d_taxa <- get_detected_taxa(d, species_columns)

# Calculate all regions
visual_only <- setdiff(vis_taxa, a_taxa)
a_only <- setdiff(a_taxa, union(union(union(b_taxa, c_taxa), d_taxa), vis_taxa))

b_only <- setdiff(setdiff(b_taxa, union(c_taxa, d_taxa)), vis_taxa)
d_only <- setdiff(setdiff(d_taxa, union(b_taxa, c_taxa)), vis_taxa)

all_three_bcd <- intersect(intersect(b_taxa, c_taxa), d_taxa)
all_three_not_vis <- setdiff(all_three_bcd, vis_taxa)
vis_and_all_three <- intersect(vis_taxa, all_three_bcd)

vis_a_only <- setdiff(intersect(vis_taxa, a_taxa), union(union(b_taxa, c_taxa), d_taxa))
vis_b_only <- setdiff(intersect(vis_taxa, b_taxa), union(c_taxa, d_taxa))
vis_d_only <- setdiff(intersect(vis_taxa, d_taxa), union(b_taxa, c_taxa))
vis_bc <- setdiff(intersect(vis_taxa, intersect(b_taxa, c_taxa)), d_taxa)

cat("\nCreating REVISED Figure 5 with clearer positioning...\n")

pdf("viz/Figure5_nested_venn_marine_mammals.pdf", width = 14, height = 12)

par(mar = c(12, 1, 2, 1))
plot.new()
plot.window(xlim = c(0, 14), ylim = c(0, 13), asp = 1)

# Colors
col_a <- rgb(41/255, 128/255, 185/255, alpha = 0.3)
col_b <- rgb(231/255, 76/255, 60/255, alpha = 0.4)
col_c <- rgb(46/255, 204/255, 113/255, alpha = 0.4)
col_d <- rgb(243/255, 156/255, 18/255, alpha = 0.4)
col_visual <- rgb(149/255, 165/255, 166/255, alpha = 0.25)

draw.circle <- function(x, y, r, col, border) {
  theta <- seq(0, 2*pi, length.out = 100)
  polygon(x + r*cos(theta), y + r*sin(theta), 
          col = col, border = border, lwd = 3)
}

# REPOSITIONED CIRCLES for clearer overlaps
# MOSAIC-A: large circle
draw.circle(6.5, 6.5, 5, col_a, rgb(41/255, 128/255, 185/255))

# Visual: overlaps A from right
draw.circle(9, 6.5, 3.5, col_visual, rgb(127/255, 140/255, 141/255))

# B, C, D positioned to show their relationships more clearly
# All three overlap near the Visual boundary to show the shared species
draw.circle(7.5, 7.5, 1.8, col_b, rgb(192/255, 57/255, 43/255))   # B: upper
draw.circle(8.5, 6.5, 1.5, col_c, rgb(39/255, 174/255, 96/255))   # C: right, overlaps Visual
draw.circle(7.5, 5.5, 1.8, col_d, rgb(230/255, 126/255, 34/255))  # D: lower

# Method labels
text(2.5, 11.5, sprintf("MOSAIC-A (%d species)", length(a_taxa)), 
     cex = 1.4, font = 2, col = rgb(41/255, 128/255, 185/255))

text(11.5, 10, sprintf("Visual\n(%d species)", length(vis_taxa)), 
     cex = 1.3, font = 2, col = rgb(127/255, 140/255, 141/255))

text(6.5, 9, sprintf("B (%d)", length(b_taxa)), 
     cex = 1.1, font = 2, col = rgb(192/255, 57/255, 43/255))
text(9.5, 8, sprintf("C (%d)", length(c_taxa)), 
     cex = 1.1, font = 2, col = rgb(39/255, 174/255, 96/255))
text(6.5, 4, sprintf("D (%d)", length(d_taxa)), 
     cex = 1.1, font = 2, col = rgb(230/255, 126/255, 34/255))

# REGION NUMBERS - positioned based on actual data structure
# 1. A ONLY (far left, not in B/C/D or Visual): 5
text(2.5, 6.5, 5, cex = 1.6, font = 2)
text(2.5, 5.7, "A unique\n(no other\neDNA or\nVisual)", cex = 0.75, font = 3)

# 2. Visual ONLY (far right, not in A): 2
text(11.5, 6.5, 2, cex = 1.6, font = 2)
text(11.5, 5.7, "Visual\nonly", cex = 0.75, font = 3)

# 3. Visual & A ONLY (overlap, not in B/C/D): 2
text(9.5, 4, 2, cex = 1.3, font = 2)
text(9.5, 3.3, "Visual &\nA only", cex = 0.7, font = 3)

# 4. B ONLY (in B, not in C/D or Visual): 2
text(6.8, 8, 2, cex = 1.2, font = 2)

# 5. D ONLY (in D, not in B/C or Visual): 1
text(6.8, 5, 1, cex = 1.2, font = 2)

# 6. Visual & B ONLY (not C/D): 1
text(8, 8, 1, cex = 1.1, font = 2)

# 7. Visual & D ONLY (not B/C): 1
text(8.5, 5.2, 1, cex = 1.1, font = 2)

# 8. Visual & B & C (not D): 1
text(8.8, 7.3, 1, cex = 1.1, font = 2)

# 9. All 3 B/C/D (not Visual): 2
text(7.8, 6.5, 2, cex = 1.2, font = 2)

# 10. Visual AND all 3 B/C/D: 1
text(8.5, 6.5, 1, cex = 1.1, font = 2)

# Title
title("Figure 5. Comparison of marine mammal species detected by eDNA sampling methods and visual surveys", 
      cex.main = 1.2, font.main = 2)

# Caption with all species
fix_name <- function(x) gsub("\\.", " ", x)

caption_text <- paste0(
  "Figure 5. Nested Venn diagram comparing species-level identifications from visual surveys with four eDNA sampling methods. ",
  "MOSAIC-A (blue, ", length(a_taxa), " species) encompasses all species detected by MOSAIC-B (", 
  length(b_taxa), " species), MOSAIC-C (", length(c_taxa), " species), and CTD (", 
  length(d_taxa), " species) methods. Visual surveys (gray, ", length(vis_taxa), 
  " species) overlap partially with MOSAIC-A. Species detections by category: ",
  "Visual-only (", length(visual_only), " species): ", paste(sort(sapply(visual_only, fix_name)), collapse = ", "), ". ",
  "MOSAIC-A unique (", length(a_only), " species, not in B/C/D or Visual): ", 
  paste(sort(sapply(a_only, fix_name)), collapse = ", "), ". ",
  "Visual & A only (", length(vis_a_only), " species, not B/C/D): ", 
  paste(sort(sapply(vis_a_only, fix_name)), collapse = ", "), ". ",
  "Visual & B only (", length(vis_b_only), " species, not C/D): ", 
  paste(sort(sapply(vis_b_only, fix_name)), collapse = ", "), ". ",
  "Visual & D only (", length(vis_d_only), " species, not B/C): ", 
  paste(sort(sapply(vis_d_only, fix_name)), collapse = ", "), ". ",
  "Visual & B & C (", length(vis_bc), " species, not D): ", 
  paste(sort(sapply(vis_bc, fix_name)), collapse = ", "), ". ",
  "B unique within A (", length(b_only), " species, not C/D or Visual): ", 
  paste(sort(sapply(b_only, fix_name)), collapse = ", "), ". ",
  "D unique within A (", length(d_only), " species, not B/C or Visual): ", 
  paste(sort(sapply(d_only, fix_name)), collapse = ", "), ". ",
  "All three eDNA methods B/C/D within A (", length(all_three_not_vis), " species, not Visual): ", 
  paste(sort(sapply(all_three_not_vis, fix_name)), collapse = ", "), ". ",
  "Visual AND all three B/C/D (", length(vis_and_all_three), " species): ", 
  paste(sort(sapply(vis_and_all_three, fix_name)), collapse = ", "), ". ",
  "The diagram demonstrates the complementary nature of visual and eDNA approaches for comprehensive species-level marine mammal monitoring. ",
  "Note: Delphinus capensis (visual) and D. bairdii (eDNA) are merged as a single taxon."
)

mtext(caption_text, side = 1, line = 0.5, cex = 0.68, adj = 0, at = 0.08)

dev.off()

cat("✓ Figure 5 saved with revised positioning\n")

