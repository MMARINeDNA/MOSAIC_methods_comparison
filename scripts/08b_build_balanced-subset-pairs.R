# =============================================================================
# v6 Stage-1 matched-subset builder (Part 2 sensitivity)
# =============================================================================
# Pairs each CTD (MOSAIC-D) cast to its nearest MOSAIC-A start, 1:1 with NO reuse,
# WITHIN the same cruise x inshore/offshore stratum (decisions 1 & 2). Optimal
# (min total distance) assignment via clue::solve_LSAP where available, else a
# deterministic greedy nearest fallback. Uses ALL 26 CTDs; no data excluded.
# Input : outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv
# Output: outputs/v6/matched_pairs_v6.csv   (26 CTD-A pairs + distances)
# =============================================================================
dir.create("outputs/v6", recursive = TRUE, showWarnings = FALSE)
md <- read.csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv",
               row.names = 1, stringsAsFactors = FALSE)
fm <- md[!duplicated(md$FilterID), ]

hav <- function(a, b, c, d) { R <- 6371; p <- pi/180
  x <- sin((c-a)*p/2)^2 + cos(a*p)*cos(c*p)*sin((d-b)*p/2)^2; 2*R*asin(pmin(1, sqrt(x))) }

# optimal 1:1 assignment (rows = CTD, cols = A; n_row <= n_col); returns col per row
assign_pairs <- function(cost) {
  if (requireNamespace("clue", quietly = TRUE)) return(as.integer(clue::solve_LSAP(cost)))
  # greedy fallback: repeatedly take the global minimum among unassigned
  out <- integer(nrow(cost)); used <- logical(ncol(cost))
  for (k in seq_len(nrow(cost))) {
    m <- cost; m[out != 0, ] <- Inf; m[, used] <- Inf
    ij <- which(m == min(m), arr.ind = TRUE)[1, ]
    out[ij[1]] <- ij[2]; used[ij[2]] <- TRUE
  }
  out
}

pairs <- list()
for (cr in c("RS2311", "SR2408")) for (zone in c("Inshore", "Offshore")) {
  D <- fm[fm$Collection.method == "MOSAIC-D" & fm$Cruise == cr & fm$Inshore.Offshore == zone, ]
  A <- fm[fm$Collection.method == "MOSAIC-A" & fm$Cruise == cr & fm$Inshore.Offshore == zone, ]
  if (!nrow(D)) next
  stopifnot(nrow(A) >= nrow(D))          # every stratum must have >= CTDs (checked: it does)
  cost <- outer(seq_len(nrow(D)), seq_len(nrow(A)),
                Vectorize(function(i, j) hav(D$Latitude.start[i], D$Longitude.start[i],
                                             A$Latitude.start[j], A$Longitude.start[j])))
  col <- assign_pairs(cost)
  for (i in seq_len(nrow(D))) pairs[[length(pairs) + 1]] <- data.frame(
    cruise = cr, zone = zone,
    CTD = D$FilterID[i], A = A$FilterID[col[i]],
    dist_km = round(cost[i, col[i]], 2),
    CTD_lat = D$Latitude.start[i], CTD_lon = D$Longitude.start[i],
    A_lat = A$Latitude.start[col[i]], A_lon = A$Longitude.start[col[i]],
    stringsAsFactors = FALSE)
}
mp <- do.call(rbind, pairs)
mp <- mp[order(mp$cruise, mp$zone, mp$dist_km), ]
write.csv(mp, "outputs/v6/matched_pairs_v6.csv", row.names = FALSE)

cat(sprintf("Matched %d CTD-A pairs (all 26 CTDs used, no A reused: %s)\n",
            nrow(mp), !any(duplicated(mp$A))))
cat("Pairs per stratum:\n"); print(table(mp$cruise, mp$zone))
cat(sprintf("Match distance (km): min %.2f  median %.2f  max %.2f\n",
            min(mp$dist_km), median(mp$dist_km), max(mp$dist_km)))
cat("\n"); print(mp[, c("cruise","zone","CTD","A","dist_km")], row.names = FALSE)
