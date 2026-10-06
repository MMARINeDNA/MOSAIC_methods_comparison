# =============================================================================
# Total distance sampled per cruise, counting overlapping samples once
# =============================================================================
# MOSAIC-B and -C filters are collected concurrently with a parent MOSAIC-A
# transect (and A filters can overlap each other in time), so summing dist_km
# over all surface filters double-counts water. Here each surface filter's
# distance is credited only for the fraction of its filtration window not already
# covered by a filter processed earlier, with MOSAIC-A taken first (it is the
# full-length transect), then B, then C. CTD casts are stationary (dist_km = 0).
# Input : outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv (retained filters)
# Output: outputs/total_distance_sampled.csv
# =============================================================================
md <- read.csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv",
               row.names = 1, stringsAsFactors = FALSE)
fm <- md[!duplicated(md$FilterID) & md$Collection.method %in% c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C"), ]
fm$t0 <- as.numeric(as.POSIXct(fm$datetime_start, tz = "UTC"))
fm$t1 <- as.numeric(as.POSIXct(fm$datetime_end,   tz = "UTC"))
stopifnot(!anyNA(fm$t0), !anyNA(fm$t1), all(fm$t1 > fm$t0), !anyNA(fm$dist_km))
fm$prio <- match(fm$Collection.method, c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C"))

# length of [a, b] not covered by the union of `cov` (a 2-column matrix of intervals)
uncovered <- function(a, b, cov) {
  if (!nrow(cov)) return(b - a)
  cov <- cov[order(cov[, 1]), , drop = FALSE]
  merged <- cov[1, , drop = FALSE]
  for (i in seq_len(nrow(cov))[-1]) {
    k <- nrow(merged)
    if (cov[i, 1] <= merged[k, 2]) merged[k, 2] <- max(merged[k, 2], cov[i, 2])
    else merged <- rbind(merged, cov[i, ])
  }
  covered <- sum(pmax(0, pmin(b, merged[, 2]) - pmax(a, merged[, 1])))
  (b - a) - covered
}

res <- do.call(rbind, lapply(split(fm, fm$Cruise), function(d) {
  d <- d[order(d$prio, d$t0), ]
  cov <- matrix(numeric(0), 0, 2); credited <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    credited[i] <- d$dist_km[i] * uncovered(d$t0[i], d$t1[i], cov) / (d$t1[i] - d$t0[i])
    cov <- rbind(cov, c(d$t0[i], d$t1[i]))
  }
  data.frame(Cruise = d$Cruise[1],
             n_filters = nrow(d),
             sum_all_surface_km  = sum(d$dist_km),
             sum_A_only_km       = sum(d$dist_km[d$Collection.method == "MOSAIC-A"]),
             nonoverlapping_km   = sum(credited))
}))
res <- rbind(res, data.frame(Cruise = "Both", n_filters = sum(res$n_filters),
                             sum_all_surface_km = sum(res$sum_all_surface_km),
                             sum_A_only_km = sum(res$sum_A_only_km),
                             nonoverlapping_km = sum(res$nonoverlapping_km)))
res[, -1:-2] <- round(res[, -1:-2], 1)
write.csv(res, "outputs/total_distance_sampled.csv", row.names = FALSE)
print(res, row.names = FALSE)
