# =============================================================================
# MOSAIC Marine Mammals — Multi-Taxon Two-Level Occupancy Model (Hierarchical)
# Version 4 (v3 + concurrent-site random intercept; distance standardized within underway filters)
# =============================================================================
#
# Changes from v1 (phyloseq RDS input, d_AB shared, dnorm(0.5) inshore prior):
#   (1) Input: loads from three flat CSVs rather than a phyloseq RDS.
#   (2) mu_b_inshore prior: dnorm(0, sd=1) — neutral inshore/offshore expectation.
#   (3) MOSAIC-A and MOSAIC-B remain merged as a single shared contrast d_AB,
#       consistent with both being direct in-line filtration methods.
#
# Two levels:
#   Level 1 — Presence/capture (phi[i, f]):
#     logit(phi[i, f]) = alpha[f]
#                      + site_term[i]                  <- NEW: concurrent-site RE
#                      + d_cruise    * cruise2[i]
#                      + b_inshore[f] * inshore[i]
#                      + d_AB        * method_AB[i]
#                      + d_C         * method_C[i]
#                      + d_vol       * vol_sc[i]
#                      + d_dist      * dist_sc[i]       <- dist active for
#                                                          underway (non-CTD) only
#     alpha[f] and b_inshore[f] are taxon-level random effects (NCP).
#     d_cruise, d_AB, d_C, d_vol, d_dist are shared across taxa.
#
#   Concurrent-site random intercept (site_term[i]):
#     Guards against pseudo-replication of filters that sample the same water
#     mass. One vessel per cruise feeds all methods, so filters whose sampling
#     time-windows OVERLAP are pumping the same water and form one "site". This
#     temporal-concurrency signal — not a shared start coordinate — is what
#     links each parallel MOSAIC-B chain to its long concurrent MOSAIC-A
#     transect, even when their start points differ by 20+ km. Sites are
#     connected components of interval overlap computed PER CRUISE (RS2311 Fall
#     vs SR2408 Summer never pool). Filters with no concurrent partner fall back
#     to a coordinate station (START coords clustered within STATION_RADIUS_KM
#     per cruise). Single membership; the site effect is a non-centered,
#     zero-mean (grand mean absorbed by alpha) random intercept shared across
#     taxa:
#       s_site[j]    <- sigma_site * z_site[j],  z_site[j] ~ N(0, 1)
#       site_term[i] <- s_site[site_idx[i]]
#
#   Distance covariate (dist_sc):
#     dist_km == 0 iff MOSAIC-D (CTD point cast); A/B/C integrate along a
#     transect. Global z-scoring mapped CTD to a constant nonzero value,
#     making d_dist re-encode the D-vs-underway contrast. dist_sc is now
#     standardized WITHIN underway filters and multiplied by an underway
#     indicator, so the term is exactly 0 for every CTD filter and d_dist is
#     identified purely by transect-length variation among A/B/C.
#
#   Level 2 — Detection per tech rep (p[i]):
#     logit(p[i]) = delta + d_time * time_sc[i]
#     Shared across all taxa.
#
# Reference levels:
#   Method:  MOSAIC-D (CTD bottle cast) — method_AB = method_C = 0
#   Cruise:  RS2311 (cruise2 = 0); SR2408 = 1
#
# NIMBLE structure (avoids NIMBLE 1.4.1 posterior_predictive bug):
#   z[i, f]     ~ dbern(phi[i, f])
#   p_eff[i, f] <- p[i] * z[i, f]
#   n_det[i, f] ~ dbinom(p_eff[i, f], n_reps[i])
#
# =============================================================================

library(nimble)
library(coda)
library(bayesplot)
library(ggplot2)

set.seed(42)
dir.create("outputs/v4/figures", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA FROM CSVs
# =============================================================================

md  <- read.csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv",
                row.names = 1, stringsAsFactors = FALSE)
tax <- read.csv("outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv",
                row.names = 1, stringsAsFactors = FALSE)
otu <- read.csv("outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv",
                row.names = 1, check.names = FALSE)

cat("Sample metadata: ", nrow(md),  "rows x", ncol(md),  "cols\n")
cat("Taxa table:      ", nrow(tax), "rows x", ncol(tax), "cols\n")
cat("OTU table:       ", nrow(otu), "rows x", ncol(otu), "cols (taxa x samples)\n\n")

stopifnot(all(colnames(otu) %in% rownames(md)))
stopifnot(all(rownames(otu) %in% rownames(tax)))

otu <- otu[, rownames(md), drop = FALSE]

# =============================================================================
# 2. DATA ALREADY FILTERED
# =============================================================================
# Input CSVs are pre-filtered: exclusions for missing coordinates, impossible
# ship speeds, and NA covariates have already been applied upstream.

stopifnot(all(!is.na(md$Volume.filtered..mL.)))
stopifnot(all(!is.na(md$dist_km)))
stopifnot(all(!is.na(md$filtration_time_min)))
stopifnot(all(!is.na(md$Inshore.Offshore)))
stopifnot(all(!is.na(md$Collection.method)))

cat("Post-filter sample count:", nrow(md), "\n\n")

# =============================================================================
# 3. DEFINE FOCAL TAXA
# =============================================================================

cat("Families present in taxonomy table:\n")
print(sort(unique(tax$family[!is.na(tax$family)])))
cat("\n")

focal_families <- c(
  "Delphinidae",
  "Phocoenidae",
  "Balaenopteridae",
  "Otariidae",
  "Phocidae",
  "Alopiidae"
)

n_families <- length(focal_families)

missing_fam <- focal_families[!focal_families %in% tax$family]
if (length(missing_fam) > 0) {
  stop(paste("Focal families not found in taxa table:",
             paste(missing_fam, collapse = ", ")))
}

# =============================================================================
# 4. BUILD FILTER-LEVEL SUFFICIENT STATISTICS
# =============================================================================
# n_det_mat[i, f] = number of positive tech reps for filter i, family f.
# n_reps_vec[i]   = total tech reps for filter i (shared across families).

filter_ids <- unique(md$FilterID)
n_filters  <- length(filter_ids)

cat("Number of filters:", n_filters, "\n")

n_reps_vec <- integer(n_filters)
for (s in seq_along(filter_ids)) {
  n_reps_vec[s] <- sum(md$FilterID == filter_ids[s])
}

cat("Tech reps per filter:\n")
print(table(n_reps_vec))
cat("\n")

n_det_mat <- matrix(0L,
                    nrow     = n_filters,
                    ncol     = n_families,
                    dimnames = list(filter_ids, focal_families))

detect_summary <- data.frame(
  family     = focal_families,
  n_asvs     = NA_integer_,
  n_det_samp = NA_integer_,
  n_det_filt = NA_integer_
)

for (f in seq_len(n_families)) {
  fam      <- focal_families[f]
  fam_taxa <- rownames(tax)[!is.na(tax$family) & tax$family == fam]
  detect_summary$n_asvs[f] <- length(fam_taxa)

  if (length(fam_taxa) == 0) {
    cat(sprintf("WARNING: No ASVs found for family '%s'. All n_det = 0.\n", fam))
    next
  }

  fam_taxa_in_otu <- fam_taxa[fam_taxa %in% rownames(otu)]
  if (length(fam_taxa_in_otu) < length(fam_taxa)) {
    cat(sprintf("WARNING: %d of %d ASVs for '%s' not found in OTU rownames.\n",
                length(fam_taxa) - length(fam_taxa_in_otu), length(fam_taxa), fam))
  }
  fam_taxa <- fam_taxa_in_otu
  if (length(fam_taxa) == 0) next

  fam_otu   <- otu[fam_taxa, rownames(md), drop = FALSE]
  fam_reads <- if (length(fam_taxa) == 1) as.numeric(fam_otu[1, ]) else colSums(fam_otu)
  names(fam_reads) <- rownames(md)

  detect_vec        <- as.integer(fam_reads > 0)
  names(detect_vec) <- names(fam_reads)
  detect_summary$n_det_samp[f] <- sum(detect_vec)

  n_det_vec <- integer(n_filters)
  for (s in seq_along(filter_ids)) {
    samp_rows    <- rownames(md)[md$FilterID == filter_ids[s]]
    n_det_vec[s] <- sum(detect_vec[samp_rows])
  }

  n_det_mat[, f]               <- n_det_vec
  detect_summary$n_det_filt[f] <- sum(n_det_vec > 0)
}

cat("=== Detection summary by family ===\n")
print(detect_summary)
cat("\n")

cat("n_det_mat column sums (total positive tech-rep observations per family):\n")
print(colSums(n_det_mat))
cat("\n")

if (all(colSums(n_det_mat) < 5)) {
  stop("CRITICAL: All families have <5 total tech-rep detections. Check OTU/taxa alignment.")
}

zero_det <- detect_summary$family[detect_summary$n_det_filt == 0]
if (length(zero_det) > 0) {
  warning(paste("Families with 0 filter-level detections (phi unidentifiable):",
                paste(zero_det, collapse = ", ")))
}

# =============================================================================
# 5. FILTER-LEVEL COVARIATES
# =============================================================================

filter_md <- md[match(filter_ids, md$FilterID), ]

# --- Inshore/Offshore ---
inshore_vec <- as.integer(filter_md$Inshore.Offshore == "Inshore")

# --- Collection method indicators (reference = MOSAIC-D) ---
# MOSAIC-A and MOSAIC-B are both direct in-line filtration and share d_AB.
# MOSAIC-C (carboy) receives a separate contrast d_C.
# MOSAIC-D (CTD cast): method_AB = method_C = 0 (reference).
method_vec <- filter_md$Collection.method
method_AB  <- as.integer(method_vec %in% c("MOSAIC-A", "MOSAIC-B"))
method_C   <- as.integer(method_vec == "MOSAIC-C")

cat("Filter counts by method (MOSAIC-D = reference):\n")
print(table(method_vec))
cat("\n")

# --- Underway indicator (non-CTD): A/B/C integrate along a transect; D is a
#     stationary point cast with dist_km == 0. ---
underway_vec <- as.integer(method_AB == 1L | method_C == 1L)
stopifnot(all((filter_md$dist_km == 0) == (underway_vec == 0L)))  # dist0 <-> CTD

# --- Continuous covariates (z-scored) ---
scale2  <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
vol_sc  <- scale2(as.numeric(filter_md$Volume.filtered..mL.))
time_sc <- scale2(as.numeric(filter_md$filtration_time_min))

# Distance: standardize WITHIN underway filters, then zero-out CTD via the
# underway indicator. This makes d_dist * dist_sc exactly 0 for every CTD
# filter, so d_dist is identified only by transect-length variation among
# A/B/C and no longer aliases the D-vs-underway contrast (see header note).
dist_km_num  <- as.numeric(filter_md$dist_km)
uw           <- underway_vec == 1L
dist_mean_uw <- mean(dist_km_num[uw])
dist_sd_uw   <- sd(dist_km_num[uw])
dist_sc      <- underway_vec * (dist_km_num - dist_mean_uw) / dist_sd_uw
cat(sprintf("dist_sc: standardized within %d underway filters (mean=%.1f km, sd=%.1f km); CTD set to 0.\n\n",
            sum(uw), dist_mean_uw, dist_sd_uw))

# =============================================================================
# 5b. CRUISE COVARIATE (fixed-effect binary contrast)
# =============================================================================
# With only two cruises a random effect is unidentifiable. d_cruise is the
# log-odds difference between SR2408 (cruise2=1) and RS2311 (reference, =0).

cruise_labels <- sort(unique(filter_md$Cruise))
n_cruises     <- length(cruise_labels)

if (n_cruises != 2) {
  warning(sprintf(
    "%d cruises found; script expects exactly 2 for a binary fixed effect.",
    n_cruises))
}

cruise2_vec <- as.integer(filter_md$Cruise == cruise_labels[2])

cat(sprintf("Cruise reference: '%s' (= 0)  |  Contrast: '%s' (= 1)\n",
            cruise_labels[1], cruise_labels[2]))
cat(sprintf("Filters: %s = %d, %s = %d\n\n",
            cruise_labels[1], sum(cruise2_vec == 0),
            cruise_labels[2], sum(cruise2_vec == 1)))

# =============================================================================
# 5c. SITE GROUPING — temporal concurrency (per cruise), coordinate fallback
# =============================================================================
# The MOSAIC methods draw from one vessel's underway supply per cruise, so
# filters whose sampling time-windows OVERLAP are pumping the same water mass
# and belong to one "site" (concurrent sampling occasion). This — not a shared
# start coordinate — is the co-location signal that matters: each parallel
# MOSAIC-B chain runs concurrently with its long parent MOSAIC-A transect (every
# B-chain in this data overlaps an A for 39-229 min), and this holds even when
# their START coordinates differ by 20+ km because the A began its long run
# upstream. Sites are built PER CRUISE (a location is not the same biodiversity
# across seasons) in three passes:
#   (a) temporal-interval overlap (concurrent sampling) -> connected components;
#   (b) CTD anchoring (NON-BRIDGING) — each temporal-singleton CTD joins the
#       site of its SINGLE nearest transect whose START or END is within
#       STATION_RADIUS_KM of the CTD point. The CTD adopts one transect's site;
#       it never merges two transects together;
#   (c) coordinate fallback — any filter still without a partner is clustered by
#       START coordinate within STATION_RADIUS_KM (per cruise).
# Single membership: each filter belongs to exactly one site.

STATION_RADIUS_KM <- 5   # CTD-anchor + coordinate-fallback radius (km)

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
       cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

t0_sec <- as.numeric(as.POSIXct(filter_md$datetime_start, tz = "UTC"))
t1_sec <- as.numeric(as.POSIXct(filter_md$datetime_end,   tz = "UTC"))
stopifnot(all(!is.na(t0_sec)), all(!is.na(t1_sec)), all(t1_sec >= t0_sec))
stopifnot(all(!is.na(filter_md$Latitude.start)), all(!is.na(filter_md$Longitude.start)))

n_f <- nrow(filter_md)

is_ctd  <- filter_md$Collection.method == "MOSAIC-D"
is_tran <- filter_md$Collection.method %in% c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C")

# --- 1) Adjacency (same cruise): temporal concurrency only ---
# Overlapping sampling intervals = same pumped water. CTD spatial anchoring is
# applied in step 2b as a NON-BRIDGING attribution (not a graph edge), so a CTD
# can never merge two otherwise-separate transects together.
adj <- matrix(FALSE, n_f, n_f)
for (i in 1:(n_f - 1)) for (j in (i + 1):n_f) {
  if (filter_md$Cruise[i] == filter_md$Cruise[j] &&
      min(t1_sec[i], t1_sec[j]) > max(t0_sec[i], t0_sec[j])) {
    adj[i, j] <- adj[j, i] <- TRUE
  }
}

# --- 2) Connected components = concurrent sampling occasions ---
site <- integer(n_f); k <- 0L
for (i in seq_len(n_f)) {
  if (site[i] != 0) next
  k <- k + 1L; stack <- i
  while (length(stack)) {
    v <- stack[1]; stack <- stack[-1]
    if (site[v] == 0) { site[v] <- k; stack <- c(stack, which(adj[v, ] & site == 0)) }
  }
}
n_temporal <- k

# --- 2b) CTD anchoring (NON-BRIDGING) ---
# Each CTD that is a temporal singleton is attributed to the site of its SINGLE
# nearest transect whose START or END lies within STATION_RADIUS_KM of the CTD
# point. The CTD adopts that one transect's existing site; no graph edge is
# added, so two transects are never merged through a shared-station CTD. CTDs
# already grouped by concurrency, or with no transect in range, are left as-is.
comp_size0 <- table(site)
for (i in which(is_ctd)) {
  if (comp_size0[as.character(site[i])] > 1) next          # already concurrent
  tr <- which(is_tran & filter_md$Cruise == filter_md$Cruise[i])
  if (length(tr) == 0) next
  d <- pmin(
    haversine_km(filter_md$Latitude.start[i], filter_md$Longitude.start[i],
                 filter_md$Latitude.start[tr], filter_md$Longitude.start[tr]),
    haversine_km(filter_md$Latitude.start[i], filter_md$Longitude.start[i],
                 filter_md$Latitude.end[tr],   filter_md$Longitude.end[tr]))
  if (min(d) < STATION_RADIUS_KM) site[i] <- site[tr[which.min(d)]]
}

# --- 3) Coordinate fallback for remaining singletons ---
# Filters still with no partner are clustered by START coordinate within each
# cruise, so repeat / near-repeat visits still share a site.
comp_size  <- table(site)
singletons <- which(comp_size[as.character(site)] == 1)
next_id    <- n_temporal
for (cr in unique(filter_md$Cruise)) {
  s <- singletons[filter_md$Cruise[singletons] == cr]
  if (length(s) < 1) next
  lab <- integer(length(s))
  for (a in seq_along(s)) {
    if (lab[a] != 0) next
    next_id <- next_id + 1L; lab[a] <- next_id
    d <- haversine_km(filter_md$Latitude.start[s[a]], filter_md$Longitude.start[s[a]],
                      filter_md$Latitude.start[s],    filter_md$Longitude.start[s])
    lab[d < STATION_RADIUS_KM & lab == 0] <- next_id
  }
  site[s] <- lab
}

# --- 4) Relabel sites to a contiguous 1..n_sites index ---
site_idx <- as.integer(factor(site))
n_sites  <- length(unique(site_idx))

cat(sprintf("Site grouping: %d temporal components + non-bridging CTD anchoring + coordinate fallback -> %d sites for %d filters\n",
            n_temporal, n_sites, n_f))
site_sz <- table(site_idx)
cat(sprintf("Sites: %d singletons | %d filters in multi-filter sites | largest site: %d filters\n",
            sum(site_sz == 1), sum(site_sz[as.character(site_idx)] > 1), max(site_sz)))
mt_chk <- table(site_idx, filter_md$Collection.method)
uw_cols <- intersect(c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C"), colnames(mt_chk))
ctd_grouped <- if ("MOSAIC-D" %in% colnames(mt_chk))
  sum(is_ctd & site_sz[as.character(site_idx)] > 1) else 0
cat(sprintf("CTD casts grouped with >=1 other filter: %d of %d | Sites mixing >=2 methods: %d | CTD+underway sites: %d\n\n",
            ctd_grouped, sum(is_ctd), sum(rowSums(mt_chk > 0) >= 2),
            if ("MOSAIC-D" %in% colnames(mt_chk))
              sum(mt_chk[, "MOSAIC-D"] > 0 & rowSums(mt_chk[, uw_cols, drop = FALSE]) > 0) else 0))

# =============================================================================
# 6. BUNDLE DATA AND CONSTANTS FOR NIMBLE
# =============================================================================

nimble_data <- list(
  n_det = n_det_mat
)

nimble_constants <- list(
  n_filters  = n_filters,
  n_families = n_families,
  n_reps     = n_reps_vec,
  cruise2    = cruise2_vec,
  inshore    = inshore_vec,
  method_AB  = method_AB,
  method_C   = method_C,
  vol_sc     = vol_sc,
  dist_sc    = dist_sc,
  time_sc    = time_sc,
  n_sites    = n_sites,
  site_idx   = site_idx,
  use_site   = 1          # 1 = site RE active; set 0 for the no-site WAIC refit
)

# =============================================================================
# 7. NIMBLE MODEL CODE
# =============================================================================

mammal_model_code <- nimbleCode({

  # ---------------------------------------------------------------------------
  # Hyperpriors — taxon-level random effects
  # ---------------------------------------------------------------------------
  mu_alpha        ~ dnorm(0, sd = 1.5)
  sigma_alpha     ~ dexp(1)

  # Neutral prior: no directional expectation for inshore/offshore effect
  mu_b_inshore    ~ dnorm(0, sd = 1)
  sigma_b_inshore ~ dexp(1)

  # ---------------------------------------------------------------------------
  # Non-centered taxon-level random effects (NCP)
  # ---------------------------------------------------------------------------
  # Sampling z ~ N(0,1) and recovering alpha/b_inshore deterministically
  # avoids Neal's funnel and improves MCMC mixing in hierarchical models.
  for (f in 1:n_families) {
    z_alpha[f]     ~ dnorm(0, sd = 1)
    alpha[f]       <- mu_alpha + sigma_alpha * z_alpha[f]

    z_b_inshore[f] ~ dnorm(0, sd = 1)
    b_inshore[f]   <- mu_b_inshore + sigma_b_inshore * z_b_inshore[f]
  }

  # ---------------------------------------------------------------------------
  # Site (concurrent-sampling-occasion) random intercept (NCP, zero-mean, shared)
  # ---------------------------------------------------------------------------
  # Zero-mean: the grand level is absorbed by alpha[f], so only the site
  # deviations are estimated. sigma_site governs how much genuine site-level
  # heterogeneity the data support (weakly identified given many small sites —
  # treat a small posterior as "little residual site structure").
  sigma_site ~ dexp(1)
  for (j in 1:n_sites) {
    z_site[j] ~ dnorm(0, sd = 1)
    s_site[j] <- sigma_site * z_site[j]
  }

  # ---------------------------------------------------------------------------
  # Shared priors — presence/capture submodel
  # ---------------------------------------------------------------------------
  # Reference method: MOSAIC-D (CTD cast).
  # d_AB: shared contrast for MOSAIC-A and MOSAIC-B (both direct in-line) vs D.
  # d_C:  contrast for MOSAIC-C (carboy) vs D.
  d_AB   ~ dnorm(0, sd = 1)
  d_C    ~ dnorm(0, sd = 1)
  d_vol  ~ dnorm(0, sd = 1)
  d_dist ~ dnorm(0, sd = 1)

  # ---------------------------------------------------------------------------
  # Cruise fixed effect — binary contrast
  # ---------------------------------------------------------------------------
  d_cruise ~ dnorm(0, sd = 1)   # SR2408 vs RS2311 (reference)

  # ---------------------------------------------------------------------------
  # Shared priors — detection submodel
  # ---------------------------------------------------------------------------
  # Prior centred at -1 (logit) -> p ~ 27% per rep at baseline.
  # Prevents the degenerate low-phi/high-p posterior corner.
  delta  ~ dnorm(-1, sd = 1)
  d_time ~ dnorm(0,  sd = 1)

  # ---------------------------------------------------------------------------
  # Likelihood
  # ---------------------------------------------------------------------------
  for (i in 1:n_filters) {

    # Per-rep detection probability — shared across taxa, varies by filter
    logit_p[i] <- delta + d_time * time_sc[i]
    p[i]       <- ilogit(logit_p[i])

    # Site effect: the filter's concurrent-sampling-occasion random intercept
    # (single membership from the temporal-overlap components + coordinate
    # fallback built in Section 5c). use_site toggles it off (0) for the
    # no-site WAIC comparison without changing the model graph.
    site_term[i] <- use_site * s_site[site_idx[i]]

    for (f in 1:n_families) {

      # Presence/capture probability
      # Reference: MOSAIC-D, RS2311, Offshore
      logit_phi[i, f] <- alpha[f]                    +
                         site_term[i]                 +
                         d_cruise  * cruise2[i]       +
                         b_inshore[f] * inshore[i]    +
                         d_AB      * method_AB[i]     +
                         d_C       * method_C[i]      +
                         d_vol     * vol_sc[i]        +
                         d_dist    * dist_sc[i]
      phi[i, f] <- ilogit(logit_phi[i, f])

      # Latent presence/capture state
      z[i, f] ~ dbern(phi[i, f])

      # Deterministic effective detection probability.
      # z=0 -> p_eff=0, forcing n_det=0 via Binomial(0, n_reps).
      # z=1 -> p_eff=p, activating the full Binomial likelihood.
      # Explicit deterministic node ensures correct NIMBLE graph traversal
      # and sampler assignment (avoids posterior_predictive bug in v1.4.1).
      p_eff[i, f] <- p[i] * z[i, f]

      n_det[i, f] ~ dbinom(p_eff[i, f], n_reps[i])
    }
  }

})

# =============================================================================
# 8. INITIAL VALUES
# =============================================================================

occ_rate        <- colMeans(n_det_mat > 0)
alpha_data_init <- qlogis(pmax(0.05, pmin(0.95, occ_rate)))

cat("Data-informed alpha[f] initialisation:\n")
print(data.frame(family     = focal_families,
                 occ_rate   = round(occ_rate, 3),
                 alpha_init = round(alpha_data_init, 2)))
cat("\n")

mu_alpha_0    <- mean(alpha_data_init)
sigma_alpha_0 <- max(0.5, sd(alpha_data_init))
z_alpha_0     <- (alpha_data_init - mu_alpha_0) / sigma_alpha_0

# z matrix: 1 wherever any detection observed, 0 elsewhere
z_init_mat <- matrix(0L, nrow = n_filters, ncol = n_families)
for (f in seq_len(n_families)) {
  z_init_mat[, f] <- as.integer(n_det_mat[, f] > 0)
}

inits <- function() list(
  mu_alpha        = mu_alpha_0,
  sigma_alpha     = sigma_alpha_0,
  mu_b_inshore    = 0,
  sigma_b_inshore = 0.5,
  z_alpha         = z_alpha_0,
  z_b_inshore     = rep(0, n_families),
  d_AB            = 0,
  d_C             = 0,
  d_vol           = 0.5,
  d_dist          = 0.5,
  d_cruise        = 0,
  delta           = -1,
  d_time          = -0.5,
  sigma_site      = 0.5,
  z_site          = rep(0, n_sites),
  z               = z_init_mat
)

# =============================================================================
# 9. PARAMETERS TO MONITOR
# =============================================================================

params_to_monitor <- c(
  # Hyperparameters
  "mu_alpha", "sigma_alpha",
  "mu_b_inshore", "sigma_b_inshore",
  # NCP standardised deviates (stochastic nodes — assess mixing on these)
  "z_alpha", "z_b_inshore",
  # Derived taxon random effects (deterministic — for inference and plots)
  "alpha", "b_inshore",
  # Shared presence/capture coefficients
  "d_cruise", "d_AB", "d_C", "d_vol", "d_dist",
  # Site random-effect scale (site deviates s_site are a high-dimensional
  # nuisance — only the scale is monitored for inference)
  "sigma_site",
  # Shared detection coefficients
  "delta", "d_time",
  # Filter-level latent probabilities
  "phi", "p"
  # Note: z (latent binary) is not monitored — high-dimensional nuisance
)

# =============================================================================
# 10. BUILD AND RUN MODEL
# =============================================================================

cat("Building NIMBLE model...\n")
mammal_model <- nimbleModel(
  code      = mammal_model_code,
  constants = nimble_constants,
  data      = nimble_data,
  inits     = inits()
)

cat("Compiling model...\n")
mammal_model_c <- compileNimble(mammal_model)

cat("Configuring MCMC...\n")
mcmc_conf <- configureMCMC(
  mammal_model,
  monitors = params_to_monitor,
  thin     = 10
)

cat("Building MCMC...\n")
mammal_mcmc   <- buildMCMC(mcmc_conf)
mammal_mcmc_c <- compileNimble(mammal_mcmc, project = mammal_model)

cat("Running MCMC (3 chains)...\n")
n_iter   <- 100000
n_burnin <- 20000
n_thin   <- 20

samples <- runMCMC(
  mammal_mcmc_c,
  niter             = n_iter + n_burnin,
  nburnin           = n_burnin,
  nchains           = 3,
  inits             = lapply(1:3, function(x) inits()),
  thin              = n_thin,
  samplesAsCodaMCMC = TRUE
)

# =============================================================================
# 11. CONVERGENCE DIAGNOSTICS
# =============================================================================

cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")

samples_combined <- do.call(rbind, lapply(samples, as.matrix))

key_params <- c(
  "mu_alpha", "sigma_alpha",
  "mu_b_inshore", "sigma_b_inshore",
  "d_cruise", "d_AB", "d_C", "d_vol", "d_dist",
  "sigma_site",
  "delta", "d_time",
  paste0("z_alpha[",     1:n_families, "]"),
  paste0("z_b_inshore[", 1:n_families, "]")
)
key_params <- intersect(key_params, colnames(samples_combined))

gelman_diag <- gelman.diag(samples[, key_params], multivariate = FALSE)
cat("Gelman-Rubin R-hat (values near 1.0 = converged):\n")
print(round(gelman_diag$psrf, 3))

ess <- effectiveSize(samples_combined[, key_params])
cat("\nEffective sample sizes:\n")
print(round(ess))

poor <- key_params[gelman_diag$psrf[, 1] > 1.05]
if (length(poor) > 0) {
  cat("\nWARNING: R-hat > 1.05 for:", paste(poor, collapse = ", "), "\n")
  cat("  Consider running longer chains.\n")
} else {
  cat("\nAll R-hat < 1.05. Convergence looks good.\n")
}

# =============================================================================
# 12. POSTERIOR SUMMARIES
# =============================================================================

cat("\n=== POSTERIOR SUMMARIES ===\n")

summarize_param <- function(mat, params) {
  params <- intersect(params, colnames(mat))
  if (length(params) == 0) return(NULL)
  t(apply(mat[, params, drop = FALSE], 2, function(x) {
    c(mean  = round(mean(x), 3),
      sd    = round(sd(x), 3),
      q2.5  = round(quantile(x, 0.025), 3),
      q50   = round(median(x), 3),
      q97.5 = round(quantile(x, 0.975), 3),
      P_pos = round(mean(x > 0), 3))
  }))
}

cat("\n--- Hyperparameters ---\n")
print(summarize_param(samples_combined,
      c("mu_alpha", "sigma_alpha", "mu_b_inshore", "sigma_b_inshore",
        "sigma_site")))

cat(sprintf(
  "\n--- Shared presence/capture covariates (method ref = MOSAIC-D | cruise ref = %s) ---\n",
  cruise_labels[1]))
print(summarize_param(samples_combined,
      c("d_cruise", "d_AB", "d_C", "d_vol", "d_dist")))

cat("\n--- Shared detection submodel ---\n")
print(summarize_param(samples_combined, c("delta", "d_time")))

cat("\n--- Taxon intercepts alpha[f] ---\n")
summ_alpha <- summarize_param(samples_combined, paste0("alpha[", 1:n_families, "]"))
if (!is.null(summ_alpha)) rownames(summ_alpha) <- focal_families
print(summ_alpha)

cat("\n--- Taxon inshore offsets b_inshore[f] ---\n")
summ_binsh <- summarize_param(samples_combined, paste0("b_inshore[", 1:n_families, "]"))
if (!is.null(summ_binsh)) rownames(summ_binsh) <- focal_families
print(summ_binsh)

# =============================================================================
# 13. POSTERIOR PLOTS
# =============================================================================

label_family_params <- function(samp_list, base, n_fam, labels) {
  old_nms <- paste0(base, "[", 1:n_fam, "]")
  new_nms <- paste0(base, "[", labels, "]")
  lapply(samp_list, function(m) {
    idx <- match(old_nms, colnames(m))
    colnames(m)[idx[!is.na(idx)]] <- new_nms[!is.na(idx)]
    m
  })
}

samples_lab <- label_family_params(samples,      "alpha",        n_families, focal_families)
samples_lab <- label_family_params(samples_lab,  "b_inshore",    n_families, focal_families)
samples_lab <- label_family_params(samples_lab,  "z_alpha",      n_families, focal_families)
samples_lab <- label_family_params(samples_lab,  "z_b_inshore",  n_families, focal_families)

samples_lab_combined <- do.call(rbind, lapply(samples_lab, as.matrix))

alpha_lab    <- paste0("alpha[",     focal_families, "]")
binsh_lab    <- paste0("b_inshore[", focal_families, "]")
subtitle_ref <- sprintf("Method ref: MOSAIC-D  |  Cruise ref: %s", cruise_labels[1])

# (a) Trace plots — scalar hyperparameters and shared coefficients
scalar_params <- intersect(
  c("mu_alpha", "sigma_alpha", "mu_b_inshore", "sigma_b_inshore",
    "d_cruise", "d_AB", "d_C", "d_vol", "d_dist", "sigma_site",
    "delta", "d_time"),
  colnames(samples_lab[[1]])
)

p_trace_scalar <- mcmc_trace(
  samples_lab,
  pars       = scalar_params,
  facet_args = list(ncol = 2, labeller = label_parsed)
) +
  labs(title    = "Trace plots — hyperparameters & shared coefficients",
       subtitle = subtitle_ref) +
  theme_bw()

ggsave("outputs/v4/figures/plot_trace_scalar.pdf", p_trace_scalar, width = 10, height = 12, dpi = 150)
cat("Saved: plot_trace_scalar.pdf\n")

# (b) Trace plots — NCP deviates
z_alpha_lab <- paste0("z_alpha[",     focal_families, "]")
z_binsh_lab <- paste0("z_b_inshore[", focal_families, "]")
re_params   <- intersect(c(z_alpha_lab, z_binsh_lab), colnames(samples_lab[[1]]))

p_trace_re <- mcmc_trace(
  samples_lab,
  pars       = re_params,
  facet_args = list(ncol = 2)
) +
  labs(title = "Trace plots — NCP deviates (z_alpha, z_b_inshore)") +
  theme_bw()

ggsave("outputs/v4/figures/plot_trace_re.pdf", p_trace_re, width = 10, height = 10, dpi = 150)
cat("Saved: plot_trace_re.pdf\n")

# (c) Interval plots — taxon intercepts
p_alpha_int <- mcmc_intervals(
  samples_lab_combined,
  pars       = intersect(alpha_lab, colnames(samples_lab_combined)),
  prob       = 0.5,
  prob_outer = 0.95
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title = "Taxon intercepts: alpha[f]  (logit scale)",
       x = "log-odds", y = NULL) +
  theme_bw()

ggsave("outputs/v4/figures/plot_alpha_intervals.pdf", p_alpha_int, width = 7, height = 4, dpi = 150)
cat("Saved: plot_alpha_intervals.pdf\n")

# (d) Interval plots — taxon inshore offsets
p_binsh_int <- mcmc_intervals(
  samples_lab_combined,
  pars       = intersect(binsh_lab, colnames(samples_lab_combined)),
  prob       = 0.5,
  prob_outer = 0.95
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title = "Taxon inshore offsets: b_inshore[f]  (logit scale)",
       x = "log-odds", y = NULL) +
  theme_bw()

ggsave("outputs/v4/figures/plot_binshore_intervals.pdf", p_binsh_int, width = 7, height = 4, dpi = 150)
cat("Saved: plot_binshore_intervals.pdf\n")

# (e) Posterior area plots — shared presence/capture coefficients
p_phi_coef <- mcmc_areas(
  samples_lab_combined,
  pars      = intersect(c("d_cruise", "d_AB", "d_C", "d_vol", "d_dist"),
                        colnames(samples_lab_combined)),
  prob      = 0.95,
  point_est = "median"
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title    = "Shared presence/capture coefficients  (logit scale)",
       subtitle = subtitle_ref,
       x        = "log-odds") +
  theme_bw()

ggsave("outputs/v4/figures/plot_phi_coef_areas.pdf", p_phi_coef, width = 7, height = 4, dpi = 150)
cat("Saved: plot_phi_coef_areas.pdf\n")

# (f) Posterior area plots — shared detection coefficients
p_det_coef <- mcmc_areas(
  samples_lab_combined,
  pars      = intersect(c("delta", "d_time"), colnames(samples_lab_combined)),
  prob      = 0.95,
  point_est = "median"
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title = "Shared detection coefficients  (logit scale)",
       x     = "log-odds") +
  theme_bw()

ggsave("outputs/v4/figures/plot_det_coef_areas.pdf", p_det_coef, width = 7, height = 3, dpi = 150)
cat("Saved: plot_det_coef_areas.pdf\n")

# (g) Combined interval plot — all parameters
all_interval_params <- intersect(
  c(scalar_params, alpha_lab, binsh_lab),
  colnames(samples_lab_combined)
)

p_all_int <- mcmc_intervals(
  samples_lab_combined,
  pars       = all_interval_params,
  prob       = 0.5,
  prob_outer = 0.95
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title    = "All parameter intervals  (50% inner, 95% outer CI)",
       subtitle = subtitle_ref,
       x        = "posterior value") +
  theme_bw()

ggsave("outputs/v4/figures/plot_all_intervals.pdf", p_all_int, width = 8, height = 10, dpi = 150)
cat("Saved: plot_all_intervals.pdf\n")

# =============================================================================
# 14. DERIVED QUANTITIES
# =============================================================================

cat("\n=== DERIVED: Mean phi by family and method ===\n")
phi_cols <- grep("^phi\\[", colnames(samples_combined), value = TRUE)
phi_samp <- samples_combined[, phi_cols, drop = FALSE]

for (f in seq_len(n_families)) {
  fam        <- focal_families[f]
  phi_f_cols <- grep(sprintf(",\\s*%d\\]$", f), phi_cols, value = TRUE)
  if (length(phi_f_cols) == 0) next
  cat(sprintf("\n  Family: %s\n", fam))
  for (m in c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C", "MOSAIC-D")) {
    idx <- which(filter_md$Collection.method == m)
    if (length(idx) == 0) next
    phi_m <- rowMeans(phi_samp[, phi_f_cols[idx], drop = FALSE])
    cat(sprintf("    %-10s  mean=%.3f  95%%CI=(%.3f, %.3f)\n",
                m, mean(phi_m), quantile(phi_m, 0.025), quantile(phi_m, 0.975)))
  }
}

cat("\n=== DERIVED: Mean p (detection probability) by method ===\n")
p_cols <- grep("^p\\[", colnames(samples_combined), value = TRUE)
p_samp <- samples_combined[, p_cols, drop = FALSE]

for (m in c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C", "MOSAIC-D")) {
  idx <- which(filter_md$Collection.method == m)
  if (length(idx) == 0) next
  p_m <- rowMeans(p_samp[, idx, drop = FALSE])
  cat(sprintf("  %-10s  mean=%.3f  95%%CI=(%.3f, %.3f)\n",
              m, mean(p_m), quantile(p_m, 0.025), quantile(p_m, 0.975)))
}

cat("\n=== DERIVED: P(>=1 detection | present+captured, 4 reps) by family+method ===\n")
for (f in seq_len(n_families)) {
  fam        <- focal_families[f]
  phi_f_cols <- grep(sprintf(",\\s*%d\\]$", f), phi_cols, value = TRUE)
  if (length(phi_f_cols) == 0) next
  cat(sprintf("\n  Family: %s\n", fam))
  for (m in c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C", "MOSAIC-D")) {
    idx <- which(filter_md$Collection.method == m)
    if (length(idx) == 0) next
    p_m  <- rowMeans(p_samp[, idx, drop = FALSE])
    p4_m <- 1 - (1 - p_m)^4
    cat(sprintf("    %-10s  mean=%.3f  95%%CI=(%.3f, %.3f)\n",
                m, mean(p4_m), quantile(p4_m, 0.025), quantile(p4_m, 0.975)))
  }
}

# =============================================================================
# 14b. WAIC MODEL COMPARISON — site random effect vs no site effect
# =============================================================================
# Marginal (over the latent z) WAIC from the monitored phi[i,f] and p[i].
# Integrating z out of  z ~ Bern(phi);  n_det ~ Binom(n_reps, p*z)  gives
#   P(y) = (1 - phi) * 1{y == 0} + phi * dbinom(y, n_reps, p),
# the "new filter observation" predictive density — the appropriate WAIC target
# for an occupancy model, and identical in form across the two fits. The reduced
# model re-runs with use_site = 0 (site term switched off; the RE nodes stay in
# the graph but do not enter phi, so its WAIC is exactly a no-site model's).
# Lower WAIC = better expected predictive accuracy.

waic_marginal <- function(samp, n_det, n_reps) {
  phi_cols <- grep("^phi\\[", colnames(samp), value = TRUE)
  m  <- regmatches(phi_cols, regexec("^phi\\[([0-9]+),\\s*([0-9]+)\\]$", phi_cols))
  ii <- as.integer(vapply(m, `[`, "", 2L))
  ff <- as.integer(vapply(m, `[`, "", 3L))
  p_cols <- paste0("p[", seq_len(nrow(n_det)), "]")
  stopifnot(all(p_cols %in% colnames(samp)), !anyNA(ii), !anyNA(ff))
  elpd <- numeric(length(phi_cols)); pw <- numeric(length(phi_cols))
  for (cc in seq_along(phi_cols)) {
    i <- ii[cc]; f <- ff[cc]; y <- n_det[i, f]; nr <- n_reps[i]
    phi_s <- samp[, phi_cols[cc]]; p_s <- samp[, p_cols[i]]
    lik_s <- (1 - phi_s) * (y == 0) + phi_s * dbinom(y, nr, p_s)
    elpd[cc] <- log(mean(lik_s))    # log pointwise predictive density
    pw[cc]   <- var(log(lik_s))     # effective-parameter penalty
  }
  elpd_point <- elpd - pw
  names(elpd_point) <- phi_cols
  list(waic = -2 * sum(elpd_point), lppd = sum(elpd), p_waic = sum(pw),
       se = 2 * sqrt(length(elpd_point) * var(elpd_point)), elpd_point = elpd_point)
}

cat("\n=== WAIC MODEL COMPARISON (site RE vs none) ===\n")
waic_site <- waic_marginal(samples_combined, n_det_mat, n_reps_vec)
cat(sprintf("Site model    : WAIC = %.1f  (p_waic = %.1f, SE = %.1f)\n",
            waic_site$waic, waic_site$p_waic, waic_site$se))

cat("Re-fitting reduced model (use_site = 0)...\n")
nimble_constants_ns <- modifyList(nimble_constants, list(use_site = 0))
model_ns    <- nimbleModel(code = mammal_model_code, constants = nimble_constants_ns,
                           data = nimble_data, inits = inits())
model_ns_c  <- compileNimble(model_ns)
conf_ns     <- configureMCMC(model_ns, monitors = c("phi", "p"), thin = 10)
mcmc_ns     <- buildMCMC(conf_ns)
mcmc_ns_c   <- compileNimble(mcmc_ns, project = model_ns)
samples_ns  <- runMCMC(mcmc_ns_c, niter = n_iter + n_burnin, nburnin = n_burnin,
                       nchains = 3, inits = lapply(1:3, function(x) inits()),
                       thin = n_thin, samplesAsCodaMCMC = TRUE)
samples_ns_combined <- do.call(rbind, lapply(samples_ns, as.matrix))
waic_nosite <- waic_marginal(samples_ns_combined, n_det_mat, n_reps_vec)
cat(sprintf("No-site model : WAIC = %.1f  (p_waic = %.1f, SE = %.1f)\n",
            waic_nosite$waic, waic_nosite$p_waic, waic_nosite$se))

# Paired difference (aligned by phi node); dWAIC > 0 favours the site model.
common     <- intersect(names(waic_site$elpd_point), names(waic_nosite$elpd_point))
diff_point <- waic_site$elpd_point[common] - waic_nosite$elpd_point[common]
dwaic      <- waic_nosite$waic - waic_site$waic
se_dwaic   <- 2 * sqrt(length(diff_point) * var(diff_point))
cat(sprintf("\ndWAIC (no-site - site) = %.1f +/- %.1f (SE of difference)\n", dwaic, se_dwaic))
cat(sprintf("=> site effect %s predictive fit; difference is %s (|dWAIC| %s 2*SE)\n",
            if (dwaic > 0) "IMPROVES" else "does NOT improve",
            if (abs(dwaic) > 2 * se_dwaic) "clear" else "within MC/statistical noise",
            if (abs(dwaic) > 2 * se_dwaic) ">" else "<="))

# =============================================================================
# 15. SAVE
# =============================================================================

save(
  samples, samples_combined,
  nimble_data, nimble_constants,
  filter_ids, filter_md, n_reps_vec, n_det_mat,
  focal_families, detect_summary,
  n_sites, site_idx, n_temporal, STATION_RADIUS_KM,
  waic_site, waic_nosite, dwaic, se_dwaic, samples_ns_combined,
  file = "outputs/v4/MOSAIC_mammal_families_twolevel_v4_results.RData"
)

cat("\nAll results saved to: MOSAIC_mammal_families_twolevel_v4_results.RData (outputs/v4/)\n")
cat("Script complete.\n")

