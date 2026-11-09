# =============================================================================
# MOSAIC Occupancy Model — Method Effectiveness: Results Extraction & Plots
# Version 3 (MOSAIC-A and MOSAIC-B merged as AB)
# =============================================================================
#
# Requires:
#   MOSAIC_mammal_families_twolevel_v3_results.RData
#   MOSAIC_taxnamed_sample_data_metadata_filtered.csv
#
# Produces six figures saved to viz/:
#   Fig 1. Method contrasts  — posterior log-odds (d_AB, d_C vs MOSAIC-D)
#   Fig 2. Marginal phi      — posterior mean ± CI by method & family
#   Fig 3. Operational P     — P(≥1 pos rep | filter) by method & family
#   Fig 4. Covariate effects — shared continuous/cruise coefficients
#   Fig 5. Volume × phi      — mechanistic scatter
#   Fig 6. Observed vs model — raw detection rate vs posterior mean phi
#
# =============================================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# AB is represented by a single colour; three-method palette
method_colours <- c(
  "MOSAIC-AB" = "#E69F00",
  "MOSAIC-C"  = "#009E73",
  "MOSAIC-D"  = "#CC79A7"
)
method_order <- c("MOSAIC-AB", "MOSAIC-C", "MOSAIC-D")

# For derived quantities we still distinguish A and B in metadata
method_colours_raw <- c(
  "MOSAIC-A"  = "#E69F00",
  "MOSAIC-B"  = "#F0C240",
  "MOSAIC-C"  = "#009E73",
  "MOSAIC-D"  = "#CC79A7"
)

family_labels <- c(
  "Delphinidae"     = "Delphinidae\n(dolphins)",
  "Phocoenidae"     = "Phocoenidae\n(porpoises)",
  "Balaenopteridae" = "Balaenopteridae\n(rorqual whales)",
  "Otariidae"       = "Otariidae\n(sea lions)",
  "Phocidae"        = "Phocidae\n(true seals)",
  "Alopiidae"       = "Alopiidae\n(mackerel sharks)"
)

post_summary <- function(x) {
  data.frame(
    mean   = mean(x),
    median = median(x),
    lo95   = quantile(x, 0.025),
    hi95   = quantile(x, 0.975),
    lo50   = quantile(x, 0.25),
    hi50   = quantile(x, 0.75),
    p_pos  = mean(x > 0)
  )
}

theme_mosaic <- function() {
  theme_bw(base_size = 11, base_family = "serif") +
    theme(
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, colour = "grey40")
    )
}

dir.create("viz/Bayesian_occupancy_modeling/additional_plots", showWarnings = FALSE)

# =============================================================================
# 1. LOAD MODEL OUTPUT AND METADATA
# =============================================================================

load("outputs/MOSAIC_mammal_families_twolevel_v3_results.RData")

md <- read.csv("outputs/MOSAIC_taxnamed_sample_data.csv",
               row.names = 1, stringsAsFactors = FALSE)

n_filters  <- length(filter_ids)
n_families <- length(focal_families)

cat("Loaded:", n_filters, "filters,", n_families, "families,",
    nrow(samples_combined), "posterior samples\n")

# =============================================================================
# 2. EXTRACT PHI AND P POSTERIOR MATRICES
# =============================================================================

phi_cols <- grep("^phi\\[", colnames(samples_combined), value = TRUE)
p_cols   <- grep("^p\\[",   colnames(samples_combined), value = TRUE)

parse_phi_idx <- function(cols) {
  m <- regmatches(cols, regexec("phi\\[(\\d+),\\s*(\\d+)\\]", cols))
  data.frame(
    col = cols,
    i   = as.integer(sapply(m, `[`, 2)),
    f   = as.integer(sapply(m, `[`, 3)),
    stringsAsFactors = FALSE
  )
}

parse_p_idx <- function(cols) {
  m <- regmatches(cols, regexec("p\\[(\\d+)\\]", cols))
  data.frame(
    col = cols,
    i   = as.integer(sapply(m, `[`, 2)),
    stringsAsFactors = FALSE
  )
}

phi_idx <- parse_phi_idx(phi_cols)
p_idx   <- parse_p_idx(p_cols)

cat("phi columns parsed:", nrow(phi_idx), "\n")
cat("p   columns parsed:", nrow(p_idx),   "\n")

# Filter-level metadata aligned to filter_ids order
filter_meta        <- filter_md[match(filter_ids, filter_md$FilterID), ]
filter_meta$n_reps <- n_reps_vec

# Assign merged AB label for grouping
filter_meta$method_group <- ifelse(
  filter_meta$Collection.method %in% c("MOSAIC-A", "MOSAIC-B"),
  "MOSAIC-AB",
  filter_meta$Collection.method
)

# =============================================================================
# 3. POSTERIOR PHI BY METHOD GROUP AND FAMILY
# =============================================================================

phi_method_family <- do.call(rbind, lapply(focal_families, function(fam) {
  f_idx <- which(focal_families == fam)
  do.call(rbind, lapply(method_order, function(meth) {
    filt_idx <- which(filter_meta$method_group == meth)
    if (length(filt_idx) == 0) {
      summ <- data.frame(mean=NA, median=NA, lo95=NA, hi95=NA,
                         lo50=NA, hi50=NA, p_pos=NA)
    } else {
      these_cols <- phi_idx$col[phi_idx$i %in% filt_idx & phi_idx$f == f_idx]
      phi_vals   <- as.vector(samples_combined[, these_cols, drop = FALSE])
      summ       <- post_summary(phi_vals)
    }
    summ$family <- fam
    summ$method <- meth
    summ
  }))
})) %>%
  mutate(
    family       = factor(family, levels = focal_families),
    method       = factor(method, levels = method_order),
    family_label = family_labels[as.character(family)]
  )

# =============================================================================
# 4. OPERATIONAL DETECTION PROBABILITY: P(≥1 pos rep | filter)
# =============================================================================

op_detect_method_family <- do.call(rbind, lapply(focal_families, function(fam) {
  f_idx <- which(focal_families == fam)
  do.call(rbind, lapply(method_order, function(meth) {
    filt_idx <- which(filter_meta$method_group == meth)
    if (length(filt_idx) == 0) {
      summ <- data.frame(mean=NA, median=NA, lo95=NA, hi95=NA,
                         lo50=NA, hi50=NA, p_pos=NA)
    } else {
      phi_these  <- phi_idx$col[phi_idx$i %in% filt_idx & phi_idx$f == f_idx]
      p_these    <- p_idx$col[p_idx$i %in% filt_idx]
      n_reps_m   <- n_reps_vec[filt_idx]
      phi_mat    <- samples_combined[, phi_these, drop = FALSE]
      p_mat      <- samples_combined[, p_these,   drop = FALSE]
      p_al1_mat  <- phi_mat * (1 - sweep(1 - p_mat, 2, n_reps_m, `^`))
      p_al1_draw <- rowMeans(p_al1_mat)
      summ       <- post_summary(p_al1_draw)
    }
    summ$family <- fam
    summ$method <- meth
    summ
  }))
})) %>%
  mutate(
    family       = factor(family, levels = focal_families),
    method       = factor(method, levels = method_order),
    family_label = family_labels[as.character(family)]
  )

# =============================================================================
# 5. OBSERVED RAW DETECTION RATES
# =============================================================================

obs_rates <- do.call(rbind, lapply(focal_families, function(fam) {
  f_idx <- which(focal_families == fam)
  do.call(rbind, lapply(method_order, function(meth) {
    filt_idx <- which(filter_meta$method_group == meth)
    n_filt   <- length(filt_idx)
    if (n_filt == 0) {
      data.frame(family=fam, method=meth, n_filt=0L,
                 n_filt_pos=NA_real_, obs_filt_rate=NA_real_)
    } else {
      n_pos <- sum(n_det_mat[filt_idx, f_idx] > 0)
      data.frame(family=fam, method=meth, n_filt=n_filt,
                 n_filt_pos=n_pos, obs_filt_rate=n_pos/n_filt)
    }
  }))
})) %>%
  mutate(
    family = factor(family, levels = focal_families),
    method = factor(method, levels = method_order)
  )

# =============================================================================
# 6. METHOD CONTRAST POSTERIORS (d_AB, d_C)
# =============================================================================

contrast_params <- c("d_AB", "d_C")
contrast_labels <- c(
  "d_AB" = "MOSAIC-AB vs MOSAIC-D",
  "d_C"  = "MOSAIC-C vs MOSAIC-D"
)
contrast_colours <- c(
  "d_AB" = method_colours["MOSAIC-AB"],
  "d_C"  = method_colours["MOSAIC-C"]
)

contrast_dens <- do.call(rbind, lapply(contrast_params, function(p) {
  if (!p %in% colnames(samples_combined)) return(NULL)
  x    <- samples_combined[, p]
  dens <- density(x, n = 512)
  data.frame(param  = p,
             label  = contrast_labels[p],
             x_val  = dens$x,
             y_val  = dens$y,
             stringsAsFactors = FALSE)
})) %>% mutate(param = factor(param, levels = contrast_params))

contrast_summ <- do.call(rbind, lapply(contrast_params, function(p) {
  if (!p %in% colnames(samples_combined)) return(NULL)
  x <- samples_combined[, p]
  cbind(data.frame(param = p, label = contrast_labels[p],
                   stringsAsFactors = FALSE),
        post_summary(x))
})) %>% mutate(param = factor(param, levels = contrast_params))

# =============================================================================
# 7. COVARIATE EFFECTS
# =============================================================================

covariate_params <- c("d_vol", "d_dist", "d_time", "d_cruise")
covariate_labels <- c(
  "d_vol"    = "Volume filtered\n(presence/capture)",
  "d_dist"   = "Distance along transect\n(presence/capture)",
  "d_time"   = "Filtration time\n(PCR detection)",
  "d_cruise" = "Cruise SR2408 vs RS2311\n(presence/capture)"
)

covariate_summ <- do.call(rbind, lapply(covariate_params, function(p) {
  if (!p %in% colnames(samples_combined)) return(NULL)
  cbind(data.frame(param = p, label = covariate_labels[p],
                   stringsAsFactors = FALSE),
        post_summary(samples_combined[, p]))
})) %>%
  mutate(
    param = factor(param, levels = rev(covariate_params)),
    label = factor(label, levels = rev(covariate_labels))
  )

# =============================================================================
# 8. VOLUME × PHI SCATTER
# =============================================================================

phi_mean_mat <- matrix(NA_real_, nrow = n_filters, ncol = n_families)
for (f in seq_len(n_families)) {
  these_cols <- phi_idx$col[phi_idx$f == f]
  fi         <- phi_idx$i[phi_idx$f == f]
  ord        <- order(fi)
  phi_mean_mat[fi[ord], f] <- colMeans(
    samples_combined[, these_cols[ord], drop = FALSE]
  )
}

vol_phi_df <- do.call(rbind, lapply(seq_len(n_families), function(f) {
  data.frame(
    family    = focal_families[f],
    filter_id = filter_ids,
    method    = filter_meta$Collection.method,
    vol_mL    = as.numeric(filter_meta$Volume.filtered..mL.),
    phi_mean  = phi_mean_mat[, f],
    stringsAsFactors = FALSE
  )
})) %>%
  filter(!is.na(phi_mean)) %>%
  mutate(
    family       = factor(family, levels = focal_families),
    method       = factor(method, levels = c("MOSAIC-A","MOSAIC-B","MOSAIC-C","MOSAIC-D")),
    family_label = family_labels[as.character(family)]
  )

# Filter-level average phi across families (for scatter)
vol_phi_avg <- vol_phi_df %>%
  group_by(filter_id, method, vol_mL) %>%
  summarise(phi_mean = mean(phi_mean, na.rm = TRUE), .groups = "drop")

# =============================================================================
# FIGURE 1: Method contrast posteriors (d_AB, d_C)
# =============================================================================

y_rug <- -0.04 * max(contrast_dens$y_val)

fig1 <- ggplot() +
  geom_ribbon(data = contrast_dens,
              aes(x = x_val, ymin = 0, ymax = y_val, fill = param),
              alpha = 0.35) +
  geom_line(data = contrast_dens,
            aes(x = x_val, y = y_val, colour = param),
            linewidth = 0.8) +
  geom_errorbarh(data = contrast_summ,
                 aes(y = y_rug, xmin = lo95, xmax = hi95, colour = param),
                 height = 0, linewidth = 1.2) +
  geom_errorbarh(data = contrast_summ,
                 aes(y = y_rug, xmin = lo50, xmax = hi50, colour = param),
                 height = 0, linewidth = 2.2) +
  geom_point(data = contrast_summ,
             aes(x = median, y = y_rug, colour = param),
             size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values   = unname(contrast_colours), guide = "none") +
  scale_colour_manual(values = unname(contrast_colours), guide = "none") +
  labs(
    title    = "Method contrasts: presence/capture log-odds vs MOSAIC-D",
    subtitle = "Posterior density with median (point), 50% CI (thick), 95% CI (thin)",
    x        = "Log-odds contrast",
    y        = "Posterior density"
  ) +
  theme_mosaic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave("viz/Bayesian_occupancy_modeling/additional_plots/fig1_method_contrasts.pdf", fig1,
       width = 6, height = 5, device = cairo_pdf)
cat("Saved: viz/Bayesian_occupancy_modeling/additional_plots/fig1_method_contrasts.pdf\n")

# =============================================================================
# FIGURE 2: Marginal phi by method and family
# =============================================================================

fig2 <- ggplot(phi_method_family %>% filter(!is.na(mean)),
               aes(x = method, y = mean, colour = method)) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.6) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 1.8) +
  geom_point(aes(shape = method), size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60") +
  facet_wrap(~ family_label, nrow = 1) +
  scale_colour_manual(values = method_colours, name = NULL) +
  scale_shape_manual(
    values = c("MOSAIC-AB" = 16, "MOSAIC-C" = 15, "MOSAIC-D" = 18),
    name   = NULL
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    title    = "Marginal presence/capture probability (\u03d5) by method and family",
    subtitle = "Posterior mean with 50% (thick) and 95% (thin) credible intervals",
    x        = NULL,
    y        = "\u03d5  (presence/capture probability)"
  ) +
  theme_mosaic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("viz/Bayesian_occupancy_modeling/additional_plots/fig2_phi_by_method_family.pdf", fig2,
       width = 11, height = 5, device = cairo_pdf)
cat("Saved: viz/Bayesian_occupancy_modeling/additional_plots/fig2_phi_by_method_family.pdf\n")

# =============================================================================
# FIGURE 3: Operational detection probability
# =============================================================================

fig3 <- ggplot(op_detect_method_family %>% filter(!is.na(mean)),
               aes(x = method, y = mean, colour = method)) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.6) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 1.8) +
  geom_point(aes(shape = method), size = 3) +
  geom_point(
    data = obs_rates %>%
      left_join(select(op_detect_method_family, family, method, family_label),
                by = c("family", "method")),
    aes(x = method, y = obs_filt_rate),
    shape = 2, size = 2.5, stroke = 0.8, colour = "black", alpha = 0.7
  ) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60") +
  facet_wrap(~ family_label, nrow = 1) +
  scale_colour_manual(values = method_colours, name = NULL) +
  scale_shape_manual(
    values = c("MOSAIC-AB" = 16, "MOSAIC-C" = 15, "MOSAIC-D" = 18),
    name   = NULL
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    title    = "Operational detection probability by method and family",
    subtitle = "P(\u22651 positive replicate | filter) — model posterior (coloured) vs raw rate (open triangles)",
    x        = NULL,
    y        = "P(\u22651 positive replicate)"
  ) +
  theme_mosaic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("viz/Bayesian_occupancy_modeling/additional_plots/fig3_operational_detection.pdf", fig3,
       width = 11, height = 5, device = cairo_pdf)
cat("Saved: viz/Bayesian_occupancy_modeling/additional_plots/fig3_operational_detection.pdf\n")

# =============================================================================
# FIGURE 4: Covariate effect sizes
# =============================================================================

fig4 <- ggplot(covariate_summ, aes(y = label, x = median)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95),
                 height = 0, linewidth = 0.7, colour = "grey30") +
  geom_errorbarh(aes(xmin = lo50, xmax = hi50),
                 height = 0, linewidth = 2, colour = "grey30") +
  geom_point(aes(fill = p_pos > 0.95 | p_pos < 0.05),
             shape = 21, size = 3.5, colour = "grey30") +
  scale_fill_manual(
    values = c("TRUE" = "#E84855", "FALSE" = "white"),
    name   = "95% CI excludes zero",
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  labs(
    title    = "Covariate effects on presence/capture and detection",
    subtitle = "Posterior median with 50% (thick) and 95% (thin) CIs; red = CI excludes zero",
    x        = "Posterior log-odds coefficient",
    y        = NULL
  ) +
  theme_mosaic() +
  theme(legend.position = "right")

ggsave("viz/Bayesian_occupancy_modeling/additional_plots/fig4_covariate_effects.pdf", fig4,
       width = 7.5, height = 4, device = cairo_pdf)
cat("Saved: viz/Bayesian_occupancy_modeling/additional_plots/fig4_covariate_effects.pdf\n")

# =============================================================================
# FIGURE 5: Volume × phi scatter
# =============================================================================

fig5 <- ggplot(vol_phi_avg,
               aes(x = vol_mL / 1000, y = phi_mean, colour = method)) +
  geom_point(aes(shape = method), alpha = 0.55, size = 2) +
  geom_smooth(
    data   = vol_phi_avg %>% filter(method %in% c("MOSAIC-A", "MOSAIC-B")),
    method = "loess", span = 0.7, se = TRUE,
    colour = method_colours["MOSAIC-AB"],
    fill   = method_colours["MOSAIC-AB"],
    alpha  = 0.15, linewidth = 0.9
  ) +
  scale_colour_manual(values = method_colours_raw, name = NULL) +
  scale_shape_manual(
    values = c("MOSAIC-A" = 16, "MOSAIC-B" = 17, "MOSAIC-C" = 15, "MOSAIC-D" = 18),
    name   = NULL
  ) +
  scale_x_log10(labels = scales::label_number(suffix = " L", scale = 1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    title    = "Filtered volume vs posterior presence/capture probability",
    subtitle = "Each point = one filter (phi averaged across families); LOESS for MOSAIC-A/B",
    x        = "Volume filtered (litres, log scale)",
    y        = expression(paste("Mean ", phi, " across families"))
  ) +
  theme_mosaic()

ggsave("viz/Bayesian_occupancy_modeling/additional_plots/fig5_volume_phi_scatter.pdf", fig5,
       width = 8, height = 5, device = cairo_pdf)
cat("Saved: viz/Bayesian_occupancy_modeling/additional_plots/fig5_volume_phi_scatter.pdf\n")

# =============================================================================
# FIGURE 6: Observed vs model-estimated filter detection rates
# =============================================================================

fig6_df <- phi_method_family %>%
  filter(!is.na(mean)) %>%
  left_join(obs_rates %>% select(family, method, obs_filt_rate),
            by = c("family", "method"))

fig6 <- ggplot(fig6_df, aes(x = obs_filt_rate, y = mean, colour = method)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.7) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.6, alpha = 0.6) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 1.8, alpha = 0.8) +
  geom_point(aes(shape = method), size = 3.5) +
  geom_text(aes(label = gsub("\n.*", "", family_label)),
            size = 2.8, family = "serif", hjust = -0.15, show.legend = FALSE) +
  scale_colour_manual(values = method_colours, name = NULL) +
  scale_shape_manual(
    values = c("MOSAIC-AB" = 16, "MOSAIC-C" = 15, "MOSAIC-D" = 18),
    name   = NULL
  ) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    title    = "Observed vs model-estimated filter-level detection rate",
    subtitle = "Points above the 1:1 line = occupancy correction lifts estimate above raw rate",
    x        = "Observed filter detection rate (raw)",
    y        = expression(paste("Model estimate: mean ", phi))
  ) +
  theme_mosaic()

ggsave("viz/Bayesian_occupancy_modeling/additional_plots/fig6_observed_vs_model.pdf", fig6,
       width = 7.5, height = 6, device = cairo_pdf)
cat("Saved: viz/Bayesian_occupancy_modeling/additional_plots/fig6_observed_vs_model.pdf\n")

# =============================================================================
# NUMERICAL SUMMARIES
# =============================================================================

cat("\n=== Operational detection probability: posterior mean [95% CI] ===\n")
op_wide <- op_detect_method_family %>%
  filter(!is.na(mean)) %>%
  mutate(fmt = sprintf("%.2f [%.2f, %.2f]", mean, lo95, hi95)) %>%
  select(family, method, fmt) %>%
  pivot_wider(names_from = method, values_from = fmt)
print(as.data.frame(op_wide))

cat("\n=== Method contrasts: log-odds vs MOSAIC-D ===\n")
print(as.data.frame(contrast_summ %>%
  select(label, median, lo95, hi95, p_pos)))

cat("\nAll figures saved to viz/Bayesian_occupancy_modeling/additional_plots/\n")

