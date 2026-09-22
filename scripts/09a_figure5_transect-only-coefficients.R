# =============================================================================
# Manuscript panels — Fig 4 (Stage 2 / v5 model), Fig-6-style coefficient panels
# =============================================================================
# Clean individual panels (NO titles/subtitles/captions) for the author to
# assemble. Posterior median with 50% + 95% credible intervals, on the log-odds
# (logit) scale, dashed reference line at 0. Reference state = SR2408 (Summer)
# surface filter at mean volume/span, offshore, average site.
# Panels:
#   panel_v5_capture_coefs   : cruise (Fall vs Summer), volume/SD, distance/SD
#   panel_v5_detection_coefs : detection intercept (delta), filtration time/SD
#   panel_v5_family_capture  : per-family baseline capture (alpha)
#   panel_v5_family_inshore  : per-family inshore offset (b_inshore)
# Output: outputs/v5/figures/panels/panel_v5_*.pdf
# =============================================================================

suppressMessages(library(ggplot2))
dir.create("outputs/v5/figures/panels", recursive = TRUE, showWarnings = FALSE)
load("outputs/v5/MOSAIC_mammal_families_twolevel_v5_results.RData")
S <- samples_combined
fam <- focal_families

# posterior interval table for a set of named columns -> tidy df (label order preserved)
ivals <- function(cols, labels) {
  do.call(rbind, Map(function(cn, lab) {
    x <- S[, cn]
    data.frame(label = lab,
               med = median(x),
               lo95 = quantile(x, .025), hi95 = quantile(x, .975),
               lo50 = quantile(x, .25),  hi50 = quantile(x, .75),
               row.names = NULL)
  }, cols, labels))
}

coef_panel <- function(df, xlab) {
  df$label <- factor(df$label, levels = rev(df$label))   # first row on top
  ggplot(df, aes(med, label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.5) +
    geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 1.4) +
    geom_point(size = 2.4, colour = "#0072B2") +
    labs(x = xlab, y = NULL) +
    theme_bw() + theme(panel.grid.minor = element_blank())
}

# Concrete meaning of "per SD": volume/distance are log10-then-z (per fold-change),
# filtration time is z-scored on raw minutes (per N minutes).
sd_vol_fold  <- 10^sd(log10(as.numeric(filter_md$Volume.filtered..mL.)))
sd_dist_fold <- 10^sd(log10(as.numeric(filter_md$dist_km)))
sd_time_min  <- sd(as.numeric(filter_md$filtration_time_min))

cap <- ivals(c("d_cruise", "d_vol", "d_dist"),
             c("Cruise (Fall vs Summer)",
               sprintf("Volume (per SD, %.1f-fold)",   sd_vol_fold),
               sprintf("Distance (per SD, %.1f-fold)", sd_dist_fold)))
det <- ivals(c("delta", "d_time"),
             c("Detection intercept",
               sprintf("Filtration time (per SD, %.0f min)", sd_time_min)))
al_cols <- sprintf("alpha[%d]", seq_along(fam))
bi_cols <- sprintf("b_inshore[%d]", seq_along(fam))
alp <- ivals(al_cols, fam)
bin <- ivals(bi_cols, fam)

specs <- list(
  panel_v5_capture_coefs   = list(cap, "log-odds of capture"),
  panel_v5_detection_coefs = list(det, "log-odds of detection"),
  panel_v5_family_capture  = list(alp, "log-odds of capture (baseline)"),
  panel_v5_family_inshore  = list(bin, "inshore offset (log-odds)"))
panels <- lapply(specs, function(s) coef_panel(s[[1]], s[[2]]))
for (nm in names(panels)) {
  ggsave(sprintf("outputs/v5/figures/panels/%s.pdf", nm), panels[[nm]],
         width = 4.6, height = 2.9, device = grDevices::cairo_pdf)
  cat("Saved:", nm, "\n")
}
saveRDS(panels, "outputs/v5/figures/panels/fig4_objs.rds")   # for assembly
# echo the numbers behind the panels
cat("\n-- capture (log-odds) --\n"); print(cap, row.names = FALSE)
cat("\n-- detection (log-odds) --\n"); print(det, row.names = FALSE)
