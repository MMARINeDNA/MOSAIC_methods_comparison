# =============================================================================
# Table 2 — MOSAIC-vs-CTD capture odds ratio: primary model + three sensitivity
# analyses (as deployed in the manuscript).
# =============================================================================
# Re-derives the MOSAIC-vs-CTD capture odds ratio for each model from its saved
# posterior and renders the table (flextable -> docx) + CSV.
# OR is the posterior median (robust for an odds ratio); 95% CrI = 2.5/97.5
# quantiles of exp(coef); P = posterior prob the contrast exceeds 0.
# =============================================================================
suppressMessages({library(flextable); library(officer)})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
dir.create("model_evolution_tracking", recursive = TRUE, showWarnings = FALSE)

models <- list(
  list(label="Primary", coef="d_method",
       file="outputs/v6/MOSAIC_mammal_families_twolevel_stage1_AvsD_v6_results.RData",
       samples="86 A, 26 CTD (112)",
       contrast="MOSAIC-A vs CTD — all data"),
  list(label="Balanced subset", coef="d_method",
       file="outputs/v6/MOSAIC_mammal_families_twolevel_stage1_AvsD_v6_matched_results.RData",
       samples="26 A, 26 CTD (52)",
       contrast="MOSAIC-A vs CTD — cruise × inshore/offshore-balanced 1:1 subset"),
  list(label="Full Dataset", coef="d_method",
       file="outputs/v6/MOSAIC_mammal_families_twolevel_stage1_allMOSAIC_v6_results.RData",
       samples="86 A, 20 B, 6 C, 26 CTD (138)",
       contrast="MOSAIC-A/B/C vs CTD — all surface methods pooled"),
  list(label="Fall 2023 pilot", coef="d_method",
       file="outputs/v5/MOSAIC_mammal_families_twolevel_stage1_AvsD_results.RData",
       samples="29 A, 20 CTD (49)",
       contrast="MOSAIC-A vs CTD — Fall 2023 pilot (Summer CTDs excluded)")
)

row_of <- function(m) {
  e <- new.env(); load(m$file, envir = e)
  sc <- if (exists("samples_combined", envir=e)) e$samples_combined else do.call(rbind, lapply(e$samples, as.matrix))
  x  <- as.numeric(sc[, m$coef]); orv <- exp(x)
  data.frame(
    Model    = m$label,
    Samples  = m$samples,
    Contrast = m$contrast,
    `Odds ratio` = sprintf("%.1f", median(orv)),
    `95% CrI` = sprintf("%.1f–%.1f", quantile(orv, .025), quantile(orv, .975)),
    `P(MOSAIC > CTD)` = sprintf("%.3f", mean(x > 0)),
    check.names = FALSE, stringsAsFactors = FALSE)
}
tab <- do.call(rbind, lapply(models, row_of))

write.csv(tab, "model_evolution_tracking/Table2_MOSAIC_vs_CTD.csv", row.names = FALSE)

foot <- paste0(
  "Odds ratios are posterior medians; 95% credible intervals are 2.5/97.5 posterior ",
  "quantiles; P is the posterior probability the contrast exceeds zero. All models are ",
  "two-level hierarchical Bayesian occupancy models (NIMBLE) with family random effects, ",
  "an inshore/offshore offset, and a concurrent-sampling site random intercept. The primary, ",
  "balanced-subset, and full-dataset models span both cruises and include a fixed ",
  "Fall-vs-Summer cruise contrast; the Fall 2023 pilot is single-cruise.")

ft <- flextable(tab)
ft <- set_caption(ft, "Table 2. MOSAIC-vs-CTD capture odds ratio across the primary model and three sensitivity analyses.")
ft <- add_footer_lines(ft, foot)
ft <- theme_booktabs(ft)
ft <- fontsize(ft, size = 9, part = "all")
ft <- fontsize(ft, size = 8, part = "footer")
ft <- align(ft, j = c("Odds ratio","95% CrI","P(MOSAIC > CTD)"), align = "center", part = "all")
ft <- autofit(ft)

save_as_docx(`Table 2` = ft, path = "model_evolution_tracking/Table2_MOSAIC_vs_CTD.docx")

cat("Wrote:\n  model_evolution_tracking/Table2_MOSAIC_vs_CTD.docx\n",
    "  model_evolution_tracking/Table2_MOSAIC_vs_CTD.csv\n", sep="")
print(tab, row.names = FALSE)
