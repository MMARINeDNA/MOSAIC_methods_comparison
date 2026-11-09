# =============================================================================
# Q2 variable-importance comparison table (volume / distance / filtration time)
# =============================================================================
# Re-derives the volume, distance, and filtration-time coefficients (as OR per SD)
# from each Q2-relevant model posterior and renders a supplement table (flextable
# -> docx) + CSV into model_evolution_tracking/.
#   v3 (no site) / v4 (+ site): single all-method models; volume/distance/time
#     estimated ALONGSIDE the method contrasts (d_AB, d_C) and the CTD; volume and
#     distance z-scored on the RAW scale.
#   v5 Stage 2: surface-only (MOSAIC-A/B/C, CTD dropped); volume and distance
#     REPLACE the method factor and are z-scored on the LOG scale.
# OR per SD = posterior median of exp(coef); 95% CrI = 2.5/97.5 quantiles;
# P = posterior prob(effect > 0). Volume/distance are capture effects; filtration
# time is a detection effect (a low P(>0) => high prob of a negative cost effect).
# =============================================================================
suppressMessages({library(flextable); library(officer)})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
dir.create("model_evolution_tracking", recursive = TRUE, showWarnings = FALSE)

models <- list(
  list(label="v3 (no site)",
       file="outputs/MOSAIC_mammal_families_twolevel_v3_results.RData"),
  list(label="v4 (+ site)",
       file="outputs/v4/MOSAIC_mammal_families_twolevel_v4_results.RData"),
  list(label="v5 Stage 2",
       file="outputs/v5/MOSAIC_mammal_families_twolevel_v5_results.RData")
)
vars <- c(d_vol="Volume", d_dist="Distance", d_time="Filtration time")

rows <- list()
for (m in models) {
  e <- new.env(); load(m$file, envir=e)
  sc <- if (exists("samples_combined",envir=e)) e$samples_combined else do.call(rbind, lapply(e$samples, as.matrix))
  for (coef in names(vars)) {
    if (!coef %in% colnames(sc)) next
    x <- as.numeric(sc[,coef]); orv <- exp(x)
    rows[[length(rows)+1]] <- data.frame(
      Variable = vars[[coef]], Model = m$label,
      `OR per SD` = sprintf("%.2f", median(orv)),
      `95% CrI`   = sprintf("%.2f–%.2f", quantile(orv,.025), quantile(orv,.975)),
      `P(effect > 0)` = sprintf("%.3f", mean(x>0)),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  rm(e); gc()
}
df <- do.call(rbind, rows)
df$Variable <- factor(df$Variable, levels=c("Volume","Distance","Filtration time"))
df <- df[order(df$Variable), ]
write.csv(df, "model_evolution_tracking/Q2_table_variable_importance.csv", row.names=FALSE)

p_time_neg <- 1 - as.numeric(df$`P(effect > 0)`[df$Variable == "Filtration time" & df$Model == "v5 Stage 2"])
foot <- c(
 paste0("OR per SD is the posterior median of exp(coefficient); 95% CrI = 2.5/97.5 posterior ",
   "quantiles; P = posterior probability the effect exceeds zero. Volume and distance are capture-",
   "submodel effects; filtration time is a detection-submodel effect, so its low P(effect > 0) ",
   sprintf("indicates a high probability of a negative (cost) effect (v5 Stage 2 P(effect < 0) = %.2f).", p_time_neg)),
 paste0("Scaling caveat: v5 Stage 2 z-scores LOG volume and LOG transect distance, so 1 SD is a ",
   "fold-change on the log scale; v3 and v4 z-score RAW volume and distance. Per-SD magnitudes are ",
   "therefore not in identical units across models, though direction and relative ranking are comparable."),
 paste0("Estimand caveat: v3 and v4 are single all-method models in which volume and distance are ",
   "estimated alongside the MOSAIC-A/B and MOSAIC-C method contrasts and the CTD; v5 Stage 2 is ",
   "surface-only (MOSAIC-A/B/C, CTD dropped) and replaces the method factor with volume and distance, ",
   "so there those predictors carry the full between-method signal. v3 has no site effect; v4 and v5 ",
   "include the concurrent-sampling site random intercept.")
)

ft <- flextable(df)
ft <- merge_v(ft, j = "Variable")
ft <- set_caption(ft, paste0("Table S. Importance of filtered volume, transect distance, and ",
      "filtration time for rare-taxon capture/detection across model specifications."))
ft <- add_footer_lines(ft, foot)
ft <- theme_booktabs(ft)
ft <- valign(ft, j="Variable", valign="top", part="body")
ft <- fontsize(ft, size=9, part="all"); ft <- fontsize(ft, size=8, part="footer")
ft <- align(ft, j=c("OR per SD","95% CrI","P(effect > 0)"), align="center", part="all")
ft <- autofit(ft)
save_as_docx(`Table S (Q2 variable importance)` = ft,
             path = "model_evolution_tracking/Q2_variable_importance_table.docx")

cat("Wrote model_evolution_tracking/Q2_variable_importance_table.docx + .csv\n\n")
print(df, row.names=FALSE)
