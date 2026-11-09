#!/usr/bin/env Rscript
# =============================================================================
# Run the whole MOSAIC analysis, in dependency order.
#   Rscript scripts/run_all.R                    # everything (hours, see scripts/README.md)
#   Rscript scripts/run_all.R --skip-models      # skip the 5 NIMBLE model scripts
#   Rscript scripts/run_all.R --from=08e         # start at the first script whose name starts with "08e"
#   Rscript scripts/run_all.R --only=08a,08c     # run just these (prefix match, comma-separated)
# Run from the repository root. Each script is run with Rscript in a fresh R
# session; stdout/stderr go to logs/<script>.log. Stops on the first error.
# With --skip-models, scripts that load a model's saved RData are skipped
# (with a message) unless that RData already exists in outputs/.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
skip_models <- "--skip-models" %in% args
from <- sub("^--from=", "", grep("^--from=", args, value = TRUE))
only <- unlist(strsplit(sub("^--only=", "", grep("^--only=", args, value = TRUE)), ","))
unknown <- setdiff(args, c("--skip-models", grep("^--(from|only)=", args, value = TRUE)))
if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = " "))

if (!file.exists("scripts/run_all.R") || !dir.exists("processed_data"))
  stop("Run from the repository root: Rscript scripts/run_all.R")

# name, is_model, RData files the script loads (checked only when models are skipped)
RD <- list(
  v6   = "outputs/v6/MOSAIC_mammal_families_twolevel_stage1_AvsD_v6_results.RData",
  v6m  = "outputs/v6/MOSAIC_mammal_families_twolevel_stage1_AvsD_v6_matched_results.RData",
  v6a  = "outputs/v6/MOSAIC_mammal_families_twolevel_stage1_allMOSAIC_v6_results.RData",
  v5p  = "outputs/v5/MOSAIC_mammal_families_twolevel_stage1_AvsD_results.RData",
  v5   = "outputs/v5/MOSAIC_mammal_families_twolevel_v5_results.RData")
steps <- list(
  list("01_csvs_to_phyloseq.R"),
  list("02_prep_detection_prob_metadata.R"),
  list("03_prep_visual_data.R"),
  list("04_MOSAIC_metadata_cleaning_for_modeling.R"),
  list("05_mosaic_figure3_familylevel.R"),
  list("06_suppfig_species-accumulation.R"),
  list("07_suppfig_nested-venn.R"),
  list("08a_model_primary_MOSAIC-A_vs_CTD.R", TRUE),
  list("08b_build_balanced-subset-pairs.R"),
  list("08b_model_A_vs_CTD_balanced-subset.R", TRUE),
  list("08c_model_all-MOSAIC_vs_CTD.R", TRUE),
  list("08d_model_A_vs_CTD_Fall-2023-pilot.R", TRUE),
  list("08e_figure4_capture-probability.R", FALSE, RD$v6),
  list("08f_table2_MOSAIC_vs_CTD.R", FALSE, unlist(RD[c("v6", "v6m", "v6a", "v5p")])),
  list("09_model_transect-only.R", TRUE),
  list("09a_figure5_transect-only-coefficients.R", FALSE, RD$v5),
  list("08.6_A_B_taxa_detection_map.R"),
  list("08.6_A_transect_source_ordination.R"),
  list("08.6_transect_PCoA_dbRDA.R"),
  list("08.6_transect_species_corroboration.R"),
  list("assemble_manuscript_figures.R", FALSE, "outputs/v5/figures/panels/fig4_objs.rds"),
  list("09.5_A-D_cetacean_detection_chisq.R"),
  list("10_rare_taxa_panels_final.R"),
  list("11_make_summary_Table_1.R"),
  list("12_total_distance_sampled.R"))
nm <- vapply(steps, `[[`, "", 1)
for (s in seq_along(steps)) if (length(steps[[s]]) < 2) steps[[s]][[2]] <- FALSE

first <- function(p) { i <- which(startsWith(nm, p)); if (!length(i)) stop("No script starts with '", p, "'"); i[1] }
sel <- seq_along(steps)
if (length(from)) sel <- sel[sel >= first(from)]
if (length(only)) sel <- intersect(sel, vapply(only, first, 1L))

for (d in c("outputs", "outputs/v5/figures/panels", "outputs/v5/figures/transect_source",
            "outputs/v5/figures/v5_surface", "outputs/v5/figures/stage1_AvsD",
            "outputs/v6/figures/stage1_v6", "outputs/v6/figures/stage1_v6_matched",
            "outputs/v6/figures/stage1_v6_allMOSAIC", "viz", "logs"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

timing <- data.frame(script = character(), status = character(), minutes = numeric())
for (i in sel) {
  f <- nm[i]; is_model <- steps[[i]][[2]]
  need <- if (length(steps[[i]]) >= 3) steps[[i]][[3]] else NULL
  if (skip_models && is_model) { message("SKIP (model) ", f); next }
  if (skip_models && !is.null(need) && !all(file.exists(need))) {
    message("SKIP (needs ", paste(need[!file.exists(need)], collapse = ", "), ") ", f); next }
  log <- file.path("logs", sub("\\.R$", ".log", f))
  message(format(Sys.time(), "%H:%M:%S"), "  RUN  ", f, "  -> ", log)
  t0 <- Sys.time()
  rc <- system2("Rscript", file.path("scripts", f), stdout = log, stderr = log)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  timing[nrow(timing) + 1, ] <- list(f, if (rc == 0) "ok" else "FAILED", mins)
  write.csv(timing, "logs/run_all_timing.csv", row.names = FALSE)
  if (rc != 0) stop(f, " failed (exit ", rc, "); see ", log)
}
print(timing)
