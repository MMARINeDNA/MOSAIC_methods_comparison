# =============================================================================
# Assemble the clean panels (saved as .rds by the panel scripts) into the
# final main-text multi-panel figures, with A/B/C/D tags. Vector cairo_pdf.
#   Fig 5 = transect-only (Stage-2 / v5) model coefficients (2 x 2)
#   Fig 7 = Q3 transect source (B-vs-A Bray-Curtis, PCoA, A/B chain map, taxon corroboration)
# Note: Fig 4 (v6 PRIMARY A-vs-CTD capture probability + OR) is produced standalone by
#       08.5_stage1_v6primary_figure.R. The species-accumulation figure was cut from the
#       final manuscript, so its assembly is removed here.
# Requires: run panels_fig4_v5_coefficients.R and the four 08.6_* panel scripts first
#           (they write the *_obj(s).rds consumed below).
# =============================================================================
suppressMessages({library(ggplot2); library(patchwork)})
d  <- "outputs/v5/figures/panels/"
rd <- function(f) readRDS(paste0(d, f))
P5 <- rd("fig4_objs.rds")   # Stage-2 / v5 coefficient panels (filename kept from panel script)
corr <- rd("corrob_obj.rds"); bc <- rd("bc_obj.rds"); pcoa <- rd("pcoa_obj.rds"); mp <- rd("abmap_obj.rds")

# --- Fig 5: transect-only (Stage-2 / v5) model coefficients (2 x 2), black text ---
fig5 <- ((P5$panel_v5_capture_coefs | P5$panel_v5_detection_coefs) /
        (P5$panel_v5_family_capture | P5$panel_v5_family_inshore)) +
        plot_annotation(tag_levels = "A") &
        theme(axis.text = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              strip.text = element_text(colour = "black"),
              legend.text = element_text(colour = "black"),
              legend.title = element_text(colour = "black"),
              plot.tag = element_text(colour = "black"),
              plot.title = element_text(colour = "black"),
              plot.subtitle = element_text(colour = "black"))
ggsave("outputs/v5/figures/v5_surface/fig5_model.pdf", fig5,
       width = 11, height = 8, device = grDevices::cairo_pdf)
ggsave("outputs/v5/figures/v5_surface/fig5_model.png", fig5,
       width = 11, height = 8, dpi = 300)

# --- Fig 7: Q3 transect source. Stacked so the map (C) and corroboration (D) each
# get a full-width row; the smaller BC (A) + PCoA (B) share the top row. The top-left
# panel is free()d on the left so panel A's y-axis is NOT forced to align with the wide
# left labels of the map/corroboration below (that alignment otherwise leaves a big gap
# in A). The map (C) is coord_fixed (aspect-locked) and therefore width-limited, so its
# "Family" key is moved to the bottom to give the map full width. Chain key lives on
# PCoA (B), suppressed on BC (A) and corroboration (D). Saved as vector PDF + 300-dpi PNG.
top7 <- (free(bc + guides(colour = "none"), side = "l") | pcoa) +
        plot_layout(widths = c(1, 1.55))   # give PCoA (B) >= BC (A) plot width
# Panel order: A/B (top), position-along-transect corroboration (C), then the map (D)
# enlarged full-width for readability. free()-ing the map on the left stops it aligning
# to panel C's family-label column (which otherwise left-pads the map with whitespace).
fig7 <- (top7 /
        (corr + guides(colour = "none")) /
        free(mp, side = "l")) +
        plot_layout(heights = c(0.8, 0.7, 2.6)) +
        plot_annotation(tag_levels = "A") &
        theme(axis.text = element_text(colour = "black"),   # black text for readability
              axis.title = element_text(colour = "black"),
              strip.text = element_text(colour = "black"),
              legend.text = element_text(colour = "black"),
              legend.title = element_text(colour = "black"),
              plot.tag = element_text(colour = "black"),
              plot.title = element_text(colour = "black"))
ggsave("outputs/v5/figures/transect_source/fig7_q3_transect_source.pdf", fig7,
       width = 11, height = 15, device = grDevices::cairo_pdf)
ggsave("outputs/v5/figures/transect_source/fig7_q3_transect_source.png", fig7,
       width = 11, height = 15, dpi = 300)

cat("Assembled: fig5_model.pdf, fig7_q3_transect_source.pdf\n")
