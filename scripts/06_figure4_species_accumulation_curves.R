# =============================================================================
# Species Accumulation Curves — MOSAIC eDNA Survey
# =============================================================================
# Four panels (2x2) of cumulative vertebrate taxa vs. effort:
#   (a) Number of PCR replicates
#   (b) Cumulative volume filtered (L)
#   (c) Cumulative filtration time (min)
#   (d) Cumulative survey distance (km)  [MOSAIC-D excluded: dist_km == 0]
#
# Two output PDFs:
#   species_accumulation_curves.pdf    — MOSAIC-A subsampled to max effort
#                                        of any other method (fair comparison)
#   species_accumulation_curves_v2.pdf — Full MOSAIC-A effort (unsubsampled)
#
# Inputs:
#   MOSAIC_taxnamed_taxa_table_metadata_filtered.csv
#   MOSAIC_taxnamed_otu_table_metadata_filtered.csv
#   MOSAIC_taxnamed_sample_data_metadata_filtered.csv
#
# Required: dplyr, tidyr, purrr, stringr, tibble, ggplot2, patchwork
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

set.seed(42)
N_PERMS <- 200
N_BINS  <- 50

# ---- 1. Load data ------------------------------------------------------------

taxa_tbl   <- read.csv("outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv",   row.names = 1)
otu_tbl    <- read.csv("outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv",    row.names = 1)
sample_tbl <- read.csv("outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv",  row.names = 1)

# ---- 2. Filter vertebrates (Chordata, excl. Appendicularia) -----------------

vert_taxa <- rownames(taxa_tbl)[
  taxa_tbl$phylum == "Chordata" & taxa_tbl$class != "Appendicularia"
]
cat(sprintf("Vertebrate taxa retained: %d\n", length(vert_taxa)))

otu_vert <- otu_tbl[rownames(otu_tbl) %in% vert_taxa, ]
pa_mat   <- t((otu_vert > 0) * 1L)
# read.csv() converts "-" to "." in colnames; restore to match sample_tbl rownames
rownames(pa_mat) <- gsub("\\.", "-", rownames(pa_mat))

# ---- 3. Join metadata -------------------------------------------------------
# sample_tbl rownames = PCR-replicate IDs matching OTU columns.
# sample_tbl$SampleID is a DIFFERENT cruise-level column — do not use
# rownames_to_column("SampleID") on sample_tbl or it will duplicate that name.
meta_join <- data.frame(
  .key                 = rownames(sample_tbl),
  Collection.method    = ifelse(sample_tbl[["Collection.method"]] == "MOSAIC-D", "CTD", sample_tbl[["Collection.method"]]),
  Volume.filtered..mL. = sample_tbl[["Volume.filtered..mL."]],
  filtration_time_min  = sample_tbl[["filtration_time_min"]],
  dist_km              = sample_tbl[["dist_km"]],
  stringsAsFactors     = FALSE
)

pa_df <- as.data.frame(pa_mat) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta_join, by = c("SampleID" = ".key")) %>%
  filter(!is.na(Collection.method))

cat(sprintf("Samples in pa_df: %d\n", nrow(pa_df)))

METHODS <- c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C", "CTD")

# ---- 4. Accumulation functions ----------------------------------------------

# Discrete: PCR replicates
accum_discrete <- function(pa_matrix, n_perms = 200, cap = NULL) {
  n     <- nrow(pa_matrix)
  use_n <- if (!is.null(cap)) min(n, cap) else n
  arr   <- as.matrix(pa_matrix)
  result <- matrix(0, nrow = n_perms, ncol = use_n)
  for (p in seq_len(n_perms)) {
    idx    <- sample(n)[seq_len(use_n)]
    cum_pa <- apply(arr[idx, , drop = FALSE], 2, cumsum)
    if (!is.matrix(cum_pa)) cum_pa <- matrix(cum_pa, ncol = 1)
    result[p, ] <- rowSums(cum_pa > 0)
  }
  tibble(x = seq_len(use_n), mean_sp = colMeans(result), sd_sp = apply(result, 2, sd))
}

# Continuous: volume, time, or distance — interpolated onto common x grid.
# Prepends (0,0) so approx() always has >= 2 points even with 1 sample.
accum_continuous <- function(pa_matrix, effort, n_perms = 200, n_bins = 50, cap = NULL) {
  effort[is.na(effort)] <- 0
  total <- sum(effort)
  if (!is.finite(total) || total == 0) return(NULL)
  use_x    <- if (!is.null(cap) && is.finite(cap)) min(total, cap) else total
  x_common <- seq(use_x / n_bins, use_x, length.out = n_bins)
  n   <- nrow(pa_matrix)
  arr <- as.matrix(pa_matrix)
  sp_at_x <- matrix(0, nrow = n_perms, ncol = n_bins)
  for (p in seq_len(n_perms)) {
    idx    <- sample(n)
    cum_ef <- cumsum(effort[idx])
    cum_pa <- apply(arr[idx, , drop = FALSE], 2, cumsum)
    if (!is.matrix(cum_pa)) cum_pa <- matrix(cum_pa, ncol = 1)
    cum_sp <- rowSums(cum_pa > 0)
    sp_at_x[p, ] <- approx(c(0, cum_ef), c(0, cum_sp), xout = x_common, rule = 2)$y
  }
  tibble(x = x_common, mean_sp = colMeans(sp_at_x), sd_sp = apply(sp_at_x, 2, sd))
}

# ---- 5. Caps and curve computation ------------------------------------------

non_A <- c("MOSAIC-B", "MOSAIC-C", "MOSAIC-D")

cap_pcr <- pa_df %>%
  filter(Collection.method %in% non_A) %>%
  count(Collection.method) %>% pull(n) %>% max()

# Continuous caps from filter-level data (collapse PCR reps by FilterID first)
pa_df_filt <- pa_df %>%
  mutate(FilterID = str_remove(SampleID, "-\\d+$")) %>%
  group_by(Collection.method, FilterID) %>%
  summarise(
    Volume.filtered..mL. = first(Volume.filtered..mL.),
    filtration_time_min  = first(filtration_time_min),
    dist_km              = first(dist_km),
    .groups = "drop"
  )

cap_vol <- pa_df_filt %>%
  filter(Collection.method %in% non_A) %>%
  group_by(Collection.method) %>%
  summarise(tot = sum(Volume.filtered..mL., na.rm = TRUE) / 1000) %>%
  pull(tot) %>% max()

cap_min <- pa_df_filt %>%
  filter(Collection.method %in% non_A) %>%
  group_by(Collection.method) %>%
  summarise(tot = sum(filtration_time_min, na.rm = TRUE)) %>%
  pull(tot) %>% max()

# km cap excludes CTD/MOSAIC-D (all dist_km == 0)
cap_km <- pa_df_filt %>%
  filter(Collection.method %in% non_A, dist_km > 0) %>%
  group_by(Collection.method) %>%
  summarise(tot = sum(dist_km, na.rm = TRUE)) %>%
  pull(tot) %>% max()

cat(sprintf("X-axis caps — PCR: %d | Vol: %.1f L | Min: %.0f min | Dist: %.1f km\n",
            cap_pcr, cap_vol, cap_min, cap_km))

compute_curves <- function(subsample_A = TRUE) {
  purrr::map_dfr(METHODS, function(method) {
    sub     <- pa_df %>% filter(Collection.method == method)
    sp_cols <- intersect(colnames(sub), vert_taxa)
    pa_m    <- sub %>% select(all_of(sp_cols))
    # (a) PCR replicates — one row per PCR rep, effort unit = 1 rep
    pcr <- accum_discrete(pa_m, N_PERMS,
                          cap = if (subsample_A && method == "MOSAIC-A") cap_pcr else NULL) %>%
      mutate(panel = "pcr", method = method)

    # (b–d) Continuous panels — collapse PCR reps to filter level first.
    # Volume, time, and distance are filter-level measurements; summing across
    # PCR reps inflates effort ~4x. Union species detections within each filter.
    sub_filt <- sub %>%
      mutate(FilterID = str_remove(SampleID, "-\\d+$")) %>%
      group_by(FilterID) %>%
      summarise(
        across(all_of(sp_cols), \(x) as.integer(max(x, na.rm = TRUE) > 0)),
        vol_L = first(Volume.filtered..mL.) / 1000,
        mins  = first(filtration_time_min),
        dist  = first(dist_km),
        .groups = "drop"
      )
    pa_filt   <- sub_filt %>% select(all_of(sp_cols))
    vol_L_f   <- sub_filt$vol_L
    mins_f    <- sub_filt$mins
    dist_f    <- sub_filt$dist

    # (b) Volume
    vol <- accum_continuous(pa_filt, vol_L_f, N_PERMS, N_BINS,
                            cap = if (subsample_A && method == "MOSAIC-A") cap_vol else NULL) %>%
      mutate(panel = "vol", method = method)

    # (c) Filtration time (min)
    tim <- accum_continuous(pa_filt, mins_f, N_PERMS, N_BINS,
                            cap = if (subsample_A && method == "MOSAIC-A") cap_min else NULL) %>%
      mutate(panel = "min", method = method)

    # (d) Survey distance (km) — CTD (MOSAIC-D) returns NULL and is omitted
    km_raw <- accum_continuous(pa_filt, dist_f, N_PERMS, N_BINS,
                               cap = if (subsample_A && method == "MOSAIC-A") cap_km else NULL)
    km <- if (!is.null(km_raw)) km_raw %>% mutate(panel = "km", method = method) else tibble()

    bind_rows(pcr, vol, tim, km)
  })
}

curves_trimmed   <- compute_curves(subsample_A = TRUE)
curves_unsampled <- compute_curves(subsample_A = FALSE)

# Fix legend order: factor levels control display sequence
method_order <- c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C", "CTD")
curves_trimmed   <- curves_trimmed   %>% mutate(method = factor(method, levels = method_order))
curves_unsampled <- curves_unsampled %>% mutate(method = factor(method, levels = method_order))

# ---- 6. Plot ----------------------------------------------------------------

METHOD_COLORS <- c(
  "MOSAIC-A" = "#1B6CA8",
  "MOSAIC-B" = "#D4752E",
  "MOSAIC-C" = "#2E9E5B",
  "CTD" = "#8B4FC7"
)

PANEL_META <- list(
  pcr = list(label = "Number of PCR replicates",        title = "(a) Per PCR Replicate"),
  vol = list(label = "Cumulative volume filtered (L)",   title = "(b) Per Liter Filtered"),
  min = list(label = "Cumulative filtration time (min)", title = "(c) Per Minute Filtered"),
  km  = list(label = "Cumulative survey distance (km)",  title = "(d) Per Kilometer")
)

make_panel <- function(data, panel_key, show_y = TRUE) {
  d    <- data %>% filter(panel == panel_key)
  meta <- PANEL_META[[panel_key]]

  p <- ggplot(d, aes(x = x, y = mean_sp, colour = method, fill = method)) +
    geom_ribbon(aes(ymin = pmax(0, mean_sp - sd_sp), ymax = mean_sp + sd_sp),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(data = d %>% group_by(method) %>% slice_tail(n = 1),
               shape = 21, size = 2.2, colour = "white", stroke = 1.2,
               aes(fill = method)) +
    scale_colour_manual(values = METHOD_COLORS, name = NULL) +
    scale_fill_manual(  values = METHOD_COLORS, name = NULL) +
    labs(title = meta$title,
         x     = meta$label,
         y     = if (show_y) "Cumulative vertebrate taxa" else NULL) +
    theme_bw(base_family = "serif", base_size = 9.5) +
    theme(
      plot.title       = element_text(size = 9.5, face = "bold", hjust = 0,
                                      margin = margin(b = 4)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e5e9ef", linewidth = 0.5,
                                      linetype = "dashed"),
      panel.border     = element_rect(colour = "#b0b8c4", linewidth = 0.7),
      axis.text        = element_text(colour = "#4a5568", size = 8),
      axis.title       = element_text(colour = "#2d3748", size = 8.5),
      axis.ticks       = element_line(colour = "#b0b8c4"),
      legend.position  = "none",
      plot.background  = element_rect(fill = "#f8f9fb", colour = NA),
      panel.background = element_rect(fill = "white")
    )

  if (max(d$x, na.rm = TRUE) >= 1000)
    p <- p + scale_x_continuous(
      labels = function(v) ifelse(v >= 1000,
                                  paste0(formatC(v / 1000, digits = 1, format = "f"), "k"),
                                  as.character(v)))
  p
}

make_figure <- function(curves, subtitle) {
  # 2x2: top=(a)PCR (b)vol  /  bottom=(c)min (d)km
  # y-axis label on left column only; MOSAIC-D absent from panel (d)
  p_pcr <- make_panel(curves, "pcr", show_y = TRUE)
  p_vol <- make_panel(curves, "vol", show_y = FALSE)
  p_min <- make_panel(curves, "min", show_y = TRUE)
  p_km  <- make_panel(curves, "km",  show_y = FALSE)

  (p_pcr | p_vol) / (p_min | p_km) +
    plot_annotation(
      title    = "MOSAIC eDNA Survey \u2014 Vertebrate Species Accumulation Curves",
      subtitle = subtitle,
      theme = theme(
        plot.title    = element_text(family = "serif", size = 11, face = "italic",
                                     hjust = 0.5, colour = "#1a2540"),
        plot.subtitle = element_text(family = "serif", size = 8,
                                     hjust = 0.5, colour = "#718096"),
        plot.background = element_rect(fill = "#f8f9fb", colour = NA)
      )
    ) &
    theme(legend.position = "bottom") &
    guides(colour = guide_legend(nrow = 1,
                                 override.aes = list(shape = 22, size = 4, stroke = 0)))
}

# ---- 7. Save ----------------------------------------------------------------

subtitle_trimmed <- paste0(
  "Mean \u00b1 1 SD \u00b7 ", N_PERMS, " permutations \u00b7 Vertebrata only \u00b7 ",
  "MOSAIC-A subsampled to max effort of other methods \u00b7 MOSAIC-D excluded from panel (d)"
)
subtitle_full <- paste0(
  "Mean \u00b1 1 SD \u00b7 ", N_PERMS, " permutations \u00b7 ",
  "Vertebrata only (Appendicularia excluded) \u00b7 MOSAIC-D excluded from panel (d)"
)

fig_trimmed   <- make_figure(curves_trimmed,   subtitle_trimmed)
fig_unsampled <- make_figure(curves_unsampled, subtitle_full)

ggsave("viz/species_accumulation_curves_zoomed.pdf",    fig_trimmed,
       width = 9, height = 8.5, units = "in", device = cairo_pdf)
ggsave("viz/species_accumulation_curves_full.pdf", fig_unsampled,
       width = 9, height = 8.5, units = "in", device = cairo_pdf)

cat("Saved species_accumulation_curves_zoomed.pdf and species_accumulation_curves_full.pdf\n")

