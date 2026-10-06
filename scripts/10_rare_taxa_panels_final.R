# ============================================================
# Species detection maps — buffered transects, one panel per species
# Buffer: opaque at centerline, fades to transparent at 5 km
# Caption: detections per collection method within 5km of positives
# CTD (MOSAIC-D): single dot if 1 species, equal-slice pie if 2+
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(maps)
  library(scales)
  library(patchwork)
})

# ── 1. Load data ──────────────────────────────────────────────
taxa <- read.csv("outputs/MOSAIC_taxnamed_taxa_table.csv",  row.names = 1)
otu  <- read.csv("outputs/MOSAIC_taxnamed_otu_table.csv",   row.names = 1)
meta <- read.csv("outputs/MOSAIC_taxnamed_sample_data.csv", row.names = 1)
colnames(otu) <- gsub("\\.", "-", colnames(otu))

# ── 2. Panel definitions (multi-species) ─────────────────────
panels_def <- list(
  list(
    panel_id     = "baleen_whales",
    title_common = "Baleen Whales",
    species = list(
      list(taxon_id = "Megaptera novaeangliae",      label = "Megaptera novaeangliae",      common = "Humpback whale",  color = "#009E73"),
      list(taxon_id = "Balaenoptera physalus",        label = "Balaenoptera physalus",        common = "Fin whale",       color = "#56B4E9"),
      list(taxon_id = "Balaenoptera musculus",        label = "Balaenoptera musculus",        common = "Blue whale",      color = "#0072B2"),
      list(taxon_id = "Balaenoptera acutorostrata",   label = "Balaenoptera acutorostrata",   common = "Minke whale",     color = "#D55E00")
    )
  ),
  list(
    panel_id     = "seals",
    title_common = "Seals",
    species = list(
      list(taxon_id = "Arctocephalus townsendi",   label = "Arctocephalus townsendi",   common = "Guadalupe fur seal", color = "#009E73"),
      list(taxon_id = "Mirounga angustirostris",   label = "Mirounga angustirostris",   common = "N. elephant seal",   color = "#CC79A7"),
      list(taxon_id = "Callorhinus ursinus",       label = "Callorhinus ursinus",       common = "N. fur seal",        color = "#E69F00")
    )
  ),
  list(
    panel_id     = "deep_divers",
    title_common = "Deep-diving cryptic cetaceans",
    species = list(
      list(taxon_id = "Berardius bairdii",  label = "Berardius bairdii",  common = "Baird's beaked whale", color = "#CC79A7"),
      list(taxon_id = "Kogia breviceps",    label = "Kogia breviceps",    common = "Pygmy sperm whale",    color = "#0072B2"),
      list(taxon_id = "Grampus griseus",    label = "Grampus griseus",    common = "Risso's dolphin",      color = "#D55E00")
    )
  ),
  list(
    panel_id     = "sharks",
    title_common = "Sharks",
    species = list(
      list(taxon_id = "Prionace glauca",        label = "Prionace glauca",        common = "Blue shark",          color = "#0072B2"),
      list(taxon_id = "Carcharodon carcharias",  label = "Carcharodon carcharias",  common = "Great white shark",   color = "#000000"),
      list(taxon_id = "Isurus oxyrinchus",       label = "Isurus oxyrinchus",       common = "Shortfin mako shark", color = "#D55E00"),
      list(taxon_id = "Alopias vulpinus",        label = "Alopias vulpinus",        common = "Thresher shark",      color = "#009E73")
    )
  )
)

# ── 3. Metadata prep ──────────────────────────────────────────
sample_totals <- colSums(otu)

meta2 <- meta %>%
  mutate(sample_id = rownames(.)) %>%
  select(sample_id, FilterID, Collection.method, Tech.rep, Cruise,
         lat_start = Latitude.start, lon_start = Longitude.start,
         lat_end   = Latitude.end,   lon_end   = Longitude.end)

all_filters <- meta2 %>%
  filter(!is.na(lat_start), !is.na(lon_start)) %>%
  mutate(
    lat_end = coalesce(lat_end, lat_start),
    lon_end = coalesce(lon_end, lon_start),
    lat_mid = (lat_start + lat_end) / 2,
    lon_mid = (lon_start + lon_end) / 2
  ) %>%
  group_by(FilterID) %>%
  slice(1) %>%
  ungroup()

cat("Total georeferenced filters:", nrow(all_filters), "\n")

# ── 4. Buffer parameters ──────────────────────────────────────
thresh_km      <- 5
km_per_deg_lat <- 111
km_per_deg_lon <- 111 * cos(32.5 * pi / 180)

buf_5km_mm  <- (thresh_km / km_per_deg_lat) * 29 * 2
lw_core     <- 0.5
lw_max      <- lw_core + buf_5km_mm
n_layers    <- 40
layer_lw    <- seq(lw_max,  lw_core, length.out = n_layers)
layer_alpha <- seq(0,       1,       length.out = n_layers)

cat(sprintf("5km buffer: %.2f mm diameter | lw range %.2f–%.2f mm\n",
            buf_5km_mm, lw_core, lw_max))

# ── 5. Segment-to-segment distance helper ────────────────────
seg_seg_dist_km <- function(ax1, ay1, ax2, ay2,
                            bx1, by1, bx2, by2) {
  dax <- ax2 - ax1; day <- ay2 - ay1
  dbx <- bx2 - bx1; dby <- by2 - by1
  wx  <- ax1 - bx1; wy  <- ay1 - by1
  a   <- dax*dax + day*day
  e   <- dbx*dbx + dby*dby
  f   <- dbx*wx  + dby*wy
  if (a < 1e-10 && e < 1e-10) return(sqrt(wx*wx + wy*wy))
  if (a < 1e-10) {
    s <- 0; t <- max(0, min(1, f / e))
  } else {
    cv <- dax*wx + day*wy
    if (e < 1e-10) {
      t <- 0; s <- max(0, min(1, -cv / a))
    } else {
      bv    <- dax*dbx + day*dby
      denom <- a*e - bv*bv
      s <- if (abs(denom) > 1e-10) max(0, min(1, (bv*f - cv*e) / denom)) else 0
      t <- (bv*s + f) / e
      if      (t < 0) { t <- 0; s <- max(0, min(1, -cv / a)) }
      else if (t > 1) { t <- 1; s <- max(0, min(1, (bv - cv) / a)) }
    }
  }
  px <- ax1 + s*dax - (bx1 + t*dbx)
  py <- ay1 + s*day - (by1 + t*dby)
  sqrt(px*px + py*py)
}

filters_to_km <- function(df) {
  df %>% mutate(
    x1 = lon_start * km_per_deg_lon, y1 = lat_start * km_per_deg_lat,
    x2 = lon_end   * km_per_deg_lon, y2 = lat_end   * km_per_deg_lat
  )
}

# ── 6. Transect membership table (built once per cruise) ──────
# Each MOSAIC-A filter defines one transect.
# Members = the A filter itself + any B/C/D within 1 km.
# Other A filters are intentionally excluded (each A is its own transect).
build_transect_table <- function(all_filters_cr) {
  af     <- all_filters_cr %>% filter(Collection.method == "MOSAIC-A")
  non_a  <- all_filters_cr %>% filter(Collection.method != "MOSAIC-A")
  af_km  <- filters_to_km(af)
  non_a_km <- filters_to_km(non_a)

  rows <- lapply(seq_len(nrow(af_km)), function(ai) {
    # Always include the A filter itself
    self_row <- data.frame(
      transect_id       = af_km$FilterID[ai],
      FilterID          = af_km$FilterID[ai],
      Collection.method = "MOSAIC-A",
      stringsAsFactors  = FALSE
    )
    if (nrow(non_a_km) == 0) return(self_row)

    dists <- sapply(seq_len(nrow(non_a_km)), function(fi)
      seg_seg_dist_km(
        af_km$x1[ai], af_km$y1[ai], af_km$x2[ai], af_km$y2[ai],
        non_a_km$x1[fi], non_a_km$y1[fi], non_a_km$x2[fi], non_a_km$y2[fi]
      ))
    close <- non_a[dists <= 1, c("FilterID", "Collection.method")]
    if (nrow(close) == 0) return(self_row)
    close$transect_id <- af_km$FilterID[ai]
    rbind(self_row, close[, c("transect_id","FilterID","Collection.method")])
  })
  do.call(rbind, rows)
}

# ── Per-species effort-rate caption ──────────────────────────
# For each method present in the cruise, report n_detecting / n_total (%)
effort_caption <- function(sp_common, pct_df, all_filters_cr) {
  det_fids <- pct_df$FilterID[pct_df$pct_mean > 0]

  method_order <- c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C", "MOSAIC-D")
  method_lbl   <- c("MOSAIC-A"="A","MOSAIC-B"="B","MOSAIC-C"="C","MOSAIC-D"="D")

  # One row per FilterID in this cruise
  cruise_fids <- all_filters_cr %>%
    distinct(FilterID, Collection.method)

  parts <- sapply(method_order, function(m) {
    mfids <- cruise_fids$FilterID[cruise_fids$Collection.method == m]
    n_tot <- length(mfids)
    if (n_tot == 0) return(NA_character_)
    n_det <- sum(mfids %in% det_fids)
    pct   <- round(n_det / n_tot * 100, 1)
    paste0(method_lbl[m], ": ", n_det, "/", n_tot, " (", pct, "%)")
  })
  parts <- parts[!is.na(parts)]

  if (all(sapply(parts, function(p) grepl("0/", p))))
    return(paste0(sp_common, ":  not detected"))

  paste0(sp_common, ":  ", paste(parts, collapse = "   "))
}

# ── 7. Per-species pct_df helper ─────────────────────────────
get_pct_df <- function(taxon_id) {
  pct_per_sample <- as.numeric(otu[taxon_id, ]) / sample_totals * 100
  pct_per_sample[is.nan(pct_per_sample)] <- 0
  names(pct_per_sample) <- colnames(otu)
  meta2 %>%
    filter(!is.na(lat_start)) %>%
    mutate(pct = pct_per_sample[sample_id]) %>%
    group_by(FilterID) %>%
    summarise(pct_mean = mean(pct, na.rm = TRUE), .groups = "drop")
}

# ── 8. Pie wedge polygon builder ─────────────────────────────
# df must have columns: FilterID, lon_mid, lat_mid, color, prop (equal slices)
# r_deg: pie radius in geographic degrees
make_ctd_pie_polys <- function(df, r_deg = 0.12, n_pts = 60) {
  do.call(rbind, lapply(split(df, df$FilterID), function(fd) {
    fd       <- fd[order(fd$species), ]          # consistent slice order
    n_sl     <- nrow(fd)
    cx       <- fd$lon_mid[1]
    cy       <- fd$lat_mid[1]
    angle_start <- pi / 2                        # 12 o'clock, clockwise
    do.call(rbind, lapply(seq_len(n_sl), function(i) {
      sweep      <- 2 * pi / n_sl               # equal slices
      angle_end  <- angle_start - sweep
      pts        <- max(3L, ceiling(n_pts / n_sl))
      theta      <- seq(angle_start, angle_end, length.out = pts)
      xs         <- c(cx, cx + r_deg * cos(theta), cx)
      ys         <- c(cy, cy + r_deg * sin(theta), cy)
      angle_start <<- angle_end
      data.frame(
        x        = xs,
        y        = ys,
        FilterID = fd$FilterID[i],
        species  = fd$species[i],
        color    = fd$color[i],
        group_id = paste0(fd$FilterID[i], "_", i),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

# ── 9. Map-building function (multi-species) ──────────────────
world_df <- map_data("world")
ca_lon   <- c(-125, -116)
ca_lat   <- c(29,    36)

method_abbrev <- c("MOSAIC-A" = "A", "MOSAIC-B" = "B",
                   "MOSAIC-C" = "C", "MOSAIC-D" = "D")

make_panel <- function(panel, fsz = 1.0, tag = NULL, axis_pt = 6 * fsz, legend_pt = 10 * fsz,
                       title_mm = 4.5 * fsz, tag_mm = 6.3 * fsz, legend_key_h = 0.55 * legend_pt / 10) {   # sizes: points (axis, legend) / mm (panel title, panel letter) / cm (legend row height)

  title_common <- panel$title_common
  cat("\n── Panel:", title_common, "──\n")

  # ── Per-species detection data ────────────────────────────────
  sp_det_list <- lapply(panel$species, function(sp) {
    cat("  Processing:", sp$label, "\n")
    pct_df <- get_pct_df(sp$taxon_id)
    filter_data <- all_filters %>%
      left_join(pct_df, by = "FilterID") %>%
      mutate(pct_mean = replace_na(pct_mean, 0))
    det <- filter_data %>%
      filter(pct_mean > 0) %>%
      mutate(method_label = method_abbrev[Collection.method],
             sp_label     = sp$label,
             sp_color     = sp$color)
    cat("    Detections:", nrow(det), "\n")
    list(sp = sp, pct_df = pct_df, det = det, filter_data = filter_data)
  })

  # ── Perpendicular offset for transect species ─────────────────
  mm_per_deg_plot <- 29
  buf_radius_deg  <- (lw_max / 2) / mm_per_deg_plot
  n_sp            <- length(sp_det_list)
  sp_offsets      <- if (n_sp == 1) 0
                    else seq(-(n_sp - 1), (n_sp - 1), by = 2) * buf_radius_deg / 2

  apply_perp_offset <- function(det, offset) {
    if (offset == 0 || nrow(det) == 0) return(det)
    det %>% mutate(
      seg_dlon  = lon_end - lon_start,
      seg_dlat  = lat_end - lat_start,
      seg_len   = sqrt(seg_dlon^2 + seg_dlat^2),
      perp_lon  = ifelse(seg_len > 1e-6, -seg_dlat / seg_len, 0),
      perp_lat  = ifelse(seg_len > 1e-6,  seg_dlon / seg_len, 1),
      lon_start = lon_start + perp_lon * offset,
      lat_start = lat_start + perp_lat * offset,
      lon_end   = lon_end   + perp_lon * offset,
      lat_end   = lat_end   + perp_lat * offset,
      lon_mid   = lon_mid   + perp_lon * offset
    ) %>% select(-seg_dlon, -seg_dlat, -seg_len, -perp_lon, -perp_lat)
  }

  sp_det_list <- lapply(seq_along(sp_det_list), function(si) {
    x <- sp_det_list[[si]]
    x$det <- apply_perp_offset(x$det, sp_offsets[si])
    x
  })

  # ── Grey effort lines ─────────────────────────────────────────
  all_det_fids <- unique(unlist(lapply(sp_det_list, function(x) x$det$FilterID)))
  grey_filters <- all_filters %>% filter(!FilterID %in% all_det_fids)

  # ── Caption text collected but not shown on map ──────────────
  cap_lines <- sapply(sp_det_list, function(x) {
    effort_caption(x$sp$common, x$pct_df, all_filters)
  })
  caption_text <- paste(cap_lines, collapse = "\n")

  # ── Assemble CTD detection table across all species ───────────
  ctd_rows <- do.call(rbind, lapply(sp_det_list, function(x) {
    pt <- x$det %>% filter(Collection.method == "MOSAIC-D")
    if (nrow(pt) == 0) return(NULL)
    pt %>% mutate(species = x$sp$label, color = x$sp$color)
  }))

  # Split into single-species dots and multi-species pies
  ctd_single <- NULL; ctd_pie_polys <- NULL
  if (!is.null(ctd_rows) && nrow(ctd_rows) > 0) {
    ctd_sp_count <- ctd_rows %>%
      group_by(FilterID) %>%
      summarise(n_sp_ctd = n_distinct(species), .groups = "drop")
    ctd_rows <- ctd_rows %>% left_join(ctd_sp_count, by = "FilterID")

    ctd_single <- ctd_rows %>%
      filter(n_sp_ctd == 1) %>%
      distinct(FilterID, .keep_all = TRUE)

    ctd_multi <- ctd_rows %>%
      filter(n_sp_ctd > 1) %>%
      mutate(prop = 1 / n_sp_ctd) %>%
      select(FilterID, lon_mid, lat_mid, species, color, prop)

    if (nrow(ctd_multi) > 0) {
      ctd_pie_polys <- make_ctd_pie_polys(ctd_multi)
      cat("  CTD pie charts:", n_distinct(ctd_multi$FilterID),
          "location(s) with", n_distinct(ctd_multi$species), "species\n")
    }
  }

  # ── Build plot ────────────────────────────────────────────────
  p <- ggplot() +
    geom_segment(
      data = grey_filters,
      aes(x = lon_start, y = lat_start, xend = lon_end, yend = lat_end),
      color = "grey65", linewidth = 0.7, alpha = 0.55,
      lineend = "round", linetype = "dashed"
    ) +
    geom_polygon(
      data = world_df, aes(x = long, y = lat, group = group),
      fill = "#c8d8b0", colour = "#7a9a5a", linewidth = 0.3
    )

  # Pass 1: buffered transect segments (non-CTD) ─────────────────
  for (x in sp_det_list) {
    det_line <- x$det %>% filter(Collection.method != "MOSAIC-D")
    col      <- x$sp$color
    if (nrow(det_line) == 0) next
    for (i in seq_len(n_layers)) {
      p <- p + geom_segment(
        data      = det_line,
        aes(x = lon_start, y = lat_start, xend = lon_end, yend = lat_end),
        color     = col, linewidth = layer_lw[i], alpha = layer_alpha[i],
        lineend   = "round"
      )
    }
  }

  # Pass 2a: single-species CTD — solid dot, no halo ─────────────
  if (!is.null(ctd_single) && nrow(ctd_single) > 0) {
    # Draw per color group so each dot gets its species color
    for (col_val in unique(ctd_single$color)) {
      sub <- ctd_single %>% filter(color == col_val)
      p <- p + geom_point(
        data  = sub,
        aes(x = lon_mid, y = lat_mid),
        color = col_val, size = 3.2, alpha = 1
      )
    }
  }

  # Pass 2b: multi-species CTD — pie chart polygons ──────────────
  if (!is.null(ctd_pie_polys) && nrow(ctd_pie_polys) > 0) {
    # Draw each wedge with its own fill color using individual calls
    for (col_val in unique(ctd_pie_polys$color)) {
      sub <- ctd_pie_polys %>% filter(color == col_val)
      p <- p + geom_polygon(
        data  = sub,
        aes(x = x, y = y, group = group_id),
        fill  = col_val, colour = "white", linewidth = 0.3, alpha = 0.95
      )
    }
  }

  # Pass 3: method labels (one layer across species so labels repel each other)
  det_labels <- bind_rows(lapply(sp_det_list, function(x) {
    det <- x$det
    if (nrow(det) == 0) return(NULL)
    det %>% mutate(
      disp_label = ifelse(Collection.method == "MOSAIC-D", "D", method_label),
      lab_col    = x$sp$color
    )
  }))
  if (nrow(det_labels) > 0) {
    p <- p + geom_label_repel(
      data          = det_labels,
      aes(x = lon_mid, y = lat_mid, label = disp_label, color = I(lab_col)),
      size = 3.0 + 1.2 * fsz, fontface = "bold", fill = "white",
      label.padding = unit(0.15, "lines"), label.size = 0.3, alpha = 0.95,
      box.padding = 0.7, point.padding = 0.5, min.segment.length = 0.2,
      segment.size = 0.3, segment.alpha = 0.7,
      force = 4, force_pull = 0.5, max.iter = 50000, max.time = 3,
      max.overlaps = Inf, seed = 42, show.legend = FALSE
    )
  }

  # Legend via off-map dummy segments ───────────────────────────
  legend_df <- do.call(rbind, lapply(panel$species, function(sp)
    data.frame(x = -999, xend = -999, y = -999, yend = -999,
               Species = sp$common, stringsAsFactors = FALSE)
  ))
  sp_col_named <- setNames(sapply(panel$species, `[[`, "color"),
                            sapply(panel$species, `[[`, "common"))
  sp_lab_named <- setNames(
    paste0("italic('", sapply(panel$species, `[[`, "label"), "')"),
    sapply(panel$species, `[[`, "common")
  )

  p <- p +
    geom_point(data = legend_df, aes(x = x, y = y, color = Species),
               shape = 16, size = 3.5 * legend_pt / 10) +
    scale_color_manual(name = "", values = sp_col_named, breaks = names(sp_col_named),
                       labels = lapply(sp_lab_named, function(l) parse(text = l))) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3.5 * legend_pt / 10))) +
    coord_fixed(ratio = 1, xlim = ca_lon, ylim = ca_lat, expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10 * fsz) +
    theme(
      legend.position      = c(0.99, 0.01),
      legend.justification = c(1, 0),
      legend.background    = element_rect(fill = alpha("white", 0.85), color = NA),
      legend.text          = element_text(size = legend_pt),
      legend.key           = element_blank(),
      legend.key.width     = unit(0.6 * legend_pt / 10, "cm"),
      legend.key.height    = unit(legend_key_h, "cm"),
      legend.spacing.y     = unit(0, "pt"),
      legend.margin        = margin(3, 6, 3, 6),
      panel.grid           = element_blank(),
      axis.text            = element_text(size = axis_pt),
      plot.margin          = margin(2, 4, 2, 4),
      panel.background     = element_rect(fill = "#d6e8f5")
    )

  # ── Title strip ───────────────────────────────────────────────
  sps     <- panel$species
  n_sps   <- length(sps)
  n_rows  <- ceiling(n_sps / 2)
  n_cols  <- min(n_sps, 2)
  row_idx <- ceiling(seq_along(sps) / 2)
  col_idx <- ((seq_along(sps) - 1) %% 2) + 1

  x_anchors <- if (n_cols == 1) 0.30 else c(0.03, 0.53)
  y_vals    <- if (n_rows == 1) 0.20  else c(0.62, 0.14)
  header_y  <- if (n_rows == 1) 0.76  else 0.90

  title_df <- data.frame(
    x      = x_anchors[col_idx],
    y      = y_vals[row_idx],
    common = sapply(sps, `[[`, "common"),
    sci    = sapply(sps, `[[`, "label"),
    color  = sapply(sps, `[[`, "color"),
    stringsAsFactors = FALSE
  )
  dot_x  <- title_df$x
  text_x <- title_df$x + 0.045

  title_p <- ggplot() +
    annotate("text", x = 0.5, y = header_y, label = panel$title_common,
             fontface = "bold", size = title_mm, hjust = 0.5, color = "grey10")
  if (!is.null(tag))   # panel letter on the same line as the panel title
    title_p <- title_p + annotate("text", x = 0.005, y = header_y, label = tag, fontface = "bold",
                                  size = tag_mm, hjust = 0, color = "grey10")
  for (si in seq_len(n_sps)) {
    title_p <- title_p +
      annotate("point", x = dot_x[si], y = title_df$y[si],
               color = title_df$color[si], size = 4.5 * fsz) +
      annotate("text",  x = text_x[si], y = title_df$y[si] + 0.12,
               label = title_df$common[si], fontface = "bold", size = 4.4 * fsz,
               hjust = 0, color = "grey10") +
      annotate("text",  x = text_x[si], y = title_df$y[si] - 0.10,
               label = title_df$sci[si], fontface = "italic", size = 4.4 * fsz,
               hjust = 0, color = "grey10")
  }
  title_p <- title_p + xlim(0, 1) + ylim(0, 1) +
    theme_void() + theme(plot.margin = margin(4, 0, 0, 0))

  strip_h <- if (n_rows == 1) 0.20 else 0.32
  cowplot::plot_grid(title_p, p, ncol = 1, rel_heights = c(strip_h, 1))
}

# ── 10. Loop over cruises ─────────────────────────────────────
cruises <- sort(unique(meta2$Cruise))
cruise_season <- c(RS2311 = "Fall 2023", SR2408 = "Summer 2024", RL2501 = "Winter 2025")
cat("\nCruises:", paste(cruises, collapse=", "), "\n")

# Accumulate caption text across cruises for separate results document
all_caption_lines <- list()

for (cruise in cruises) {
  cat("\n====", cruise, "====\n")
  all_filters_cr   <<- all_filters %>% filter(Cruise == cruise)
  all_filters_orig  <- all_filters
  cat("Filters in cruise:", nrow(all_filters_cr), "\n")

  # Collect captions for this cruise (computed inside make_panel but not shown)
  cruise_captions <- list()

  panels <- lapply(seq_along(panels_def), function(i) {
    panel <- panels_def[[i]]
    all_filters <<- all_filters_cr
    p <- make_panel(panel, tag = LETTERS[i])
    all_filters <<- all_filters_orig
    p
  })

  # Re-extract captions by re-running effort_caption outside make_panel
  for (panel in panels_def) {
    cap_lines <- sapply(panel$species, function(sp) {
      pct_df <- get_pct_df(sp$taxon_id)
      all_filters_tmp <- all_filters_cr
      effort_caption(sp$common, pct_df, all_filters_tmp)
    })
    cruise_captions[[panel$title_common]] <- cap_lines
  }
  all_caption_lines[[cruise]] <- cruise_captions

  combined <- (panels[[1]] | panels[[2]]) /
              (panels[[3]] | panels[[4]]) +
    plot_annotation(
      title    = paste0("Rare Taxa Detections – ", cruise_season[[cruise]]),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5)
      )
    )

  panel_file <- paste0("viz/rare_taxa_panels_", cruise, "-pie.png")
  ggsave(panel_file, combined, width = 14, height = 12, dpi = 300, bg = "white")
  cat("Panel map saved:", panel_file, "\n")

  for (i in seq_along(panels)) {
    pid      <- panels_def[[i]]$panel_id
    ind_file <- paste0("viz/rare_taxa_", cruise, "_", pid, "-pie.png")
    ggsave(ind_file, panels[[i]], width = 7, height = 6.5, dpi = 300, bg = "white")
  }
  cat("Individual maps saved for", cruise, "\n")

  # ── Journal-quality PDF ───────────────────────────────────────
  cat("  Rendering journal PDF...\n")
  jpanels <- lapply(seq_along(panels_def), function(i) {
    panel <- panels_def[[i]]
    all_filters <<- all_filters_cr
    p <- make_panel(panel, fsz = 0.58, tag = LETTERS[i], axis_pt = 6, legend_pt = 8, title_mm = 3.7, tag_mm = 4.6, legend_key_h = 0.34)
    all_filters <<- all_filters_orig
    p
  })

  jcombined <- (jpanels[[1]] | jpanels[[2]]) /
               (jpanels[[3]] | jpanels[[4]]) +
    plot_annotation(
      title    = paste0("Rare Taxa Detections – ", cruise_season[[cruise]]),
      theme = theme(
        text          = element_text(family = "serif"),
        plot.title    = element_text(face = "bold", size = 8, hjust = 0.5, family = "serif"),
        plot.subtitle = element_text(size = 5.5, color = "grey40", hjust = 0.5,
                                     family = "serif", lineheight = 1.3)
      )
    )

  pdf_file <- paste0("viz/rare_taxa_panels_", cruise, ".pdf")
  ggsave(pdf_file, jcombined, width = 7.2, height = 7.5, device = cairo_pdf, bg = "white")
  cat("Journal PDF saved:", pdf_file, "\n")
}

# ── 11. Write caption document ────────────────────────────────
cap_doc_lines <- c(
  "MOSAIC Survey — Figure Captions (Detection Rates by Method)",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d")),
  "",
  "For each species: detections / total cruise filters per MOSAIC method (% of effort).",
  "Methods: A = MOSAIC-A (underway continuous), B = MOSAIC-B (underway discrete),",
  "         C = MOSAIC-C (underway discrete, alternate), D = MOSAIC-D (CTD station)",
  ""
)

for (cruise in names(all_caption_lines)) {
  cap_doc_lines <- c(cap_doc_lines,
    paste0("══════════════════════════════════════"),
    paste0("Cruise: ", cruise),
    paste0("══════════════════════════════════════"),
    ""
  )
  for (panel_title in names(all_caption_lines[[cruise]])) {
    cap_doc_lines <- c(cap_doc_lines,
      paste0("  ── ", panel_title, " ──"),
      paste0("  ", all_caption_lines[[cruise]][[panel_title]]),
      ""
    )
  }
}

#cap_file <- "viz/MOSAIC_figure_captions.txt"
#writeLines(cap_doc_lines, cap_file)
#cat("\nCaption document saved:", cap_file, "\n")

cat("\nAll done!\n")

