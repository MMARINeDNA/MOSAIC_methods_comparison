# MOSAIC sampling effort summary table
# Rows: Cruise x Method (CTD split by depth category)
# Cols: n filters, Volume (L), Distance (km), Depth
# Outputs: mosaic_tbl data frame, flextable (.docx + .html), kableExtra (.html)

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(kableExtra)
  library(officer)
})

# ── Load data ─────────────────────────────────────────────────
meta <- read.csv(
  "outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv",
  row.names = 1
) %>%
  filter(!is.na(FilterID)) %>%
  mutate(sample_id = rownames(.))

taxa <- read.csv(
  "outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv",
  row.names = 1
)
otu <- read.csv(
  "outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv",
  row.names = 1
)
colnames(otu) <- gsub("[.]", "-", colnames(otu))

# ── Target taxa ────────────────────────────────────────────────
cet_fams   <- c("Delphinidae","Balaenopteridae","Phocoenidae",
                "Ziphiidae","Physeteridae","Artiodactyla order")
cet_taxa   <- rownames(taxa)[taxa[,"family"] %in% cet_fams]
shark_taxa <- rownames(taxa)[taxa[,"class"] == "Chondrichthyes"]

# ── Filter-level detections ────────────────────────────────────
# For each sample (row), check if any target taxon has reads
det_vec <- function(sids, tgt) {
  sids <- sids[sids %in% colnames(otu)]
  tgt  <- tgt[tgt %in% rownames(otu)]
  if (!length(sids) || !length(tgt)) return(setNames(rep(FALSE, length(sids)), sids))
  setNames(colSums(otu[tgt, sids, drop = FALSE]) > 0, sids)
}

meta <- meta %>%
  mutate(
    det_cet   = det_vec(sample_id, cet_taxa)[sample_id],
    det_shark = det_vec(sample_id, shark_taxa)[sample_id]
  )

# ── Helper: filter-level detection rates ──────────────────────
# Collapse reps to filter, then compute % filters detecting
det_rates <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(c(group_vars, "FilterID")))) %>%
    summarise(
      filt_cet   = any(det_cet),
      filt_shark = any(det_shark),
      .groups    = "drop"
    ) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      pct_cet   = round(mean(filt_cet)   * 100, 1),
      pct_shark = round(mean(filt_shark) * 100, 1),
      .groups   = "drop"
    )
}

# ── Classify CTD depth category ───────────────────────────────
# RS2311: <250 / >250 already in metadata
# SR2408: all CTD samples are mesopelagic (>250 m)
meta <- meta %>%
  mutate(ctd_depth_cat = case_when(
    Collection.method != "MOSAIC-D" ~ NA_character_,
    Cruise == "SR2408"              ~ "CTD >250 m",
    Depth..m. == "<250"            ~ "CTD <250 m",
    TRUE                           ~ "CTD >250 m"
  ))

# ── Helper: summarise to filter level ─────────────────────────
summarise_effort <- function(df, group_vars) {
  df %>%
    distinct(across(all_of(c(group_vars, "FilterID",
                              "dist_km", "Volume.filtered..mL.")))) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_filters = n(),
      vol_mean  = round(mean(Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
      vol_min   = round(min( Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
      vol_max   = round(max( Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
      dist_mean = round(mean(dist_km, na.rm = TRUE), 1),
      dist_min  = round(min( dist_km, na.rm = TRUE), 1),
      dist_max  = round(max( dist_km, na.rm = TRUE), 1),
      .groups   = "drop"
    )
}

# ── Per-cruise stats ──────────────────────────────────────────
uw_stats <- meta %>%
  filter(Collection.method != "MOSAIC-D") %>%
  summarise_effort(c("Cruise", "Collection.method")) %>%
  left_join(
    meta %>% filter(Collection.method != "MOSAIC-D") %>%
      det_rates(c("Cruise", "Collection.method")),
    by = c("Cruise", "Collection.method")
  ) %>%
  rename(Method = Collection.method)

ctd_stats <- meta %>%
  filter(Collection.method == "MOSAIC-D") %>%
  summarise_effort(c("Cruise", "ctd_depth_cat")) %>%
  left_join(
    meta %>% filter(Collection.method == "MOSAIC-D") %>%
      det_rates(c("Cruise", "ctd_depth_cat")),
    by = c("Cruise", "ctd_depth_cat")
  ) %>%
  rename(Method = ctd_depth_cat)

# ── All-cruises summary rows ───────────────────────────────────
uw_det_all <- meta %>%
  filter(Collection.method != "MOSAIC-D") %>%
  group_by(FilterID) %>%
  summarise(filt_cet=any(det_cet), filt_shark=any(det_shark), .groups="drop") %>%
  summarise(pct_cet=round(mean(filt_cet)*100,1), pct_shark=round(mean(filt_shark)*100,1))

uw_all <- meta %>%
  filter(Collection.method != "MOSAIC-D") %>%
  distinct(FilterID, dist_km, Volume.filtered..mL.) %>%
  summarise(
    n_filters = n(),
    vol_mean  = round(mean(Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
    vol_min   = round(min( Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
    vol_max   = round(max( Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
    dist_mean = round(mean(dist_km, na.rm = TRUE), 1),
    dist_min  = round(min( dist_km, na.rm = TRUE), 1),
    dist_max  = round(max( dist_km, na.rm = TRUE), 1)
  ) %>%
  bind_cols(uw_det_all) %>%
  mutate(Cruise = "All cruises", Method = "Underway")

ctd_det_all <- meta %>%
  filter(Collection.method == "MOSAIC-D") %>%
  group_by(FilterID) %>%
  summarise(filt_cet=any(det_cet), filt_shark=any(det_shark), .groups="drop") %>%
  summarise(pct_cet=round(mean(filt_cet)*100,1), pct_shark=round(mean(filt_shark)*100,1))

ctd_all <- meta %>%
  filter(Collection.method == "MOSAIC-D") %>%
  distinct(FilterID, dist_km, Volume.filtered..mL.) %>%
  summarise(
    n_filters = n(),
    vol_mean  = round(mean(Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
    vol_min   = round(min( Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
    vol_max   = round(max( Volume.filtered..mL., na.rm = TRUE) / 1000, 1),
    dist_mean = round(mean(dist_km, na.rm = TRUE), 1),
    dist_min  = round(min( dist_km, na.rm = TRUE), 1),
    dist_max  = round(max( dist_km, na.rm = TRUE), 1)
  ) %>%
  bind_cols(ctd_det_all) %>%
  mutate(Cruise = "All cruises", Method = "CTD")

# ── Depth strings from metadata ────────────────────────────────
depth_per_group <- meta %>%
  filter(Collection.method == "MOSAIC-D") %>%
  distinct(Cruise, ctd_depth_cat, FilterID, Depth..m.) %>%
  group_by(Cruise, ctd_depth_cat) %>%
  count(Depth..m.) %>%
  mutate(entry = paste0(Depth..m., " (n=", n, ")")) %>%
  summarise(Depth = paste(entry, collapse = "; "), .groups = "drop") %>%
  rename(Method = ctd_depth_cat)

depth_all_ctd <- meta %>%
  filter(Collection.method == "MOSAIC-D") %>%
  distinct(FilterID, Depth..m.) %>%
  count(Depth..m.) %>%
  mutate(entry = paste0(Depth..m., " (n=", n, ")")) %>%
  summarise(Depth = paste(entry, collapse = "; ")) %>%
  pull(Depth)

# ── Format helper ─────────────────────────────────────────────
fmt_range <- function(mean, lo, hi, na_if_zero = FALSE) {
  ifelse(
    na_if_zero & mean == 0, "NA",
    paste0(mean, " (", lo, "-", hi, ")")
  )
}

# ── Assemble mosaic_tbl ───────────────────────────────────────
method_order <- c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C", "CTD <250 m", "CTD >250 m")

per_cruise <- bind_rows(uw_stats, ctd_stats) %>%
  left_join(depth_per_group, by = c("Cruise", "Method")) %>%
  mutate(
    Depth = case_when(
      grepl("MOSAIC", Method) ~ "Surface (~5 m)",
      TRUE                    ~ Depth
    ),
    Method = factor(Method, levels = method_order)
  ) %>%
  arrange(Cruise, Method) %>%
  mutate(Method = as.character(Method))

all_rows <- bind_rows(uw_all, ctd_all) %>%
  mutate(
    Depth = ifelse(Method == "Underway", "Surface (~5 m)", depth_all_ctd)
  )

mosaic_tbl <- bind_rows(per_cruise, all_rows) %>%
  transmute(
    Cruise          = Cruise,
    Method          = Method,
    `n filters`     = n_filters,
    `Volume (L)`    = fmt_range(vol_mean, vol_min, vol_max),
    `Distance (km)` = fmt_range(dist_mean, dist_min, dist_max, na_if_zero = TRUE),
    Depth           = Depth,
    `% filters detecting cetaceans` = pct_cet,
    `% filters detecting sharks`    = pct_shark
  )

cat("── mosaic_tbl ──────────────────────────────────────────────\n")
print(mosaic_tbl, n = 20)

# ── Flextable ─────────────────────────────────────────────────
grey_vline <- fp_border(width = 0.5, color = "#cccccc")
thick      <- fp_border(width = 1.5, color = "black")
thin       <- fp_border(width = 0.5, color = "black")
mid        <- fp_border(width = 1.0, color = "black")

cruise_breaks <- mosaic_tbl %>%
  mutate(row = row_number()) %>%
  group_by(Cruise) %>%
  summarise(last_row = max(row), .groups = "drop") %>%
  filter(last_row < nrow(mosaic_tbl)) %>%
  pull(last_row)

mosaic_flextable <- mosaic_tbl %>%
  flextable() %>%
  merge_v(j = "Cruise") %>%
  valign(j = "Cruise", valign = "top") %>%
  bold(part = "header") %>%
  # text cols left in both header and body; numeric cols centred
  align(align = "center", part = "all") %>%
  align(j = c("Cruise", "Method", "Depth"), align = "left", part = "all") %>%
  # column widths
  width(j = "Cruise",         width = 0.85) %>%
  width(j = "Method",         width = 0.90) %>%
  width(j = "n filters",      width = 0.55) %>%
  width(j = "Volume (L)",     width = 0.95) %>%
  width(j = "Distance (km)",  width = 1.00) %>%
  width(j = "Depth",          width = 1.90) %>%
  width(j = "% filters detecting cetaceans", width = 1.2) %>%
  width(j = "% filters detecting sharks",    width = 1.1) %>%
  # borders: remove all, then add back selectively
  border_remove() %>%
  hline_top(border = thick, part = "header") %>%
  hline_bottom(border = thick, part = "header") %>%
  hline_bottom(border = thick, part = "body") %>%
  hline(i = cruise_breaks, border = mid, part = "body") %>%
  # grey vertical dividers between columns (inner only)
  vline(border = grey_vline, part = "all") %>%
  # font
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  set_table_properties(layout = "autofit")

# ── kableExtra ────────────────────────────────────────────────
cruise_row_counts <- mosaic_tbl %>%
  mutate(Cruise = factor(Cruise, levels = unique(Cruise))) %>%
  count(Cruise, .drop = FALSE) %>%
  { setNames(.$n, as.character(.$Cruise)) }

mosaic_kable <- mosaic_tbl %>%
  select(-Cruise) %>%
  kbl(
    booktabs = TRUE,
    linesep  = "",
    align    = c("l", "r", "r", "r", "l", "r", "r")
  ) %>%
  kable_styling(
    latex_options     = c("hold_position", "scale_down"),
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width        = FALSE,
    font_size         = 10
  ) %>%
  pack_rows(index = cruise_row_counts, bold = TRUE, italic = FALSE) %>%
  footnote(
    general           = "Volume and distance reported as mean (min-max). Distance not applicable (NA) for CTD (stationary casts). CTD depth from metadata.",
    general_title     = "Note: ",
    footnote_as_chunk = TRUE,
    threeparttable    = TRUE
  )

# ── Save ──────────────────────────────────────────────────────
save_as_docx(mosaic_flextable,
             path = "outputs/summary_table1.docx")
save_as_html(mosaic_flextable,
             path = "outputs/summary_table1.html")
writeLines(as.character(mosaic_kable),
           "outputs/summary_table1_kable.html")
