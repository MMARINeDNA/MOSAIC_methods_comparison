library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

OUT_PATH <- "viz/mosaic_figure3_familylevel_metadata_filtered.pdf"
OUT_W    <- 7.5
OUT_H    <- 11

# ── Colours (Wong 2011) ───────────────────────────────────────────────────────
COLORS <- c(
  Actinopteri    = "#0072B2",
  Mammalia       = "#D55E00",
  Chondrichthyes = "#009E73",
  Aves           = "#E69F00",
  Unassigned     = "#BBBBBB"
)
COLORS_LIGHT <- c(
  Actinopteri    = "#b3d4ea",
  Mammalia       = "#f0c0a8",
  Chondrichthyes = "#a8ddd0",
  Aves           = "#f5dfa0"
)

METHOD_ORDER  <- c("MOSAIC-A", "MOSAIC-B", "MOSAIC-C", "MOSAIC-D")
METHOD_LABELS <- c("MOSAIC-A"="MOSAIC-A","MOSAIC-B"="MOSAIC-B",
                   "MOSAIC-C"="MOSAIC-C","MOSAIC-D"="CTD")
VERT_CLASSES  <- c("Actinopteri","Mammalia","Chondrichthyes","Aves")
GROUP_ORDER   <- c("Actinopteri","Mammalia","Chondrichthyes","Aves","Unassigned")

# ── Load and filter data from CSVs ────────────────────────────────────────────
OTU_PATH <- "outputs/MOSAIC_taxnamed_otu_table_metadata_filtered.csv"
TAX_PATH <- "outputs/MOSAIC_taxnamed_taxa_table_metadata_filtered.csv"
SAM_PATH <- "outputs/MOSAIC_taxnamed_sample_data_metadata_filtered.csv"

otu <- read.csv(OTU_PATH, row.names = 1, check.names = FALSE)
tax <- read.csv(TAX_PATH, row.names = 1, check.names = FALSE)
sam <- read.csv(SAM_PATH, row.names = 1, check.names = FALSE)

# Restore hyphens in column names (R converts them to dots on read)
colnames(otu) <- gsub("\\.", "-", colnames(otu))

# Align samples
shared <- intersect(colnames(otu), rownames(sam))
otu <- otu[, shared]; sam <- sam[shared, ]

# Keep only the four focal methods
sam <- sam[sam$Collection.method %in% METHOD_ORDER & !is.na(sam$Collection.method), ]
otu <- otu[, rownames(sam)]

# Family-level filter: exclude taxa whose family field is a higher-rank placeholder
fam_ok <- !grepl("\\border\\b|\\bclass\\b|\\bphylum\\b",
                 tax$family, ignore.case = TRUE, perl = TRUE)
message("Excluded (below family): ", sum(!fam_ok), " — ",
        paste(rownames(tax)[!fam_ok], collapse = ", "))
otu <- otu[fam_ok, ]; tax <- tax[fam_ok, ]

# Helper: logical mask for a taxonomic group
grp_mask <- function(g) {
  if (g == "Unassigned") !(tax$class %in% VERT_CLASSES)
  else !is.na(tax$class) & tax$class == g
}

method_v <- sam$Collection.method

# ── Panel A/B: reads and taxa richness ────────────────────────────────────────
read_rows <- taxa_rows <- list()
for (m in METHOD_ORDER) {
  sids <- rownames(sam)[method_v == m]
  sub  <- otu[, sids, drop = FALSE]
  for (g in GROUP_ORDER) {
    mask <- grp_mask(g)
    read_rows[[length(read_rows)+1]] <- data.frame(method=m, group=g,
      reads=sum(sub[mask, ]), stringsAsFactors=FALSE)
    taxa_rows[[length(taxa_rows)+1]] <- data.frame(method=m, group=g,
      n_taxa=sum(rowSums(sub[mask, , drop=FALSE]) > 0), stringsAsFactors=FALSE)
  }
}

df_reads <- do.call(rbind, read_rows) %>%
  group_by(method) %>% mutate(pct = reads / sum(reads) * 100) %>% ungroup() %>%
  mutate(group  = factor(group,  levels = rev(GROUP_ORDER)),
         method = factor(METHOD_LABELS[method], levels = rev(METHOD_LABELS)))

df_taxa <- do.call(rbind, taxa_rows) %>%
  group_by(method) %>% mutate(pct = n_taxa / sum(n_taxa) * 100) %>% ungroup() %>%
  mutate(group  = factor(group,  levels = rev(GROUP_ORDER)),
         method = factor(METHOD_LABELS[method], levels = rev(METHOD_LABELS)))

# ── Panel C: exclusive vertebrate taxa ────────────────────────────────────────
unique_rows <- list()
for (m in METHOD_ORDER) {
  m_sids <- rownames(sam)[method_v == m]
  o_sids <- rownames(sam)[method_v != m & method_v %in% METHOD_ORDER]
  m_pres <- rowSums(otu[, m_sids, drop=FALSE]) > 0
  o_pres <- if (length(o_sids) > 0) rowSums(otu[, o_sids, drop=FALSE]) > 0
            else rep(FALSE, nrow(otu))
  for (g in VERT_CLASSES) {
    mask <- grp_mask(g)
    unique_rows[[length(unique_rows)+1]] <- data.frame(method=m, group=g,
      n_unique=sum(mask & m_pres & !o_pres), stringsAsFactors=FALSE)
  }
}

df_unique <- do.call(rbind, unique_rows) %>%
  mutate(group  = factor(group,  levels = rev(VERT_CLASSES)),
         method = factor(METHOD_LABELS[method], levels = rev(METHOD_LABELS)))

# ── Panel D: per-PCR-replicate detection rate ─────────────────────────────────
detect_rows <- list()
for (m in METHOD_ORDER) {
  sids <- rownames(sam)[method_v == m]
  n    <- length(sids)
  for (g in VERT_CLASSES) {
    mask <- grp_mask(g)
    det  <- sum(colSums(otu[mask, sids, drop=FALSE]) > 0)
    detect_rows[[length(detect_rows)+1]] <- data.frame(method=m, group=g,
      n_total=n, detected=det, not_detected=n-det, stringsAsFactors=FALSE)
  }
}

df_detect <- do.call(rbind, detect_rows) %>%
  pivot_longer(c(detected, not_detected), names_to="status", values_to="count") %>%
  group_by(method, group) %>% mutate(pct = count / n_total * 100) %>% ungroup() %>%
  mutate(method   = factor(METHOD_LABELS[method], levels = rev(METHOD_LABELS)),
         group    = factor(group, levels = VERT_CLASSES),
         status   = factor(status, levels = c("not_detected","detected")),
         fill_key = paste0(status,"_",group))

det_colors <- c(
  setNames(COLORS[VERT_CLASSES], paste0("detected_",    VERT_CLASSES)),
  setNames(COLORS_LIGHT,         paste0("not_detected_", VERT_CLASSES))
)

# ── Theme ─────────────────────────────────────────────────────────────────────
th <- function(bs = 9) {
  theme_classic(base_family = "serif", base_size = bs) +
  theme(
    axis.line          = element_line(linewidth = 0.4),
    axis.ticks         = element_line(linewidth = 0.4),
    axis.ticks.length  = unit(2.5, "pt"),
    panel.grid.major.x = element_line(colour = "#eeeeee", linewidth = 0.3),
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = bs, hjust = 0,
                                      margin = margin(b = 4)),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "italic", size = bs - 0.5),
    plot.margin        = margin(4, 8, 4, 4)
  )
}

# ── End-of-bar total labels ───────────────────────────────────────────────────
reads_totals <- df_reads %>%
  group_by(method) %>%
  summarise(total_reads = sum(reads), .groups = "drop") %>%
  mutate(label = ifelse(total_reads >= 1e6,
                        paste0("n = ", sprintf("%.1f", total_reads/1e6), "M reads"),
                        paste0("n = ", formatC(total_reads, format="d", big.mark=","), " reads")))

taxa_totals <- df_taxa %>%
  group_by(method) %>%
  summarise(total_taxa = sum(n_taxa), .groups = "drop") %>%
  mutate(label = paste0("n = ", total_taxa, " taxa"))

# ── Panels ────────────────────────────────────────────────────────────────────
pA <- ggplot(df_reads,  aes(pct, method, fill = group)) +
  geom_col(width = 0.65) +
  geom_text(data = reads_totals, aes(x = 101, y = method, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.8, family = "serif",
            colour = "#666666", fontface = "italic") +
  scale_fill_manual(values = COLORS, drop = FALSE) +
  scale_x_continuous(labels = \(x) paste0(x, "%"), expand = c(0,0),
                     limits = c(0, 145)) +
  labs(title = "(A)  Proportion of sequencing reads", x = NULL, y = NULL) + th()

pB <- ggplot(df_taxa,   aes(pct, method, fill = group)) +
  geom_col(width = 0.65) +
  geom_text(data = taxa_totals, aes(x = 101, y = method, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.8, family = "serif",
            colour = "#666666", fontface = "italic") +
  scale_fill_manual(values = COLORS, drop = FALSE) +
  scale_x_continuous(labels = \(x) paste0(x, "%"), expand = c(0,0),
                     limits = c(0, 145)) +
  labs(title = "(B)  Proportion of taxa richness", x = NULL, y = NULL) + th()

unique_totals <- df_unique %>%
  group_by(method) %>%
  summarise(total = sum(n_unique), .groups = "drop")
max_unique <- max(unique_totals$total)

pC <- ggplot(df_unique, aes(n_unique, method, fill = group)) +
  geom_col(width = 0.65) +
  geom_text(data = unique_totals,
            aes(x = max_unique * 1.04, y = method, label = total),
            inherit.aes = FALSE, hjust = 0, size = 2.8, family = "serif",
            colour = "#666666", fontface = "italic") +
  scale_fill_manual(values = COLORS, drop = FALSE) +
  scale_x_continuous(expand = c(0,0), limits = c(0, max_unique * 1.25),
                     breaks = seq(0, 80, 20)) +
  labs(title = "(C)  Exclusive vertebrate detections",
       x = "Number of vertebrate taxa", y = NULL) + th()

detect_labels <- df_detect %>%
  filter(status == "detected") %>%
  select(method, group, count, n_total) %>%
  mutate(label = paste0(count, " / ", n_total))

pD <- ggplot(df_detect, aes(pct, method, fill = fill_key)) +
  geom_col(width = 0.65) +
  geom_text(data = detect_labels,
            aes(x = 101, y = method, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.5, family = "serif",
            colour = "#666666") +
  facet_wrap(~ group, ncol = 1, strip.position = "left") +
  scale_fill_manual(values = det_colors, drop = FALSE) +
  scale_x_continuous(labels = \(x) paste0(x, "%"), expand = c(0,0), limits = c(0, 145)) +
  labs(title = "(D)  Proportion of PCR replicates detecting each vertebrate class",
       x = "Proportion of filters (%)", y = NULL) +
  th() + theme(strip.text      = element_text(face = "italic", size = 8),
               strip.placement = "outside",
               axis.text.y     = element_text(size = 8),
               panel.spacing   = unit(3, "pt"))

# ── Legends ───────────────────────────────────────────────────────────────────
leg1 <- get_legend(
  ggplot(data.frame(group = factor(GROUP_ORDER, levels = GROUP_ORDER), x=1, y=1),
         aes(x, y, fill = group)) + geom_col() +
  scale_fill_manual(values = COLORS,
    labels = c(expression(italic("Actinopteri")), expression(italic("Mammalia")),
               expression(italic("Chondrichthyes")), expression(italic("Aves")),
               "Unassigned"),
    guide = guide_legend(nrow = 1, title = NULL)) +
  theme_void() +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 8.5, family = "serif"),
        legend.key.size  = unit(10, "pt"),
        legend.spacing.x = unit(4, "pt"))
)

leg2 <- get_legend(
  ggplot(data.frame(
    status = factor(c("Detected (>=1 read)","Not detected"),
                    levels = c("Detected (>=1 read)","Not detected")),
    fill = c(COLORS["Mammalia"], COLORS_LIGHT["Mammalia"]), x=1, y=1),
    aes(x, y, fill = status)) + geom_col() +
  scale_fill_manual(
    values = setNames(c(COLORS["Mammalia"], COLORS_LIGHT["Mammalia"]),
                      c("Detected (>=1 read)", "Not detected")),
    guide = guide_legend(nrow = 1, title = NULL)) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 8, family = "serif"),
        legend.key.size = unit(9, "pt"))
)

# ── Caption ───────────────────────────────────────────────────────────────────
n_per_method <- table(method_v)
n_excl       <- sum(!fam_ok)
n_total_taxa <- nrow(tax) + n_excl
cap <- paste0(
  "Cruises RS2311 (Fall 2023, November) and SR2408 (Summer 2024, July-August).  ",
  "Sample sizes: ",
  paste(paste0(METHOD_LABELS[names(n_per_method)], " n = ", as.integer(n_per_method)),
        collapse = "; "), ".  ",
  "MOSAIC-B and MOSAIC-C: SR2408 only.  ",
  "Taxonomic filter: ", nrow(tax), " of ", n_total_taxa,
  " taxa resolved to >=family level retained; ",
  n_excl, " taxa identified only to order or above excluded.  ",
  "Panel D: a filter is counted as detected if it yielded >=1 read for that class."
)
cap_grob <- ggdraw() +
  draw_label(cap, x = 0.02, y = 0.5, hjust = 0, vjust = 0.5,
             fontfamily = "serif", size = 6.5, colour = "#555555", lineheight = 1.4)

# ── Assemble ──────────────────────────────────────────────────────────────────
top3 <- plot_grid(pA, pB, pC, ncol = 1, rel_heights = c(1,1,1), align = "v", axis = "lr")
fig  <- plot_grid(top3, pD, ncol = 1, rel_heights = c(3,4), align = "v", axis = "lr")

final <- plot_grid(fig, leg1, leg2, cap_grob,
                   ncol = 1, rel_heights = c(1, 0.04, 0.035, 0.1))

save_plot(OUT_PATH, final, base_width = OUT_W, base_height = OUT_H,
          dpi = 300, device = cairo_pdf)
message("Saved: ", OUT_PATH)

