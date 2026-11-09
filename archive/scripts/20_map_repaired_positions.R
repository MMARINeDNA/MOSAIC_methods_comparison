# Map of the sample positions written by script 1.5 (ship-track repaired metadata): ship track in grey, each filter's
# start-to-end segment coloured by method, repaired transects (start or end from the ship track) drawn thick,
# CTD casts (on their CalCOFI stations) as points. Diagnostic only.
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(data.table); library(maps)})
qc <- read.csv("outputs/geodata_QC_report.csv")
trk <- list(RS2311 = "ship_tracks/SR2311_combined_gps_coordinates.csv", SR2408 = "ship_tracks/SR2408_combined_gps_coordinates.csv")
ttl <- c(RS2311 = "Fall 2023 (RS2311)", SR2408 = "Summer 2024 (SR2408)")
cols <- c("MOSAIC-A" = "#0072B2", "MOSAIC-B" = "#E69F00", "MOSAIC-C" = "#009E73", "MOSAIC-D (CTD)" = "#D55E00")
qc$method <- ifelse(qc$Collection.method == "MOSAIC-D", "MOSAIC-D (CTD)", qc$Collection.method)
qc$repaired <- qc$src_start == "ship track" | qc$src_end == "ship track"
world <- map_data("world")
mk <- function(cr) {
  g <- fread(trk[[cr]], select = c("latitude", "longitude"))[seq(1, .N, by = 120)]
  d <- filter(qc, Cruise == cr, !is.na(Latitude.start))
  ggplot() + geom_polygon(data = world, aes(long, lat, group = group), fill = "#e3ead6", colour = "#9aab85", linewidth = 0.2) +
    geom_path(data = g, aes(longitude, latitude), colour = "grey70", linewidth = 0.3) +
    geom_segment(data = filter(d, method != "MOSAIC-D (CTD)"), aes(x = Longitude.start, y = Latitude.start, xend = Longitude.end, yend = Latitude.end,
                 colour = method, linewidth = repaired)) +
    geom_point(data = filter(d, method == "MOSAIC-D (CTD)"), aes(Longitude.start, Latitude.start, colour = method), size = 1.8) +
    scale_colour_manual(values = cols, name = NULL) + scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.6), name = "Repaired from\nship track") +
    coord_quickmap(xlim = c(-125, -116.8), ylim = c(29.6, 35.2)) + labs(title = ttl[[cr]], x = NULL, y = NULL) + theme_bw(base_size = 11)
}
library(patchwork)
p <- (mk("RS2311") | mk("SR2408")) + plot_layout(guides = "collect") +
  plot_annotation(title = "Sample positions after repair from the ship GPS track (script 1.5; no workbook)",
                  subtitle = "Grey = ship track; thick = start and/or end replaced from the track; points = CTD casts on their labelled CalCOFI station")
ggsave("outputs/repaired_positions_map.png", p, width = 14, height = 7, dpi = 170)
