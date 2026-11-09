# Fall34 (MOSAIC-A, Fall 2023): original field coordinates vs workbook-corrected coordinates, with neighbouring
# Fall transects and CalCOFI Line 86.7 stations for context. Diagnostic only.
suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(ggplot2); library(ggrepel); library(maps)})
orig <- data.frame(sample_data(readRDS("outputs/MOSAIC.ps.RDS")), check.names = FALSE)
corr <- read.csv("outputs/geodata_corrected_MOSAIC_metadata.csv", check.names = FALSE, row.names = 1)
pick <- function(d, src) d %>% filter(Cruise == "RS2311", Collection.method == "MOSAIC-A") %>%
  distinct(FilterID, Latitude.start, Longitude.start, Latitude.end, Longitude.end) %>%
  transmute(FilterID, src, lat1 = as.numeric(Latitude.start), lon1 = as.numeric(Longitude.start),
            lat2 = as.numeric(Latitude.end), lon2 = as.numeric(Longitude.end))
tr <- bind_rows(pick(orig, "Original"), pick(corr, "Workbook-corrected"))
f34 <- filter(tr, FilterID == "Fall34")
oth <- filter(tr, FilterID != "Fall34", src == "Workbook-corrected",
              pmax(lat1, lat2) < 34.2, pmin(lat1, lat2) > 33.4, pmax(lon1, lon2) < -118, pmin(lon1, lon2) > -119)
st <- read.csv("data/reference/CalCOFIStationOrder.csv", check.names = FALSE) %>%
  filter(Line == 86.7 | (Line > 86 & Line < 88), `Lat (dec)` > 33.4, `Lat (dec)` < 34.2, `Lon (dec)` > -119, `Lon (dec)` < -118) %>%
  transmute(label = paste(Line, Sta), lat = `Lat (dec)`, lon = `Lon (dec)`)
cols <- c("Original" = "#D55E00", "Workbook-corrected" = "#0072B2")
p <- ggplot() +
  geom_polygon(data = map_data("world"), aes(long, lat, group = group), fill = "#c8d8b0", colour = "#7a9a5a") +
  geom_segment(data = oth, aes(x = lon1, y = lat1, xend = lon2, yend = lat2), colour = "grey55", linewidth = 0.6) +
  geom_point(data = st, aes(lon, lat), shape = 4, colour = "grey30") +
  geom_text_repel(data = st, aes(lon, lat, label = label), size = 3, colour = "grey30", max.overlaps = 20) +
  geom_segment(data = f34, aes(x = lon1, y = lat1, xend = lon2, yend = lat2, colour = src), linewidth = 1.4,
               arrow = arrow(length = unit(0.18, "cm"))) +
  geom_point(data = f34, aes(lon1, lat1, colour = src), size = 3) +
  scale_colour_manual(values = cols, name = "Fall34") +
  coord_quickmap(xlim = c(-118.9, -118.0), ylim = c(33.55, 34.1)) +
  labs(title = "Fall34 (MOSAIC-A, 15:34-16:14): original vs workbook-corrected",
       subtitle = "Dot = start, arrow = end; grey = other corrected Fall MOSAIC-A transects; x = CalCOFI stations",
       x = NULL, y = NULL) + theme_bw()
ggsave("outputs/Fall34_original_vs_corrected.png", p, width = 9, height = 6.5, dpi = 200)
print(f34)
