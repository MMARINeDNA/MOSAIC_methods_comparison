# For each timed MOSAIC-A/B filter, compare original and workbook-corrected start/end positions with the ship GPS
# track at the log's start/end clock times under two local->UTC offsets (PDT = +7 h, PST = +8 h). Diagnostic only.
suppressPackageStartupMessages({library(phyloseq); library(data.table); library(geodist)})
tracks <- list(RS2311 = "ship_tracks/SR2311_combined_gps_coordinates.csv", SR2408 = "ship_tracks/SR2408_combined_gps_coordinates.csv")
tracks <- lapply(tracks, function(f) { g <- fread(file.path("/Users/Eldridge/ClaudeCode/research/MOSAIC1", f), select = 1:4)
  g[, t := as.numeric(as.POSIXct(timestamp, tz = "UTC"))]; setkey(g, t); g })
orig <- data.frame(sample_data(readRDS("outputs/MOSAIC.ps.RDS")), check.names = FALSE)
corr <- read.csv("outputs/geodata_corrected_MOSAIC_metadata.csv", check.names = FALSE, row.names = 1)
f <- unique(corr[corr$Collection.method %in% c("MOSAIC-A", "MOSAIC-B") & grepl(":", corr$Filtration.time.start) & grepl(":", corr$Filtration.time.end),
                 c("Cruise", "FilterID", "Date..YYYYMMDD.", "Filtration.time.start", "Filtration.time.end", "Latitude.start", "Longitude.start", "Latitude.end", "Longitude.end")])
o <- unique(orig[, c("FilterID", "Latitude.start", "Longitude.start", "Latitude.end", "Longitude.end")]); names(o)[-1] <- paste0("o_", names(o)[-1])
f <- merge(f, o, by = "FilterID")
at <- function(g, utc) { i <- g[, findInterval(utc, t)]; i <- pmax(i, 1); g[i, .(latitude, longitude)] }
km <- function(la1, lo1, la2, lo2) mapply(function(a, b, c, d) geodist(c(lat = a, lon = b), c(lat = c, lon = d), measure = "geodesic") / 1000, la1, lo1, la2, lo2)
res <- do.call(rbind, lapply(seq_len(nrow(f)), function(k) {
  r <- f[k, ]; g <- tracks[[r$Cruise]]; if (is.null(g)) return(NULL)
  d <- as.Date(as.character(r$Date..YYYYMMDD.), "%Y%m%d")
  ut <- function(hm, off) as.numeric(as.POSIXct(paste(d, hm), tz = "UTC")) + off * 3600
  out <- data.frame(Cruise = r$Cruise, FilterID = r$FilterID, date = d)
  for (off in c(7, 8)) {
    s <- at(g, ut(r$Filtration.time.start, off)); e <- at(g, ut(r$Filtration.time.end, off))
    out[[paste0("corr_vs_track_km_off", off)]] <- max(km(r$Latitude.start, r$Longitude.start, s$latitude, s$longitude), km(r$Latitude.end, r$Longitude.end, e$latitude, e$longitude))
    out[[paste0("orig_vs_track_km_off", off)]] <- max(km(r$o_Latitude.start, r$o_Longitude.start, s$latitude, s$longitude), km(r$o_Latitude.end, r$o_Longitude.end, e$latitude, e$longitude))
  }
  out }))
write.csv(res, "outputs/workbook_vs_ship_track_timezone_check.csv", row.names = FALSE)
res$workbook_offset <- ifelse(res$corr_vs_track_km_off7 < 0.5, "+7 (PDT)", ifelse(res$corr_vs_track_km_off8 < 0.5, "+8 (PST)", "neither"))
print(table(res$date, res$workbook_offset))
