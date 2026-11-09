# Evidence on CTD station labels: for each CTD cast, (a) the CalCOFI stations nearest the end of the previous and the start of
# the next MOSAIC-A transect (positions from script 1.5), and (b) every place the ship was stationary (< 0.5 kt for >= 8 min)
# from 5 h before to 1 h after the filtration start, with the nearest CalCOFI station. Diagnostic only; changes nothing.
suppressPackageStartupMessages({library(dplyr); library(data.table)})
options(width = 250)
md  <- read.csv("outputs/geodata_corrected_MOSAIC_metadata.csv", row.names = 1, check.names = FALSE)
qc  <- read.csv("outputs/geodata_QC_report.csv")
stn <- read.csv("data/reference/CalCOFIStationOrder.csv", check.names = FALSE)
stn <- data.frame(lab = paste(as.numeric(stn$Line), as.numeric(stn$Sta)), lat = stn$`Lat (dec)`, lon = stn$`Lon (dec)`)
hav <- function(a, b, c, e) { p <- pi/180; 2*6371*asin(pmin(1, sqrt(sin((c-a)*p/2)^2 + cos(a*p)*cos(c*p)*sin((e-b)*p/2)^2))) }
near <- function(la, lo) { k <- hav(la, lo, stn$lat, stn$lon); i <- which.min(k); c(lab = stn$lab[i], km = round(k[i], 1)) }
trk <- list(RS2311 = "ship_tracks/SR2311_combined_gps_coordinates.csv", SR2408 = "ship_tracks/SR2408_combined_gps_coordinates.csv")
G <- lapply(trk, function(f) { g <- fread(f, select = c("timestamp", "latitude", "longitude"))[seq(1, .N, by = 120)]   # ~ 1 fix per minute
  g[, t := as.numeric(as.POSIXct(substr(timestamp, 1, 19), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))]
  g[, sp := c(NA, hav(latitude[-.N], longitude[-.N], latitude[-1], longitude[-1]) / (diff(t) / 3600) / 1.852)]; g })
hm <- function(v) { v <- gsub("^([0-9]+)\\.([0-9]+)$", "\\1:\\2", trimws(as.character(v))); p <- strsplit(v, ":"); sapply(p, function(z) if (length(z) == 2) as.numeric(z[1]) * 60 + as.numeric(z[2]) else NA_real_) }
f <- md[!duplicated(md$FilterID), ]
f$t0 <- as.numeric(as.POSIXct(sprintf("%s %02d:%02d", f$Date..YYYYMMDD., hm(f$Filtration.time.start) %/% 60, hm(f$Filtration.time.start) %% 60), format = "%Y%m%d %H:%M", tz = "America/Los_Angeles"))
f$t1 <- f$t0 + ((hm(f$Filtration.time.end) - hm(f$Filtration.time.start)) %% 1440) * 60
tr <- f[f$Collection.method == "MOSAIC-A" & !is.na(f$t0) & !is.na(f$Latitude.end), ]
ctd <- f[f$Collection.method == "MOSAIC-D", ]
ctd <- ctd[!duplicated(paste(ctd$Cruise, ifelse(is.na(ctd$CTD.Cast), ctd$FilterID, ctd$CTD.Cast))), ]
res <- do.call(rbind, lapply(seq_len(nrow(ctd)), function(i) {
  r <- ctd[i, ]; a <- tr[tr$Cruise == r$Cruise, ]
  pv <- a[a$t1 <= r$t0 + 60, ]; pv <- pv[which.max(pv$t1), ]; nx <- a[a$t0 >= r$t0 - 60, ]; nx <- nx[which.min(nx$t0), ]
  pe <- if (nrow(pv)) near(pv$Latitude.end, pv$Longitude.end) else c(lab = NA, km = NA)
  ns <- if (nrow(nx)) near(nx$Latitude.start, nx$Longitude.start) else c(lab = NA, km = NA)
  g <- G[[r$Cruise]][t >= r$t0 - 5 * 3600 & t <= r$t0 + 3600]; g[, stop := !is.na(sp) & sp < 0.5]; g[, run := rleid(stop)]
  st <- g[stop == TRUE, .(n = .N, t = t[1], lat = mean(latitude), lon = mean(longitude)), by = run][n >= 8]
  stops <- if (nrow(st)) paste(sapply(seq_len(nrow(st)), function(k) { nn <- near(st$lat[k], st$lon[k]); sprintf("%s(%skm,%dmin,%+.0fmin)", nn["lab"], nn["km"], st$n[k], (st$t[k] - r$t0) / 60) }), collapse = " ") else ""
  data.frame(Cruise = r$Cruise, Cast = r$CTD.Cast, FilterID = r$FilterID, label = qc$label_station[match(r$FilterID, qc$FilterID)],
             prev_transect_end = paste0(pe["lab"], " (", pe["km"], "km)"), next_transect_start = paste0(ns["lab"], " (", ns["km"], "km)"), ship_stops = stops)
}))
res$label_agrees <- mapply(function(l, p, n, s) !is.na(l) && (grepl(paste0("^", l, " "), p) || grepl(paste0("^", l, " "), n) || grepl(paste0("(^| )", l, "\\("), s)), res$label, res$prev_transect_end, res$next_transect_start, res$ship_stops)
write.csv(res, "outputs/ctd_label_evidence_ship_stops_neighbours.csv", row.names = FALSE)
print(res[order(res$Cruise, res$Cast), ], row.names = FALSE, right = FALSE)
