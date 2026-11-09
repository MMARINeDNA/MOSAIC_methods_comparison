# =============================================================================
# CTD station from the transects either side of the cast (time-wise)
# =============================================================================
# The ship stops on a CalCOFI station for the CTD: the preceding MOSAIC-A/B transect ends
# there and the next one starts there. For each CTD filter, find the A/B transect that
# ended most recently before its filtration start and the one that starts next, take the
# official CalCOFI station nearest each point, and compare with the position used for the
# CTD and with the station named in the sample label.
# Input : outputs/MOSAIC_taxnamed_sample_data.csv (positions as used by the pipeline),
#         data/reference/CalCOFIStationOrder.csv
# Output: outputs/ctd_station_from_neighbours.csv
# =============================================================================
suppressPackageStartupMessages(library(dplyr)); options(width = 250)
num <- function(v) suppressWarnings(as.numeric(v))
hav <- function(a, b, c, e) { p <- pi/180; 2*6371*asin(pmin(1, sqrt(sin((c-a)*p/2)^2 + cos(a*p)*cos(c*p)*sin((e-b)*p/2)^2))) }
st <- read.csv("data/reference/CalCOFIStationOrder.csv", check.names = FALSE)
st <- data.frame(lab = paste(as.numeric(st$Line), as.numeric(st$Sta)), lat = st$`Lat (dec)`, lon = st$`Lon (dec)`)
nst <- function(la, lo) { if (is.na(la)) return(c(NA, NA)); k <- hav(la, lo, st$lat, st$lon); i <- which.min(k); c(i, k[i]) }
lab_of <- function(la, lo) { r <- nst(la, lo); if (is.na(r[1])) NA_character_ else sprintf("%s (%.1f km)", st$lab[r[1]], r[2]) }

md <- read.csv("outputs/MOSAIC_taxnamed_sample_data.csv", row.names = 1, stringsAsFactors = FALSE)
f <- md[!duplicated(md$FilterID), ]
f$t0 <- as.POSIXct(f$datetime_start, tz = "UTC"); f$dur <- num(f$filtration_time_min); f$t1 <- f$t0 + f$dur * 60
ab <- f[f$Collection.method %in% c("MOSAIC-A", "MOSAIC-B") & !is.na(f$t1), ]
ctd <- f[f$Collection.method == "MOSAIC-D", ]

res <- lapply(seq_len(nrow(ctd)), function(i) {
  x <- ctd[i, ]; a <- ab[ab$Cruise == x$Cruise, ]
  prev <- a[a$t1 <= x$t0 + 1, ]; prev <- prev[which.max(prev$t1), ]
  nxt  <- a[a$t0 >= x$t0 - 1, ]; nxt  <- nxt[which.min(nxt$t0), ]
  data.frame(FilterID = x$FilterID, Cruise = x$Cruise, filt_start = format(x$t0, "%m-%d %H:%M"),
    used_position_station = lab_of(x$Latitude.start, x$Longitude.start),
    label_station = ifelse(is.na(x$Station) | x$Station == "NA", NA, x$Station),
    prev_transect = prev$FilterID, prev_end_station = lab_of(prev$Latitude.end, prev$Longitude.end),
    prev_gap_h = round(as.numeric(difftime(x$t0, prev$t1, units = "hours")), 1),
    next_transect = nxt$FilterID, next_start_station = lab_of(nxt$Latitude.start, nxt$Longitude.start),
    next_gap_h = round(as.numeric(difftime(nxt$t0, x$t0, units = "hours")), 1), stringsAsFactors = FALSE)
})
r <- bind_rows(res)
stn_only <- function(s) sub(" \\(.*", "", s)
r$prev_matches_used <- stn_only(r$prev_end_station) == stn_only(r$used_position_station)
r$next_matches_used <- stn_only(r$next_start_station) == stn_only(r$used_position_station)
r$prev_next_agree   <- stn_only(r$prev_end_station) == stn_only(r$next_start_station)
write.csv(r, "outputs/ctd_station_from_neighbours.csv", row.names = FALSE)
cat("CTD filters:", nrow(r), "\n")
cat("previous transect ends at the used station:", sum(r$prev_matches_used, na.rm = TRUE),
    "| next transect starts there:", sum(r$next_matches_used, na.rm = TRUE),
    "| either:", sum(r$prev_matches_used | r$next_matches_used, na.rm = TRUE), "\n")
cat("previous end station == next start station:", sum(r$prev_next_agree, na.rm = TRUE), "\n\n")
print(r %>% select(FilterID, filt_start, used_position_station, label_station, prev_end_station, prev_gap_h, next_start_station, next_gap_h, prev_matches_used, next_matches_used), row.names = FALSE)
