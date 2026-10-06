# =============================================================================
# 1.5  Sample-location validation and repair from the ship GPS track
# =============================================================================
# Run after 01_csvs_to_phyloseq.R and before 02_prep_detection_prob_metadata.R.
#
# What it does
#   1. reads the original phyloseq object (script 01) and keeps the real samples (not blanks);
#   2. adds the Inshore.Offshore labels;
#   3. validates each transect's field-entered start and end position against the ship's GPS track at the
#      logged filtration start and end clock times, and repairs the positions that disagree or are missing;
#   4. sets stationary (point) samples to end = start, places CTD casts where the ship was stationary
#      before filtering began (ship track), falling back to the CalCOFI station label, and recomputes transect distance (dist_km);
#   5. runs QC checks (hard checks stop the script; soft flags go to a report);
#   6. writes the sample table that script 02 reads. The control samples (field blanks) are kept in it, with their
#      collection method removed (the Fall 2023 blanks were entered as MOSAIC-A; the Summer 2024 control has no
#      method) and no position; script 02 separates them from the real samples.
#
# Rules for which positions are used (all configurable in the SETTINGS block below)
#   transect methods (track_methods): each end is checked separately. The log clock time (local time, converted to
#       UTC with `local_tz`, so daylight-saving changes are handled) is looked up in the ship track. If the field
#       position is within `track_tolerance_km` of the track it is kept; if it is further away, or missing, it is replaced
#       by the track position. A transect whose cruise has no track file, or whose time is missing, keeps its field position.
#   end time implausible: an end time that makes the filtration longer than `max_filtration_min` is tried 12 h later
#       (an AM/PM typo); if that gives a plausible duration it is used for the track lookup, otherwise the end is
#       not repaired.
#   point methods (point_methods; stationary CTD casts): the filtration time is not the cast time (water is filtered
#       afterwards while the ship moves on), so the clock time cannot be looked up directly. The ship was stationary
#       (< stop_speed_kt for >= stop_min_min) during the cast: the cast is placed at the mean ship position of the LATEST
#       such stop that began no later than stop_tol_min after the filtration start (looking back stop_window_h hours).
#       Without a stop: a field position within max_station_km of the CalCOFI station in the sample's station label is kept;
#       if it is missing or further away, the station's position is used (needs paths$station_file); with no usable station
#       label the field position is kept. End = start, distance 0. The label is then only a check (soft flag).
#   every other method (e.g. a carboy filled along the transect and filtered later): field coordinates are kept.
#   For every filter the final positions are also compared with the ship track (is the position ever on the track
#   that day?) and the result is reported as a soft flag.
# =============================================================================
suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(tibble); library(geodist); library(data.table)})

# =============================== SETTINGS =====================================
paths <- list(
  phyloseq           = "outputs/MOSAIC.ps.RDS",                                   # original phyloseq object (script 01)
  inshore_sheet      = "data/MOSAIC_sample_data_inshore.offshore_added.csv",      # hand-made; needs an Inshore.Offshore column
  station_file       = "data/reference/CalCOFIStationOrder.csv",                  # OPTIONAL: set to NULL to skip station fill/checks
  out_metadata       = "outputs/geodata_corrected_MOSAIC_metadata.csv",           # read by script 02
  out_qc             = "outputs/geodata_QC_report.csv",
  out_inshore_export = "outputs/MOSAIC_sample_data_to_add_inshore-offshore.csv")
# station_file columns: Line, Sta, `Lat (dec)`, `Lon (dec)` (the CalCOFI "Station Order" table)

# Ship GPS tracks, one file per cruise, named by the sample table's Cruise value. Each needs a UTC time column,
# a latitude column and a longitude column (names set in track_columns). A cruise with no entry here is not repaired.
ship_tracks   <- list(RS2311 = "ship_tracks/SR2311_combined_gps_coordinates.csv",
                      SR2408 = "ship_tracks/SR2408_combined_gps_coordinates.csv")
track_columns <- c(time = "timestamp", lat = "latitude", lon = "longitude")   # time must be UTC, e.g. "2023-11-04 15:44:49.65+00:00"
local_tz      <- "America/Los_Angeles"       # time zone of the Filtration.time.* clock entries (handles daylight saving)

# Study area: any start/end position outside this box stops the script. Set to NULL to skip the check.
study_area <- list(lat = c(29, 36),            # degrees north: c(min, max)
                   lon = c(-125, -116.5))      # degrees east (west is negative): c(min, max)

keep_devices        <- c("CTD", "Underway system")   # Collection.device values that are samples (the rest are blanks)
control_sample_type <- "Control"                     # Sample.type of field blanks / controls: they have no collection method
point_device        <- "CTD"                         # Collection.device of stationary samples: end = start, distance 0
point_methods       <- "MOSAIC-D"                    # Collection.method of stationary samples
track_methods       <- c("MOSAIC-A", "MOSAIC-B")     # transect methods validated/repaired from the ship track
track_tolerance_km  <- 2                             # a field position further than this from the track at its logged time is replaced
max_gap_min         <- 5                             # nearest GPS fix must be within this many minutes of the logged time
stop_speed_kt       <- 0.5                           # a ship speed below this is 'stationary' (CTD cast)
stop_min_min        <- 8                             # ...for at least this many minutes
stop_window_h       <- 5                             # a cast is looked for this many hours before the filtration start
stop_tol_min        <- 5                             # ...and a stop may begin up to this many minutes after the filtration start (clock rounding)
max_filtration_min  <- 480                           # a longer filtration means the end time is a typo
max_speed_kt        <- 14.8                          # QC: implied ship speed above this is flagged
min_transect_km     <- 0.5                           # QC: a transect shorter than this is flagged
max_station_km      <- 2                             # QC: a point sample further than this from a station is flagged
# ==============================================================================

# ---- 1. original phyloseq object -> sample table ----------------------------
ps  <- readRDS(paths$phyloseq)
sd0 <- data.frame(sample_data(ps), check.names = FALSE)
MOSAIC_sample_data <- sd0[sd0$Collection.device %in% keep_devices, ]      # same samples script 02 keeps

# Controls (field blanks): every sample that is not one of the real devices. They have no collection method; the
# Fall 2023 blanks were entered as MOSAIC-A, which made them count as MOSAIC-A samples. Remove the method for
# controls (the Summer 2024 control already has none). The raw sheet and script 01 are left unchanged.
controls <- sd0[!(sd0$Collection.device %in% keep_devices), ]
is_control <- controls$Sample.type %in% control_sample_type
stopifnot("a sample outside keep_devices is not a control" = all(is_control))
had_method <- controls$Collection.method[!is.na(controls$Collection.method)]
controls$Collection.method[is_control] <- NA
if (length(had_method)) message("Collection method removed from ", length(had_method), " control rows (",
  paste(names(table(had_method)), table(had_method), sep = ": ", collapse = ", "), "): ",
  paste(unique(controls$FilterID[!is.na(sd0[rownames(controls), "Collection.method"])]), collapse = ", "))

# export used to make the hand-edited Inshore/Offshore sheet (station 60 = cutoff)
write.csv(MOSAIC_sample_data, paths$out_inshore_export)

# Inshore.Offshore labels: taken from the hand-made sheet (by sample name). No other column of
# that sheet is used (its dist_km / depth.rank were pasted from an earlier run, its CTD end
# positions were filled in by hand and 5 CTD bottle-number labels were altered by Excel).
io <- read.csv(paths$inshore_sheet, row.names = 1, check.names = FALSE)
stopifnot(setequal(rownames(io), rownames(MOSAIC_sample_data)), !anyNA(io$Inshore.Offshore))
MOSAIC_sample_data$Inshore.Offshore <- io[rownames(MOSAIC_sample_data), "Inshore.Offshore"]

# coordinates and volume are character in the phyloseq object ("NC"/"NA" for missing)
coord_orig <- c("Latitude.start", "Longitude.start", "Latitude.end", "Longitude.end")
for (k in c(coord_orig, "Volume.filtered..mL.")) MOSAIC_sample_data[[k]] <- suppressWarnings(as.numeric(MOSAIC_sample_data[[k]]))
orig_coords <- MOSAIC_sample_data[, coord_orig]                       # field entries, kept for the report

hav <- function(a, b, c, e) { p <- pi/180; 2*6371*asin(pmin(1, sqrt(sin((c-a)*p/2)^2 + cos(a*p)*cos(c*p)*sin((e-b)*p/2)^2))) }

# ---- 2. ship tracks ------------------------------------------------------------
tracks <- lapply(ship_tracks, function(f) {
  stopifnot("ship track file not found" = file.exists(f))
  g <- fread(f, select = unname(track_columns)); setnames(g, unname(track_columns), c("time", "lat", "lon"))
  g[, t := as.numeric(as.POSIXct(substr(time, 1, 19), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))]   # UTC epoch seconds
  g <- g[!is.na(t) & !is.na(lat) & !is.na(lon), .(t, lat, lon)]; setkey(g, t); g })
no_track <- setdiff(unique(MOSAIC_sample_data$Cruise), names(tracks))
if (length(no_track)) message("No ship track for cruise(s): ", paste(no_track, collapse = ", "), " - their field positions are kept")

# position of the ship at each time in `tt` (epoch s): nearest GPS fix, NA if none within max_gap_min
track_at <- function(cruise, tt) {
  out <- data.frame(lat = rep(NA_real_, length(tt)), lon = NA_real_, gap_min = NA_real_)
  for (cr in intersect(unique(cruise[!is.na(tt)]), names(tracks))) {
    g <- tracks[[cr]]; i <- which(cruise == cr & !is.na(tt))
    j <- pmin(pmax(findInterval(tt[i], g$t), 1L), nrow(g) - 1L)
    k <- ifelse(abs(tt[i] - g$t[j]) <= abs(g$t[j + 1L] - tt[i]), j, j + 1L)
    out$lat[i] <- g$lat[k]; out$lon[i] <- g$lon[k]; out$gap_min[i] <- abs(tt[i] - g$t[k]) / 60 }
  out[!is.na(out$gap_min) & out$gap_min > max_gap_min, c("lat", "lon")] <- NA
  out }
# how far is a position from the track within +-window_h of t_ref, and when (minutes from t_ref) is the track nearest to it?
near_track <- function(cruise, lat, lon, t_ref, window_h = 36) {
  out <- data.frame(km = rep(NA_real_, length(lat)), offset_min = NA_real_)
  for (i in which(!is.na(lat) & !is.na(lon) & !is.na(t_ref) & cruise %in% names(tracks))) {
    g <- tracks[[cruise[i]]]; r <- findInterval(c(t_ref[i] - window_h * 3600, t_ref[i] + window_h * 3600), g$t)
    if (r[2] <= r[1]) next
    ix <- (r[1] + 1L):r[2]; d <- hav(lat[i], lon[i], g$lat[ix], g$lon[ix]); m <- which.min(d)
    out$km[i] <- d[m]; out$offset_min[i] <- (g$t[ix[m]] - t_ref[i]) / 60 }
  out }

# every stationary period of the ship (one row per stop: cruise, start time, minutes, mean position), from 1-minute fixes
stops <- rbindlist(lapply(names(tracks), function(cr) {
  g <- tracks[[cr]][!duplicated(floor(t / 60))]
  g[, sp := c(NA, hav(lat[-.N], lon[-.N], lat[-1], lon[-1]) / (diff(t) / 3600) / 1.852)]
  g[, slow := !is.na(sp) & sp < stop_speed_kt][, run := rleid(slow)]
  st <- g[slow == TRUE, .(n = .N, t = t[1], lat = mean(lat), lon = mean(lon)), by = run][n >= stop_min_min]
  st[, cruise := cr][, run := NULL] }))

# ---- 3. filtration clock times -------------------------------------------------
# Known typo "6.05" -> "6:05" (as in script 02). Times are local clock; the end is start + duration (may pass midnight).
hm <- function(v) { v <- gsub("^([0-9]+)\\.([0-9]+)$", "\\1:\\2", trimws(as.character(v)))
  ok <- grepl("^[0-9]{1,2}:[0-9]{2}$", v); out <- rep(NA_real_, length(v))
  p <- strsplit(v[ok], ":"); out[ok] <- sapply(p, function(z) as.numeric(z[1]) * 60 + as.numeric(z[2])); out }
s_min <- hm(MOSAIC_sample_data$Filtration.time.start); e_min <- hm(MOSAIC_sample_data$Filtration.time.end)
dur <- (e_min - s_min) %% 1440
# an end time that makes the filtration implausibly long is tried 12 h later (AM/PM typo); otherwise it is not used
too_long <- !is.na(dur) & dur > max_filtration_min
dur12    <- (e_min + 720 - s_min) %% 1440
end_12h  <- too_long & dur12 <= max_filtration_min
dur[end_12h] <- dur12[end_12h]
dur[too_long & !end_12h] <- NA
t_start <- as.numeric(as.POSIXct(sprintf("%s %02d:%02d", MOSAIC_sample_data$Date..YYYYMMDD., s_min %/% 60, s_min %% 60),
                                 format = "%Y%m%d %H:%M", tz = local_tz))                      # UTC epoch seconds
t_end   <- t_start + dur * 60
MOSAIC_sample_data$dur_min <- dur

# ---- 4. transects: validate against the ship track, repair ------------------------
is_trk <- MOSAIC_sample_data$Collection.method %in% track_methods
cr     <- MOSAIC_sample_data$Cruise
ts <- track_at(cr, ifelse(is_trk, t_start, NA_real_))
te <- track_at(cr, ifelse(is_trk, t_end,   NA_real_))
off_s <- hav(orig_coords$Latitude.start, orig_coords$Longitude.start, ts$lat, ts$lon)    # field position vs track at logged time
off_e <- hav(orig_coords$Latitude.end,   orig_coords$Longitude.end,   te$lat, te$lon)
rep_s <- !is.na(ts$lat) & (is.na(off_s) | off_s > track_tolerance_km)
rep_e <- !is.na(te$lat) & (is.na(off_e) | off_e > track_tolerance_km)
MOSAIC_sample_data[rep_s, c("Latitude.start", "Longitude.start")] <- ts[rep_s, c("lat", "lon")]
MOSAIC_sample_data[rep_e, c("Latitude.end",   "Longitude.end")]   <- te[rep_e, c("lat", "lon")]
src_start <- ifelse(rep_s, "ship track", ifelse(is_trk & !is.na(ts$lat), "original (matches track)", "original"))
src_end   <- ifelse(rep_e, "ship track", ifelse(is_trk & !is.na(te$lat), "original (matches track)", "original"))
# evidence for the repaired ends: where on the track is the field position? (a position that sits on the track a few
# minutes earlier/later than logged points to a clock-time error rather than a coordinate error)
ev_s <- near_track(cr, ifelse(rep_s, orig_coords$Latitude.start, NA), orig_coords$Longitude.start, t_start)
ev_e <- near_track(cr, ifelse(rep_e, orig_coords$Latitude.end,   NA), orig_coords$Longitude.end,   t_end)

# ---- 5. point samples and distances ---------------------------------------------
use_stations <- !is.null(paths$station_file) && file.exists(paths$station_file)
if (use_stations) {
  stn <- read.csv(paths$station_file, check.names = FALSE)
  stn <- data.frame(lab = paste(as.numeric(stn$Line), as.numeric(stn$Sta)), lat = stn$`Lat (dec)`, lon = stn$`Lon (dec)`)
} else message("Station fill and checks skipped (paths$station_file is NULL or missing).")
# station label "86.7 80" / "93.3_80" -> "86.7 80"; anything else -> NA
lab_key <- function(x) { x <- trimws(as.character(x)); x[x %in% c("NA", "")] <- NA
  vapply(strsplit(x, "[ _]+"), function(z) if (length(z) == 2 && !anyNA(suppressWarnings(as.numeric(z))))
    paste(as.numeric(z[1]), as.numeric(z[2])) else NA_character_, "") }
is_point <- MOSAIC_sample_data$Collection.device == point_device
is_pm    <- MOSAIC_sample_data$Collection.method %in% point_methods
key <- lab_key(MOSAIC_sample_data$Station)
# cast position from the ship's stationary period before filtering began
stop_lat <- stop_lon <- stop_offset <- stop_dur <- rep(NA_real_, nrow(MOSAIC_sample_data))
for (i in which(is_pm & !is.na(t_start) & cr %in% stops$cruise)) {
  a <- stops[cruise == cr[i] & t <= t_start[i] + stop_tol_min * 60 & t >= t_start[i] - stop_window_h * 3600]
  if (nrow(a)) { a <- a[which.max(t)]; stop_lat[i] <- a$lat; stop_lon[i] <- a$lon; stop_offset[i] <- (a$t - t_start[i]) / 60; stop_dur[i] <- a$n }
}
has_stop <- !is.na(stop_lat)
off_stop_km <- ifelse(has_stop, hav(orig_coords$Latitude.start, orig_coords$Longitude.start, stop_lat, stop_lon), NA_real_)
MOSAIC_sample_data$Latitude.start[has_stop]  <- stop_lat[has_stop]
MOSAIC_sample_data$Longitude.start[has_stop] <- stop_lon[has_stop]
src_start[has_stop] <- "ship stop"
if (use_stations) {
  labelled <- is_pm & !has_stop & !is.na(key) & key %in% stn$lab
  off_label_km <- ifelse(labelled, hav(orig_coords$Latitude.start, orig_coords$Longitude.start, stn$lat[match(key, stn$lab)], stn$lon[match(key, stn$lab)]), NA_real_)
  fill <- labelled & (is.na(off_label_km) | off_label_km > max_station_km)    # field position missing or not within max_station_km of its labelled station
  MOSAIC_sample_data$Latitude.start[fill]  <- stn$lat[match(key[fill], stn$lab)]
  MOSAIC_sample_data$Longitude.start[fill] <- stn$lon[match(key[fill], stn$lab)]
  src_start[fill] <- "CalCOFI station label"
}
MOSAIC_sample_data$Longitude.end[is_point] <- MOSAIC_sample_data$Longitude.start[is_point]
MOSAIC_sample_data$Latitude.end[is_point]  <- MOSAIC_sample_data$Latitude.start[is_point]
src_end[is_point] <- src_start[is_point]
start <- cbind(MOSAIC_sample_data$Longitude.start, MOSAIC_sample_data$Latitude.start)
end   <- cbind(MOSAIC_sample_data$Longitude.end,   MOSAIC_sample_data$Latitude.end)
colnames(start) <- colnames(end) <- c("x", "y")
ok <- stats::complete.cases(start, end)
MOSAIC_sample_data$dist_km <- NA_real_
MOSAIC_sample_data$dist_km[ok] <- geodist(start[ok, , drop = FALSE], end[ok, , drop = FALSE], paired = TRUE, measure = "geodesic") / 1000
MOSAIC_sample_data$dist_km[is_point & !is.na(MOSAIC_sample_data$Latitude.start)] <- 0

# ---- 6. QC: does the corrected sheet make sense? (checked before script 02 reads it) ----
# Hard checks stop the script. Soft flags are written to the QC report and printed; they are things
# to look at, not errors. Station checks run only if paths$station_file is set and the file exists.
nearest <- function(la, lo) { if (!use_stations || is.na(la)) return(list(lab = NA_character_, km = NA_real_))
  k <- hav(la, lo, stn$lat, stn$lon); i <- which.min(k); list(lab = stn$lab[i], km = k[i]) }

q <- MOSAIC_sample_data
q$row <- rownames(q)
q$t0 <- as.POSIXct(t_start, origin = "1970-01-01", tz = "UTC"); q$t1 <- as.POSIXct(t_end, origin = "1970-01-01", tz = "UTC")
q$src_start <- src_start; q$src_end <- src_end; q$end_time_12h <- end_12h
q$orig_start_vs_track_km <- off_s; q$orig_end_vs_track_km <- off_e
q$ev_start_offset_min <- ev_s$offset_min; q$ev_end_offset_min <- ev_e$offset_min
if (!exists("off_label_km")) off_label_km <- rep(NA_real_, nrow(q))
q$off_label_km <- off_label_km; q$off_stop_km <- off_stop_km; q$stop_offset_min <- stop_offset; q$stop_dur_min <- stop_dur
q$ev_start_km <- ev_s$km; q$ev_end_km <- ev_e$km

# --- hard checks (every row) ---
rep_ok <- tapply(seq_len(nrow(q)), q$FilterID, function(i)
  all(sapply(q[i, c(coord_orig, "dist_km")], function(v) all(is.na(v)) || (!anyNA(v) && diff(range(v)) < 1e-6))))
stopifnot("replicates of one filter have different positions" = all(rep_ok))
stopifnot("a coordinate is not a valid latitude/longitude" =
          all(abs(c(q$Latitude.start, q$Latitude.end)) <= 90, na.rm = TRUE) &
          all(abs(c(q$Longitude.start, q$Longitude.end)) <= 180, na.rm = TRUE))
if (!is.null(study_area)) {
  lat_all <- c(q$Latitude.start, q$Latitude.end); lon_all <- c(q$Longitude.start, q$Longitude.end)
  outside <- !is.na(lat_all) & (lat_all < study_area$lat[1] | lat_all > study_area$lat[2] |
                                lon_all < study_area$lon[1] | lon_all > study_area$lon[2])
  stopifnot("a start or end position is outside study_area (see SETTINGS)" = !any(outside))
}
stopifnot("a point sample has different start/end or a distance other than 0" =
          all(q$Latitude.start[is_point] == q$Latitude.end[is_point] & q$dist_km[is_point] == 0, na.rm = TRUE))
stopifnot("negative distance" = all(q$dist_km >= 0, na.rm = TRUE))

# --- soft flags (one row per filter) ---
f <- q[!duplicated(q$FilterID), ]
f$speed_kt <- f$dist_km / (f$dur_min / 60) / 1.852
nr <- lapply(seq_len(nrow(f)), function(i) nearest(f$Latitude.start[i], f$Longitude.start[i]))
f$nearest_station <- vapply(nr, `[[`, "", "lab"); f$km_to_nearest_station <- round(vapply(nr, `[[`, 0, "km"), 2)
f$label_station <- lab_key(f$Station)
# is the final position ever on the ship track that day? (t_ref = filtration start, or midday local if no time)
t_ref <- ifelse(is.na(as.numeric(f$t0)), as.numeric(as.POSIXct(paste(f$Date..YYYYMMDD., "12:00"), format = "%Y%m%d %H:%M", tz = local_tz)), as.numeric(f$t0))
on_s <- near_track(f$Cruise, f$Latitude.start, f$Longitude.start, t_ref)
on_e <- near_track(f$Cruise, f$Latitude.end,   f$Longitude.end,   t_ref)
f$start_to_track_km <- round(on_s$km, 2); f$end_to_track_km <- round(on_e$km, 2)
# stations at the end of the previous / start of the next timed transect (time-wise neighbours of each point sample)
timed <- f[!(f$Collection.method %in% point_methods) & !is.na(f$t1) & !is.na(f$Latitude.end), ]
f$prev_end_station <- f$next_start_station <- NA_character_
if (use_stations) for (i in which(f$Collection.method %in% point_methods & !is.na(f$t0))) {
  a <- timed[timed$Cruise == f$Cruise[i], ]
  pv <- a[a$t1 <= f$t0[i] + 60, ]; pv <- pv[which.max(pv$t1), ]
  nx <- a[a$t0 >= f$t0[i] - 60, ]; nx <- nx[which.min(nx$t0), ]
  if (nrow(pv)) f$prev_end_station[i]   <- nearest(pv$Latitude.end,   pv$Longitude.end)$lab
  if (nrow(nx)) f$next_start_station[i] <- nearest(nx$Latitude.start, nx$Longitude.start)$lab
}
flag <- function(cond, txt) ifelse(!is.na(cond) & cond, paste0(txt, "; "), "")
is_tr <- !(f$Collection.method %in% point_methods)       # transects of any method
is_pt <- f$Collection.method %in% point_methods
f$flags <- paste0(
  flag(f$src_start == "ship track" & !is.na(f$orig_start_vs_track_km), sprintf("start repaired from ship track (field position %.0f km off)", f$orig_start_vs_track_km)),
  flag(f$src_start == "ship track" & is.na(f$orig_start_vs_track_km),  "start was missing; taken from ship track"),
  flag(f$src_end   == "ship track" & !is.na(f$orig_end_vs_track_km) & is_tr, sprintf("end repaired from ship track (field position %.0f km off)", f$orig_end_vs_track_km)),
  flag(f$src_end   == "ship track" & is.na(f$orig_end_vs_track_km) & is_tr,  "end was missing; taken from ship track"),
  flag(f$end_time_12h, "end time looked 12 h off; +12 h used for the track lookup"),
  flag(f$src_start == "ship stop" & !is.na(f$off_stop_km) & f$off_stop_km > max_station_km, sprintf("field position was %.0f km from the ship stop; ship stop used", f$off_stop_km)),
  flag(f$src_start == "ship stop" & is.na(f$off_stop_km), "position taken from the ship stop (no field position)"),
  flag(is_pt & f$src_start != "ship stop" & !is.na(f$Latitude.start), "no ship stop found before filtering; station label / field position used"),
  flag(f$src_start == "CalCOFI station label" & is.na(f$off_label_km), "position taken from the CalCOFI station label (no field position)"),
  flag(f$src_start == "CalCOFI station label" & !is.na(f$off_label_km) & f$off_label_km > max_station_km, sprintf("field position was %.0f km from the labelled station; station used", f$off_label_km)),
  flag(is.na(f$Latitude.start), "no start position"),
  flag(is_tr & is.na(f$Latitude.end), "no end position"),
  flag(is_tr & is.na(f$dur_min), "no usable filtration end time"),
  flag(is_tr & f$speed_kt > max_speed_kt & f$Collection.method %in% track_methods, sprintf("implied speed > %g kt", max_speed_kt)),
  flag(is_tr & f$dist_km < min_transect_km, sprintf("transect < %g km", min_transect_km)),
  flag(f$src_start != "ship track" & !is.na(f$start_to_track_km) & f$start_to_track_km > track_tolerance_km, "start position never on the ship track that day"),
  flag(f$src_end   != "ship track" & is_tr & !is.na(f$end_to_track_km) & f$end_to_track_km > track_tolerance_km, "end position never on the ship track that day"),
  flag(is_pt & f$km_to_nearest_station > max_station_km, sprintf("point sample > %g km from a station", max_station_km)),
  flag(is_pt & !is.na(f$label_station) & f$label_station != f$nearest_station, "point sample is not on its labelled station"),
  flag(is_pt & !is.na(f$prev_end_station) & f$nearest_station != f$prev_end_station & f$nearest_station != f$next_start_station,
       "neither neighbouring transect meets this point sample's station"))
qc <- f[, c("FilterID", "Cruise", "Collection.method", "src_start", "src_end",
            "Latitude.start", "Longitude.start", "Latitude.end", "Longitude.end", "dist_km", "dur_min", "speed_kt",
            "orig_start_vs_track_km", "orig_end_vs_track_km", "off_stop_km", "stop_offset_min", "stop_dur_min", "ev_start_offset_min", "ev_end_offset_min",
            "start_to_track_km", "end_to_track_km",
            "nearest_station", "km_to_nearest_station", "label_station", "prev_end_station", "next_start_station", "flags")]
orig_f <- orig_coords[match(qc$FilterID, MOSAIC_sample_data$FilterID), ]
qc <- cbind(qc[, 1:5], orig_Latitude.start = orig_f$Latitude.start, orig_Longitude.start = orig_f$Longitude.start,
            orig_Latitude.end = orig_f$Latitude.end, orig_Longitude.end = orig_f$Longitude.end, qc[, -(1:5)])
write.csv(qc, paths$out_qc, row.names = FALSE)
cat(sprintf("\nQC: %d filters checked, hard checks passed; %d filters carry soft flags (%s)\n", nrow(qc), sum(qc$flags != ""), paths$out_qc))
print(qc[qc$flags != "", c("FilterID", "Collection.method", "src_start", "src_end", "dist_km", "speed_kt", "flags")], row.names = FALSE)

# ---- 7. report and write ----------------------------------------------------
cat("\nposition source, start / end (rows):\n"); print(table(MOSAIC_sample_data$Collection.method, src_start)); print(table(MOSAIC_sample_data$Collection.method, src_end))
cat("rows with no start position:", sum(is.na(MOSAIC_sample_data$Latitude.start)),
    " | rows with no end position:", sum(is.na(MOSAIC_sample_data$Latitude.end)), "\n")
out <- MOSAIC_sample_data[, setdiff(names(MOSAIC_sample_data), "dur_min")]                       # real samples: same columns as before (+ Inshore.Offshore, dist_km)
ctrl <- controls; for (k in c(coord_orig, "Volume.filtered..mL.")) ctrl[[k]] <- suppressWarnings(as.numeric(ctrl[[k]]))
ctrl$Inshore.Offshore <- NA_character_; ctrl$dist_km <- NA_real_                                  # controls have no position
all_samples <- rbind(out, ctrl[, names(out)])[rownames(sd0), ]                                     # original row order, controls included
stopifnot(nrow(all_samples) == nrow(sd0), !anyNA(rownames(all_samples)))
write.csv(all_samples, paths$out_metadata)
cat("Saved:", paths$out_metadata, "\n")
