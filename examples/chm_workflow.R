# =============================================================================
# copc4R -- Example: CHM Workflow with lidR
# =============================================================================
#
# Demonstrates a complete canopy height model (CHM) pipeline using
# copc4R + lidR, pulling data directly from Microsoft Planetary Computer's
# USGS 3DEP COPC collection via HTTP range reads (no bulk downloads).
#
#   1. Define a WGS 84 AOI rectangle
#   2. STAC search + SAS token from Planetary Computer
#   3. Loop over ALL overlapping tiles -- copc_info(), read_copc(),
#      as_las() -> merge into one LAS object
#   4–6. lidR pipeline -- classify_ground -> normalize_height ->
#                        rasterize_canopy
#
# Required packages:
#   install.packages(c("httr", "jsonlite", "sf", "lidR", "terra"))
#
# Data: https://planetarycomputer.microsoft.com/dataset/3dep-lidar-copc
# =============================================================================

library(copc4R)
library(lidR)
library(sf)
library(httr)
library(jsonlite)

cat("=== copc4R CHM Workflow (Planetary Computer 3DEP) ===\n\n")

# -- 1. Define a rectangular AOI in WGS 84 ----------------------------------
# Vegetated area near the NAU School of Forestry, Flagstaff AZ.
xmin <- -111.654;  xmax <- -111.649
ymin <-   35.1682; ymax <-   35.1707

aoi_coords <- matrix(c(
  xmin, ymin,
  xmax, ymin,
  xmax, ymax,
  xmin, ymax,
  xmin, ymin   # close the ring
), ncol = 2, byrow = TRUE)

aoi <- st_sf(
  id       = 1,
  geometry = st_sfc(st_polygon(list(aoi_coords)), crs = 4326)
)

cat("AOI (WGS 84):\n")
print(st_bbox(aoi))
cat("\n")

# -- 2. Search Planetary Computer STAC for overlapping 3DEP tiles -----------
bbox_vec <- c(xmin, ymin, xmax, ymax)

resp <- POST(
  "https://planetarycomputer.microsoft.com/api/stac/v1/search",
  content_type_json(),
  body = toJSON(list(
    collections = list("3dep-lidar-copc"),
    bbox         = bbox_vec,
    limit        = 50
  ), auto_unbox = TRUE)
)
stop_for_status(resp)
items <- fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                  simplifyVector = FALSE)$features
cat("STAC search returned", length(items), "tile(s)\n")
if (length(items) == 0L) stop("No 3DEP COPC tiles found for this AOI.")

# Get a SAS token to sign tile URLs
sas <- fromJSON(content(
  GET("https://planetarycomputer.microsoft.com/api/sas/v1/token/3dep-lidar-copc"),
  as = "text", encoding = "UTF-8"
))$token

# -- 3. Read ALL points inside AOI from every overlapping tile --------------
all_las <- list()

for (i in seq_along(items)) {
  item <- items[[i]]
  href <- item$assets$data$href
  url  <- paste0(href, "?", sas)

  cat(sprintf("-- Tile %d/%d: %s --\n", i, length(items), item$id))

  # Quick metadata peek
  info <- copc_info(url)
  cat(sprintf("  Points in tile: %s   Spacing: %.2f\n",
              format(info$point_count, big.mark = ","), info$spacing))

  # Discover the tile's native CRS from its WKT VLR
  hdr <- read_copc_header(url)
  tile_crs <- NULL
  for (v in hdr[["Variable Length Records"]]) {
    wkt <- v[["WKT OGC COORDINATE SYSTEM"]]
    if (!is.null(wkt) && nzchar(wkt)) { tile_crs <- st_crs(wkt); break }
  }
  if (is.null(tile_crs)) {
    cat("  ⚠ Cannot determine CRS -- skipping\n")
    next
  }

  # Reproject AOI to the tile's CRS and range-read all points inside
  aoi_proj <- st_transform(aoi, tile_crs)
  result   <- read_copc(url,
                        aoi      = aoi_proj,
                        select   = "xyzicrnap",
                        filter   = "-drop_withheld",
                        progress = TRUE)

  if (nrow(result$data) > 0L) {
    cat(sprintf("  Points read: %s\n", format(nrow(result$data), big.mark = ",")))
    all_las[[length(all_las) + 1L]] <- as_las(result)
  } else {
    cat("  0 points inside AOI -- skipping\n")
  }
}

if (length(all_las) == 0L)
  stop("No points were read from any tile.")

# Merge tiles into a single LAS object
las <- do.call(rbind, all_las)
cat(sprintf("\nTotal points from %d tile(s): %s\n",
            length(all_las), format(npoints(las), big.mark = ",")))
cat("Columns:", paste(names(las@data), collapse = ", "), "\n")
cat("\nClassification counts:\n")
print(table(las$Classification))
cat("\n")

# -- 4. Classify ground with CSF --------------------------------------------
cat("Classifying ground (CSF)...\n")
las <- classify_ground(las, algorithm = csf(
  sloop_smooth     = FALSE,
  class_threshold  = 0.5,
  cloth_resolution = 0.5,
  rigidness        = 1L,
  time_step        = 0.65
))
cat("Ground points:", sum(las$Classification == 2L), "\n\n")

# -- 5. Normalize heights with TIN ------------------------------------------
cat("Normalizing heights (TIN)...\n")
nlas <- normalize_height(las, algorithm = tin())
cat("Height range:",
    round(min(nlas$Z), 2), "–", round(max(nlas$Z), 2), "\n\n")

# -- 6. Rasterize canopy height model --------------------------------------
cat("Building CHM (pitfree)...\n")
chm <- rasterize_canopy(nlas, res = 0.5,
                        algorithm = pitfree(
                          thresholds = c(0, 2, 5, 10, 15),
                          max_edge   = c(0, 1.5)
                        ))
cat("CHM raster:", terra::nrow(chm), "x", terra::ncol(chm),
    "cells at", terra::res(chm)[1], "m\n")
cat("CHM range: ", round(terra::minmax(chm)[1], 2), "–",
    round(terra::minmax(chm)[2], 2), "m\n\n")

# -- 7. Display ------------------------------------------------------------
terra::plot(chm, col = lidR::height.colors(50),
            main = "Canopy Height Model -- NAU / Flagstaff (copc4R + 3DEP)")

cat("=== CHM workflow complete ===\n")
