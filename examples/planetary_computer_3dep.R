# =============================================================================
# copc4R -- Example: Microsoft Planetary Computer -- USGS 3DEP COPC
# =============================================================================
#
# Demonstrates two approaches for fetching USGS 3DEP LiDAR data from
# Microsoft's Planetary Computer:
#
#   A. The easy way -- download_3dep_copc() wraps everything in one call
#   B. The manual way -- STAC query + per-tile COPC range reads
#
# Both approaches use COPC range reads under the hood: only the octree
# nodes overlapping the AOI are fetched, not entire tiles.
#
# Data: https://planetarycomputer.microsoft.com/dataset/3dep-lidar-copc
# STAC: https://planetarycomputer.microsoft.com/api/stac/v1
#
# Required packages:
#   install.packages(c("httr", "jsonlite", "sf", "lidR"))
# =============================================================================

library(sf)
library(copc4R)

# =============================================================================
# 1. Define the AOI as an sf polygon
# =============================================================================
# Small polygon near Flagstaff, AZ (NAU campus area).
# Change this aoi to your area of interest.

aoi_file <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82.gpkg",
                        package = "copc4R")
aoi <- read_sf(aoi_file)

cat("Area of interest (WGS 84):\n")
print(st_bbox(aoi))
cat("\n")

dest_dir <- path.expand("~/3dep_tiles")


# =============================================================================
# Approach A -- One-liner with download_3dep_copc()
# =============================================================================
# Handles STAC search, SAS token, CRS reprojection, range reads,
# polygon clipping, and .laz writing automatically.

cat("===============================================================\n")
cat("  Approach A: download_3dep_copc()  (recommended)\n")
cat("===============================================================\n\n")

result <- download_3dep_copc(
  aoi       = aoi,
  dest_dir  = dest_dir,
  select    = "*",
  filter    = "-keep_voxel 2",     # 2-m voxel thinning
  merge     = TRUE,                # single merged output
  overwrite = TRUE,
  progress  = TRUE,
  verbose   = TRUE
)

cat("\nOutput files:\n")
print(result)

# Read the merged file back in with copc4R or lidR
if (file.exists(result$file[1])) {
  cat(sprintf("\nMerged file: %s (%.1f MB, %s pts)\n",
              basename(result$file[1]),
              file.size(result$file[1]) / 1e6,
              format(result$n_points[1], big.mark = ",")))
}


# =============================================================================
# Approach B -- Manual STAC query + per-tile COPC range reads
# =============================================================================
# For advanced users who want full control over tile iteration.

cat("\n===============================================================\n")
cat("  Approach B: manual STAC + COPC range reads\n")
cat("===============================================================\n\n")

library(httr)
library(jsonlite)
library(lidR)

# 2. Search the Planetary Computer STAC API
bbox_vec <- as.numeric(st_bbox(aoi))

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

cat("STAC search returned", length(items), "tile(s)\n\n")
if (length(items) == 0)
  stop("No 3DEP COPC tiles found for this AOI.")

# 3. Get a SAS token to sign the URLs
sas <- fromJSON(content(
  GET("https://planetarycomputer.microsoft.com/api/sas/v1/token/3dep-lidar-copc"),
  as = "text", encoding = "UTF-8"
))$token

# 4. Read ONLY the AOI portion from each tile via COPC range reads
all_las <- list()

for (i in seq_along(items)) {
  item <- items[[i]]
  href <- item$assets$data$href
  signed_url <- paste0(href, "?", sas)

  cat(sprintf("-- Tile %d/%d: %s --\n", i, length(items), item$id))

  # Read header to discover the tile's native CRS
  hdr <- read_copc_header(signed_url)

  tile_crs <- NULL
  for (v in hdr[["Variable Length Records"]]) {
    wkt <- v[["WKT OGC COORDINATE SYSTEM"]]
    if (!is.null(wkt) && nzchar(wkt)) { tile_crs <- st_crs(wkt); break }
  }
  if (is.null(tile_crs))
    stop("Cannot determine CRS for tile ", item$id)

  # Reproject AOI into the tile's native CRS and range-read
  aoi_proj <- st_transform(aoi, tile_crs)
  res <- read_copc(signed_url, aoi = aoi_proj,
                   filter = "-keep_voxel 2",
                   progress = TRUE)

  if (nrow(res$data) > 0L) {
    cat(sprintf("  Points: %s\n", format(nrow(res$data), big.mark = ",")))
    all_las[[length(all_las) + 1L]] <- as_las(res)
  } else {
    cat("  0 points in AOI -- skipping\n")
  }
}

if (length(all_las) == 0L)
  stop("No points were read from any tile.")

# 5. Merge and write
las <- do.call(rbind, all_las)
cat(sprintf("\nTotal points from %d tile(s): %s\n",
            length(all_las), format(npoints(las), big.mark = ",")))

out_laz <- file.path(dest_dir, "flagstaff_aoi_manual.laz")
writeLAS(las, out_laz)
cat("Wrote:", out_laz, "\n")

# 6. Visualize
cat("\nClassification counts:\n")
print(table(las$Classification))

plot(las)

cat("\n=== Workflow complete ===\n")
