# =============================================================================
# copc4R  –  Example: Microsoft Planetary Computer – USGS 3DEP COPC
# =============================================================================
#
# This example demonstrates how to:
#   1. Define an AOI polygon in WGS 84 (lon/lat)
#   2. Query the Planetary Computer STAC API for overlapping 3DEP COPC tiles
#   3. Let the user choose which dataset to load (interactive)
#   4. Read COPC point cloud data within the AOI bounding box
#   5. Clip the points to the polygon boundary
#   6. Visualize the result
#
# The Planetary Computer hosts nationwide USGS 3DEP LiDAR data as
# Cloud Optimized Point Clouds (COPC). Each tile is a .copc.laz file
# served from Azure Blob Storage.
#
# Data info: https://planetarycomputer.microsoft.com/dataset/3dep-lidar-copc
# STAC API:  https://planetarycomputer.microsoft.com/api/stac/v1
#
# Required packages:
#   install.packages(c("httr", "jsonlite", "sf", "lidR"))
# =============================================================================

library(copc4R)
library(lidR)
library(lidRviewer)

# ── Helper functions ──────────────────────────────────────────────────────────

#' Search the Planetary Computer for 3DEP COPC tiles intersecting a bbox
#'
#' @param bbox  Numeric vector c(xmin, ymin, xmax, ymax) in WGS 84 (EPSG:4326)
#' @param limit Max number of items to return (max 1000)
#' @return A list of STAC items (features)
search_3dep_copc <- function(bbox, limit = 50) {
  if (!requireNamespace("httr", quietly = TRUE))
    stop("Package 'httr' is required. Install with: install.packages('httr')")
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')")

  stac_url <- "https://planetarycomputer.microsoft.com/api/stac/v1/search"

  body <- list(
    collections = list("3dep-lidar-copc"),
    bbox         = bbox,
    limit        = limit
  )

  resp <- httr::POST(
    stac_url,
    httr::content_type_json(),
    body = jsonlite::toJSON(body, auto_unbox = TRUE)
  )
  httr::stop_for_status(resp)

  result <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"),
                               simplifyVector = FALSE)

  cat("STAC search returned", length(result$features), "item(s)\n")
  result$features
}

#' Get a signed COPC URL for a STAC item (appends SAS token)
#'
#' @param item A single STAC item (list)
#' @return Character string: the signed URL
sign_planetary_url <- function(item) {
  if (!requireNamespace("httr", quietly = TRUE))
    stop("Package 'httr' is required.")
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required.")

  href <- item$assets$data$href
  if (is.null(href))
    stop("No 'data' asset found in this STAC item.")

  sas_url <- "https://planetarycomputer.microsoft.com/api/sas/v1/token/3dep-lidar-copc"
  resp <- httr::GET(sas_url)
  httr::stop_for_status(resp)
  token_data <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))

  sep <- if (grepl("\\?", href)) "&" else "?"
  paste0(href, sep, token_data$token)
}

#' Extract the CRS from a STAC item's projection metadata
#'
#' The 3DEP COPC items have proj:epsg = NULL but include a full
#' proj:projjson object. This helper extracts the CRS reliably.
#'
#' @param item A single STAC item (list)
#' @return An sf CRS object
get_item_crs <- function(item) {
  props <- item$properties

  # 1. Try proj:epsg (often NULL for 3DEP, but check anyway)
  epsg <- props[["proj:epsg"]]
  if (!is.null(epsg) && !is.na(epsg)) {
    crs <- sf::st_crs(as.integer(epsg))
    if (!is.na(crs)) return(crs)
  }

  # 2. Try proj:projjson — the most reliable source for 3DEP COPC
  projjson <- props[["proj:projjson"]]
  if (!is.null(projjson)) {
    if (!is.null(projjson$components)) {
      horiz <- projjson$components[[1]]
      horiz_epsg <- horiz$id$code
      if (!is.null(horiz_epsg)) {
        crs <- sf::st_crs(as.integer(horiz_epsg))
        if (!is.na(crs)) {
          cat("  CRS from STAC proj:projjson -> EPSG:", horiz_epsg,
              "(", horiz$name, ")\n")
          return(crs)
        }
      }
    }
    json_str <- jsonlite::toJSON(projjson, auto_unbox = TRUE)
    crs <- tryCatch(sf::st_crs(json_str), error = function(e) sf::st_crs(NA))
    if (!is.na(crs)) return(crs)
  }

  # 3. Try proj:wkt2
  wkt2 <- props[["proj:wkt2"]]
  if (!is.null(wkt2) && nzchar(wkt2)) {
    crs <- tryCatch(sf::st_crs(wkt2), error = function(e) sf::st_crs(NA))
    if (!is.na(crs)) return(crs)
  }

  stop("Could not determine CRS from STAC item. ",
       "Check item$properties for proj:epsg, proj:projjson, or proj:wkt2.")
}

#' Extract the USGS project name from a STAC item ID
#'
#' @param item A single STAC item (list)
#' @return Character: the 3dep:usgs_id or a cleaned-up project name
get_item_project <- function(item) {
  # The 3dep:usgs_id groups tiles by acquisition/project
  usgs_id <- item$properties[["3dep:usgs_id"]]
  if (!is.null(usgs_id)) return(usgs_id)
  # Fallback: extract from the item ID (e.g., "USGS_LPC_AZ_Coconino_2019_...")
  sub("^USGS_LPC_", "", sub("_[^_]+$", "", item$id))
}

#' Summarise STAC items into a table of unique acquisitions/projects
#'
#' @param items List of STAC items
#' @return data.frame with project name, tile count, date range, point count
summarise_acquisitions <- function(items) {
  info <- lapply(items, function(it) {
    data.frame(
      project    = get_item_project(it),
      tile_id    = it$id,
      start_date = it$properties[["start_datetime"]] %||% NA_character_,
      end_date   = it$properties[["end_datetime"]]   %||% NA_character_,
      n_points   = it$properties[["pc:count"]]        %||% NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, info)
}

# Null-coalescing operator (if not already defined)
`%||%` <- function(a, b) if (is.null(a)) b else a


# =============================================================================
# 1. Define the AOI as an sf polygon
# =============================================================================
# Create a polygon in WGS 84 (lon/lat). Change these coordinates to your
# area of interest. This is a roughly 700 m × 230 m area near Flagstaff, AZ.

if (!requireNamespace("sf", quietly = TRUE))
  stop("Package 'sf' is required. Install with: install.packages('sf')")

library(sf)

aoi_coords <- matrix(c(
  -111.6720, 35.1700,
  -111.6640, 35.1700,
  -111.6640, 35.1720,
  -111.6720, 35.1720,
  -111.6720, 35.1700    # close the ring
), ncol = 2, byrow = TRUE)

aoi_sf <- st_sf(
  id       = 1,
  geometry = st_sfc(st_polygon(list(aoi_coords)), crs = 4326)
)

cat("Area of interest (WGS 84):\n")
print(st_bbox(aoi_sf))
cat("\n")


# =============================================================================
# 2. Search the Planetary Computer STAC API
# =============================================================================
aoi_bbox <- as.numeric(st_bbox(aoi_sf))
items    <- search_3dep_copc(aoi_bbox, limit = 50)

if (length(items) == 0)
  stop("No 3DEP COPC tiles found for this AOI. ",
       "Try a different location or a larger polygon.")


# =============================================================================
# 3. Show available acquisitions and let the user choose
# =============================================================================
acq_df <- summarise_acquisitions(items)

# Group by project (acquisition)
projects <- unique(acq_df$project)

cat("\n──────────────────────────────────────────────────────────────────\n")
cat("  Available LiDAR acquisitions overlapping the AOI\n")
cat("──────────────────────────────────────────────────────────────────\n\n")

for (i in seq_along(projects)) {
  proj  <- projects[i]
  tiles <- acq_df[acq_df$project == proj, ]
  dates <- unique(na.omit(c(tiles$start_date, tiles$end_date)))
  date_range <- if (length(dates) > 0) {
    paste(substr(min(dates), 1, 10), "to", substr(max(dates), 1, 10))
  } else {
    "unknown"
  }
  total_pts <- sum(tiles$n_points, na.rm = TRUE)
  cat(sprintf("  [%d] %s\n", i, proj))
  cat(sprintf("      Tiles : %d\n", nrow(tiles)))
  cat(sprintf("      Dates : %s\n", date_range))
  cat(sprintf("      Points: %s (across all tiles)\n\n",
              format(total_pts, big.mark = ",")))
}

# ── Interactive selection ─────────────────────────────────────────────────────
if (length(projects) == 1L) {
  choice <- 1L
  cat("Only one acquisition available — selecting it automatically.\n")
} else if (interactive()) {
  cat("Enter the number of the acquisition to load (1-", length(projects), "): ", sep = "")
  choice <- as.integer(readline())
  if (is.na(choice) || choice < 1L || choice > length(projects))
    stop("Invalid selection.")
} else {
  choice <- 1L
  cat("Non-interactive session — defaulting to [1].\n")
}

selected_project <- projects[choice]
selected_tiles   <- acq_df$tile_id[acq_df$project == selected_project]
selected_items   <- Filter(function(it) it$id %in% selected_tiles, items)

cat("\nSelected acquisition:", selected_project, "\n")
cat("Tiles to process   :", length(selected_items), "\n")


# =============================================================================
# 4. Get the CRS and reproject the AOI
# =============================================================================
# All tiles in one acquisition share the same CRS.
tile_crs <- get_item_crs(selected_items[[1]])
aoi_proj <- st_transform(aoi_sf, tile_crs)

cat("\nAOI reprojected to tile CRS:\n")
print(st_bbox(aoi_proj))

aoi_bbox_proj <- as.numeric(st_bbox(aoi_proj))


# =============================================================================
# 5. Read COPC data from each tile and merge
# =============================================================================
# For each tile that overlaps the AOI, read up to 1,000,000 points within
# the AOI bounding box, then combine them.

all_data    <- list()
all_headers <- list()

for (i in seq_along(selected_items)) {
  it  <- selected_items[[i]]
  cat(sprintf("\n── Tile %d/%d: %s ──\n", i, length(selected_items), it$id))

  url <- sign_planetary_url(it)
  hdr <- read_copc_header(url)
  cat(sprintf("  Points in tile: %s\n",
              format(hdr[["Number of point records"]], big.mark = ",")))

  # Check that the AOI bbox actually overlaps the tile extent
  tile_xmin <- hdr[["Min X"]]; tile_xmax <- hdr[["Max X"]]
  tile_ymin <- hdr[["Min Y"]]; tile_ymax <- hdr[["Max Y"]]

  if (aoi_bbox_proj[1] > tile_xmax || aoi_bbox_proj[3] < tile_xmin ||
      aoi_bbox_proj[2] > tile_ymax || aoi_bbox_proj[4] < tile_ymin) {
    cat("  Skipping — AOI does not overlap this tile's extent.\n")
    next
  }

  res <- read_copc(
    url,
    bbox       = aoi_bbox_proj,
    select     = "xyzicrnR",
    max_points = 1000000,
    progress   = TRUE
  )

  if (nrow(res$data) > 0L) {
    cat(sprintf("  Points read: %s\n", format(nrow(res$data), big.mark = ",")))
    all_data[[length(all_data) + 1L]]    <- res$data
    all_headers[[length(all_headers) + 1L]] <- res$header
  } else {
    cat("  No points returned for this tile's overlap region.\n")
  }
}

if (length(all_data) == 0L)
  stop("No points were read from any tile. Check the AOI and CRS.")

# Merge data from all tiles
merged_data <- data.table::rbindlist(all_data)

# Use the first tile's header as a base, update extents
result <- list(data = merged_data, header = all_headers[[1]])

cat(sprintf("\nTotal points from %d tile(s): %s\n",
            length(all_data), format(nrow(merged_data), big.mark = ",")))


# =============================================================================
# 6. Clip points to the polygon
# =============================================================================
pts_sf <- st_as_sf(merged_data, coords = c("X", "Y"), crs = tile_crs)
inside <- st_intersects(pts_sf, aoi_proj, sparse = FALSE)[, 1]

n_bbox   <- nrow(merged_data)
n_inside <- sum(inside)
cat(sprintf("\nPoints inside bounding box : %s\n",
            format(n_bbox, big.mark = ",")))
cat(sprintf("Points inside polygon      : %s\n",
            format(n_inside, big.mark = ",")))
cat(sprintf("Points removed by clipping : %s\n",
            format(n_bbox - n_inside, big.mark = ",")))

result$data <- merged_data[inside, ]

if (nrow(result$data) == 0L)
  stop("All points were outside the polygon after clipping.")

# Update header extents
result$header[["Min X"]] <- min(result$data$X)
result$header[["Max X"]] <- max(result$data$X)
result$header[["Min Y"]] <- min(result$data$Y)
result$header[["Max Y"]] <- max(result$data$Y)
result$header[["Min Z"]] <- min(result$data$Z)
result$header[["Max Z"]] <- max(result$data$Z)
result$header[["Number of point records"]] <- nrow(result$data)


# =============================================================================
# 7. Convert to lidR LAS object and visualize
# =============================================================================
if (requireNamespace("lidR", quietly = TRUE)) {
  library(lidR)

  las <- as_las(result)
  cat("\nLAS object after polygon clipping:\n")
  print(las)

  cat("\nClassification counts:\n")
  print(table(las$Classification))

  # Plot las
  view(las)
}


# =============================================================================
# 8. Summary
# =============================================================================
cat("\n=== Workflow complete ===\n")
cat("  Source      : Microsoft Planetary Computer – USGS 3DEP COPC\n")
cat("  Acquisition:", selected_project, "\n")
cat("  Tiles used  :", length(all_data), "\n")
cat("  Points kept :", format(nrow(result$data), big.mark = ","), "\n")
cat("  Max points  : 1,000,000 per tile\n")
