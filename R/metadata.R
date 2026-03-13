# -- Metadata convenience functions -------------------------------------------
#
# "Header-first" lightweight calls that don't pull point data unless needed.
# ---------------------------------------------------------------------------

#' Get COPC octree metadata
#'
#' Returns the COPC Info VLR contents: octree center, halfsize, spacing,
#' hierarchy root location, and GPS-time range.  Only the file header is
#' read -- no point data is fetched.
#'
#' @param path_or_url Character. Path or URL to a `.copc.laz` file.
#'
#' @return A named list with elements \code{center_x}, \code{center_y},
#'   \code{center_z}, \code{halfsize}, \code{spacing},
#'   \code{root_hier_offset}, \code{root_hier_size},
#'   \code{gpstime_minimum}, \code{gpstime_maximum}, \code{crs}
#'   (CRS WKT string if available), and \code{point_count}.
#'
#' @examples
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#' info <- copc_info(f)
#' info$spacing
#' info$point_count
#' info$point_format
#'
#' @export
copc_info <- function(path_or_url) {
  hdr <- read_copc_header(path_or_url)
  copc <- hdr[["COPC Info"]]

  # Extract CRS from WKT VLR if present
  crs_wkt <- NULL
  vlrs <- hdr[["Variable Length Records"]]
  for (v in vlrs) {
    if (!is.null(v[["WKT OGC COORDINATE SYSTEM"]])) {
      crs_wkt <- v[["WKT OGC COORDINATE SYSTEM"]]
      break
    }
  }

  list(
    center_x          = copc$center_x,
    center_y          = copc$center_y,
    center_z          = copc$center_z,
    halfsize          = copc$halfsize,
    spacing           = copc$spacing,
    root_hier_offset  = copc$root_hier_offset,
    root_hier_size    = copc$root_hier_size,
    gpstime_minimum   = copc$gpstime_minimum,
    gpstime_maximum   = copc$gpstime_maximum,
    point_count       = hdr[["Number of point records"]],
    point_format      = hdr[["Point Data Format ID"]],
    crs_wkt           = crs_wkt
  )
}


#' Get COPC file spatial bounds
#'
#' Returns the 3-D bounding box from the LAS header.  No points are read.
#'
#' @param path_or_url Character. Path or URL to a `.copc.laz` file.
#' @param as_sf Logical.  If \code{TRUE} and \pkg{sf} is installed,
#'   return an \code{sf} bounding-box polygon with CRS.  Default \code{FALSE}.
#'
#' @return If \code{as_sf = FALSE}, a named numeric vector
#'   \code{c(xmin, ymin, zmin, xmax, ymax, zmax)}.
#'   If \code{as_sf = TRUE}, an \code{sf::st_bbox()} object.
#'
#' @examples
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#' copc_bounds(f)
#'
#' \donttest{
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   copc_bounds(f, as_sf = TRUE)
#' }
#' }
#'
#' @export
copc_bounds <- function(path_or_url, as_sf = FALSE) {
  hdr <- read_copc_header(path_or_url)

  bounds <- c(
    xmin = hdr[["Min X"]], ymin = hdr[["Min Y"]], zmin = hdr[["Min Z"]],
    xmax = hdr[["Max X"]], ymax = hdr[["Max Y"]], zmax = hdr[["Max Z"]]
  )

  if (as_sf) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required for as_sf = TRUE.", call. = FALSE)

    bb <- sf::st_bbox(c(
      xmin = bounds[["xmin"]], ymin = bounds[["ymin"]],
      xmax = bounds[["xmax"]], ymax = bounds[["ymax"]]
    ))

    # Try to attach CRS
    vlrs <- hdr[["Variable Length Records"]]
    for (v in vlrs) {
      if (!is.null(v[["WKT OGC COORDINATE SYSTEM"]])) {
        attr(bb, "crs") <- sf::st_crs(v[["WKT OGC COORDINATE SYSTEM"]])
        break
      }
    }
    return(bb)
  }

  bounds
}


#' Estimate point density for an AOI
#'
#' Quickly estimates the point density (points per square unit) without
#' reading actual point data.  Uses the COPC hierarchy to count points
#' in nodes overlapping the query, then divides by the AOI area.
#'
#' @param path_or_url Character. Path or URL to a `.copc.laz` file.
#' @param bbox Numeric vector \code{c(xmin, ymin, xmax, ymax)}, or
#'   \code{NULL} to use the full file extent.
#' @param aoi An \code{sf}/\code{sfc} polygon object.  Used instead of
#'   \code{bbox} if provided (its bounding box is used for the query).
#'
#' @return A named list with \code{estimated_points}, \code{area},
#'   \code{density} (points per unit area), and \code{nodes_checked}.
#'
#' @examples
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#' dens <- copc_density(f)
#' dens$density
#'
#' @export
copc_density <- function(path_or_url, bbox = NULL, aoi = NULL) {
  if (!is.null(aoi)) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required when using 'aoi' parameter.", call. = FALSE)
    bbox <- as.numeric(sf::st_bbox(aoi))
  }

  hdr <- read_copc_header(path_or_url)

  if (is.null(bbox)) {
    bbox <- c(hdr[["Min X"]], hdr[["Min Y"]], hdr[["Max X"]], hdr[["Max Y"]])
  }
  stopifnot(is.numeric(bbox), length(bbox) == 4L)

  area <- (bbox[3] - bbox[1]) * (bbox[4] - bbox[2])

  # Use C++ node selection to count points without decompressing
  path <- .resolve_path(path_or_url, progress = FALSE)
  node_info <- cpp_count_nodes(path, as.numeric(bbox))

  list(
    estimated_points = node_info$total_points,
    area             = area,
    density          = if (area > 0) node_info$total_points / area else Inf,
    nodes_checked    = node_info$num_nodes
  )
}
