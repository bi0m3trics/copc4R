# -- Octree-based point count estimation & adaptive splitting ---------------
#
# Implements the pattern Howard Butler describes: "query the octree to
# get estimated maximum point counts for a query window, and then splitting
# that down to some threshold."
#
# This enables adaptive chunking where work units are balanced by estimated
# point count rather than fixed spatial area.
# ---------------------------------------------------------------------------

#' Estimate point counts for spatial windows via the COPC octree
#'
#' Queries the COPC octree hierarchy to estimate point counts for one or
#' more bounding boxes **without decompressing any point data**.  This is
#' very fast even for remote files because it only reads the hierarchy
#' pages.
#'
#' Use this to pre-flight a processing job: estimate the cost of each tile,
#' identify dense areas that need smaller work units, or filter out empty
#' regions before committing to expensive reads.
#'
#' @param path_or_url Character. Path or URL to a `.copc.laz` file.
#' @param bboxes A list of numeric(4) vectors `c(xmin, ymin, xmax, ymax)`,
#'   or a single numeric(4) vector for one window. If `NULL`, uses the
#'   full file extent.
#' @param max_depth Integer. Maximum octree depth to descend.
#'   `-1L` (default) uses all levels for the most accurate estimate.
#'
#' @return A data.frame with one row per bbox:
#'   \describe{
#'     \item{xmin, ymin, xmax, ymax}{The query window.}
#'     \item{estimated_points}{Estimated point count from the hierarchy.}
#'     \item{num_nodes}{Number of octree nodes that overlap the window.}
#'     \item{area}{Area of the query window in CRS units squared.}
#'     \item{density}{Estimated points per unit area.}
#'   }
#'
#' @examples
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#'
#' # Single window
#' copc_estimate_points(f)
#'
#' # Multiple windows
#' hdr <- read_copc_header(f)
#' xmid <- (hdr[["Min X"]] + hdr[["Max X"]]) / 2
#' ymid <- (hdr[["Min Y"]] + hdr[["Max Y"]]) / 2
#' bboxes <- list(
#'   c(hdr[["Min X"]], hdr[["Min Y"]], xmid, ymid),
#'   c(xmid, ymid, hdr[["Max X"]], hdr[["Max Y"]])
#' )
#' copc_estimate_points(f, bboxes)
#'
#' @export
copc_estimate_points <- function(path_or_url,
                                 bboxes    = NULL,
                                 max_depth = -1L) {

  path <- .resolve_path(path_or_url, progress = FALSE)

  if (is.null(bboxes)) {
    hdr <- read_copc_header(path_or_url)
    bboxes <- list(c(hdr[["Min X"]], hdr[["Min Y"]],
                     hdr[["Max X"]], hdr[["Max Y"]]))
  }

  # Accept a single numeric(4) as a convenience

  if (is.numeric(bboxes) && length(bboxes) == 4L)
    bboxes <- list(bboxes)

  stopifnot(is.list(bboxes), length(bboxes) >= 1L)
  max_depth <- as.integer(max_depth)

  results <- vector("list", length(bboxes))
  for (i in seq_along(bboxes)) {
    bb <- as.numeric(bboxes[[i]])
    stopifnot(length(bb) == 4L)

    node_info <- cpp_count_nodes(path, bb)
    area <- (bb[3] - bb[1]) * (bb[4] - bb[2])

    results[[i]] <- data.frame(
      xmin             = bb[1],
      ymin             = bb[2],
      xmax             = bb[3],
      ymax             = bb[4],
      estimated_points = node_info$total_points,
      num_nodes        = node_info$num_nodes,
      area             = area,
      density          = if (area > 0) node_info$total_points / area else Inf,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Adaptively split a bbox until estimated points fall below a threshold
#'
#' Recursively subdivides a spatial window using the COPC octree to
#' estimate point counts.
#' Splitting continues (quadtree-style) until every sub-window
#' has an estimated point count at or below `max_points`, or the minimum
#' tile side length (`min_size`) is reached.
#'
#' This implements the pattern described by Howard Butler: *"query the
#' octree to get estimated maximum point counts for a query window, and
#' then splitting that down to some threshold"* -- useful for balancing
#' "work units" of some size/memory/complexity.
#'
#' @param path_or_url Character. Path or URL to a `.copc.laz` file.
#' @param bbox Numeric(4) `c(xmin, ymin, xmax, ymax)`. The initial window
#'   to split. `NULL` (default) uses the full file extent.
#' @param max_points Numeric. Target maximum estimated point count per
#'   sub-window. Default `1e6` (1 million).
#' @param min_size Numeric. Minimum side length of a sub-window (in CRS
#'   units). Splitting stops even if the count threshold is exceeded.
#'   Default `50`.
#' @param collar Numeric. Buffer distance to add around each resulting
#'   tile for processing overlap (trimmed off after processing).
#'   Default `0` (no collar).
#'
#' @return A data.frame with one row per resulting tile:
#'   \describe{
#'     \item{xmin, ymin, xmax, ymax}{Core tile extent (before collar).}
#'     \item{read_xmin, read_ymin, read_xmax, read_ymax}{Extent
#'       including the collar -- use this for the actual read.}
#'     \item{estimated_points}{Estimate from the octree.}
#'     \item{depth}{Recursion depth (0 = original window).}
#'   }
#'
#' @seealso [copc_estimate_points()], [copc_apply()]
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#'
#' # Split file extent into tiles of <= 5000 points each
#' tiles <- copc_adaptive_split(f, max_points = 5000)
#' print(tiles)
#'
#' # With 10m collar for spatial algorithms
#' tiles <- copc_adaptive_split(f, max_points = 5000, collar = 10)
#' }
#'
#' @export
copc_adaptive_split <- function(path_or_url,
                                bbox       = NULL,
                                max_points = 1e6,
                                min_size   = 50,
                                collar     = 0) {

  path <- .resolve_path(path_or_url, progress = FALSE)

  if (is.null(bbox)) {
    hdr <- read_copc_header(path_or_url)
    bbox <- c(hdr[["Min X"]], hdr[["Min Y"]],
              hdr[["Max X"]], hdr[["Max Y"]])
  }

  stopifnot(is.numeric(bbox), length(bbox) == 4L)
  stopifnot(is.numeric(max_points), length(max_points) == 1L, max_points > 0)
  stopifnot(is.numeric(min_size), length(min_size) == 1L, min_size > 0)
  stopifnot(is.numeric(collar), length(collar) == 1L, collar >= 0)

  # Recursive quadtree split
  tiles <- list()

  .split <- function(bb, depth) {
    node_info <- cpp_count_nodes(path, as.numeric(bb))
    est <- node_info$total_points

    width  <- bb[3] - bb[1]
    height <- bb[4] - bb[2]

    # Stop splitting if: under threshold, or tile is too small
    if (est <= max_points || width <= min_size || height <= min_size) {
      tiles[[length(tiles) + 1L]] <<- data.frame(
        xmin = bb[1], ymin = bb[2], xmax = bb[3], ymax = bb[4],
        read_xmin = bb[1] - collar, read_ymin = bb[2] - collar,
        read_xmax = bb[3] + collar, read_ymax = bb[4] + collar,
        estimated_points = est,
        depth = depth,
        stringsAsFactors = FALSE
      )
      return()
    }

    # Quadrant split
    xmid <- (bb[1] + bb[3]) / 2
    ymid <- (bb[2] + bb[4]) / 2

    .split(c(bb[1], bb[2], xmid, ymid),  depth + 1L)  # SW
    .split(c(xmid,  bb[2], bb[3], ymid),  depth + 1L)  # SE
    .split(c(bb[1], ymid,  xmid,  bb[4]), depth + 1L)  # NW
    .split(c(xmid,  ymid,  bb[3], bb[4]), depth + 1L)  # NE
  }

  .split(bbox, 0L)

  if (length(tiles) == 0L)
    return(data.frame(xmin = numeric(0), ymin = numeric(0),
                      xmax = numeric(0), ymax = numeric(0),
                      read_xmin = numeric(0), read_ymin = numeric(0),
                      read_xmax = numeric(0), read_ymax = numeric(0),
                      estimated_points = numeric(0), depth = integer(0)))

  do.call(rbind, tiles)
}
