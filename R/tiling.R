# -- Deterministic tiling: read_copc_tiles() ------------------------------
#
# Splits an extent into a regular grid and reads one tile at a time.
# Returns a list of read_copc() results, one per tile.
# ---------------------------------------------------------------------------

#' Read COPC data as a grid of tiles
#'
#' Divides an area of interest into a regular grid of rectangular tiles
#' and reads each tile independently.
#'
#' @param path_or_url Character. Path or URL to a \code{.copc.laz} file.
#' @param extent Numeric vector \code{c(xmin, ymin, xmax, ymax)}, an
#'   \code{sf::st_bbox()} object, or an \code{sf} geometry.
#'   Default \code{NULL} uses the full file extent.
#' @param tile_size Numeric scalar.  Side length of each square tile
#'   (in CRS units).
#' @param overlap Numeric scalar.  Buffer distance added around each
#'   tile to ensure seamless coverage (default 0).
#' @param select Character.  Column selection string.
#' @param filter Character.  Attribute filter string.
#' @param max_depth Integer.  Maximum octree depth (\code{-1} = all).
#' @param resolution Numeric.  Target point spacing (in CRS units);
#'   overrides \code{max_depth} if both are given.
#' @param max_points Numeric.  Per-tile max point cap.
#' @param progress Logical.
#' @param as_las Logical.  If \code{TRUE} and \pkg{lidR} is available,
#'   each tile is returned as a \code{LAS} object instead of a raw list.
#'
#' @return A list of tile results.  Each element is the output of
#'   \code{read_copc()} (or a \code{LAS} object if \code{as_las = TRUE}).
#'   The list carries an attribute \code{"tile_bboxes"} with the bbox
#'   of each tile as a 4-column matrix.
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#' tiles <- read_copc_tiles(f, tile_size = 200, progress = FALSE)
#' length(tiles)  # number of tiles
#' nrow(tiles[[1]]$data)  # points in first tile
#' }
#'
#' @export
read_copc_tiles <- function(path_or_url,
                            extent     = NULL,
                            tile_size,
                            overlap    = 0,
                            select     = NULL,
                            filter     = NULL,
                            max_depth  = -1L,
                            resolution = NULL,
                            max_points = Inf,
                            progress   = TRUE,
                            as_las     = FALSE) {
  stopifnot(is.numeric(tile_size), length(tile_size) == 1L, tile_size > 0)
  stopifnot(is.numeric(overlap), length(overlap) == 1L, overlap >= 0)

  # -- Determine extent ----------------------------------------------
  if (is.null(extent)) {
    hdr <- read_copc_header(path_or_url)
    extent <- c(
      xmin = hdr[["Min X"]], ymin = hdr[["Min Y"]],
      xmax = hdr[["Max X"]], ymax = hdr[["Max Y"]]
    )
  } else if (inherits(extent, "bbox")) {
    extent <- as.numeric(extent)[1:4]
  } else if (inherits(extent, "sf") || inherits(extent, "sfc")) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' required.", call. = FALSE)
    extent <- as.numeric(sf::st_bbox(extent))[1:4]
  }
  stopifnot(is.numeric(extent), length(extent) == 4L)
  names(extent) <- c("xmin", "ymin", "xmax", "ymax")

  # -- Build tile grid -----------------------------------------------
  xs <- seq(extent["xmin"], extent["xmax"], by = tile_size)
  ys <- seq(extent["ymin"], extent["ymax"], by = tile_size)

  tiles <- expand.grid(col = seq_along(xs), row = seq_along(ys))
  tile_bboxes <- matrix(NA_real_, nrow = nrow(tiles), ncol = 4,
                        dimnames = list(NULL, c("xmin","ymin","xmax","ymax")))
  for (i in seq_len(nrow(tiles))) {
    ci <- tiles$col[i]
    ri <- tiles$row[i]
    tile_bboxes[i, ] <- c(
      xs[ci]     - overlap,
      ys[ri]     - overlap,
      min(xs[ci] + tile_size + overlap, extent["xmax"] + overlap),
      min(ys[ri] + tile_size + overlap, extent["ymax"] + overlap)
    )
  }

  # -- Read each tile ------------------------------------------------
  n_tiles <- nrow(tile_bboxes)
  results <- vector("list", n_tiles)

  for (i in seq_len(n_tiles)) {
    if (progress) {
      message(sprintf("Tile %d/%d", i, n_tiles))
    }
    bb <- tile_bboxes[i, ]
    result <- read_copc(
      path_or_url = path_or_url,
      bbox        = bb,
      select      = select,
      filter      = filter,
      max_depth   = max_depth,
      resolution  = resolution,
      max_points  = max_points,
      progress    = FALSE
    )
    if (as_las) {
      result <- as_las(result)
    }
    results[[i]] <- result
  }

  attr(results, "tile_bboxes") <- tile_bboxes
  results
}
