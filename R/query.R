# -- Query helpers: polygon AOI, buffered corridors, thin/sample -----------
#
# These are convenience wrappers around read_copc() that handle sf geometry
# pre-processing (polygon/line buffering, CRS transforms, etc.).
# ---------------------------------------------------------------------------

# Suppress R CMD check NOTEs for data.table NSE variables
utils::globalVariables(c(".grp", ".idx", ".d"))

# read_copc_polygon() -- back-compat wrapper; prefer read_copc(aoi = polygon).
read_copc_polygon <- function(path_or_url,
                              polygon,
                              zrange     = NULL,
                              select     = NULL,
                              filter     = NULL,
                              max_depth  = -1L,
                              resolution = NULL,
                              max_points = Inf,
                              progress   = TRUE) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for read_copc_polygon().", call. = FALSE)

  read_copc(
    path_or_url = path_or_url,
    aoi         = polygon,
    zrange      = zrange,
    select      = select,
    filter      = filter,
    max_depth   = max_depth,
    resolution  = resolution,
    max_points  = max_points,
    progress    = progress
  )
}


# read_copc_corridor() -- back-compat wrapper; prefer read_copc(aoi = line, buffer = width).
read_copc_corridor <- function(path_or_url,
                               line,
                               width,
                               zrange     = NULL,
                               select     = NULL,
                               filter     = NULL,
                               max_depth  = -1L,
                               resolution = NULL,
                               max_points = Inf,
                               progress   = TRUE) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for read_copc_corridor().", call. = FALSE)

  stopifnot(is.numeric(width), length(width) == 1L, width > 0)

  corridor <- sf::st_buffer(line, dist = width)
  read_copc(
    path_or_url = path_or_url,
    aoi         = corridor,
    zrange      = zrange,
    select      = select,
    filter      = filter,
    max_depth   = max_depth,
    resolution  = resolution,
    max_points  = max_points,
    progress    = progress
  )
}


#' Thin / sample a COPC read to N points
#'
#' Reads from COPC and then thins the result.  Supports two modes:
#' \enumerate{
#'   \item \strong{Random sampling} (default): randomly keeps at most
#'     \code{n} points.
#'   \item \strong{Voxel-based thinning}: if \code{voxel_size} is given,
#'     overlays a 3-D voxel grid and keeps one point per occupied cell.
#'     This is equivalent to PDAL's \code{filters.voxeldownsize}.
#' }
#' Combine with \code{max_depth} or \code{resolution} for a fast
#' LOD preview.
#'
#' @inheritParams read_copc
#' @param n Integer.  Maximum number of points to keep after sampling.
#'   Used when \code{voxel_size} is \code{NULL}.
#' @param voxel_size Numeric scalar.  Side length of cubic voxel cells
#'   (in CRS units).  When supplied, voxel-based thinning is performed
#'   instead of random sampling.  \code{n} is then used as an optional
#'   secondary cap (thin to voxels first, then randomly thin to \code{n}
#'   if the result is still larger).
#' @param mode Character. Voxel thinning mode when \code{voxel_size} is
#'   given.
#'   \describe{
#'     \item{\code{"first"}}{Keep the first point encountered in each
#'       voxel (fastest).}
#'     \item{\code{"random"}}{Keep one random point per voxel (default).}
#'     \item{\code{"center"}}{Keep the point nearest the voxel center
#'       (most spatially regular).}
#'   }
#'   Ignored when \code{voxel_size} is \code{NULL}.
#'
#' @return Same as \code{\link{read_copc}} but with thinned
#'   \code{$data}.
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#'
#' # Random sample of up to 500 points
#' result <- read_copc_sample(f, n = 500, progress = FALSE)
#' nrow(result$data)
#'
#' # Voxel-based thinning (1-meter cells, keep nearest to center)
#' result_vox <- read_copc_sample(f, voxel_size = 1.0,
#'                                mode = "center",
#'                                progress = FALSE)
#' nrow(result_vox$data)
#' }
#'
#' @export
read_copc_sample <- function(path_or_url,
                             n          = Inf,
                             voxel_size = NULL,
                             mode       = c("random", "first", "center"),
                             bbox       = NULL,
                             aoi        = NULL,
                             zrange     = NULL,
                             select     = NULL,
                             filter     = NULL,
                             max_depth  = -1L,
                             resolution = NULL,
                             progress   = TRUE) {
  mode <- match.arg(mode)
  stopifnot(is.numeric(n), length(n) == 1L, n > 0)

  result <- read_copc(
    path_or_url = path_or_url,
    bbox        = bbox,
    aoi         = aoi,
    zrange      = zrange,
    select      = select,
    filter      = filter,
    max_depth   = max_depth,
    resolution  = resolution,
    max_points  = Inf,
    progress    = progress
  )

  nr <- nrow(result$data)
  if (nr == 0L) return(result)

  # -- Voxel-based thinning -------------------------------------------
  if (!is.null(voxel_size)) {
    stopifnot(is.numeric(voxel_size), length(voxel_size) == 1L, voxel_size > 0)

    dt <- result$data

    # Compute integer voxel indices
    vx <- as.integer(floor(dt$X / voxel_size))
    vy <- as.integer(floor(dt$Y / voxel_size))
    vz <- as.integer(floor(dt$Z / voxel_size))

    if (mode == "first") {
      # Keep first occurrence in each voxel: use duplicated on the key
      key <- paste(vx, vy, vz, sep = "_")
      keep_idx <- which(!duplicated(key))

    } else if (mode == "random") {
      # Shuffle rows, then keep first per voxel (= random per voxel)
      perm <- sample.int(nr)
      key <- paste(vx[perm], vy[perm], vz[perm], sep = "_")
      keep_perm <- perm[!duplicated(key)]
      keep_idx <- sort(keep_perm)

    } else if (mode == "center") {
      # Keep the point nearest to the voxel center
      cx <- (vx + 0.5) * voxel_size
      cy <- (vy + 0.5) * voxel_size
      cz <- (vz + 0.5) * voxel_size
      dist_sq <- (dt$X - cx)^2 + (dt$Y - cy)^2 + (dt$Z - cz)^2

      grp <- paste(vx, vy, vz, sep = "_")
      dt_tmp <- data.table::data.table(
        .idx = seq_len(nr), .grp = grp, .d = dist_sq
      )
      keep_idx <- dt_tmp[, list(.idx = .idx[which.min(.d)]), by = .grp]$.idx
      keep_idx <- sort(keep_idx)
    }

    result$data <- dt[keep_idx, ]

    if (progress) {
      message(sprintf("Voxel thin (%.2f, %s): %d -> %d points",
                      voxel_size, mode, nr, length(keep_idx)))
    }

    nr <- nrow(result$data)
  }

  # -- Secondary random cap to n --------------------------------------
  if (nr > n) {
    idx <- sort(sample.int(nr, n))
    result$data <- result$data[idx, ]
  }

  result
}
