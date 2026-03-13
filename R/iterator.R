# -- Streaming / iterator API: read_copc_iter() ---------------------------
#
# Closure-based iterator over read_copc() results.
# Not exported -- reads all points into memory on first $yield().
# ---------------------------------------------------------------------------

# read_copc_iter -- back-compat; loads all points into memory, then pages them.
read_copc_iter <- function(path_or_url,
                           bbox        = NULL,
                           aoi         = NULL,
                           zrange      = NULL,
                           select      = NULL,
                           filter      = NULL,
                           max_depth   = -1L,
                           resolution  = NULL,
                           batch_nodes = 10L,
                           progress    = TRUE) {

  stopifnot(is.character(path_or_url), length(path_or_url) == 1L)
  batch_nodes <- as.integer(batch_nodes)

  # -- AOI -> bbox ------------------------------------------------------
  clip_geom <- NULL
  if (!is.null(aoi)) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required for the 'aoi' parameter.", call. = FALSE)
    aoi_bbox <- as.numeric(sf::st_bbox(aoi))
    if (is.null(bbox)) bbox <- aoi_bbox
    clip_geom <- sf::st_geometry(aoi)
    if (!inherits(clip_geom, "sfc"))
      clip_geom <- sf::st_sfc(clip_geom, crs = sf::st_crs(aoi))
  }

  # -- Resolve path ----------------------------------------------------
  path <- .resolve_path(path_or_url, progress = progress)
  .check_copc_ext(path)

  if (is.null(bbox))   bbox   <- numeric(0)
  if (is.null(zrange)) zrange <- numeric(0)
  if (is.null(select)) select <- "*"
  if (is.null(filter)) filter <- ""

  # -- Read header once ------------------------------------------------
  hdr_list <- cpp_read_copc_header(path)

  # -- State: batch index ----------------------------------------------
  # We use read_copc() for each batch, controlling chunks via
  # max_points (approximate).  For a true chunk-level iterator we'd
  # need C++ support for returning partial results; as a practical
  # first implementation, we use max_points per yield.
  #
  # Strategy: read all points with a batch-sized cap, then use an

  # offset-based approach via the data already read.

  # Read everything once (streamed, not full download) and store
  all_result <- NULL
  current_offset <- 0L
  batch_size <- NULL  # determined from first full read

  .full_read <- function() {
    if (!is.null(all_result)) return()
    all_result <<- read_copc(
      path_or_url = path_or_url,
      bbox        = if (length(bbox) == 4L) bbox else NULL,
      aoi         = if (!is.null(clip_geom)) aoi else NULL,
      zrange      = if (length(zrange) == 2L) zrange else NULL,
      select      = select,
      filter      = filter,
      max_depth   = max_depth,
      resolution  = resolution,
      progress    = progress
    )
    total_rows <- nrow(all_result$data)
    # Approximate batch_size from batch_nodes and average node density
    # Default: split into chunks of ~batch_nodes * avg_pts_per_node
    if (total_rows > 0L) {
      batch_size <<- max(1000L, as.integer(total_rows / max(1L, total_rows %/% (batch_nodes * 100L))))
      # Cap at reasonable size
      batch_size <<- min(batch_size, total_rows)
    } else {
      batch_size <<- 0L
    }
  }

  # -- Build iterator environment --------------------------------------
  iter <- new.env(parent = emptyenv())

  iter$yield <- function() {
    .full_read()
    if (is.null(all_result) || current_offset >= nrow(all_result$data)) {
      return(NULL)
    }
    end_row <- min(current_offset + batch_size, nrow(all_result$data))
    batch <- all_result$data[(current_offset + 1L):end_row, ]
    current_offset <<- end_row
    batch
  }

  iter$has_next <- function() {
    .full_read()
    !is.null(all_result) && current_offset < nrow(all_result$data)
  }

  iter$header <- function() {
    .full_read()
    all_result$header
  }

  iter$reset <- function() {
    current_offset <<- 0L
    invisible(NULL)
  }

  iter$collect <- function(max_points = Inf) {
    .full_read()
    if (is.null(all_result)) return(data.table::data.table())
    if (is.finite(max_points) && nrow(all_result$data) > max_points) {
      return(all_result$data[seq_len(max_points), ])
    }
    all_result$data
  }

  iter$progress_info <- function() {
    .full_read()
    total <- nrow(all_result$data)
    list(
      total_points = total,
      yielded      = current_offset,
      remaining    = total - current_offset,
      batch_size   = batch_size
    )
  }

  class(iter) <- c("copc_iterator", "environment")
  iter
}


#' @exportS3Method
print.copc_iterator <- function(x, ...) {
  info <- tryCatch(x$progress_info(), error = function(e)
    list(total_points = "?", yielded = 0, remaining = "?", batch_size = "?"))
  cat("copc_iterator\n")
  cat("  Total points:", format(info$total_points, big.mark = ","), "\n")
  cat("  Yielded:     ", format(info$yielded, big.mark = ","), "\n")
  cat("  Remaining:   ", format(info$remaining, big.mark = ","), "\n")
  cat("  Batch size:  ", format(info$batch_size, big.mark = ","), "\n")
  invisible(x)
}
