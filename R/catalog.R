# -- lidR catalog ergonomics -----------------------------------------------
#
# Builds a LAScatalog-like object from one or more COPC URLs/paths.
# If lidR >= 4.0 is installed we create a real LAScatalog; otherwise we
# return a lightweight "copc_catalog" S3 object that can be chunked.
# ---------------------------------------------------------------------------

#' Build a LAScatalog from COPC files
#'
#' Creates a \code{lidR::LAScatalog} (when \pkg{lidR} is installed) from
#' one or more COPC paths/URLs.  Headers are read via the COPC hierarchy
#' without downloading full files, and the resulting catalog inherits
#' lidR processing options (\code{opt_chunk_size()}, etc.).
#'
#' @param paths Character vector.  Paths or URLs to \code{.copc.laz} files.
#' @param chunk_size Numeric.  Default chunk size in CRS units
#'   (0 = file-level processing).  Sets \code{opt_chunk_size()}.
#' @param chunk_buffer Numeric.  Buffer around each chunk.
#'   Sets \code{opt_chunk_buffer()}.
#' @param select Character.  Default column selection for processing.
#' @param filter Character.  Default attribute filter string.
#' @param progress Logical.
#'
#' @return A \code{lidR::LAScatalog} object if \pkg{lidR} is installed,
#'   otherwise a \code{copc_catalog} S3 list with file metadata.
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#' ctg <- read_copc_catalog(f, progress = FALSE)
#' print(ctg)
#' }
#'
#' @export
read_copc_catalog <- function(paths,
                              chunk_size   = 0,
                              chunk_buffer = 0,
                              select       = "*",
                              filter       = "",
                              progress     = TRUE) {
  stopifnot(is.character(paths), length(paths) >= 1L)

  # -- Collect headers -------------------------------------------------
  hdrs <- lapply(paths, function(p) {
    if (progress) message("Reading header: ", basename(sub("[?#].*$", "", p)))
    tryCatch(
      read_copc_header(p),
      error = function(e) {
        warning("Skipping ", p, ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  })

  # Remove failed reads
  ok <- !vapply(hdrs, is.null, FALSE)
  paths <- paths[ok]
  hdrs  <- hdrs[ok]
  if (length(hdrs) == 0L) stop("No valid COPC headers found.", call. = FALSE)

  # -- Build spatial data.frame ----------------------------------------
  file_info <- data.frame(
    filename = paths,
    Min.X    = vapply(hdrs, `[[`, 0.0, "Min X"),
    Max.X    = vapply(hdrs, `[[`, 0.0, "Max X"),
    Min.Y    = vapply(hdrs, `[[`, 0.0, "Min Y"),
    Max.Y    = vapply(hdrs, `[[`, 0.0, "Max Y"),
    Min.Z    = vapply(hdrs, `[[`, 0.0, "Min Z"),
    Max.Z    = vapply(hdrs, `[[`, 0.0, "Max Z"),
    Number.of.point.records = vapply(hdrs, `[[`, 0.0, "Number of point records"),
    Point.Data.Format.ID    = vapply(hdrs, `[[`, 0L,  "Point Data Format ID"),
    stringsAsFactors = FALSE
  )

  # -- Try to create a lidR LAScatalog ---------------------------------
  if (requireNamespace("lidR", quietly = TRUE)) {
    # lidR::readLAScatalog() expects actual files, so we build
    # a LAScatalog from scratch using the header information.
    # We leverage lidR >= 4.0 which accepts data.frames.
    ctg <- tryCatch({
      .build_lidR_catalog(file_info, hdrs, chunk_size, chunk_buffer,
                          select, filter)
    }, error = function(e) {
      warning("Could not build lidR::LAScatalog: ", conditionMessage(e),
              "\nFalling back to copc_catalog S3.", call. = FALSE)
      NULL
    })
    if (!is.null(ctg)) return(ctg)
  }

  # -- Fallback: lightweight S3 catalog --------------------------------
  catalog <- structure(
    list(
      files   = file_info,
      headers = hdrs,
      options = list(
        chunk_size   = chunk_size,
        chunk_buffer = chunk_buffer,
        select       = select,
        filter       = filter
      )
    ),
    class = "copc_catalog"
  )
  catalog
}


#' @export
print.copc_catalog <- function(x, ...) {
  cat("copc_catalog with", nrow(x$files), "file(s)\n")
  cat("  Total points:", format(sum(x$files$Number.of.point.records),
                                big.mark = ","), "\n")
  cat("  Extent X: [",
      round(min(x$files$Min.X), 2), ",",
      round(max(x$files$Max.X), 2), "]\n")
  cat("  Extent Y: [",
      round(min(x$files$Min.Y), 2), ",",
      round(max(x$files$Max.Y), 2), "]\n")
  cat("  Options: chunk_size =", x$options$chunk_size,
      ", filter = '", x$options$filter, "'\n")
  invisible(x)
}


#' Process a copc_catalog or LAScatalog tile-by-tile
#'
#' Iterates over chunks of a \code{copc_catalog} or \code{lidR::LAScatalog},
#' reading and processing each chunk via a user-supplied function.
#'
#' When \code{chunk_size > 0} (set via \code{\link{read_copc_catalog}}),
#' each file's extent is divided into a regular tile grid.  Only the COPC
#' octree nodes overlapping each tile's bounding box are fetched -- so
#' memory usage is proportional to a single chunk, not the entire dataset.
#'
#' @param catalog A \code{copc_catalog} or \code{lidR::LAScatalog} object
#'   (as returned by \code{\link{read_copc_catalog}}).
#' @param fun A function that takes a single \code{read_copc()} result
#'   (list with \code{$data} and \code{$header}) and returns a value.
#' @param ... Additional arguments passed to \code{fun}.
#' @param cores Integer.  Number of parallel workers (default 1 =
#'   sequential).  When \code{cores > 1}, chunks are processed in
#'   parallel via \code{parallel::parLapply()} with a PSOCK cluster.
#'   Each worker loads \pkg{copc4R} and independently fetches its
#'   chunk from the COPC file(s) -- ideal for remote data where
#'   network I/O is the bottleneck.
#' @param progress Logical.
#'
#' @return A list of results, one per chunk.
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#' ctg <- read_copc_catalog(f, chunk_size = 100, progress = FALSE)
#' results <- copc_catalog_apply(ctg, function(x) nrow(x$data), progress = FALSE)
#' }
#'
#' \donttest{
#' # Parallel processing (4 cores)
#' results <- copc_catalog_apply(ctg, function(x) {
#'   data.frame(n = nrow(x$data), z_mean = mean(x$data$Z))
#' }, cores = 4L, progress = FALSE)
#' }
#'
#' @export
copc_catalog_apply <- function(catalog, fun, ..., cores = 1L, progress = TRUE) {
  UseMethod("copc_catalog_apply")
}

#' @export
copc_catalog_apply.copc_catalog <- function(catalog, fun, ...,
                                       cores = 1L,
                                       progress = TRUE) {
  opts   <- catalog$options
  extras <- list(...)

  # -- Build chunk specifications ----------------------------------
  # Each chunk = list(file, bbox).  Only the COPC octree nodes that
  # overlap the bbox are fetched, keeping peak memory bounded to a
  # single chunk regardless of overall dataset size.
  chunks <- list()

  for (i in seq_len(nrow(catalog$files))) {
    fi <- catalog$files[i, ]
    f  <- fi$filename

    if (opts$chunk_size > 0) {
      xs  <- seq(fi$Min.X, fi$Max.X, by = opts$chunk_size)
      ys  <- seq(fi$Min.Y, fi$Max.Y, by = opts$chunk_size)
      buf <- opts$chunk_buffer

      for (xi in xs) {
        for (yi in ys) {
          chunks[[length(chunks) + 1L]] <- list(
            file = f,
            bbox = c(xi - buf, yi - buf,
                     min(xi + opts$chunk_size + buf, fi$Max.X + buf),
                     min(yi + opts$chunk_size + buf, fi$Max.Y + buf))
          )
        }
      }
    } else {
      chunks[[length(chunks) + 1L]] <- list(file = f, bbox = NULL)
    }
  }

  n_chunks <- length(chunks)
  if (progress)
    message(sprintf("copc_catalog_apply: %d chunk(s) from %d file(s)%s",
                    n_chunks, nrow(catalog$files),
                    if (as.integer(cores) > 1L)
                      sprintf(", %d cores", as.integer(cores)) else ""))

  # -- Chunk worker ------------------------------------------------
  .process_chunk <- function(chunk, .fun, .opts, .extras) {
    result <- copc4R::read_copc(
      path_or_url = chunk$file,
      bbox        = chunk$bbox,
      select      = .opts$select,
      filter      = .opts$filter,
      progress    = FALSE
    )
    if (nrow(result$data) == 0L) return(NULL)
    do.call(.fun, c(list(result), .extras))
  }

  # -- Dispatch: parallel or sequential ----------------------------
  cores <- min(as.integer(cores), n_chunks)

  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(copc4R))
    results <- parallel::parLapply(
      cl, chunks, .process_chunk,
      .fun = fun, .opts = opts, .extras = extras
    )
  } else {
    results <- vector("list", n_chunks)
    for (i in seq_along(chunks)) {
      if (progress)
        message(sprintf("  chunk %d/%d", i, n_chunks))
      results[[i]] <- .process_chunk(chunks[[i]], fun, opts, extras)
    }
  }

  results
}


#' @export
copc_catalog_apply.LAScatalog <- function(catalog, fun, ...,
                                     cores = 1L,
                                     progress = TRUE) {
  files  <- catalog@data$filename
  sel    <- lidR::opt_select(catalog)
  flt    <- lidR::opt_filter(catalog)
  csz    <- lidR::opt_chunk_size(catalog)
  cbuf   <- lidR::opt_chunk_buffer(catalog)
  extras <- list(...)
  opts   <- list(select = sel, filter = flt,
                 chunk_size = csz, chunk_buffer = cbuf)

  # Build chunks
  chunks <- list()
  for (i in seq_along(files)) {
    f <- files[i]
    if (csz > 0) {
      hdr <- read_copc_header(f)
      xs  <- seq(hdr[["Min X"]], hdr[["Max X"]], by = csz)
      ys  <- seq(hdr[["Min Y"]], hdr[["Max Y"]], by = csz)
      for (xi in xs) {
        for (yi in ys) {
          chunks[[length(chunks) + 1L]] <- list(
            file = f,
            bbox = c(xi - cbuf, yi - cbuf,
                     min(xi + csz + cbuf, hdr[["Max X"]] + cbuf),
                     min(yi + csz + cbuf, hdr[["Max Y"]] + cbuf))
          )
        }
      }
    } else {
      chunks[[length(chunks) + 1L]] <- list(file = f, bbox = NULL)
    }
  }

  n_chunks <- length(chunks)
  if (progress)
    message(sprintf("copc_catalog_apply: %d chunk(s) from %d file(s)%s",
                    n_chunks, length(files),
                    if (as.integer(cores) > 1L)
                      sprintf(", %d cores", as.integer(cores)) else ""))

  .process_chunk <- function(chunk, .fun, .opts, .extras) {
    result <- copc4R::read_copc(
      path_or_url = chunk$file,
      bbox        = chunk$bbox,
      select      = .opts$select,
      filter      = .opts$filter,
      progress    = FALSE
    )
    if (nrow(result$data) == 0L) return(NULL)
    do.call(.fun, c(list(result), .extras))
  }

  cores <- min(as.integer(cores), n_chunks)

  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(copc4R))
    results <- parallel::parLapply(
      cl, chunks, .process_chunk,
      .fun = fun, .opts = opts, .extras = extras
    )
  } else {
    results <- vector("list", n_chunks)
    for (i in seq_along(chunks)) {
      if (progress)
        message(sprintf("  chunk %d/%d", i, n_chunks))
      results[[i]] <- .process_chunk(chunks[[i]], fun, opts, extras)
    }
  }

  results
}


# -- Internal: build a real lidR LAScatalog --------------------------------

.build_lidR_catalog <- function(file_info, hdrs, chunk_size, chunk_buffer,
                                select, filter) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for lidR catalog creation.", call. = FALSE)

  # Create sf polygons from file extents
  polys <- lapply(seq_len(nrow(file_info)), function(i) {
    fi <- file_info[i, ]
    sf::st_polygon(list(matrix(c(
      fi$Min.X, fi$Min.Y,
      fi$Max.X, fi$Min.Y,
      fi$Max.X, fi$Max.Y,
      fi$Min.X, fi$Max.Y,
      fi$Min.X, fi$Min.Y
    ), ncol = 2, byrow = TRUE)))
  })
  geom <- sf::st_sfc(polys)

  # Try to get CRS from first header
  crs <- sf::NA_crs_
  for (v in hdrs[[1]][["Variable Length Records"]]) {
    if (!is.null(v[["WKT OGC COORDINATE SYSTEM"]])) {
      crs <- sf::st_crs(v[["WKT OGC COORDINATE SYSTEM"]])
      break
    }
  }
  sf::st_crs(geom) <- crs

  file_sf <- sf::st_sf(file_info, geometry = geom)

  # Build LAScatalog
  ctg <- new("LAScatalog")
  ctg@data <- file_sf

  # Set processing options
  lidR::opt_chunk_size(ctg)   <- chunk_size
  lidR::opt_chunk_buffer(ctg) <- chunk_buffer
  lidR::opt_select(ctg)       <- select
  lidR::opt_filter(ctg)       <- filter

  ctg
}


# ========================================================================
# copc_apply() -- COPC-native parallel chunk processing engine
# ========================================================================
#
# The spiritual successor to copc_catalog_apply(): splits a spatial extent
# into a regular grid of chunks, reads each chunk (plus buffer) from
# one or more COPC files via range queries, converts to lidR::LAS,
# applies a user function, clips buffer from the output, and optionally
# auto-merges all chunk results.
#
# Parallelism is handled by the future framework so users can choose
# plan(multisession), plan(multicore), plan(cluster), etc.
# ========================================================================

#' Apply a function over COPC data in parallel spatial chunks
#'
#' Partitions one or more COPC point cloud files into a regular grid of
#' spatial chunks, reads each chunk (plus a buffer zone) via COPC range
#' queries, converts to a \code{lidR::LAS} object, and applies a
#' user-supplied function.
#'
#' This is the COPC-native equivalent of \code{lidR::catalog_apply()}
#' and the successor to \code{\link{copc_catalog_apply}()}.
#' Unlike a traditional \code{LAScatalog} workflow that reads from local
#' \code{.las}/\code{.laz} tiles, \code{copc_apply()} queries the COPC
#' octree spatial index and streams only the needed data for each chunk
#' -- via HTTP range reads for remote files, or direct seeks for local
#' files.
#'
#' @section Chunk grid:
#' The processing area is divided into a regular grid of squares with
#' side length \code{chunk_size} (in the COPC files' native CRS units,
#' typically metres for projected CRS).  The grid origin is snapped to
#' \code{chunk_alignment} so that the tiling is reproducible across
#' different AOIs.
#'
#' @section Buffer handling:
#' Each chunk is read with an additional \code{chunk_buffer} margin on
#' all four sides.  This is critical for algorithms that rely on spatial
#' context beyond the chunk boundary (ground classification, height
#' normalisation, canopy models, tree segmentation, etc.).
#'
#' After \code{FUN} returns, the buffer is \strong{automatically clipped}
#' from the output so that the final product covers exactly the chunk
#' footprint -- no overlaps and no edge artifacts when results are merged.
#'
#' Clipping is type-aware:
#' \itemize{
#'   \item \code{sf} / \code{sfc} -- \code{sf::st_crop()}
#'   \item \code{SpatRaster}      -- \code{terra::crop()}
#'   \item \code{LAS}             -- \code{lidR::clip_rectangle()}
#'   \item \code{data.frame} / \code{data.table} with X, Y columns --
#'     row-filtered to chunk bbox
#'   \item All other types         -- returned unchanged (no clipping)
#' }
#'
#' @section Parallelism:
#' Parallelism is controlled by the \pkg{future} framework.
#' Set your plan \emph{before} calling \code{copc_apply()}:
#' \preformatted{
#'   library(future)
#'   plan(multisession, workers = 4)
#'   # ... copc_apply(...) ...
#'   plan(sequential)
#' }
#' Each worker reads its chunks independently from the COPC source(s),
#' so network I/O is the primary bottleneck for remote files and the
#' work scales well.  If no parallel plan is set (or if \pkg{future.apply}
#' is not installed), chunks are processed sequentially.
#'
#' @param source Character vector of COPC file paths or HTTP(S) URLs,
#'   or a \code{copc_catalog} object as returned by
#'   \code{\link{read_copc_catalog}}.
#' @param FUN Function applied to each chunk.
#'   Receives a single \code{lidR::LAS} object (chunk + buffer points)
#'   and must return a result or \code{NULL} for empty / skipped chunks.
#' @param ... Additional arguments forwarded to \code{FUN}.
#' @param aoi Optional \code{sf}/\code{sfc} polygon restricting the
#'   processing extent.  Automatically reprojected to the COPC native
#'   CRS if needed.
#' @param chunk_size Numeric. Side length of each square chunk in CRS
#'   units.  Default \code{200}.
#' @param chunk_buffer Numeric. Buffer distance added around each chunk
#'   for spatial context.  Points in the buffer are available to
#'   \code{FUN} but are clipped from the output.  Default \code{15}.
#' @param chunk_alignment Numeric length-2 vector \code{c(x, y)} that
#'   shifts the grid origin for reproducible tiling.
#'   Default \code{c(0, 0)}.
#' @param select Character. Column selection string forwarded to
#'   \code{\link{read_copc}()}.  Default \code{"*"} (all columns).
#' @param filter Character. Attribute filter string forwarded to
#'   \code{\link{read_copc}()}.
#' @param packages Character vector of additional package names to load
#'   in parallel workers.  Core packages (\pkg{copc4R}, \pkg{lidR},
#'   \pkg{sf}, \pkg{terra}, \pkg{data.table}) are auto-detected and
#'   always included when available.
#' @param automerge Logical.  When \code{TRUE} (default), results are
#'   merged based on type: \code{sf} via \code{rbind()},
#'   \code{SpatRaster} via \code{terra::mosaic()},
#'   \code{data.frame} via \code{rbind()}, \code{LAS} via
#'   \code{rbind()}.  Unrecognised types are returned as a list.
#' @param progress Logical.  When \code{TRUE} (default), displays a
#'   text progress bar during sequential processing and prints
#'   summary messages.  For parallel (\pkg{future}) execution the
#'   \pkg{progressr} package is supported: wrap the call in
#'   \code{progressr::with_progress()} to see updates.
#' @param max_points_per_chunk Numeric or \code{NULL}.  When set,
#'   after building the regular chunk grid each chunk is checked
#'   against the COPC octree to estimate its point count.  Chunks
#'   exceeding this threshold are recursively split (quadtree-style)
#'   into smaller sub-chunks until the estimate drops below the limit
#'   (or \code{min_split_size} is reached).  This implements the
#'   adaptive work-unit pattern: balance memory/complexity per chunk
#'   rather than using a fixed spatial grid.
#'   Default \code{NULL} (disabled -- use uniform spatial grid only).
#' @param min_split_size Numeric.  Minimum side length (CRS units) for
#'   adaptive sub-chunking.  Below this threshold splitting stops even
#'   if \code{max_points_per_chunk} is exceeded.  Default \code{25}.
#' @param plot Logical.  When \code{TRUE}, opens a graphics window
#'   showing the chunk grid and colours each tile as it is processed
#'   (green = success, orange = empty/skipped).  This mimics the
#'   spatial progress display used by \pkg{lidR}'s catalog engine.
#'   In sequential mode the plot updates in real time; in parallel
#'   mode the final status is drawn after all chunks complete.
#'   Default \code{FALSE}.
#'
#' @return If \code{automerge = TRUE} and all non-NULL results share a
#'   recognised type, the merged object.
#'   Otherwise, a list of per-chunk results (\code{NULL} entries for
#'   empty or failed chunks).
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#'
#' # Simple per-chunk summary (sequential)
#' results <- copc_apply(f, function(las) {
#'   data.frame(n = lidR::npoints(las), zmean = mean(las$Z))
#' }, chunk_size = 100, chunk_buffer = 0, progress = FALSE)
#' print(results)
#' }
#'
#' \donttest{
#' # Parallel processing with the future framework
#' library(future)
#' plan(multisession, workers = 4)
#'
#' results <- copc_apply(
#'   source = "https://example.com/data.copc.laz",
#'   FUN = function(las) {
#'     las <- lidR::classify_ground(las, lidR::csf())
#'     nlas <- lidR::normalize_height(las, lidR::tin())
#'     lidR::rasterize_canopy(nlas, res = 1, lidR::pitfree())
#'   },
#'   chunk_size = 500,
#'   chunk_buffer = 30,
#'   automerge = TRUE
#' )
#' plan(sequential)
#' }
#'
#' @export
copc_apply <- function(source,
                       FUN,
                       ...,
                       aoi                  = NULL,
                       chunk_size           = 200,
                       chunk_buffer         = 15,
                       chunk_alignment      = c(0, 0),
                       max_points_per_chunk = NULL,
                       min_split_size       = 25,
                       select               = "*",
                       filter               = "",
                       packages             = NULL,
                       automerge            = TRUE,
                       progress             = TRUE,
                       plot                 = FALSE) {

  # -- validation --------------------------------------------------------
  if (!requireNamespace("lidR", quietly = TRUE))
    stop("Package 'lidR' is required for copc_apply().", call. = FALSE)
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for copc_apply().", call. = FALSE)

  FUN <- match.fun(FUN)
  stopifnot(is.numeric(chunk_size), length(chunk_size) == 1L, chunk_size > 0)
  stopifnot(is.numeric(chunk_buffer), length(chunk_buffer) == 1L,
            chunk_buffer >= 0)
  stopifnot(is.numeric(chunk_alignment), length(chunk_alignment) == 2L)

  extras <- list(...)

  # -- resolve source to character paths ---------------------------------
  if (inherits(source, "copc_catalog")) {
    paths <- source$files$filename
  } else if (inherits(source, "LAScatalog")) {
    paths <- source@data$filename
  } else {
    paths <- as.character(source)
  }
  stopifnot(length(paths) >= 1L)

  # -- read headers & extract CRS / extents ------------------------------
  if (progress) message("copc_apply: reading ", length(paths), " header(s)...")
  hdrs <- lapply(paths, read_copc_header)

  file_crs <- .detect_file_crs(hdrs[[1L]])

  file_bboxes <- do.call(rbind, lapply(hdrs, function(h) {
    c(xmin = h[["Min X"]], ymin = h[["Min Y"]],
      xmax = h[["Max X"]], ymax = h[["Max Y"]])
  }))
  if (!is.matrix(file_bboxes))
    file_bboxes <- matrix(file_bboxes, nrow = 1,
                          dimnames = list(NULL, c("xmin","ymin","xmax","ymax")))

  # -- determine processing extent ---------------------------------------
  if (!is.null(aoi)) {
    if (!is.null(file_crs)) aoi <- sf::st_transform(aoi, file_crs)
    bb <- as.numeric(sf::st_bbox(aoi))
    proc_ext <- c(xmin = bb[1], ymin = bb[2], xmax = bb[3], ymax = bb[4])
  } else {
    proc_ext <- c(
      xmin = min(file_bboxes[, "xmin"]),
      ymin = min(file_bboxes[, "ymin"]),
      xmax = max(file_bboxes[, "xmax"]),
      ymax = max(file_bboxes[, "ymax"])
    )
  }

  # -- build aligned chunk grid ------------------------------------------
  grid <- .copc_chunk_grid(proc_ext, chunk_size, chunk_alignment)
  n_chunks <- nrow(grid)

  if (n_chunks == 0L) {
    warning("Degenerate processing extent -- no chunks created.", call. = FALSE)
    return(if (automerge) NULL else list())
  }

  # -- adaptive splitting (optional) -------------------------------------
  # If max_points_per_chunk is set, subdivide chunks that exceed the

  # threshold using the octree hierarchy (no decompression needed).
  if (!is.null(max_points_per_chunk)) {
    stopifnot(is.numeric(max_points_per_chunk),
              length(max_points_per_chunk) == 1L,
              max_points_per_chunk > 0)
    stopifnot(is.numeric(min_split_size),
              length(min_split_size) == 1L,
              min_split_size > 0)

    if (progress)
      message(sprintf("  Adaptive splitting: max %s pts/chunk, min size %.0fm",
                      format(max_points_per_chunk, big.mark = ","),
                      min_split_size))

    new_grid <- list()

    for (i in seq_len(nrow(grid))) {
      bb <- c(grid$xmin[i], grid$ymin[i], grid$xmax[i], grid$ymax[i])
      sub_tiles <- .adaptive_quadtree_split(
        paths, bb, max_points_per_chunk, min_split_size
      )
      for (st in sub_tiles) {
        new_grid[[length(new_grid) + 1L]] <- st
      }
    }

    grid <- do.call(rbind, new_grid)
    n_chunks <- nrow(grid)

    if (progress)
      message(sprintf("  After adaptive split: %d chunks", n_chunks))
  }

  # -- find which source files overlap each chunk (+buffer) --------------
  chunk_files <- vector("list", n_chunks)
  for (i in seq_len(n_chunks)) {
    bx <- c(grid$xmin[i] - chunk_buffer, grid$ymin[i] - chunk_buffer,
            grid$xmax[i] + chunk_buffer, grid$ymax[i] + chunk_buffer)
    hit <- which(
      file_bboxes[, "xmin"] < bx[3] & file_bboxes[, "xmax"] > bx[1] &
      file_bboxes[, "ymin"] < bx[4] & file_bboxes[, "ymax"] > bx[2]
    )
    chunk_files[[i]] <- paths[hit]
  }

  # drop chunks with no overlapping files
  keep <- vapply(chunk_files, function(x) length(x) > 0L, FALSE)
  grid        <- grid[keep, , drop = FALSE]
  chunk_files <- chunk_files[keep]
  n_active    <- nrow(grid)

  if (n_active == 0L) {
    warning("No chunks overlap the source file(s).", call. = FALSE)
    return(if (automerge) NULL else list())
  }

  if (progress) {
    message(sprintf(
      "copc_apply: %d chunks (%.0fm grid, %.0fm buffer) from %d file(s)",
      n_active, chunk_size, chunk_buffer, length(paths)))
  }

  # -- initialise plot progress ------------------------------------------
  if (plot) {
    .plot_progress_init(grid, file_bboxes, chunk_buffer)
  }

  # -- build task list ---------------------------------------------------
  tasks <- lapply(seq_len(n_active), function(i) {
    list(
      id          = i,
      chunk_bbox  = c(grid$xmin[i], grid$ymin[i],
                      grid$xmax[i], grid$ymax[i]),
      buffer_bbox = c(grid$xmin[i] - chunk_buffer,
                      grid$ymin[i] - chunk_buffer,
                      grid$xmax[i] + chunk_buffer,
                      grid$ymax[i] + chunk_buffer),
      files       = chunk_files[[i]]
    )
  })

  # -- chunk worker (closure over select, filter, FUN, extras) -----------
  .worker <- function(task) {
    tryCatch({
      # Read from every overlapping file with buffer
      parts <- list()
      for (f in task$files) {
        res <- copc4R::read_copc(
          path_or_url = f,
          bbox        = task$buffer_bbox,
          select      = select,
          filter      = filter,
          progress    = FALSE
        )
        if (nrow(res$data) > 0L)
          parts[[length(parts) + 1L]] <- copc4R::as_las(res)
      }
      if (length(parts) == 0L) return(NULL)

      # Merge tiles that contribute to this chunk
      las <- if (length(parts) == 1L) parts[[1L]] else do.call(rbind, parts)

      # Apply user function
      out <- do.call(FUN, c(list(las), extras))

      # Clip buffer from result
      out <- .clip_chunk_buffer(out, task$chunk_bbox)

      # Wrap terra objects so they survive serialisation back to the
      # main process (external pointers are not valid across processes).
      .wrap_for_transport(out)

    }, error = function(e) {
      warning(sprintf("copc_apply chunk %d: %s",
                      task$id, conditionMessage(e)), call. = FALSE)
      NULL
    })
  }

  # -- dispatch: future or sequential ------------------------------------
  use_future <- requireNamespace("future.apply", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE) &&
    !inherits(future::plan(), "sequential")

  if (use_future) {
    # Auto-detect packages for workers
    worker_pkgs <- "copc4R"
    for (.pkg in c("lidR", "sf", "terra", "data.table", "rlas"))
      if (requireNamespace(.pkg, quietly = TRUE))
        worker_pkgs <- c(worker_pkgs, .pkg)
    worker_pkgs <- unique(c(worker_pkgs, packages))

    plan_label <- class(future::plan())[1L]
    n_workers  <- tryCatch(future::nbrOfWorkers(),
                           error = function(e) NA_integer_)
    if (progress)
      message(sprintf("  Dispatching to %s (%s workers)...",
                      plan_label,
                      if (is.na(n_workers)) "?" else n_workers))

    if (progress &&
        requireNamespace("progressr", quietly = TRUE)) {
      p <- progressr::progressor(along = tasks)
      results <- future.apply::future_lapply(
        tasks, function(t) { out <- .worker(t); p(); out },
        future.seed     = TRUE,
        future.packages = worker_pkgs
      )
    } else {
      results <- future.apply::future_lapply(
        tasks, .worker,
        future.seed     = TRUE,
        future.packages = worker_pkgs
      )
    }
  } else {
    results <- vector("list", n_active)
    if (progress && n_active > 1L) {
      pb <- utils::txtProgressBar(min = 0, max = n_active, style = 3)
      on.exit(close(pb), add = TRUE)
      for (i in seq_len(n_active)) {
        results[[i]] <- .worker(tasks[[i]])
        utils::setTxtProgressBar(pb, i)
        if (plot) .plot_progress_update(grid, i, !is.null(results[[i]]))
      }
    } else {
      for (i in seq_len(n_active)) {
        results[[i]] <- .worker(tasks[[i]])
        if (plot) .plot_progress_update(grid, i, !is.null(results[[i]]))
      }
    }
  }

  n_ok <- sum(!vapply(results, is.null, FALSE))

  # -- parallel: draw final status on plot ------------------------------
  if (plot && use_future) {
    for (i in seq_len(n_active)) {
      .plot_progress_update(grid, i, !is.null(results[[i]]))
    }
  }

  if (progress)
    message(sprintf("copc_apply: %d/%d chunks returned results", n_ok, n_active))

  # Unwrap any packed terra objects before merging / returning
  results <- lapply(results, .unwrap_from_transport)

  if (automerge) return(.copc_automerge(results))
  results
}


# -- copc_apply helpers ----------------------------------------------------

#' @noRd
.detect_file_crs <- function(hdr) {
  for (v in hdr[["Variable Length Records"]]) {
    wkt <- v[["WKT OGC COORDINATE SYSTEM"]]
    if (!is.null(wkt) && nzchar(wkt))
      return(sf::st_crs(wkt))
  }
  NULL
}


#' Build an aligned chunk grid covering an extent
#' @noRd
.copc_chunk_grid <- function(extent, size, alignment) {
  ext <- unname(extent)
  if (ext[3] <= ext[1] || ext[4] <= ext[2])
    return(data.frame(xmin = numeric(0), ymin = numeric(0),
                      xmax = numeric(0), ymax = numeric(0)))

  x0 <- alignment[1] + floor((ext[1] - alignment[1]) / size) * size
  y0 <- alignment[2] + floor((ext[2] - alignment[2]) / size) * size
  x1 <- alignment[1] + ceiling((ext[3] - alignment[1]) / size) * size
  y1 <- alignment[2] + ceiling((ext[4] - alignment[2]) / size) * size

  # Ensure at least one chunk even for tiny extents
  if (x1 <= x0) x1 <- x0 + size
  if (y1 <= y0) y1 <- y0 + size

  xs <- seq(x0, x1 - size, by = size)
  ys <- seq(y0, y1 - size, by = size)

  g <- expand.grid(xi = seq_along(xs), yi = seq_along(ys))

  data.frame(
    xmin = xs[g$xi],
    ymin = ys[g$yi],
    xmax = xs[g$xi] + size,
    ymax = ys[g$yi] + size
  )
}


#' Clip a chunk result to remove the buffer region
#'
#' Dispatches on the result type: sf, SpatRaster, LAS, data.frame.
#' Unknown types are returned unchanged.
#' @noRd
.clip_chunk_buffer <- function(result, bbox) {
  if (is.null(result)) return(NULL)

  # sf / sfc
  if (requireNamespace("sf", quietly = TRUE) &&
      (inherits(result, "sf") || inherits(result, "sfc"))) {
    box <- sf::st_bbox(c(xmin = bbox[1], ymin = bbox[2],
                         xmax = bbox[3], ymax = bbox[4]),
                       crs = sf::st_crs(result))
    return(suppressWarnings(sf::st_crop(result, box)))
  }

  # SpatRaster
  if (requireNamespace("terra", quietly = TRUE) &&
      inherits(result, "SpatRaster")) {
    ext <- terra::ext(bbox[1], bbox[3], bbox[2], bbox[4])
    return(terra::crop(result, ext))
  }

  # LAS
  if (inherits(result, "LAS")) {
    return(lidR::clip_rectangle(result,
                                bbox[1], bbox[2], bbox[3], bbox[4]))
  }

  # data.frame / data.table with X, Y columns
  if (is.data.frame(result) && all(c("X", "Y") %in% names(result))) {
    inside <- result$X >= bbox[1] & result$X <= bbox[3] &
              result$Y >= bbox[2] & result$Y <= bbox[4]
    return(result[inside, , drop = FALSE])
  }

  result
}


#' Wrap terra objects for safe transport across process boundaries
#'
#' SpatRaster and SpatVector use C++ external pointers that become
#' invalid when serialised to a parallel worker and back.  This wraps
#' them into inert R objects that can be safely transferred, mirroring
#' what lidR's engine does internally.
#' @noRd
.wrap_for_transport <- function(x) {
  if (is.null(x)) return(NULL)
  if (requireNamespace("terra", quietly = TRUE)) {
    if (inherits(x, "SpatRaster") || inherits(x, "SpatVector"))
      return(terra::wrap(x))
  }
  x
}

#' Restore packed terra objects after transport
#' @noRd
.unwrap_from_transport <- function(x) {
  if (is.null(x)) return(NULL)
  if (requireNamespace("terra", quietly = TRUE)) {
    if (inherits(x, "PackedSpatRaster") || inherits(x, "PackedSpatVector"))
      return(terra::unwrap(x))
  }
  x
}


#' Auto-merge a list of chunk results
#'
#' Dispatches on the first non-NULL result's type.
#' Falls back to returning the list if merging fails.
#' @noRd
.copc_automerge <- function(results) {
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) return(NULL)
  if (length(results) == 1L) return(results[[1L]])

  first <- results[[1L]]

  # sf
  if (inherits(first, "sf"))
    return(tryCatch(do.call(rbind, results),
                    error = function(e) results))

  # SpatRaster -- mosaic handles non-overlapping tiles cleanly
  if (inherits(first, "SpatRaster"))
    return(tryCatch({
      src <- terra::sprc(results)
      terra::mosaic(src)
    }, error = function(e) results))

  # data.frame / data.table
  if (is.data.frame(first))
    return(tryCatch(do.call(rbind, results),
                    error = function(e) results))

  # LAS
  if (inherits(first, "LAS"))
    return(tryCatch(do.call(rbind, results),
                    error = function(e) results))

  results
}


# -- Plot progress helpers -------------------------------------------------

#' Initialise the chunk-grid progress plot
#'
#' Draws the processing extent, source file outlines, and the chunk grid.
#' Subsequent calls to \code{.plot_progress_update()} colour individual
#' tiles as they complete.
#' @noRd
.plot_progress_init <- function(grid, file_bboxes, buffer) {
  # overall extent (with a small margin for labels)
  xr <- range(grid$xmin, grid$xmax)
  yr <- range(grid$ymin, grid$ymax)
  dx <- diff(xr) * 0.04
  dy <- diff(yr) * 0.04

  graphics::par(mar = c(2.5, 2.5, 1.5, 0.5))

  graphics::plot.new()
  graphics::plot.window(xlim = c(xr[1] - dx, xr[2] + dx),
                        ylim = c(yr[1] - dy, yr[2] + dy),
                        asp  = 1)
  graphics::title(main = sprintf("copc_apply: %d chunks", nrow(grid)),
                  cex.main = 0.9)
  graphics::axis(1, cex.axis = 0.7)
  graphics::axis(2, cex.axis = 0.7)

  # Draw source-file extents
  if (is.matrix(file_bboxes)) {
    for (fi in seq_len(nrow(file_bboxes))) {
      graphics::rect(file_bboxes[fi, "xmin"], file_bboxes[fi, "ymin"],
                     file_bboxes[fi, "xmax"], file_bboxes[fi, "ymax"],
                     border = "grey40", lwd = 1.5, lty = 2, col = NA)
    }
  }

  # Draw chunk grid (unfilled)
  graphics::rect(grid$xmin, grid$ymin, grid$xmax, grid$ymax,
                 col = NA, border = "grey70", lwd = 0.4)

  # Chunk number labels (skip if too many)
  n <- nrow(grid)
  if (n <= 200) {
    cx <- (grid$xmin + grid$xmax) / 2
    cy <- (grid$ymin + grid$ymax) / 2
    cex_label <- if (n > 100) 0.4 else if (n > 50) 0.5 else 0.6
    graphics::text(cx, cy, labels = seq_len(n), cex = cex_label,
                   col = "grey50")
  }
}


#' Colour a single chunk on the progress plot
#'
#' @param grid  The chunk grid data.frame.
#' @param i     The chunk index to update (1-based).
#' @param ok    Logical; \code{TRUE} = success (green),
#'              \code{FALSE} = empty/skipped (orange).
#' @noRd
.plot_progress_update <- function(grid, i, ok) {
  fill <- if (ok) grDevices::adjustcolor("forestgreen", alpha.f = 0.45)
          else    grDevices::adjustcolor("orange",      alpha.f = 0.35)
  graphics::rect(grid$xmin[i], grid$ymin[i],
                 grid$xmax[i], grid$ymax[i],
                 col = fill, border = "grey70", lwd = 0.4)
}


# -- Adaptive quadtree split (used by copc_apply max_points_per_chunk) -----

#' Recursively split a bbox until estimated points are below threshold
#'
#' Queries the COPC octree hierarchy across all `paths` to sum point
#' estimates.  Splits quadtree-style when the sum exceeds `max_pts`.
#' @noRd
.adaptive_quadtree_split <- function(paths, bbox, max_pts, min_size) {
  # Sum estimates across all overlapping files
  total_est <- 0L
  for (p in paths) {
    path <- .resolve_path(p, progress = FALSE)
    node_info <- tryCatch(
      cpp_count_nodes(path, as.numeric(bbox)),
      error = function(e) list(total_points = 0L)
    )
    total_est <- total_est + node_info$total_points
  }

  width  <- bbox[3] - bbox[1]
  height <- bbox[4] - bbox[2]

  if (total_est <= max_pts || width <= min_size || height <= min_size) {
    return(list(data.frame(
      xmin = bbox[1], ymin = bbox[2],
      xmax = bbox[3], ymax = bbox[4],
      stringsAsFactors = FALSE
    )))
  }

  # Split into 4 quadrants
  xmid <- (bbox[1] + bbox[3]) / 2
  ymid <- (bbox[2] + bbox[4]) / 2

  quads <- list(
    c(bbox[1], bbox[2], xmid,    ymid),     # SW
    c(xmid,    bbox[2], bbox[3], ymid),     # SE
    c(bbox[1], ymid,    xmid,    bbox[4]),  # NW
    c(xmid,    ymid,    bbox[3], bbox[4])   # NE
  )

  result <- list()
  for (q in quads) {
    sub <- .adaptive_quadtree_split(paths, q, max_pts, min_size)
    result <- c(result, sub)
  }
  result
}
