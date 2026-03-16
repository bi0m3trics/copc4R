#' Read a COPC (.copc.laz) point cloud
#'
#' Reads point records from a Cloud Optimized Point Cloud file, applying
#' optional spatial queries via the COPC hierarchy.
#'
#' @param path_or_url Character. One of:
#'   \itemize{
#'     \item A local path to a \code{.copc.laz} file.
#'     \item A direct HTTP(S) URL to a \code{.copc.laz} file (streamed via
#'       HTTP range reads).
#'     \item A \strong{STAC search endpoint URL} (e.g.
#'       \code{"https://planetarycomputer.microsoft.com/api/stac/v1/search"}).
#'       When a STAC URL is detected, the function queries the catalog using
#'       \code{aoi}, discovers all intersecting COPC tiles, streams only the
#'       AOI portion from each tile, and returns a single merged result.
#'       Tile signing (e.g. Planetary Computer SAS tokens) is handled
#'       automatically.
#'   }
#' @param bbox Numeric vector of length 4: \code{c(xmin, ymin, xmax, ymax)}.
#'   Only points within this 2-D bounding box are returned.
#'   \code{NULL} (default) reads all points.
#'   You may also pass a numeric(4) directly to \code{aoi} as a shorthand.
#' @param aoi Spatial query region.  Accepts any of the following -- the
#'   function figures out the right thing to do automatically:
#'   \itemize{
#'     \item An \code{sf} / \code{sfc} \strong{POLYGON or MULTIPOLYGON}: the
#'       COPC octree is queried via the polygon bounding box and points are
#'       clipped to the exact polygon boundary.
#'     \item An \code{sf} / \code{sfc} \strong{LINESTRING or MULTILINESTRING}:
#'       buffered by \code{buffer} (required) to form a corridor polygon, then
#'       handled as above.
#'     \item An \code{sf} / \code{sfc} \strong{POINT or MULTIPOINT}: buffered
#'       by \code{buffer} (required) to form a circle/union polygon, then
#'       handled as above.
#'     \item A \strong{numeric vector of length 4}
#'       \code{c(xmin, ymin, xmax, ymax)}: equivalent to supplying \code{bbox}.
#'   }
#'   The AOI is automatically reprojected to the file's native CRS when
#'   the CRS of \code{aoi} differs from the file.  Requires \pkg{sf}.
#' @param buffer Numeric. Buffer distance (in the AOI's CRS units) applied
#'   when \code{aoi} is a point or line geometry.  Required for those
#'   geometry types; ignored for polygons.  Default \code{NULL}.
#' @param zrange Numeric vector of length 2: \code{c(zmin, zmax)}.
#'   Additional Z-range filter. \code{NULL} (default) = no Z filter.
#' @param select Character. Column selection using lidR-style single-letter codes:
#'   \code{"x"} X, \code{"y"} Y, \code{"z"} Z, \code{"i"} Intensity,
#'   \code{"t"} gpstime, \code{"r"} ReturnNumber, \code{"n"} NumberOfReturns,
#'   \code{"c"} Classification, \code{"a"} ScanAngleRank,
#'   \code{"s"} PointSourceID, \code{"u"} UserData, \code{"p"} flags,
#'   \code{"R"} RGB (PDRF\eqn{\ge}7), \code{"N"} NIR (PDRF 8).
#'   \code{"*"} or \code{NULL} (default) returns all columns.
#' @param filter Character.  LAStools/lidR-style filter string.  Predicates
#'   can be combined, e.g.
#'   \code{"-keep_class 2 6 -drop_withheld -keep_voxel 0.5"}.
#'   See the \sQuote{Point filters} section for all supported predicates.
#' @param max_depth Integer.  Maximum octree depth.  Lower values yield
#'   sparser (coarser) clouds, useful for previews or LOD sampling.
#'   \code{-1L} (default) fetches all levels.
#' @param resolution Numeric.  Target point spacing (CRS units).  When
#'   supplied the appropriate \code{max_depth} is computed automatically.
#'   Overrides \code{max_depth}.  Equivalent to PDAL \code{readers.copc}
#'   \code{resolution}.
#' @param max_points Numeric.  Cap on total points returned.
#'   \code{Inf} (default) = no cap.
#' @param threads Integer.  Parallel download threads for remote files.
#'   Ignored for local files.  Default \code{4L}.
#' @param collections Character.  STAC collection ID(s) to search when
#'   \code{path_or_url} is a STAC endpoint.  Automatically set to
#'   \code{"3dep-lidar-copc"} for the Microsoft Planetary Computer.
#'   Only needed for other STAC catalogs.  Default \code{NULL}.
#' @param progress Logical. Print progress messages?  Default \code{TRUE}.
#'
#' @section Point filters:
#' The \code{filter} parameter accepts LAStools-style predicates; combine
#' any number in a single string.
#' \describe{
#'   \item{Classification}{\code{-keep_class <c1> <c2> ...},
#'     \code{-drop_class <c1> <c2> ...},
#'     \code{-keep_ground} (class 2),
#'     \code{-drop_noise} (classes 7 & 18)}
#'   \item{Return number}{\code{-keep_first}, \code{-keep_last},
#'     \code{-keep_single}, \code{-keep_return <r1> ...},
#'     \code{-drop_return <r1> ...}}
#'   \item{Flags}{\code{-drop_withheld}, \code{-drop_overlap}}
#'   \item{Intensity}{\code{-keep_intensity_above <val>},
#'     \code{-keep_intensity_below <val>}}
#'   \item{Z-value}{\code{-keep_z_above <val>}, \code{-keep_z_below <val>}}
#'   \item{Scan angle}{\code{-keep_scan_angle_above <deg>},
#'     \code{-keep_scan_angle_below <deg>}}
#'   \item{Thinning}{\code{-keep_random_fraction <0-1>},
#'     \code{-keep_every_nth <n>},
#'     \code{-keep_voxel <size>} (3-D voxel grid; keeps the point nearest
#'     the voxel center per occupied cell; applied after other predicates)}
#' }
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{$data}}{A \code{data.table} of point attributes.
#'       Extra Bytes dimensions (if present in the file) are included
#'       as additional numeric columns.}
#'     \item{\code{$header}}{A named list of header fields following
#'       rlas naming conventions.}
#'   }
#'
#' @examples
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#'
#' # Read all points
#' result <- read_copc(f, progress = FALSE)
#' head(result$data)
#' nrow(result$data)
#'
#' \donttest{
#' # Bounding box as aoi shorthand (numeric(4))
#' hdr <- read_copc_header(f)
#' xmid <- (hdr[["Min X"]] + hdr[["Max X"]]) / 2
#' ymid <- (hdr[["Min Y"]] + hdr[["Max Y"]]) / 2
#' ground <- read_copc(f, aoi = c(xmid-50, ymid-50, xmid+50, ymid+50),
#'                     filter = "-keep_class 2", progress = FALSE)
#'
#' # Polygon AOI -- auto-clips to boundary
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   gpkg <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82.gpkg",
#'                       package = "copc4R")
#'   bldg <- sf::st_read(gpkg, quiet = TRUE)
#'   result <- read_copc(f, aoi = bldg, progress = FALSE)
#'
#'   # Point AOI with buffer -- auto-buffers and clips
#'   pt <- sf::st_sf(geometry = sf::st_sfc(
#'     sf::st_point(c(xmid, ymid)), crs = sf::st_crs(bldg)))
#'   result_pt <- read_copc(f, aoi = pt, buffer = 20, progress = FALSE)
#' }
#'
#' # Coarse LOD preview (fewer points)
#' preview <- read_copc(f, max_depth = 2L, progress = FALSE)
#'
#' # Voxel thinning via the filter string
#' thinned <- read_copc(f, filter = "-keep_voxel 1.0", progress = FALSE)
#' }
#'
#' @export
#' @importFrom data.table data.table setDT
read_copc <- function(path_or_url,
                      bbox        = NULL,
                      aoi         = NULL,
                      buffer      = NULL,
                      zrange      = NULL,
                      select      = NULL,
                      filter      = NULL,
                      max_depth   = -1L,
                      resolution  = NULL,
                      max_points  = Inf,
                      threads     = 4L,
                      collections = NULL,
                      progress    = TRUE) {

  # --- Input validation -----------------------------------------------
  if (!is.character(path_or_url) || length(path_or_url) != 1L ||
      is.na(path_or_url) || !nzchar(trimws(path_or_url)))
    stop("'path_or_url' must be a single non-empty string.", call. = FALSE)
  path_or_url <- trimws(path_or_url)

  # -- STAC endpoint detection --------------------------------------------
  # Delegate to the multi-tile STAC reader when path_or_url looks like a
  # STAC search API endpoint.  Heuristics (in order of confidence):
  #   1. Contains /stac/v*/search or /api/stac
  #   2. Ends with /search (common STAC API pattern)
  #   3. Does NOT end in .laz / .copc.laz (i.e. not a direct file URL)
  # Static catalog roots (catalog.json) are NOT search APIs -- supply
  # the specific /search endpoint for those.
  is_stac_search <- grepl("^https?://", path_or_url, ignore.case = TRUE) &&
    !grepl("\\.copc\\.laz|\\.laz$", path_or_url, ignore.case = TRUE) &&
    (
      grepl("/stac/v[0-9]/search|/api/stac", path_or_url, ignore.case = TRUE) ||
      grepl("/search$|/search\\?",           path_or_url, ignore.case = TRUE) ||
      grepl("/stac/",                        path_or_url, ignore.case = TRUE)
    )

  if (is_stac_search) {
    if (is.null(aoi))
      stop("'aoi' is required when 'path_or_url' is a STAC search endpoint.\n",
           "Supply a polygon sf/sfc, point/line with 'buffer', or numeric(4) bbox.",
           call. = FALSE)
    return(.read_copc_from_stac(
      stac_url    = path_or_url,
      aoi         = aoi,
      buffer      = buffer,
      collections = collections,
      zrange      = zrange,
      select      = select,
      filter      = filter,
      max_depth   = max_depth,
      resolution  = resolution,
      max_points  = max_points,
      threads     = threads,
      progress    = progress
    ))
  }

  # Handle HTTP(S) URLs: download to a session-cached tempfile
  path <- .resolve_path(path_or_url, progress = progress)

  # Warn early if file extension doesn't look like COPC
  .check_copc_ext(path)

  # -- AOI handling ---------------------------------------------------
  # Use .normalise_aoi() to robustly handle sf/sfc/sfg/numeric inputs,
  # empty geometries, missing CRS, geometry collections, etc.
  clip_geom <- NULL
  if (!is.null(aoi)) {
    aoi_norm <- .normalise_aoi(aoi, buffer = buffer, caller = "read_copc")

    if (aoi_norm$type == "bbox") {
      # numeric(4) shorthand -- bbox only, no polygon clipping
      if (is.null(bbox)) bbox <- aoi_norm$bbox
    } else {
      aoi_sfc <- aoi_norm$sfc

      # Auto-reproject to the file's native CRS when they differ
      aoi_crs  <- sf::st_crs(aoi_sfc)
      file_crs <- tryCatch(
        .crs_from_copc_header(read_copc_header(path_or_url)),
        error = function(e) sf::st_crs(NA)
      )
      if (!is.na(file_crs) && !is.na(aoi_crs) &&
          !identical(aoi_crs, file_crs)) {
        if (progress)
          message("Auto-reprojecting AOI from ", aoi_crs$input,
                  " to file CRS (", file_crs$input, ")")
        aoi_sfc <- tryCatch(
          sf::st_transform(aoi_sfc, file_crs),
          error = function(e)
            stop("CRS reprojection failed: ", conditionMessage(e), call. = FALSE)
        )
      }

      if (is.null(bbox)) bbox <- as.numeric(sf::st_bbox(aoi_sfc))
      clip_geom <- aoi_sfc
    }
  }

  if (is.null(bbox))   bbox   <- numeric(0)
  if (is.null(zrange)) zrange <- numeric(0)
  if (is.null(select)) select <- "*"
  if (is.null(filter)) filter <- ""

  # -- Extract -keep_voxel from filter (handled in R, not C++) ------
  voxel_size <- NULL
  if (nzchar(filter)) {
    m <- regmatches(filter, regexec("-keep_voxel\\s+([0-9.eE+\\-]+)", filter))
    if (length(m[[1]]) == 2L) {
      voxel_size <- as.numeric(m[[1]][2L])
      if (is.na(voxel_size) || voxel_size <= 0)
        stop("-keep_voxel requires a positive numeric value.", call. = FALSE)
      # Strip the token from the filter string before passing to C++
      filter <- trimws(sub("-keep_voxel\\s+[0-9.eE+\\-]+", "", filter))
    }
  }

  stopifnot(is.numeric(bbox) && (length(bbox) == 0L || length(bbox) == 4L))
  stopifnot(is.numeric(zrange) && (length(zrange) == 0L || length(zrange) == 2L))
  stopifnot(is.character(select), length(select) == 1L)
  stopifnot(is.character(filter), length(filter) == 1L)
  stopifnot(is.numeric(max_points), length(max_points) == 1L)
  max_depth <- as.integer(max_depth)
  threads   <- as.integer(threads)

  # -- Resolution -> max_depth conversion -----------------------------
  # COPC spacing at octree level d = spacing / 2^d.
  # To achieve target resolution r: d = floor(log2(spacing / r))
  if (!is.null(resolution)) {
    stopifnot(is.numeric(resolution), length(resolution) == 1L, resolution > 0)
    info <- copc_info(path_or_url)
    spacing <- info$spacing
    if (!is.null(spacing) && spacing > 0) {
      computed_depth <- as.integer(floor(log2(spacing / resolution)))
      computed_depth <- max(0L, computed_depth)
      if (progress) {
        message(sprintf(
          "Resolution %.2f -> max_depth %d (root spacing=%.2f)",
          resolution, computed_depth, spacing))
      }
      max_depth <- computed_depth
    } else {
      warning("Could not determine COPC spacing; ignoring 'resolution' parameter.")
    }
  }

  # -- Cache-aware orchestration for HTTP URLs ------------------------
  is_http    <- grepl("^https?://", path, ignore.case = TRUE)
  use_cache  <- .cache_enabled() && is_http
  prefetched <- NULL

  if (use_cache) {
    # Select nodes via C++, then check chunk cache
    node_info <- cpp_select_nodes(path, as.numeric(bbox),
                                  as.numeric(zrange), max_depth)
    n_nodes <- length(node_info$offset)

    if (n_nodes > 0L) {
      prefetched  <- vector("list", n_nodes)
      uncached    <- logical(n_nodes)

      for (i in seq_len(n_nodes)) {
        cached <- .cache_get(path, node_info$offset[i], node_info$byte_size[i])
        if (!is.null(cached)) {
          prefetched[[i]] <- cached
        } else {
          uncached[i] <- TRUE
        }
      }

      cache_hits  <- sum(!uncached)
      uncached_idx <- which(uncached)

      if (progress && cache_hits > 0L)
        message(sprintf("Cache: %d/%d chunks hit", cache_hits, n_nodes))

      # Parallel-fetch uncached chunks, then store them in cache
      if (length(uncached_idx) > 0L) {
        fetched <- cpp_fetch_raw_chunks_parallel(
          path,
          node_info$offset[uncached_idx],
          node_info$byte_size[uncached_idx],
          threads
        )
        for (j in seq_along(uncached_idx)) {
          i <- uncached_idx[j]
          prefetched[[i]] <- fetched[[j]]
          .cache_put(path, node_info$offset[i],
                     node_info$byte_size[i], fetched[[j]])
        }
      }
    }
  }

  result <- .Call("_copc4R_cpp_read_copc",
                  path,
                  as.numeric(bbox),
                  as.numeric(zrange),
                  as.character(select),
                  as.double(max_points),
                  as.logical(progress),
                  max_depth,
                  as.character(filter),
                  if (use_cache) 1L else threads,
                  prefetched,
                  PACKAGE = "copc4R")

  # Ensure $data is a proper data.table
  if (!inherits(result$data, "data.table")) {
    data.table::setDT(result$data)
  }
  # Assign key alloc.col for data.table
  data.table::alloc.col(result$data)

  # -- Post-clip to AOI polygon -------------------------------------
  if (!is.null(clip_geom) && nrow(result$data) > 0L) {
    pts_sf <- sf::st_as_sf(result$data, coords = c("X", "Y"),
                           remove = FALSE,
                           crs = sf::st_crs(clip_geom))
    inside <- sf::st_intersects(pts_sf, sf::st_union(clip_geom),
                                sparse = FALSE)[, 1]
    result$data <- result$data[inside, ]
  }

  # -- Post-filter: voxel thinning (keep point nearest voxel center) -
  if (!is.null(voxel_size) && nrow(result$data) > 0L) {
    n_before <- nrow(result$data)
    result$data <- .voxel_center_thin(result$data, voxel_size)
    if (progress)
      message(sprintf("Voxel thin (%.2f, center): %s -> %s points",
                      voxel_size,
                      format(n_before, big.mark = ","),
                      format(nrow(result$data), big.mark = ",")))
  }

  result
}


#' Read COPC header only
#'
#' Reads and parses the LAS 1.4 header, VLRs, EVLRs, and COPC Info VLR
#' from a .copc.laz file without reading any point data.
#'
#' @param path_or_url Character. Path or URL to a .copc.laz file.
#'
#' @return A named list of header fields following rlas naming conventions,
#'   including a \code{$`COPC Info`} sub-list with octree metadata.
#'
#' @examples
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#' hdr <- read_copc_header(f)
#' hdr[["File Signature"]]
#' hdr[["Number of point records"]]
#' hdr[["Point Data Format ID"]]
#' hdr[["COPC Info"]]$spacing
#'
#' @export
read_copc_header <- function(path_or_url) {
  stopifnot(is.character(path_or_url), length(path_or_url) == 1L)

  # Check header cache first (avoids a round-trip for repeated calls)
  cached <- .header_cache_get(path_or_url)
  if (!is.null(cached)) return(cached)

  path <- .resolve_path(path_or_url, progress = FALSE)
  .check_copc_ext(path)

  result <- .Call("_copc4R_cpp_read_copc_header",
                  path,
                  PACKAGE = "copc4R")

  # Cache the result (keyed by original path_or_url AND resolved path)
  .header_cache_put(path_or_url, result)
  if (!identical(path_or_url, path))
    .header_cache_put(path, result)

  result
}

# -- Internal: resolve path/URL -----------------------------------------
# Session-level cache for downloaded remote files.
.url_cache <- new.env(parent = emptyenv())

.resolve_path <- function(path_or_url, progress = FALSE) {
  # Check for HTTP(S) URL
  if (grepl("^https?://", path_or_url, ignore.case = TRUE)) {
    # Check URL extension before downloading
    .check_copc_ext_url(path_or_url)

    # If C++ HTTP range-read support is compiled in, pass the URL directly
    # so the COPC hierarchy can be traversed via range requests (only the
    # header + relevant octree nodes are fetched, not the whole file).
    if (cpp_has_http_support()) {
      if (progress) message("Using HTTP range reads (streaming)")
      return(path_or_url)
    }

    # Fallback: C++ does not have curl -- download the entire file.
    if (progress)
      message("HTTP range-read support not available; downloading full file")

    # Return cached tempfile if already downloaded this session
    cached <- .url_cache[[path_or_url]]
    if (!is.null(cached) && file.exists(cached)) {
      if (progress) message("Using cached download: ", cached)
      return(cached)
    }

    # Download to a tempfile
    # Temporarily increase timeout -- large COPC tiles (100+ MB) can
    # exceed R's default 60-second limit, especially from Azure Europe.
    tmp <- tempfile(fileext = ".copc.laz")
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = max(old_timeout, 600L))

    if (progress) {
      # Attempt to get file size via a HEAD request for progress reporting
      file_size <- NA_real_
      tryCatch({
        if (requireNamespace("curl", quietly = TRUE)) {
          h <- curl::new_handle()
          curl::handle_setopt(h, nobody = TRUE)
          head_resp <- curl::curl_fetch_memory(path_or_url, handle = h)
          cl <- head_resp$headers
          # Parse Content-Length from raw headers
          hdr_text <- rawToChar(cl)
          m <- regmatches(hdr_text,
                          regexpr("Content-Length:\\s*[0-9]+", hdr_text,
                                  ignore.case = TRUE))
          if (length(m) == 1L)
            file_size <- as.numeric(sub(".*:\\s*", "", m))
        }
      }, error = function(e) NULL)

      if (!is.na(file_size) && file_size > 0) {
        message(sprintf("Downloading: %s (%.1f MB)",
                        basename(sub("[?#].*$", "", path_or_url)),
                        file_size / 1e6))
      } else {
        message(sprintf("Downloading: %s",
                        basename(sub("[?#].*$", "", path_or_url))))
      }

      # Use curl for download with a real progress bar if available
      if (requireNamespace("curl", quietly = TRUE)) {
        h <- curl::new_handle()
        curl::handle_setopt(h, connecttimeout = 60L)
        tryCatch(
          curl::curl_download(path_or_url, tmp, handle = h, quiet = FALSE),
          error = function(e) {
            # Fall back to base download.file on curl error
            utils::download.file(path_or_url, tmp, mode = "wb", quiet = FALSE)
          }
        )
      } else {
        utils::download.file(path_or_url, tmp, mode = "wb", quiet = FALSE)
      }
    } else {
      utils::download.file(path_or_url, tmp, mode = "wb", quiet = TRUE)
    }

    if (!file.exists(tmp) || file.size(tmp) == 0L)
      stop("Download failed for: ", path_or_url, call. = FALSE)
    if (progress) message("Downloaded ", round(file.size(tmp) / 1e6, 1), " MB")

    .url_cache[[path_or_url]] <- tmp
    return(tmp)
  }

  # Local file: verify existence
  if (!file.exists(path_or_url))
    stop("File not found: ", path_or_url, call. = FALSE)
  path_or_url
}

.check_copc_ext <- function(path) {
  bn <- tolower(basename(path))
  if (grepl("\\.copc\\.laz$", bn)) return(invisible(NULL))
  if (grepl("\\.la[sz]$", bn)) {
    stop(
      "The file '", basename(path), "' appears to be a standard LAS/LAZ ",
      "file, not a COPC file. copc4R only reads Cloud Optimized Point ",
      "Cloud (.copc.laz) files.\n\n",
      "To convert LAS/LAZ to COPC, use one of:\n",
      "  - copc4R:  copc4R::write_copc(lidR::readLAS('input.laz'), 'output.copc.laz')\n",
      "  - lasR:    lasR::exec(lasR::write_copc('output.copc.laz'), on = 'input.laz')\n",
      "  - pdal:    pdal translate input.laz output.copc.laz\n",
      "  - untwine: untwine --files input.laz -o output_dir\n",
      "  - See https://copc.io for more information.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.check_copc_ext_url <- function(url) {
  # Strip query string / fragment before checking extension
  clean <- sub("[?#].*$", "", url)
  bn <- tolower(basename(clean))
  if (grepl("\\.copc\\.laz$", bn)) return(invisible(NULL))
  if (grepl("\\.la[sz]$", bn)) {
    stop(
      "The URL appears to point to a standard LAS/LAZ file ('", bn, "'), ",
      "not a COPC file. copc4R only reads Cloud Optimized Point Cloud ",
      "(.copc.laz) files.\n\n",
      "To convert LAS/LAZ to COPC, use one of:\n",
      "  - copc4R:  copc4R::write_copc(lidR::readLAS('input.laz'), 'output.copc.laz')\n",
      "  - lasR:    lasR::exec(lasR::write_copc('output.copc.laz'), on = 'input.laz')\n",
      "  - pdal:    pdal translate input.laz output.copc.laz\n",
      "  - untwine: untwine --files input.laz -o output_dir\n",
      "  - See https://copc.io for more information.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


# -- Internal: voxel center thinning -----------------------------------------
#
# Overlays a 3-D voxel grid of cell size `voxel_size` on the point cloud
# and keeps, for each occupied voxel, the single point closest to the
# voxel center (barycenter).
#
# @param dt  A data.table with at least X, Y, Z columns.
# @param voxel_size  Positive numeric, cell size in CRS units.
# @return  A filtered data.table (same columns as input).
# @noRd

.voxel_center_thin <- function(dt, voxel_size) {
  # Compute integer voxel indices
  vx <- as.integer(floor(dt$X / voxel_size))
  vy <- as.integer(floor(dt$Y / voxel_size))
  vz <- as.integer(floor(dt$Z / voxel_size))

  # Voxel center coordinates
  cx <- (vx + 0.5) * voxel_size
  cy <- (vy + 0.5) * voxel_size
  cz <- (vz + 0.5) * voxel_size

  # Squared distance from each point to its voxel center
  dist_sq <- (dt$X - cx)^2 + (dt$Y - cy)^2 + (dt$Z - cz)^2

  # Group by voxel, keep the index of the nearest point
  grp <- paste(vx, vy, vz, sep = "_")
  tmp <- data.table::data.table(
    .idx = seq_len(nrow(dt)),
    .grp = grp,
    .d   = dist_sq
  )
  keep_idx <- sort(tmp[, list(.idx = .idx[which.min(.d)]), by = .grp]$.idx)

  dt[keep_idx, ]
}
