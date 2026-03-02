#' Read a COPC (.copc.laz) point cloud
#'
#' Reads point records from a Cloud Optimized Point Cloud file, applying
#' optional spatial queries via the COPC hierarchy.
#'
#' @param path_or_url Character. Path to a local .copc.laz file, or an
#'   HTTP(S) URL.  Remote files are transparently downloaded to a
#'   temporary file the first time they are accessed in a session.
#' @param bbox Numeric vector of length 4: \code{c(xmin, ymin, xmax, ymax)}.
#'   Only points within this 2-D bounding box are returned.
#'   \code{NULL} (default) reads all points.
#' @param zrange Numeric vector of length 2: \code{c(zmin, zmax)}.
#'   Additional Z-range filter. \code{NULL} (default) = no Z filter.
#' @param select Character. Column selection string following lidR-style
#'   semantics:
#'   \itemize{
#'     \item \code{"xyz"} = X, Y, Z
#'     \item \code{"i"} = Intensity
#'     \item \code{"t"} = gpstime
#'     \item \code{"r"} = ReturnNumber
#'     \item \code{"n"} = NumberOfReturns
#'     \item \code{"c"} = Classification
#'     \item \code{"a"} = ScanAngleRank
#'     \item \code{"s"} = PointSourceID
#'     \item \code{"u"} = UserData
#'     \item \code{"p"} = flags (ScanDirectionFlag, EdgeOfFlightline, Synthetic/Keypoint/Withheld/Overlap, ScannerChannel)
#'     \item \code{"R"} = RGB (if PDRF >= 7)
#'     \item \code{"N"} = NIR (if PDRF 8)
#'     \item \code{"*"} or \code{NULL} = all available columns
#'   }
#' @param max_points Numeric. Maximum number of points to return.
#'   \code{Inf} (default) returns all selected points.
#' @param progress Logical. Print progress messages? Default \code{TRUE}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{$data}}{A \code{data.table} of point attributes.}
#'     \item{\code{$header}}{A named list of header fields following
#'       rlas naming conventions.}
#'   }
#'
#' @export
#' @importFrom data.table data.table setDT
read_copc <- function(path_or_url,
                      bbox       = NULL,
                      zrange     = NULL,
                      select     = NULL,
                      max_points = Inf,
                      progress   = TRUE) {

  stopifnot(is.character(path_or_url), length(path_or_url) == 1L)

  # Handle HTTP(S) URLs: download to a session-cached tempfile
  path <- .resolve_path(path_or_url, progress = progress)

  # Warn early if file extension doesn't look like COPC
  .check_copc_ext(path)

  if (is.null(bbox))   bbox   <- numeric(0)
  if (is.null(zrange)) zrange <- numeric(0)
  if (is.null(select)) select <- "*"

  stopifnot(is.numeric(bbox) && (length(bbox) == 0L || length(bbox) == 4L))
  stopifnot(is.numeric(zrange) && (length(zrange) == 0L || length(zrange) == 2L))
  stopifnot(is.character(select), length(select) == 1L)
  stopifnot(is.numeric(max_points), length(max_points) == 1L)

  result <- .Call("_copc4R_cpp_read_copc",
                  path,
                  as.numeric(bbox),
                  as.numeric(zrange),
                  as.character(select),
                  as.double(max_points),
                  as.logical(progress),
                  PACKAGE = "copc4R")

  # Ensure $data is a proper data.table
  if (!inherits(result$data, "data.table")) {
    data.table::setDT(result$data)
  }
  # Assign key alloc.col for data.table
  data.table::alloc.col(result$data)

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
#' @export
read_copc_header <- function(path_or_url) {
  stopifnot(is.character(path_or_url), length(path_or_url) == 1L)

  path <- .resolve_path(path_or_url, progress = FALSE)
  .check_copc_ext(path)

  .Call("_copc4R_cpp_read_copc_header",
        path,
        PACKAGE = "copc4R")
}

# ── Internal: resolve path/URL ─────────────────────────────────────────
# Session-level cache for downloaded remote files.
.url_cache <- new.env(parent = emptyenv())

.resolve_path <- function(path_or_url, progress = FALSE) {
  # Check for HTTP(S) URL
  if (grepl("^https?://", path_or_url, ignore.case = TRUE)) {
    # Check URL extension before downloading
    .check_copc_ext_url(path_or_url)

    # Return cached tempfile if already downloaded this session
    cached <- .url_cache[[path_or_url]]
    if (!is.null(cached) && file.exists(cached)) {
      if (progress) message("Using cached download: ", cached)
      return(cached)
    }

    # Download to a tempfile
    # Temporarily increase timeout — large COPC tiles (100+ MB) can
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
      "  - pdal:    pdal translate input.laz output.copc.laz\n",
      "  - untwine: untwine --files input.laz -o output_dir\n",
      "  - See https://copc.io for more information.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
