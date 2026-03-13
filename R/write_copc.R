#' Write a COPC LAZ File
#'
#' Writes point-cloud data to a Cloud Optimized Point Cloud (COPC) LAZ
#' file.
#'
#' The output is a valid `.copc.laz` file readable by any COPC-aware
#' software (PDAL, QGIS, CloudCompare, copc4R, \dots) **and** any
#' standard LAZ reader (the point data is stored as a normal LAZ stream
#' with a chunk table).
#'
#' @param x Point-cloud data in one of four forms (see *Input types*).
#' @param filename Character.
#'   Output file path (should end in `.copc.laz`).
#' @param max_depth Integer.
#'   Maximum octree depth.
#'     \code{-1} (default) = auto (targets ~10 000 pts/leaf, max 8).
#'     \code{0} = all points in a single root node.
#' @param filter Character.
#'   Optional LAStools-style attribute filter applied **before** writing,
#'   e.g. `"-keep_class 2 6"`.  Only matching points are written.
#' @param progress Logical.  Print progress messages?
#'
#' @return Invisibly returns `filename`.
#'
#' @section Input types:
#' `write_copc()` is designed to "just work" with whatever you have:
#'
#' \describe{
#'   \item{`read_copc()` result}{A list with `$data` + `$header`.
#'     Header metadata (CRS, scale, offset, point format) is preserved.}
#'   \item{`lidR::LAS` object}{CRS, scale factors, offsets, and point
#'     format are extracted from the LAS header automatically.}
#'   \item{`data.table` / `data.frame`}{Must contain at least `X`, `Y`,
#'     `Z`.  Smart defaults are used: scale = 0.001, offsets computed
#'     from data extent, point format inferred from columns present.}
#'   \item{`sf` object}{Point or polygon geometry with LiDAR attributes
#'     -- coordinates extracted, CRS preserved.}
#' }
#'
#' @section Point format:
#' COPC requires LAS 1.4, point data record format 6, 7, or 8.
#' When no header is provided, the format is inferred from columns:
#'
#' \describe{
#'   \item{Format 6}{XYZ + intensity, classification, return info, \dots}
#'   \item{Format 7}{Format 6 + R, G, B}
#'   \item{Format 8}{Format 7 + NIR}
#' }
#'
#' @section Smart defaults (bare data.frame / data.table):
#' When writing a bare table without header metadata, `write_copc()`
#' mirrors the conventions used by lidR and rlas:
#'
#' \itemize{
#'   \item **Scale factors**: 0.001 (millimetre precision)
#'   \item **Offsets**: `floor(min(X))`, `floor(min(Y))`, 0
#'   \item **Bbox / point count**: computed from data
#'   \item **Point format**: inferred from column names
#'   \item **CRS**: none (unless an `sf` object with CRS is passed)
#' }
#'
#' @section Column name normalisation:
#' Common variants are silently mapped to rlas conventions:
#' `x`->`X`, `y`->`Y`, `z`->`Z`, `gps_time`->`gpstime`,
#' `intensity`->`Intensity`, `classification`->`Classification`,
#' `return_number`->`ReturnNumber`, etc.
#'
#' @examples
#' f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                  package = "copc4R")
#'
#' # Round-trip from read_copc()
#' result <- read_copc(f, progress = FALSE)
#' tmp <- tempfile(fileext = ".copc.laz")
#' write_copc(result, tmp, progress = FALSE)
#' back <- read_copc(tmp, progress = FALSE)
#' nrow(back$data)
#' unlink(tmp)
#'
#' \donttest{
#' # Write a lidR LAS object
#' if (requireNamespace("lidR", quietly = TRUE)) {
#'   las <- as_las(result)
#'   write_copc(las, tempfile(fileext = ".copc.laz"), progress = FALSE)
#' }
#'
#' # Write a bare data.frame -- smart defaults kick in
#' df <- data.frame(X = runif(100), Y = runif(100), Z = runif(100))
#' write_copc(df, tempfile(fileext = ".copc.laz"), progress = FALSE)
#' }
#'
#' @export
write_copc <- function(x,
                       filename,
                       max_depth = -1L,
                       filter    = NULL,
                       progress  = TRUE) {

    # -- 1. Unpack input into (dt, hdr) --------------------------------
    unpacked <- .unpack_input(x)
    dt  <- unpacked$dt
    hdr <- unpacked$hdr

    # -- 2. Normalise column names -------------------------------------
    dt <- .normalise_colnames(dt)

    if (!all(c("X", "Y", "Z") %in% names(dt)))
        stop("Input must contain at least X, Y, Z columns.", call. = FALSE)

    if (nrow(dt) == 0L)
        stop("Cannot write a file with 0 points.", call. = FALSE)

    # -- 3. Apply filter (if requested) --------------------------------
    if (!is.null(filter) && nzchar(filter)) {
        dt <- .apply_write_filter(dt, filter)
        if (nrow(dt) == 0L)
            stop("Filter removed all points -- nothing to write.",
                 call. = FALSE)
        if (progress)
            message(sprintf("write_copc: %s points after filter",
                            format(nrow(dt), big.mark = ",")))
    }

    # -- 4. Fill in smart defaults for missing header fields -----------
    hdr <- .fill_header_defaults(dt, hdr)

    # -- 5. Validate path ----------------------------------------------
    filename <- normalizePath(filename, mustWork = FALSE)
    if (!grepl("\\.copc\\.laz$", filename, ignore.case = TRUE))
        warning("Output filename does not end in '.copc.laz'.",
                call. = FALSE)

    if (!dir.exists(dirname(filename)))
        dir.create(dirname(filename), recursive = TRUE)

    # -- 6. Write ------------------------------------------------------
    cpp_write_copc(
        filename        = filename,
        data            = dt,
        header          = hdr,
        max_depth_input = as.integer(max_depth),
        progress        = progress
    )

    invisible(filename)
}


# =======================================================================
# Internal helpers
# =======================================================================

#' Unpack various input types into a data.table + header list
#' @noRd
.unpack_input <- function(x) {
    # -- read_copc() result (list with $data + $header) -------------
    if (is.list(x) && !is.data.frame(x) &&
        all(c("data", "header") %in% names(x))) {
        dt  <- data.table::as.data.table(x$data)
        hdr <- x$header
        return(list(dt = dt, hdr = hdr))
    }

    # -- lidR LAS object -----------------------------------------------
    if (inherits(x, "LAS")) {
        dt  <- data.table::as.data.table(x@data)
        hdr <- as.list(x@header@PHB)

        # Extract VLRs
        if (length(x@header@VLR) > 0)
            hdr[["Variable Length Records"]] <- x@header@VLR

        # Extract CRS as WKT and inject as a pseudo-VLR for the C++ writer
        crs <- sf::st_crs(x)
        if (!is.na(crs)) {
            hdr <- .inject_wkt_vlr(hdr, crs$wkt)
        }

        return(list(dt = dt, hdr = hdr))
    }

    # -- sf object (point cloud stored as sf points) ---------------
    if (inherits(x, "sf")) {
        crs <- sf::st_crs(x)
        coords <- sf::st_coordinates(x)
        dt <- data.table::as.data.table(sf::st_drop_geometry(x))
        data.table::set(dt, j = "X", value = coords[, 1])
        data.table::set(dt, j = "Y", value = coords[, 2])
        if (ncol(coords) >= 3)
            data.table::set(dt, j = "Z", value = coords[, 3])

        hdr <- list()
        if (!is.na(crs))
            hdr <- .inject_wkt_vlr(hdr, crs$wkt)

        return(list(dt = dt, hdr = hdr))
    }

    # -- bare data.frame / data.table ----------------------------------
    if (inherits(x, "data.frame") || inherits(x, "data.table")) {
        dt <- data.table::as.data.table(x)
        return(list(dt = dt, hdr = list()))
    }

    stop(paste0(
        "`x` must be one of:\n",
        "  - a read_copc() result (list with $data + $header)\n",
        "  - a lidR::LAS object\n",
        "  - a data.frame or data.table with X, Y, Z columns\n",
        "  - an sf object with point geometry"
    ), call. = FALSE)
}


#' Normalise common column name variants to rlas conventions
#' @noRd
.normalise_colnames <- function(dt) {
    name_map <- c(
        # Case variants
        "x" = "X", "y" = "Y", "z" = "Z",
        "intensity" = "Intensity",
        "classification" = "Classification",
        "returnnumber" = "ReturnNumber",
        "numberofreturns" = "NumberOfReturns",
        "userdata" = "UserData",
        "pointsourceid" = "PointSourceID",
        "scandirectionflag" = "ScanDirectionFlag",
        "edgeofflightline" = "EdgeOfFlightline",
        "scannerchannel" = "ScannerChannel",
        "scanangle" = "ScanAngleRank",
        "scananglerank" = "ScanAngleRank",
        "scan_angle" = "ScanAngleRank",
        "scan_angle_rank" = "ScanAngleRank",
        "gps_time" = "gpstime",
        "gpstime" = "gpstime",
        "r" = "R", "g" = "G", "b" = "B",
        "red" = "R", "green" = "G", "blue" = "B",
        "nir" = "NIR",
        # Underscore/snake_case variants
        "return_number" = "ReturnNumber",
        "number_of_returns" = "NumberOfReturns",
        "user_data" = "UserData",
        "point_source_id" = "PointSourceID",
        "scan_direction_flag" = "ScanDirectionFlag",
        "edge_of_flight_line" = "EdgeOfFlightline",
        "scanner_channel" = "ScannerChannel",
        "synthetic_flag" = "Synthetic_flag",
        "keypoint_flag" = "Keypoint_flag",
        "withheld_flag" = "Withheld_flag",
        "overlap_flag" = "Overlap_flag"
    )

    nms <- names(dt)
    lower <- tolower(nms)
    for (i in seq_along(nms)) {
        # Try exact match first, then lowercase
        if (nms[i] %in% names(name_map)) {
            data.table::setnames(dt, nms[i], name_map[nms[i]])
        } else if (lower[i] %in% names(name_map)) {
            data.table::setnames(dt, nms[i], name_map[lower[i]])
        }
    }

    dt
}


#' Inject a WKT string as a VLR entry in the header
#' @noRd
.inject_wkt_vlr <- function(hdr, wkt_string) {
    if (is.null(wkt_string) || !nzchar(wkt_string)) return(hdr)

    wkt_vlr <- list(
        `user ID`                 = "LASF_Projection",
        `record ID`               = 2112L,
        `description`             = "OGC WKT",
        `WKT OGC COORDINATE SYSTEM` = wkt_string
    )

    if (is.null(hdr[["Variable Length Records"]])) {
        hdr[["Variable Length Records"]] <- list(wkt_vlr)
    } else {
        # Replace existing WKT VLR if present
        vlrs <- hdr[["Variable Length Records"]]
        replaced <- FALSE
        for (i in seq_along(vlrs)) {
            if (!is.null(vlrs[[i]][["WKT OGC COORDINATE SYSTEM"]])) {
                vlrs[[i]] <- wkt_vlr
                replaced <- TRUE
                break
            }
        }
        if (!replaced) vlrs[[length(vlrs) + 1L]] <- wkt_vlr
        hdr[["Variable Length Records"]] <- vlrs
    }

    hdr
}


#' Fill in smart header defaults for any missing fields
#'
#' Mirrors lidR/rlas conventions:
#'   - Scale = 0.001 (mm precision)
#'   - X/Y offset = floor(min(X/Y)); Z offset = 0
#'   - Point format inferred from columns
#'   - Bbox + point count from data
#' @noRd
.fill_header_defaults <- function(dt, hdr) {

    # Scale factors
    if (is.null(hdr[["X scale factor"]]))
        hdr[["X scale factor"]] <- 0.001
    if (is.null(hdr[["Y scale factor"]]))
        hdr[["Y scale factor"]] <- 0.001
    if (is.null(hdr[["Z scale factor"]]))
        hdr[["Z scale factor"]] <- 0.001

    # Offsets -- compute from data extent (lidR convention)
    if (is.null(hdr[["X offset"]]))
        hdr[["X offset"]] <- floor(min(dt$X, na.rm = TRUE))
    if (is.null(hdr[["Y offset"]]))
        hdr[["Y offset"]] <- floor(min(dt$Y, na.rm = TRUE))
    if (is.null(hdr[["Z offset"]]))
        hdr[["Z offset"]] <- 0

    # Point format -- infer from columns present
    if (is.null(hdr[["Point Data Format ID"]])) {
        has_rgb <- all(c("R", "G", "B") %in% names(dt))
        has_nir <- "NIR" %in% names(dt)
        hdr[["Point Data Format ID"]] <-
            if (has_nir) 8L else if (has_rgb) 7L else 6L
    }

    hdr
}


#' Apply a LAStools-style filter string to a data.table (in R)
#'
#' Supports the most common filters used with copc4R:
#'   -keep_class, -drop_class, -keep_first, -keep_last,
#'   -keep_single, -drop_withheld, -drop_overlap,
#'   -keep_z_above, -keep_z_below, -drop_z_above, -drop_z_below
#' @noRd
.apply_write_filter <- function(dt, filter_str) {
    tokens <- strsplit(trimws(filter_str), "\\s+")[[1L]]
    keep <- rep(TRUE, nrow(dt))
    i <- 1L

    while (i <= length(tokens)) {
        tok <- tokens[i]

        if (tok == "-keep_class") {
            classes <- integer(0)
            while (i + 1L <= length(tokens) &&
                   !startsWith(tokens[i + 1L], "-")) {
                i <- i + 1L
                classes <- c(classes, as.integer(tokens[i]))
            }
            if ("Classification" %in% names(dt))
                keep <- keep & dt$Classification %in% classes

        } else if (tok == "-drop_class") {
            classes <- integer(0)
            while (i + 1L <= length(tokens) &&
                   !startsWith(tokens[i + 1L], "-")) {
                i <- i + 1L
                classes <- c(classes, as.integer(tokens[i]))
            }
            if ("Classification" %in% names(dt))
                keep <- keep & !(dt$Classification %in% classes)

        } else if (tok == "-keep_first") {
            if ("ReturnNumber" %in% names(dt))
                keep <- keep & dt$ReturnNumber == 1L

        } else if (tok == "-keep_last") {
            if (all(c("ReturnNumber", "NumberOfReturns") %in% names(dt)))
                keep <- keep & dt$ReturnNumber == dt$NumberOfReturns

        } else if (tok == "-keep_single") {
            if ("NumberOfReturns" %in% names(dt))
                keep <- keep & dt$NumberOfReturns == 1L

        } else if (tok == "-drop_withheld") {
            if ("Withheld_flag" %in% names(dt))
                keep <- keep & (dt$Withheld_flag == 0L |
                                dt$Withheld_flag == FALSE)

        } else if (tok == "-drop_overlap") {
            if ("Overlap_flag" %in% names(dt))
                keep <- keep & (dt$Overlap_flag == 0L |
                                dt$Overlap_flag == FALSE)

        } else if (tok %in% c("-keep_z_above", "-drop_z_below")) {
            i <- i + 1L
            val <- as.numeric(tokens[i])
            keep <- keep & dt$Z >= val

        } else if (tok %in% c("-keep_z_below", "-drop_z_above")) {
            i <- i + 1L
            val <- as.numeric(tokens[i])
            keep <- keep & dt$Z <= val

        } else if (tok == "-keep_return") {
            returns <- integer(0)
            while (i + 1L <= length(tokens) &&
                   !startsWith(tokens[i + 1L], "-")) {
                i <- i + 1L
                returns <- c(returns, as.integer(tokens[i]))
            }
            if ("ReturnNumber" %in% names(dt))
                keep <- keep & dt$ReturnNumber %in% returns

        } else if (tok == "-drop_return") {
            returns <- integer(0)
            while (i + 1L <= length(tokens) &&
                   !startsWith(tokens[i + 1L], "-")) {
                i <- i + 1L
                returns <- c(returns, as.integer(tokens[i]))
            }
            if ("ReturnNumber" %in% names(dt))
                keep <- keep & !(dt$ReturnNumber %in% returns)
        }
        # Unknown filters are silently ignored (they may be read-time only)

        i <- i + 1L
    }

    dt[keep, ]
}
