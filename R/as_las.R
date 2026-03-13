#' Convert copc4R output to a lidR LAS object
#'
#' If the \pkg{lidR} package is installed, \code{as_las()} converts the
#' output of \code{read_copc()} to a \code{lidR::LAS} object.
#'
#' @param x A list as returned by \code{read_copc()}, with elements
#'   \code{$data} (a \code{data.table}) and \code{$header} (a named list
#'   following rlas header conventions).
#'
#' @return A \code{lidR::LAS} object.
#'
#' @details
#' This function requires \pkg{lidR} to be installed (it is listed in
#' Suggests, not Imports).
#'
#' The header list must contain at minimum the rlas-standard fields
#' (Version Major/Minor, Point Data Format ID, scale factors, offsets,
#' extents).  \code{read_copc()} already produces a compatible header.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("lidR", quietly = TRUE)) {
#'   f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
#'                    package = "copc4R")
#'   result <- read_copc(f, max_depth = 2L, progress = FALSE)
#'   las <- as_las(result)
#'   print(las)
#' }
#' }
#'
#' @export
as_las <- function(x) {
  if (!requireNamespace("lidR", quietly = TRUE)) {
    stop("Package 'lidR' is required for as_las(). ",
         "Install it with: install.packages('lidR')",
         call. = FALSE)
  }

  if (!is.list(x) || is.null(x$data) || is.null(x$header)) {
    stop("'x' must be a list with $data and $header, as returned by read_copc().",
         call. = FALSE)
  }

  dt <- x$data
  hdr_list <- x$header

  # -- Ensure column types match lidR expectations ----------------------
  # lidR expects X/Y/Z as numeric; Intensity/ReturnNumber/etc. as integer;
  # flag fields as logical.
  num_cols <- c("X", "Y", "Z", "gpstime")
  int_cols <- c("Intensity", "ReturnNumber", "NumberOfReturns",
                "Classification", "UserData", "PointSourceID",
                "ScanDirectionFlag", "EdgeOfFlightline",
                "ScanAngleRank", "ScanAngle",
                "ScannerChannel",
                "R", "G", "B", "NIR")
  lgl_cols <- c("Synthetic_flag", "Keypoint_flag", "Withheld_flag",
                "Overlap_flag")

  for (col in num_cols) {
    if (col %in% names(dt) && !is.numeric(dt[[col]])) {
      data.table::set(dt, j = col, value = as.numeric(dt[[col]]))
    }
  }
  for (col in int_cols) {
    if (col %in% names(dt) && !is.integer(dt[[col]])) {
      data.table::set(dt, j = col, value = as.integer(dt[[col]]))
    }
  }
  for (col in lgl_cols) {
    if (col %in% names(dt) && !is.logical(dt[[col]])) {
      data.table::set(dt, j = col, value = as.logical(dt[[col]]))
    }
  }

  # -- Build LASheader from the rlas-style list -------------------------
  las_header <- lidR::LASheader(hdr_list)

  # -- Construct LAS object ---------------------------------------------
  las <- lidR::LAS(dt, las_header, check = FALSE)

  las
}
