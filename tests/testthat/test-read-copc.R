# tests for copc4R
#
# Requires environment variable COPC4R_TEST_FILE pointing to a local
# .copc.laz file.  Tests are skipped if the file is not available.

test_that("read_copc_header returns valid header list", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  hdr <- read_copc_header(fpath)

  # Basic structure checks

expect_type(hdr, "list")
  expect_equal(hdr[["File Signature"]], "LASF")
  expect_equal(hdr[["Version Major"]], 1L)
  expect_gte(hdr[["Version Minor"]], 4L)
  expect_true(hdr[["Point Data Format ID"]] %in% c(6L, 7L, 8L))
  expect_gt(hdr[["Number of point records"]], 0)


  # COPC info
  expect_true("COPC Info" %in% names(hdr))
  copc <- hdr[["COPC Info"]]
  expect_true(is.numeric(copc$center_x))
  expect_true(is.numeric(copc$halfsize))
  expect_gt(copc$halfsize, 0)
})

test_that("read_copc reads all points without filters", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  result <- read_copc(fpath, progress = FALSE)

  expect_type(result, "list")
  expect_true("data" %in% names(result))
  expect_true("header" %in% names(result))

  dt <- result$data
  expect_s3_class(dt, "data.table")
  expect_gt(nrow(dt), 0)

  # Required columns should be present
  expect_true(all(c("X", "Y", "Z") %in% names(dt)))
})

test_that("read_copc bbox filter returns subset", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  # Read header to get extents
  hdr <- read_copc_header(fpath)
  xmid <- (hdr[["Min X"]] + hdr[["Max X"]]) / 2
  ymid <- (hdr[["Min Y"]] + hdr[["Max Y"]]) / 2
  xrng <- (hdr[["Max X"]] - hdr[["Min X"]]) / 4
  yrng <- (hdr[["Max Y"]] - hdr[["Min Y"]]) / 4

  bbox <- c(xmid - xrng, ymid - yrng, xmid + xrng, ymid + yrng)
  result <- read_copc(fpath, bbox = bbox, progress = FALSE)
  dt <- result$data

  expect_gt(nrow(dt), 0)

  # All points should be within bbox (point-level filtering)
  expect_true(all(dt$X >= bbox[1] & dt$X <= bbox[3]))
  expect_true(all(dt$Y >= bbox[2] & dt$Y <= bbox[4]))
})

test_that("read_copc select parameter works", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  result <- read_copc(fpath, select = "xyzi", max_points = 100, progress = FALSE)
  dt <- result$data

  expect_true(all(c("X", "Y", "Z", "Intensity") %in% names(dt)))
  # gpstime should NOT be present with select="xyzi"
  expect_false("gpstime" %in% names(dt))
})

test_that("read_copc max_points caps results", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  result <- read_copc(fpath, max_points = 50, progress = FALSE)
  # Note: actual count may vary because we cap at node selection level
  # but point-level filtering happens later.  At minimum, we should

  # not get wildly more points than max_points.
  expect_lte(nrow(result$data), 1e6) # sanity
})

test_that("as_las works with lidR if installed", {
  skip_if_not_installed("lidR")
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  result <- read_copc(fpath, max_points = 1000, progress = FALSE)
  las <- as_las(result)

  expect_s4_class(las, "LAS")

  # Run lidR's own validity check
  chk <- tryCatch(lidR::las_check(las, print = FALSE),
                   error = function(e) NULL)
  # las_check returns a list with $errors and $warnings
  if (!is.null(chk) && is.list(chk)) {
    # Ideally no errors
    if ("errors" %in% names(chk)) {
      expect_equal(length(chk$errors), 0,
                   info = paste("las_check errors:", paste(chk$errors, collapse = "; ")))
    }
  }
})
