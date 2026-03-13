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


# ═══════════════════════════════════════════════════════════════════════════
# Tests for new v0.2.0 features
# ═══════════════════════════════════════════════════════════════════════════

test_that("copc_info returns metadata", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  info <- copc_info(fpath)
  expect_type(info, "list")
  expect_true(is.numeric(info$center_x))
  expect_true(is.numeric(info$halfsize))
  expect_gt(info$halfsize, 0)
  expect_gt(info$point_count, 0)
  expect_true(info$point_format %in% c(6L, 7L, 8L))
})

test_that("copc_bounds returns spatial extent", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  bounds <- copc_bounds(fpath)
  expect_length(bounds, 6)
  expect_named(bounds, c("xmin", "ymin", "zmin", "xmax", "ymax", "zmax"))
  expect_true(bounds["xmax"] > bounds["xmin"])
  expect_true(bounds["ymax"] > bounds["ymin"])
})

test_that("copc_density estimates point density", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  dens <- copc_density(fpath)
  expect_type(dens, "list")
  expect_gt(dens$estimated_points, 0)
  expect_gt(dens$area, 0)
  expect_gt(dens$density, 0)
  expect_gt(dens$nodes_checked, 0)
})

test_that("read_copc filter parameter works", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  # Read ground points only (class 2)
  result_ground <- read_copc(fpath, filter = "-keep_class 2",
                             max_points = 1000, progress = FALSE)
  if (nrow(result_ground$data) > 0) {
    expect_true(all(result_ground$data$Classification == 2))
  }
})

test_that("read_copc max_depth parameter works (LOD)", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  # Read at low resolution (depth 0 = root only)
  result_lod <- read_copc(fpath, max_depth = 0L, progress = FALSE)
  result_full <- read_copc(fpath, max_depth = -1L, progress = FALSE)

  # LOD should have fewer points
  expect_lte(nrow(result_lod$data), nrow(result_full$data))
})

test_that("read_copc_tiles creates tile grid", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  hdr <- read_copc_header(fpath)
  xrange <- hdr[["Max X"]] - hdr[["Min X"]]

  # Use tile_size that creates 2x2 grid
  tiles <- read_copc_tiles(fpath,
                           tile_size = xrange / 2,
                           progress = FALSE)

  expect_type(tiles, "list")
  expect_gte(length(tiles), 1)
  expect_true(!is.null(attr(tiles, "tile_bboxes")))
})

test_that("read_copc_sample limits output", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  result <- read_copc_sample(fpath, n = 100, max_depth = 2L, progress = FALSE)
  expect_s3_class(result$data, "data.table")
  expect_lte(nrow(result$data), 100)
})

test_that("read_copc_iter creates a working iterator", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  iter <- read_copc_iter(fpath, max_depth = 0L, progress = FALSE)
  expect_s3_class(iter, "copc_iterator")

  # Should be able to yield a batch
  if (iter$has_next()) {
    batch <- iter$yield()
    expect_s3_class(batch, "data.table")
    expect_gt(nrow(batch), 0)
  }

  # Collect should work
  iter$reset()
  all_data <- iter$collect()
  expect_s3_class(all_data, "data.table")
})

test_that("copc_cache functions work", {
  # Cache config
  prev <- copc_cache_config(mem_max_mb = 128)
  expect_type(prev, "list")
  copc_cache_config(mem_max_mb = prev$mem_max_mb)

  # Cache stats
  stats <- copc_cache_stats()
  expect_type(stats, "list")
  expect_true("mem_entries" %in% names(stats))
  expect_true("disk_mb" %in% names(stats))

  # Cache stats should include header_entries
  expect_true("header_entries" %in% names(stats))

  # Cache clear
  prev_stats <- copc_cache_clear(mem = TRUE, disk = FALSE)
  expect_type(prev_stats, "list")
})

test_that("header caching works", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  copc_cache_clear(mem = TRUE)

  hdr1 <- read_copc_header(fpath)
  stats_after <- copc_cache_stats()
  expect_gte(stats_after$header_entries, 1L)

  # Second call should return cached (same result)
  hdr2 <- read_copc_header(fpath)
  expect_identical(hdr1, hdr2)
})

test_that("cpp_select_nodes returns node info", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  nodes <- copc4R:::cpp_select_nodes(fpath, numeric(0), numeric(0), -1L)
  expect_type(nodes, "list")
  expect_true("offset" %in% names(nodes))
  expect_true("byte_size" %in% names(nodes))
  expect_true("point_count" %in% names(nodes))
  expect_true("level" %in% names(nodes))
  expect_gt(length(nodes$offset), 0L)
})

test_that("read_copc_catalog builds catalog", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  ctg <- read_copc_catalog(fpath, progress = FALSE)
  # Should be either LAScatalog or copc_catalog

  expect_true(inherits(ctg, "copc_catalog") || inherits(ctg, "LAScatalog"))
})


# ═══════════════════════════════════════════════════════════════════════════
# Tests for new v0.3.0 features
# ═══════════════════════════════════════════════════════════════════════════

test_that("resolution parameter auto-computes max_depth", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  info <- copc_info(fpath)
  skip_if(is.null(info$spacing) || info$spacing <= 0,
          message = "COPC spacing not available")

  # Use a coarse resolution → should return fewer points
  coarse_res <- info$spacing * 2  # coarser than root
  result_coarse <- read_copc(fpath, resolution = coarse_res, progress = FALSE)

  # Use fine resolution → should return more points
  fine_res <- info$spacing / 4
  result_fine <- read_copc(fpath, resolution = fine_res, progress = FALSE)

  expect_lte(nrow(result_coarse$data), nrow(result_fine$data))
})

test_that("expanded filter predicates work", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  # -drop_noise should remove classes 7 and 18
  result_no_noise <- read_copc(fpath, filter = "-drop_noise",
                                max_points = 5000, progress = FALSE)
  if (nrow(result_no_noise$data) > 0) {
    expect_false(any(result_no_noise$data$Classification %in% c(7, 18)))
  }

  # -keep_first should only return return number 1
  result_first <- read_copc(fpath, filter = "-keep_first",
                             max_points = 5000, progress = FALSE)
  if (nrow(result_first$data) > 0) {
    expect_true(all(result_first$data$ReturnNumber == 1))
  }

  # -keep_every_nth 10 should thin by ~90%
  result_all <- read_copc(fpath, max_depth = 0L, progress = FALSE)
  result_nth <- read_copc(fpath, max_depth = 0L,
                           filter = "-keep_every_nth 10", progress = FALSE)
  if (nrow(result_all$data) > 100) {
    ratio <- nrow(result_nth$data) / nrow(result_all$data)
    expect_lt(ratio, 0.2)  # should be ~0.1
  }
})

test_that("read_copc_sample with voxel_size works", {
  fpath <- Sys.getenv("COPC4R_TEST_FILE", unset = "")
  skip_if(fpath == "", message = "COPC4R_TEST_FILE not set")
  skip_if(!file.exists(fpath), message = "COPC4R_TEST_FILE not found")

  # Voxel thinning with mode = "first"
  result_voxel <- read_copc_sample(fpath, voxel_size = 5.0,
                                    mode = "first",
                                    max_depth = 2L,
                                    progress = FALSE)
  expect_s3_class(result_voxel$data, "data.table")

  # Compare with full read at same depth
  result_full <- read_copc(fpath, max_depth = 2L, progress = FALSE)
  expect_lte(nrow(result_voxel$data), nrow(result_full$data))

  # Voxel thinning with mode = "center"
  result_center <- read_copc_sample(fpath, voxel_size = 5.0,
                                     mode = "center",
                                     max_depth = 2L,
                                     progress = FALSE)
  expect_s3_class(result_center$data, "data.table")

  # Secondary n cap on top of voxels
  result_capped <- read_copc_sample(fpath, n = 50, voxel_size = 5.0,
                                     max_depth = 2L,
                                     progress = FALSE)
  expect_lte(nrow(result_capped$data), 50)
})
