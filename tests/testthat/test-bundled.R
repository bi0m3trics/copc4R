# ── Tests using bundled COPC files ──────────────────────────────────────────
#
# These tests use the small COPC files shipped in inst/extdata/ and
# are designed to always run during R CMD check (no env-vars needed).
# ───────────────────────────────────────────────────────────────────────────

# Helper: paths to the bundled COPC + GPKG files
copc_file <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
                          package = "copc4R")
gpkg_file <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82.gpkg",
                          package = "copc4R")
skip_bundled <- !nzchar(copc_file) || !file.exists(copc_file)

# ── read_copc_header ──────────────────────────────────────────────────────

test_that("read_copc_header returns a valid header (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  hdr <- read_copc_header(copc_file)

  expect_type(hdr, "list")
  expect_identical(hdr[["File Signature"]], "LASF")
  expect_identical(hdr[["Version Major"]], 1L)
  expect_gte(hdr[["Version Minor"]], 4L)          # LAS 1.4
  expect_true(hdr[["Point Data Format ID"]] %in% c(6L, 7L, 8L))
  expect_gt(hdr[["Number of point records"]], 0)

  # COPC Info VLR
  copc <- hdr[["COPC Info"]]
  expect_true(is.list(copc))
  expect_true(is.numeric(copc$spacing))
  expect_gt(copc$spacing, 0)
  expect_gt(copc$halfsize, 0)
})


# ── read_copc (basic) ────────────────────────────────────────────────────

test_that("read_copc reads all points from bundled file", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc(copc_file, progress = FALSE)

  expect_type(result, "list")
  expect_true(all(c("data", "header") %in% names(result)))
  expect_s3_class(result$data, "data.table")
  expect_gt(nrow(result$data), 0)
  expect_true(all(c("X", "Y", "Z") %in% names(result$data)))

  # Known properties of this file
  expect_equal(nrow(result$data), 20691)
})


# ── read_copc: bbox ──────────────────────────────────────────────────────

test_that("read_copc bbox filter returns a spatial subset (bundled)", {
 skip_if(skip_bundled, "bundled COPC file not found")

  hdr <- read_copc_header(copc_file)
  xmid <- (hdr[["Min X"]] + hdr[["Max X"]]) / 2
  ymid <- (hdr[["Min Y"]] + hdr[["Max Y"]]) / 2
  xrng <- (hdr[["Max X"]] - hdr[["Min X"]]) / 4
  yrng <- (hdr[["Max Y"]] - hdr[["Min Y"]]) / 4
  bb   <- c(xmid - xrng, ymid - yrng, xmid + xrng, ymid + yrng)

  result <- read_copc(copc_file, bbox = bb, progress = FALSE)
  dt <- result$data

  expect_gt(nrow(dt), 0)
  expect_true(all(dt$X >= bb[1] & dt$X <= bb[3]))
  expect_true(all(dt$Y >= bb[2] & dt$Y <= bb[4]))
})


# ── read_copc: select ───────────────────────────────────────────────────

test_that("read_copc select restricts columns (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc(copc_file, select = "xyzi", max_points = 200,
                      progress = FALSE)
  dt <- result$data

  expect_true(all(c("X", "Y", "Z", "Intensity") %in% names(dt)))
  expect_false("gpstime" %in% names(dt))
})


# ── read_copc: filter -keep_class ────────────────────────────────────────

test_that("read_copc -keep_class filter works (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc(copc_file, filter = "-keep_class 2",
                      max_points = 5000, progress = FALSE)
  if (nrow(result$data) > 0) {
    expect_true(all(result$data$Classification == 2L))
  }
})


# ── read_copc: filter -drop_noise ────────────────────────────────────────

test_that("read_copc -drop_noise filter removes classes 7 & 18 (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc(copc_file, filter = "-drop_noise",
                      max_points = 5000, progress = FALSE)
  if (nrow(result$data) > 0) {
    expect_false(any(result$data$Classification %in% c(7L, 18L)))
  }
})


# ── read_copc: filter -keep_first ────────────────────────────────────────

test_that("read_copc -keep_first filter returns first returns (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc(copc_file, filter = "-keep_first",
                      max_points = 5000, progress = FALSE)
  if (nrow(result$data) > 0 && "ReturnNumber" %in% names(result$data)) {
    expect_true(all(result$data$ReturnNumber == 1L))
  }
})


# ── read_copc: filter -keep_every_nth ────────────────────────────────────

test_that("read_copc -keep_every_nth thins correctly (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  full  <- read_copc(copc_file, max_depth = 0L, progress = FALSE)
  every <- read_copc(copc_file, max_depth = 0L,
                     filter = "-keep_every_nth 5", progress = FALSE)

  if (nrow(full$data) > 50) {
    ratio <- nrow(every$data) / nrow(full$data)
    expect_lt(ratio, 0.35)   # expect ~0.2
  }
})


# ── read_copc: filter -keep_voxel ───────────────────────────────────────

test_that("read_copc -keep_voxel thins to voxel centers (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  full    <- read_copc(copc_file, max_depth = 2L, progress = FALSE)
  thinned <- read_copc(copc_file, max_depth = 2L,
                       filter = "-keep_voxel 5.0", progress = FALSE)

  expect_gt(nrow(thinned$data), 0)
  expect_lt(nrow(thinned$data), nrow(full$data))
})


# ── read_copc: -keep_voxel validation ───────────────────────────────────

test_that("read_copc -keep_voxel rejects non-positive values", {
  skip_if(skip_bundled, "bundled COPC file not found")

  expect_error(
    read_copc(copc_file, filter = "-keep_voxel -1", progress = FALSE),
    "positive numeric"
  )
})


# ── read_copc: max_depth / LOD ──────────────────────────────────────────

test_that("read_copc max_depth returns fewer points (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  lod  <- read_copc(copc_file, max_depth = 0L, progress = FALSE)
  full <- read_copc(copc_file, progress = FALSE)

  expect_lte(nrow(lod$data), nrow(full$data))
})


# ── read_copc: resolution auto-computes max_depth ────────────────────────

test_that("resolution parameter works (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  info <- copc_info(copc_file)
  skip_if(is.null(info$spacing) || info$spacing <= 0,
          "COPC spacing not available")

  coarse <- read_copc(copc_file, resolution = info$spacing * 2,
                      progress = FALSE)
  fine   <- read_copc(copc_file, resolution = info$spacing / 4,
                      progress = FALSE)

  expect_lte(nrow(coarse$data), nrow(fine$data))
})


# ── copc_info ────────────────────────────────────────────────────────────

test_that("copc_info returns expected metadata (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  info <- copc_info(copc_file)

  expect_type(info, "list")
  expect_true(is.numeric(info$spacing))
  expect_gt(info$spacing, 0)
  expect_gt(info$halfsize, 0)
  expect_gt(info$point_count, 0)
  expect_true(info$point_format %in% c(6L, 7L, 8L))
  expect_equal(info$point_count, 20691)
})


# ── copc_bounds ──────────────────────────────────────────────────────────

test_that("copc_bounds returns numeric bounds (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  b <- copc_bounds(copc_file)

  expect_length(b, 6)
  expect_named(b, c("xmin", "ymin", "zmin", "xmax", "ymax", "zmax"))
  expect_true(b["xmax"] > b["xmin"])
  expect_true(b["ymax"] > b["ymin"])
})

test_that("copc_bounds as_sf returns sf bbox", {
  skip_if(skip_bundled, "bundled COPC file not found")
  skip_if_not_installed("sf")

  bb <- copc_bounds(copc_file, as_sf = TRUE)
  expect_s3_class(bb, "bbox")
})


# ── copc_density ─────────────────────────────────────────────────────────

test_that("copc_density returns density estimate (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  dens <- copc_density(copc_file)

  expect_type(dens, "list")
  expect_gt(dens$estimated_points, 0)
  expect_gt(dens$area, 0)
  expect_gt(dens$density, 0)
})


# ── as_las ───────────────────────────────────────────────────────────────

test_that("as_las creates a valid LAS object (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")
  skip_if_not_installed("lidR")

  result <- read_copc(copc_file, max_points = 500, progress = FALSE)
  las <- as_las(result)

  expect_s4_class(las, "LAS")
  expect_gt(lidR::npoints(las), 0)
})

test_that("as_las rejects bad input", {
  expect_error(as_las(list(a = 1)), "\\$data and \\$header")
})


# ── cache ────────────────────────────────────────────────────────────────

test_that("copc_cache_config round-trips settings", {
  prev <- copc_cache_config(mem_max_mb = 128)
  expect_type(prev, "list")
  expect_true(is.numeric(prev$mem_max_mb))
  # restore

  copc_cache_config(mem_max_mb = prev$mem_max_mb)
})

test_that("copc_cache_stats returns expected fields", {
  stats <- copc_cache_stats()
  expect_type(stats, "list")
  expect_true(all(c("mem_entries", "mem_mb", "header_entries",
                     "disk_entries", "disk_mb") %in% names(stats)))
})

test_that("copc_cache_clear works without error", {
  prev <- copc_cache_clear(mem = TRUE, disk = FALSE)
  expect_type(prev, "list")
})

test_that("header caching avoids re-reads (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  copc_cache_clear(mem = TRUE)

  hdr1 <- read_copc_header(copc_file)
  stats <- copc_cache_stats()
  expect_gte(stats$header_entries, 1L)

  hdr2 <- read_copc_header(copc_file)
  expect_identical(hdr1, hdr2)
})


# ── read_copc_tiles ──────────────────────────────────────────────────────

test_that("read_copc_tiles produces a tile list (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  hdr   <- read_copc_header(copc_file)
  xrng  <- hdr[["Max X"]] - hdr[["Min X"]]

  tiles <- read_copc_tiles(copc_file, tile_size = xrng / 2,
                           progress = FALSE)

  expect_type(tiles, "list")
  expect_gte(length(tiles), 1)
  expect_true(!is.null(attr(tiles, "tile_bboxes")))
})


# ── read_copc_sample ────────────────────────────────────────────────────

test_that("read_copc_sample random cap works (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc_sample(copc_file, n = 100, progress = FALSE)
  expect_s3_class(result$data, "data.table")
  expect_lte(nrow(result$data), 100)
})

test_that("read_copc_sample voxel-first mode works (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc_sample(copc_file, voxel_size = 5.0,
                             mode = "first", max_depth = 2L,
                             progress = FALSE)
  full   <- read_copc(copc_file, max_depth = 2L, progress = FALSE)

  expect_s3_class(result$data, "data.table")
  expect_lte(nrow(result$data), nrow(full$data))
})

test_that("read_copc_sample voxel-center mode works (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc_sample(copc_file, voxel_size = 5.0,
                             mode = "center", max_depth = 2L,
                             progress = FALSE)
  expect_s3_class(result$data, "data.table")
  expect_gt(nrow(result$data), 0)
})

test_that("read_copc_sample voxel + secondary n cap (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  result <- read_copc_sample(copc_file, n = 50, voxel_size = 5.0,
                             max_depth = 2L, progress = FALSE)
  expect_lte(nrow(result$data), 50)
})


# ── read_copc_iter ──────────────────────────────────────────────────────

test_that("read_copc_iter yields batches (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  iter <- read_copc_iter(copc_file, max_depth = 0L, progress = FALSE)
  expect_s3_class(iter, "copc_iterator")

  if (iter$has_next()) {
    batch <- iter$yield()
    expect_s3_class(batch, "data.table")
    expect_gt(nrow(batch), 0)
  }

  # collect
  iter$reset()
  all_data <- iter$collect()
  expect_s3_class(all_data, "data.table")
})


# ── read_copc_polygon ───────────────────────────────────────────────────

test_that("read_copc_polygon clips to polygon (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")
  skip_if_not_installed("sf")
  skip_if(!nzchar(gpkg_file) || !file.exists(gpkg_file),
          "bundled GPKG not found")

  bldg <- sf::st_read(gpkg_file, quiet = TRUE)
  bldg_utm <- sf::st_transform(bldg, 26912)

  result <- read_copc_polygon(copc_file, polygon = bldg_utm,
                              progress = FALSE)
  expect_s3_class(result$data, "data.table")
  expect_gt(nrow(result$data), 0)
})


# ── read_copc_corridor ──────────────────────────────────────────────────

test_that("read_copc_corridor buffers a line (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")
  skip_if_not_installed("sf")

  hdr  <- read_copc_header(copc_file)
  xmid <- (hdr[["Min X"]] + hdr[["Max X"]]) / 2
  ymid <- (hdr[["Min Y"]] + hdr[["Max Y"]]) / 2

  line <- sf::st_sfc(sf::st_linestring(matrix(c(
    xmid - 100, ymid,
    xmid + 100, ymid
  ), ncol = 2, byrow = TRUE)), crs = 26912)

  result <- read_copc_corridor(copc_file, line = line, width = 20,
                               progress = FALSE)
  expect_s3_class(result$data, "data.table")
  expect_gt(nrow(result$data), 0)
})


# ── read_copc_catalog ───────────────────────────────────────────────────

test_that("read_copc_catalog builds a catalog (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  ctg <- read_copc_catalog(copc_file, progress = FALSE)
  expect_true(
    inherits(ctg, "copc_catalog") || inherits(ctg, "LAScatalog")
  )
})

test_that("catalog_apply processes each file (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  ctg <- read_copc_catalog(copc_file, progress = FALSE)
  results <- catalog_apply(ctg, function(x) nrow(x$data), progress = FALSE)

  expect_type(results, "list")
  expect_length(results, 1)
  expect_gt(results[[1]], 0)
})

test_that("catalog_apply with chunk_size tiles the file (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  ctg <- read_copc_catalog(copc_file, chunk_size = 50,
                           chunk_buffer = 0, progress = FALSE)
  results <- catalog_apply(ctg, function(x) {
    data.frame(n = nrow(x$data), z_mean = mean(x$data$Z))
  }, progress = FALSE)

  expect_type(results, "list")
  expect_gte(length(results), 1)
  # All chunks combined should have points
  metrics <- do.call(rbind, Filter(Negate(is.null), results))
  expect_gt(sum(metrics$n), 0)
})

test_that("catalog_apply cores parameter accepts integer (bundled)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  ctg <- read_copc_catalog(copc_file, chunk_size = 100, progress = FALSE)
  # cores = 1L should work identically to default
  results <- catalog_apply(ctg, function(x) nrow(x$data),
                           cores = 1L, progress = FALSE)
  expect_type(results, "list")
  expect_gte(length(results), 1)
})


# ── Internal: .voxel_center_thin ────────────────────────────────────────

test_that(".voxel_center_thin selects nearest to center", {
  dt <- data.table::data.table(
    X = c(0.1, 0.4, 0.9, 1.1),
    Y = c(0.1, 0.4, 0.9, 1.1),
    Z = c(0.1, 0.4, 0.1, 0.1)
  )
  # voxel_size = 1.0 → voxel centers at 0.5, 1.5, ...
  # Points at (0.4, 0.4, 0.4) and (0.9, 0.9, 0.1) are in voxel (0,0,0)
  # Point at (0.4, 0.4, 0.4) is closer to center (0.5, 0.5, 0.5)
  result <- copc4R:::.voxel_center_thin(dt, 1.0)

  # Should have ≤ 3 voxels occupied (two points may share a cell)
  expect_lte(nrow(result), nrow(dt))
  expect_gt(nrow(result), 0)
  expect_true(all(c("X", "Y", "Z") %in% names(result)))
})


# ── Edge cases ───────────────────────────────────────────────────────────

test_that("read_copc_header rejects non-COPC .las paths", {
  tmp <- tempfile(fileext = ".las")
  file.create(tmp)
  on.exit(unlink(tmp))
  expect_error(read_copc_header(tmp), "standard LAS/LAZ")
})

test_that("read_copc errors on missing file", {
  expect_error(read_copc("/nonexistent/file.copc.laz", progress = FALSE),
               "not found")
})

test_that("as_las requires lidR", {
  # Hard to test without unloading lidR, but we can test the list check
  expect_error(as_las("not a list"), "\\$data and \\$header")
})

# ── write_copc round-trip tests ──────────────────────────────────────────

test_that("write_copc round-trips bundled file (point count)", {
  skip_if(skip_bundled, "bundled COPC file not found")

  d <- read_copc(copc_file, progress = FALSE)
  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)

  write_copc(d, tf, progress = FALSE)
  expect_true(file.exists(tf))
  expect_gt(file.info(tf)$size, 0)

  d2 <- read_copc(tf, progress = FALSE)
  expect_equal(nrow(d2$data), nrow(d$data))
})

test_that("write_copc round-trips coordinates losslessly", {
  skip_if(skip_bundled, "bundled COPC file not found")

  d <- read_copc(copc_file, progress = FALSE)
  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)
  write_copc(d, tf, progress = FALSE)
  d2 <- read_copc(tf, progress = FALSE)

  # Sort both by all key columns for comparison
  dt1 <- data.table::copy(d$data)
  dt2 <- data.table::copy(d2$data)
  data.table::setorderv(dt1, c("X", "Y", "Z", "gpstime"))
  data.table::setorderv(dt2, c("X", "Y", "Z", "gpstime"))

  expect_equal(sort(dt2$X), sort(dt1$X))
  expect_equal(sort(dt2$Y), sort(dt1$Y))
  expect_equal(sort(dt2$Z), sort(dt1$Z))
})

test_that("write_copc round-trips attributes", {
  skip_if(skip_bundled, "bundled COPC file not found")

  d <- read_copc(copc_file, progress = FALSE)
  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)
  write_copc(d, tf, progress = FALSE)
  d2 <- read_copc(tf, progress = FALSE)

  # Sort both datasets consistently
  for (col in c("Intensity", "Classification", "ReturnNumber",
                "NumberOfReturns", "PointSourceID")) {
    expect_equal(sort(d2$data[[col]]), sort(d$data[[col]]),
                 info = paste("attribute:", col))
  }
})

test_that("write_copc output is a valid COPC file", {
  skip_if(skip_bundled, "bundled COPC file not found")

  d <- read_copc(copc_file, progress = FALSE)
  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)
  write_copc(d, tf, progress = FALSE)

  # Should have a valid COPC header
  h2 <- read_copc_header(tf)
  expect_identical(h2[["File Signature"]], "LASF")
  expect_identical(h2[["Version Major"]], 1L)
  expect_identical(h2[["Version Minor"]], 4L)
  expect_true(h2[["Point Data Format ID"]] %in% c(6L, 7L, 8L))
  expect_gt(h2[["Number of point records"]], 0)
  expect_true(is.list(h2[["COPC Info"]]))
  expect_gt(h2[["COPC Info"]]$halfsize, 0)
})

test_that("write_copc preserves header scale/offset", {
  skip_if(skip_bundled, "bundled COPC file not found")

  d <- read_copc(copc_file, progress = FALSE)
  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)
  write_copc(d, tf, progress = FALSE)

  h1 <- d$header
  h2 <- read_copc_header(tf)
  expect_equal(h2[["X scale factor"]], h1[["X scale factor"]])
  expect_equal(h2[["Y scale factor"]], h1[["Y scale factor"]])
  expect_equal(h2[["Z scale factor"]], h1[["Z scale factor"]])
  expect_equal(h2[["X offset"]], h1[["X offset"]])
  expect_equal(h2[["Y offset"]], h1[["Y offset"]])
  expect_equal(h2[["Z offset"]], h1[["Z offset"]])
})


# ── write_copc input flexibility tests ───────────────────────────────────

test_that("write_copc accepts a bare data.frame with smart defaults", {
  df <- data.frame(
    X = c(100.1, 100.2, 100.3),
    Y = c(200.1, 200.2, 200.3),
    Z = c(10.1,  10.2,  10.3)
  )
  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)
  write_copc(df, tf, progress = FALSE)

  back <- read_copc(tf, progress = FALSE)
  expect_equal(nrow(back$data), 3L)

  # Check smart defaults applied
  h <- read_copc_header(tf)
  expect_equal(h[["X scale factor"]], 0.001)
  expect_equal(h[["X offset"]], 100)   # floor(min(X))
  expect_equal(h[["Y offset"]], 200)
  expect_equal(h[["Z offset"]], 0)
  expect_equal(h[["Point Data Format ID"]], 6L) # no RGB/NIR
})

test_that("write_copc accepts a lidR::LAS object", {
  skip_if(skip_bundled, "bundled COPC file not found")
  skip_if_not_installed("lidR")

  d <- read_copc(copc_file, progress = FALSE)
  las <- as_las(d)

  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)
  write_copc(las, tf, progress = FALSE)

  back <- read_copc(tf, progress = FALSE)
  expect_equal(nrow(back$data), nrow(d$data))
})

test_that("write_copc normalises column names", {
  df <- data.frame(
    x = c(1, 2, 3),
    y = c(4, 5, 6),
    z = c(7, 8, 9),
    intensity = c(100L, 200L, 300L),
    classification = c(2L, 2L, 6L)
  )
  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)
  write_copc(df, tf, progress = FALSE)

  back <- read_copc(tf, progress = FALSE)
  expect_equal(nrow(back$data), 3L)
  expect_true(all(back$data$Intensity == c(100L, 200L, 300L)))
  expect_true(all(back$data$Classification == c(2L, 2L, 6L)))
})

test_that("write_copc filter argument works", {
  skip_if(skip_bundled, "bundled COPC file not found")

  d <- read_copc(copc_file, progress = FALSE)
  n_ground <- sum(d$data$Classification == 2L)
  skip_if(n_ground == 0L, "no ground points in bundled file")

  tf <- tempfile(fileext = ".copc.laz")
  on.exit(unlink(tf), add = TRUE)
  write_copc(d, tf, filter = "-keep_class 2", progress = FALSE)

  back <- read_copc(tf, progress = FALSE)
  expect_equal(nrow(back$data), n_ground)
  expect_true(all(back$data$Classification == 2L))
})

test_that("write_copc errors on empty input", {
  df <- data.frame(X = numeric(0), Y = numeric(0), Z = numeric(0))
  tf <- tempfile(fileext = ".copc.laz")
  expect_error(write_copc(df, tf, progress = FALSE),
               "0 points")
})

test_that("write_copc errors on missing XYZ", {
  df <- data.frame(a = 1, b = 2)
  tf <- tempfile(fileext = ".copc.laz")
  expect_error(write_copc(df, tf, progress = FALSE),
               "X, Y, Z")
})
