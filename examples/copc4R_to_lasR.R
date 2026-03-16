# =============================================================================
# copc4R_to_lasR.R -- copc4R ➜ lasR / lidR processing handoff
# =============================================================================
#
# Demonstrates how to use copc4R for cloud-native data discovery and I/O,
# then hand off to lasR (or lidR) for point-cloud processing.
#
# The VPC format is the glue:
#   copc4R writes spec-conformant VPCs that lasR (and QGIS) can consume.
#
# Workflow overview:
#   1. Discover and fetch remote COPC tiles (copc4R)
#   2. Write a VPC index (copc4R)
#   3. Process with lasR pipelines -or- use as_las() + lidR
#
# Required packages:
#   install.packages("copc4R")
#   install.packages("lasR", repos = "https://r-lidar.r-universe.dev")
#   install.packages(c("lidR", "sf", "terra"))
# =============================================================================

library(copc4R)

cat("===================================================================\n")
cat(" copc4R ➜ lasR / lidR handoff demo\n")
cat("===================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Data discovery: use copc4R to find and inspect remote COPC tiles
# ---------------------------------------------------------------------------
# copc4R is designed for cloud-native COPC access -- it reads only the
# octree hierarchy pages via HTTP range requests, so you never download
# full files just to peek at the metadata.
#
# For this example we use the bundled local files; in practice you'd
# point at remote URLs on S3, Azure Blob Storage, or Planetary Computer.
# ---------------------------------------------------------------------------

f <- system.file("extdata",
                 "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
                 package = "copc4R")
cat("Source file:", f, "\n\n")

# Quick metadata peek (no point decompression)
info <- copc_info(f)
cat("Points:", format(info$point_count, big.mark = ","), "\n")
cat("CRS:   ", if (nzchar(info$crs_wkt)) sf::st_crs(info$crs_wkt)$Name else "(none)", "\n\n")

# Octree-based point estimation (header + hierarchy only)
est <- copc_estimate_points(f)
cat(sprintf("Estimated points: %s  |  density: %.1f pts/m²\n\n",
            format(est$estimated_points, big.mark = ","),
            est$density))

# ---------------------------------------------------------------------------
# 2. Write a VPC index
# ---------------------------------------------------------------------------
# The VPC is a STAC ItemCollection (GeoJSON FeatureCollection) that
# every tool in the ecosystem understands: lasR, QGIS, PDAL wrench.
#
# copc4R's write_copc_vpc() is fully spec-conformant:
#   • 3-D bounding boxes
#   • proj:wkt2 + proj:projjson for robust CRS round-tripping
#   • pc:indexed = TRUE (COPC octree = built-in spatial index)
#   • Compound CRS handling
#   • Optional hierarchical overviews for QGIS seamless zoom
# ---------------------------------------------------------------------------

vpc_path <- file.path(tempdir(), "demo_handoff.vpc")
write_copc_vpc(f, vpc_path, progress = TRUE)

cat("\nVPC written:", vpc_path, "\n")
cat("Contents (first 30 lines):\n")
cat(readLines(vpc_path, n = 30), sep = "\n")
cat("\n...\n\n")

# ---------------------------------------------------------------------------
# 3a. Read the VPC back -- fast path (no header re-reads)
# ---------------------------------------------------------------------------
# With trust_vpc = TRUE (the default), read_copc_vpc() builds the catalog
# entirely from VPC metadata.  For large projects (1000s of tiles) this
# skips thousands of HTTP round-trips.

ctg <- read_copc_vpc(vpc_path, trust_vpc = TRUE, progress = TRUE)
cat("\nCatalog class:", class(ctg)[1], "\n")
print(ctg)

# ---------------------------------------------------------------------------
# 3b. Option A: hand off to lasR for pipeline processing
# ---------------------------------------------------------------------------
# lasR can consume VPC files directly as input to exec().
# This is the recommended path for heavy processing (DTM, CHM,
# classification, tree segmentation, ...).

if (requireNamespace("lasR", quietly = TRUE)) {
  cat("\n--- lasR pipeline processing via VPC ---\n")

  # lasR reads the .vpc file directly -- no need to go through copc4R
  pipeline <- lasR::reader() +
    lasR::rasterize(2, "zmax")

  # exec() on the VPC file
  result <- lasR::exec(pipeline, on = vpc_path)
  cat("lasR result class:", class(result), "\n")

  if (requireNamespace("terra", quietly = TRUE)) {
    cat("Raster dimensions:", paste(dim(result), collapse = " x "), "\n")
  }

  # More complex pipeline example:
  # pipeline <- lasR::classify_with_csf() +
  #   lasR::triangulate(filter = lasR::keep_ground()) +
  #   lasR::dtm(1) +
  #   lasR::normalize() +
  #   lasR::chm(0.5) +
  #   lasR::write_las(paste0(tempdir(), "/*_normalized.laz"))
  # lasR::exec(pipeline, on = vpc_path, ncores = 4)

} else {
  cat("\nlasR not installed -- skipping lasR pipeline demo.\n")
  cat("Install with:\n")
  cat("  install.packages('lasR', repos = 'https://r-lidar.r-universe.dev')\n\n")
}

# ---------------------------------------------------------------------------
# 3c. Option B: copc4R read + as_las() + lidR processing
# ---------------------------------------------------------------------------
# When you need fine-grained spatial queries (e.g., small AOI from a huge
# remote dataset), copc4R's range-read capability is unmatched.
# Read just the points you need, then convert to LAS for lidR.

if (requireNamespace("lidR", quietly = TRUE)) {
  cat("\n--- copc4R ➜ lidR processing ---\n")

  # Read a spatial subset via COPC octree query
  hdr <- read_copc_header(f)
  xmid <- (hdr[["Min X"]] + hdr[["Max X"]]) / 2
  ymid <- (hdr[["Min Y"]] + hdr[["Max Y"]]) / 2
  small_bbox <- c(xmid - 25, ymid - 25, xmid + 25, ymid + 25)

  result <- read_copc(f, bbox = small_bbox, progress = FALSE)
  cat(sprintf("Read %s points from bbox query\n",
              format(nrow(result$data), big.mark = ",")))

  # Convert to lidR LAS object
  las <- as_las(result)
  cat("LAS class:", class(las)[1], "\n")
  cat("LAS points:", lidR::npoints(las), "\n")

  # Now use any lidR function
  # las_ground <- lidR::classify_ground(las, lidR::csf())
  # dtm <- lidR::rasterize_terrain(las_ground, 1, lidR::tin())
  # chm <- lidR::rasterize_canopy(las, 0.5, lidR::p2r())

} else {
  cat("\nlidR not installed -- skipping lidR demo.\n")
}

# ---------------------------------------------------------------------------
# 4. Summary: when to use which tool
# ---------------------------------------------------------------------------
cat("\n===================================================================\n")
cat(" When to use which tool:\n")
cat("===================================================================\n\n")
cat(" copc4R  → Cloud-native data discovery, spatial queries, VPC I/O\n")
cat("            • HTTP range reads (fetch only what you need)\n")
cat("            • Octree estimation (no decompression)\n")
cat("            • Adaptive splitting for balanced workloads\n")
cat("            • Spec-conformant VPC writing\n\n")
cat(" lasR    → High-performance local/batch processing\n")
cat("            • Pipeline architecture (single-pass, C++ speed)\n")
cat("            • DTM, CHM, classification, segmentation\n")
cat("            • Reads VPC natively via exec(pipeline, on = 'file.vpc')\n\n")
cat(" lidR    → Interactive / exploratory analysis in R\n")
cat("            • Rich R-native API (LAS objects, ggplot, etc.)\n")
cat("            • Works with as_las() output from copc4R\n\n")
cat(" VPC     → The interoperability glue\n")
cat("            • STAC ItemCollection format\n")
cat("            • Understood by copc4R, lasR, QGIS, PDAL wrench\n")
cat("            • copc4R writes it; everyone else reads it\n")
cat("===================================================================\n")
