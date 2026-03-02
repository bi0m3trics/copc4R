# =============================================================================
# copc4R  –  Example: CHM workflow with lidR
# =============================================================================
#
# This example demonstrates:
#   1. read_copc_header()  – inspect a COPC file's metadata
#   2. read_copc()         – read point data with spatial filtering
#   3. as_las()            – convert to a lidR::LAS object
#   4. lidR processing     – classify_ground, normalize_height, rasterize_canopy
#
# NOTE: The USGS file at
#   https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/...
# is a regular LAZ, NOT a COPC file. copc4R only reads COPC files that
# contain the COPC octree hierarchy VLR (typically named *.copc.laz).
#
# We use the Autzen Stadium COPC from PDAL's public S3 bucket — it has
# classified vegetation and ground points ideal for a CHM workflow.
# =============================================================================

library(copc4R)
library(lidR)

# ── 1. Source: Autzen Stadium COPC from public S3 ─────────────────────────────
# read_copc() and read_copc_header() accept both local paths and URLs.
# Remote files are downloaded to a session-cached tempfile automatically.
#
# To read directly from S3:
  copc_file <- "https://s3.amazonaws.com/hobu-lidar/autzen-classified.copc.laz"
#
# For a local copy (if you've placed it in inst/extdata):
# copc_file <- system.file("extdata", "autzen-classified.copc.laz",
                         # package = "copc4R")
if (copc_file == "")
  copc_file <- "inst/extdata/autzen-classified.copc.laz"

# ── 2. Inspect the header ─────────────────────────────────────────────────────
hdr <- read_copc_header(copc_file)
cat("LAS Version :", hdr[["Version Major"]], ".", hdr[["Version Minor"]], "\n")
cat("Point Format:", hdr[["Point Data Format ID"]], "\n")
cat("Total Points:", format(hdr[["Number of point records"]], big.mark = ","), "\n")
cat("X extent    :", hdr[["Min X"]], "–", hdr[["Max X"]], "\n")
cat("Y extent    :", hdr[["Min Y"]], "–", hdr[["Max Y"]], "\n")
cat("Z extent    :", hdr[["Min Z"]], "–", hdr[["Max Z"]], "\n")

copc <- hdr[["COPC Info"]]
cat("\nCOPC Info:\n")
cat("  Center   :", copc$center_x, ",", copc$center_y, "\n")
cat("  Halfsize :", copc$halfsize, "m\n")
cat("  GPS range:", copc$gpstime_minimum, "–", copc$gpstime_maximum, "\n")

# ── 3. Read a spatial subset of points ─────────────────────────────────────────
# Pull a 500 ft × 500 ft window (the CRS units are US survey feet).
# select = "xyzicrnR" gives us XYZ, Intensity, Classification,
# ReturnNumber, NumberOfReturns, and RGB.
bbox <- c(636500, 849000, 637000, 849500)  # xmin, ymin, xmax, ymax
cat("\nReading points inside bbox:", bbox, "...\n")
result <- read_copc(copc_file,
                    bbox       = bbox,
                    select     = "xyzicrnR",
                    max_points = 500000,
                    progress   = TRUE)

cat("Points read:", format(nrow(result$data), big.mark = ","), "\n")
cat("Columns    :", paste(names(result$data), collapse = ", "), "\n")

# Quick summary
cat("\nClassification counts:\n")
print(table(result$data$Classification))

# ── 4. Convert to lidR::LAS ───────────────────────────────────────────────────
las <- as_las(result)
cat("\nLAS object:\n")
print(las)

# ── 5. Classify ground with CSF (Cloth Simulation Filter) ─────────────────────
cat("\nClassifying ground with CSF...\n")
las <- classify_ground(las, algorithm = csf(sloop_smooth    = FALSE,
                                            class_threshold = 0.5,
                                            cloth_resolution = 0.5,
                                            rigidness       = 1L,
                                            time_step       = 0.65))
cat("Ground points:", sum(las$Classification == 2L), "\n")

# ── 6. Normalize heights with TIN ─────────────────────────────────────────────
cat("Normalizing heights (TIN)...\n")
nlas <- normalize_height(las, algorithm = tin())
cat("Height range after normalisation:",
    round(min(nlas$Z), 2), "–", round(max(nlas$Z), 2), "m\n")

# ── 7. Rasterize canopy height model with pitfree() ──────────────────────────
cat("Building CHM (pitfree)...\n")
chm <- rasterize_canopy(nlas, res = 0.5, algorithm = pitfree(thresholds = c(0, 2, 5, 10, 15),
                                                             max_edge   = c(0, 1.5)))
cat("CHM raster :", terra::nrow(chm), "x", terra::ncol(chm),
    "cells at", terra::res(chm)[1], "m resolution\n")
cat("CHM range  :", round(terra::minmax(chm)[1], 2), "–",
    round(terra::minmax(chm)[2], 2), "m\n")

# ── 8. Display ────────────────────────────────────────────────────────────────
terra::plot(chm, col  = lidR::height.colors(50),
            main = "Canopy Height Model - Autzen Stadium (copc4R + lidR)")
