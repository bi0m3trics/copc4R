# =============================================================================
# vpc_demo.R -- Virtual Point Cloud (VPC) workflow demo
# =============================================================================
#
# Demonstrates the VPC (Virtual Point Cloud) support in copc4R.
#
# A VPC is a lightweight JSON file (STAC ItemCollection) that references
# multiple COPC tiles. QGIS and PDAL wrench use this as a "virtual layer"
# so the whole tiled dataset is treated as a single point cloud.
#
# This workflow:
#   1. Reads local COPC tile(s)
#   2. Writes a .vpc file (valid STAC ItemCollection)
#   3. Reads the .vpc back as a catalog
#   4. Processes the catalog with copc_apply()
#
# The resulting .vpc can be loaded in QGIS 3.32+ as a native point cloud
# layer -- zoom out to see the overview, zoom in for full-res tiles.
#
# Required packages:
#   install.packages(c("copc4R", "jsonlite", "lidR", "sf"))
# =============================================================================

library(copc4R)

cat("===================================================================\n")
cat(" copc4R -- Virtual Point Cloud (VPC) Workflow Demo\n")
cat("===================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Get example COPC file(s)
# ---------------------------------------------------------------------------
f <- system.file("extdata", "NoAZCampus_SWFSC_Bldg82_2mVoxel.copc.laz",
                 package = "copc4R")
cat("Source file:", f, "\n")

# In real-world usage you'd have many tiles:
# files <- list.files("data/tiles/", "\\.copc\\.laz$", full.names = TRUE)
files <- f

# ---------------------------------------------------------------------------
# 2. Write a VPC index
# ---------------------------------------------------------------------------
vpc_path <- file.path(tempdir(), "demo.vpc")

cat("\nWriting VPC to:", vpc_path, "\n")
write_copc_vpc(files, vpc_path, progress = TRUE)

# The .vpc file is valid GeoJSON -- inspect it
cat("\nVPC file contents:\n")
cat(readLines(vpc_path, n = 20), sep = "\n")
cat("\n...\n")

# ---------------------------------------------------------------------------
# 3. Read the VPC back as a catalog
# ---------------------------------------------------------------------------
cat("\nReading VPC back as catalog:\n")
ctg <- read_copc_vpc(vpc_path, progress = TRUE)
print(ctg)

# The VPC metadata is available as an attribute:
vpc_meta <- attr(ctg, "vpc_metadata")
cat("\nVPC source:", vpc_meta$source, "\n")
cat("Number of tiles:", vpc_meta$n_features, "\n")

# ---------------------------------------------------------------------------
# 4. Use octree-based point estimation for pre-flight planning
# ---------------------------------------------------------------------------
cat("\n--- Octree-based point count estimation ---\n")
est <- copc_estimate_points(f)
cat(sprintf("Full extent: ~%s points, density %.1f pts/m²\n",
            format(est$estimated_points, big.mark = ","),
            est$density))

# Multiple windows
hdr <- read_copc_header(f)
xmid <- (hdr[["Min X"]] + hdr[["Max X"]]) / 2
ymid <- (hdr[["Min Y"]] + hdr[["Max Y"]]) / 2
quads <- list(
  SW = c(hdr[["Min X"]], hdr[["Min Y"]], xmid, ymid),
  NE = c(xmid, ymid, hdr[["Max X"]], hdr[["Max Y"]])
)
est_q <- copc_estimate_points(f, quads)
cat("\nPer-quadrant estimates:\n")
print(est_q)

# ---------------------------------------------------------------------------
# 5. Adaptive splitting -- density-aware work units
# ---------------------------------------------------------------------------
cat("\n--- Adaptive splitting (a la Howard Butler's pattern) ---\n")
cat("Splitting file extent into tiles of <= 5000 estimated points...\n")
tiles <- copc_adaptive_split(f, max_points = 5000, collar = 5)
cat(sprintf("Result: %d adaptive tiles\n", nrow(tiles)))
print(tiles[, c("xmin", "ymin", "xmax", "ymax",
                "estimated_points", "depth")])

# ---------------------------------------------------------------------------
# 6. Process adaptive tiles (manual loop with collar)
# ---------------------------------------------------------------------------
cat("\n--- Processing adaptive tiles with collar ---\n")
results <- lapply(seq_len(nrow(tiles)), function(i) {
  # Read with the collar (expanded bbox)
  result <- read_copc(f,
    bbox = c(tiles$read_xmin[i], tiles$read_ymin[i],
             tiles$read_xmax[i], tiles$read_ymax[i]),
    progress = FALSE)

  if (nrow(result$data) == 0L) return(NULL)

  # Process (e.g., compute mean Z)
  z_mean <- mean(result$data$Z)

  # Trim collar: keep only points in the core tile
  core <- result$data[
    result$data$X >= tiles$xmin[i] & result$data$X <= tiles$xmax[i] &
    result$data$Y >= tiles$ymin[i] & result$data$Y <= tiles$ymax[i], ]

  data.frame(
    tile = i,
    n_with_collar = nrow(result$data),
    n_core = nrow(core),
    z_mean = z_mean
  )
})
results <- do.call(rbind, Filter(Negate(is.null), results))
cat("\nTile results:\n")
print(results)

cat("\n===================================================================\n")
cat(" VPC demo complete! The .vpc file can be opened directly in QGIS.\n")
cat("===================================================================\n")
