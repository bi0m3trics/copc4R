# =============================================================================
# catalog_demo.R -- COPC-Native Parallel Chunk Processing with copc_apply()
# =============================================================================
#
# Demonstrates the full power of copc4R's processing engine:
#   1. Define a WGS 84 AOI rectangle
#   2. STAC search for USGS 3DEP COPC tiles on Planetary Computer
#   3. copc_apply() -- partition into 200 m chunks with 15 m buffer,
#      process each chunk in parallel:
#        ground classification -> height normalisation -> CHM ->
#        tree detection -> crown segmentation -> crown metrics
#   4. Auto-merged sf result of crown-level forest metrics
#
# This workflow:
#   - Never downloads full point cloud tiles
#   - Reads only the COPC octree nodes overlapping each chunk (range reads)
#   - Bounds peak memory to one chunk per worker
#   - Parallelises via the future framework (plan(multisession))
#   - Automatically clips buffer from results -- no edge artifacts
#
# Required packages:
#   install.packages(c("copc4R", "lidR", "sf", "terra", "httr",
#                      "jsonlite", "future", "future.apply"))
#
# Data: https://planetarycomputer.microsoft.com/dataset/3dep-lidar-copc
# =============================================================================

library(copc4R)
library(lidR)
library(sf)
library(terra)
library(httr)
library(jsonlite)
library(future)

cat("===================================================================\n")
cat(" copc4R -- Parallel chunk processing demo (copc_apply)\n")
cat("===================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Define the AOI rectangle (WGS 84)
# ---------------------------------------------------------------------------
# Forested area near Bryce Canyon, UT -- approx 1 km x 1 km
coords <- matrix(
  c(
    -111.5649374, 35.2452949,
    -111.5697439, 35.2407034,
    -111.5633925, 35.2356911,
    -111.5637787, 35.2292061,
    -111.5557964, 35.2291009,
    -111.5525778, 35.2319754,
    -111.5506895, 35.2324311,
    -111.5485008, 35.2335178,
    -111.5448101, 35.2358313,
    -111.5447672, 35.2447692,
    -111.5532215, 35.2460309,
    -111.5556248, 35.2463464,
    -111.5564831, 35.2461361,
    -111.5624054, 35.2447692,
    -111.5638216, 35.2447692,
    -111.5649374, 35.2452949
  ),
  ncol = 2, byrow = TRUE
)

aoi <- st_sf(
  id       = 1,
  geometry = st_sfc(st_polygon(list(coords)), crs = 4326)
)

cat("AOI (WGS 84):\n")
print(st_bbox(aoi))
cat("\n")

# ---------------------------------------------------------------------------
# 2. Search Planetary Computer STAC for overlapping 3DEP COPC tiles
# ---------------------------------------------------------------------------
cat("Searching Planetary Computer STAC for 3DEP tiles...\n")

bbox_vec <- as.numeric(st_bbox(aoi))

resp <- POST(
  "https://planetarycomputer.microsoft.com/api/stac/v1/search",
  content_type_json(),
  body = toJSON(list(
    collections = list("3dep-lidar-copc"),
    bbox         = bbox_vec,
    limit        = 100
  ), auto_unbox = TRUE)
)
stop_for_status(resp)
items <- fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                  simplifyVector = FALSE)$features

cat(sprintf("  Found %d tile(s)\n", length(items)))
if (length(items) == 0L) stop("No 3DEP COPC tiles found for this AOI.")

# Sign tile URLs with a SAS token
sas <- fromJSON(content(
  GET("https://planetarycomputer.microsoft.com/api/sas/v1/token/3dep-lidar-copc"),
  as = "text", encoding = "UTF-8"
))$token

tile_urls <- vapply(items, function(item) {
  paste0(item$assets$data$href, "?", sas)
}, character(1))

for (i in seq_along(items))
  cat(sprintf("  [%d] %s\n", i, items[[i]]$id))
cat("\n")

# ---------------------------------------------------------------------------
# 3. Set up parallel plan
# ---------------------------------------------------------------------------
n_workers <- min(parallelly::availableCores(logical = FALSE), 6L)
plan(multisession, workers = n_workers)
cat(sprintf("Parallel plan: multisession with %d workers\n\n", n_workers))

# ---------------------------------------------------------------------------
# 4. Define the per-chunk processing function
# ---------------------------------------------------------------------------
# FUN receives a LAS object (chunk + 15 m buffer) and must return a result.
# copc_apply() clips the buffer from the output automatically after FUN
# returns, so no manual clipping is required.

process_chunk <- function(las) {

  # Skip chunks with too few points for a meaningful analysis
  if (lidR::npoints(las) < 100L) return(NULL)

  # -- Ground classification (CSF) ------------------------------------
  # Buffer ensures the ground surface is correctly estimated at the
  # chunk edges -- TIN / CSF needs neighbour points beyond the boundary.
  las <- lidR::classify_ground(las, algorithm = lidR::csf())

  # -- Normalise heights (TIN) ----------------------------------------
  nlas <- lidR::normalize_height(las, algorithm = lidR::tin())

  # -- Canopy height model (pit-free) ---------------------------------
  chm <- lidR::rasterize_canopy(nlas, res = 0.5,
                                algorithm = lidR::pitfree(
                                  thresholds = c(0, 10, 20),
                                  max_edge   = c(0, 1.5)
                                ))

  # -- Individual tree detection (local maximum filter) ---------------
  ttops <- lidR::locate_trees(chm, algorithm = lidR::lmf(ws = 2.5))
  if (is.null(ttops) || nrow(ttops) == 0L) return(NULL)

  # -- Crown segmentation (Dalponte 2016 watershed) -------------------
  algo <- lidR::dalponte2016(chm, ttops)
  nlas <- lidR::segment_trees(nlas, algo)

  # -- Crown-level metrics --------------------------------------------
  # Returns an sf object with one row per crown.  Buffer-region crowns
  # are present here but will be clipped by copc_apply() automatically.
  metrics <- lidR::crown_metrics(nlas, ~list(z_max = max(Z), z_mean = mean(Z)))
  return(metrics)
}

# ---------------------------------------------------------------------------
# 5. Run copc_apply()
# ---------------------------------------------------------------------------
cat("Starting copc_apply()...\n")
cat(sprintf("  Chunk size : 200 m\n"))
cat(sprintf("  Buffer     : 15 m\n"))
cat(sprintf("  Filter     : -drop_withheld\n\n"))

t0 <- proc.time()

results <- copc_apply(
  source       = tile_urls,
  FUN          = process_chunk,
  aoi          = aoi,
  chunk_size   = 200,
  chunk_buffer = 15,
  select       = "*",
  filter       = "-drop_withheld",
  automerge    = TRUE,
  progress     = TRUE,
  plot = TRUE
)

elapsed <- (proc.time() - t0)["elapsed"]

# Reset to sequential -- always clean up parallel plans
plan(sequential)

# ---------------------------------------------------------------------------
# 6. Results
# ---------------------------------------------------------------------------
cat("\n===================================================================\n")
cat(" RESULTS\n")
cat("===================================================================\n")
cat(sprintf("  Elapsed    : %.1f seconds\n", elapsed))

if (inherits(results, "sf") && nrow(results) > 0L) {
  cat(sprintf("  Crowns     : %s\n", format(nrow(results), big.mark = ",")))
  cat(sprintf("  z_max range: %.1f -- %.1f m\n",
              min(results$z_max, na.rm = TRUE),
              max(results$z_max, na.rm = TRUE)))
  cat(sprintf("  z_mean     : %.1f -- %.1f m\n",
              min(results$z_mean, na.rm = TRUE),
              max(results$z_mean, na.rm = TRUE)))

  cat("\n  Plotting crown z_max...\n")
  plot(results["z_max"],
       main = "Crown z_max (m) -- copc4R + 3DEP + copc_apply()",
       key.pos = 4)
} else {
  cat("  Result type: ", class(results), "\n")
  cat("  (No crowns detected -- check AOI and tile coverage)\n")
}

cat("===================================================================\n")
