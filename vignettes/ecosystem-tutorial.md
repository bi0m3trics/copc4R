# copc4R in the R LiDAR Ecosystem — lidR, lasR, and Cloud-Native Workflows

> Andrew Sánchez Meador · Northern Arizona University
> copc4R 0.3.0 · March 2026

---

## Contents

1. [The big picture](#the-big-picture)
2. [Installation](#installation)
3. [What each package does](#what-each-package-does)
4. [Tutorial 1 — Quick start: point + buffer → lidR](#tutorial-1--quick-start)
5. [Tutorial 2 — Canopy Height Model from Planetary Computer](#tutorial-2--chm-from-planetary-computer)
6. [Tutorial 3 — Download and catalog with VPC](#tutorial-3--download-and-catalog-with-vpc)
7. [Tutorial 4 — Large-area processing with lasR](#tutorial-4--large-area-processing-with-lasr)
8. [Deep dive — lasR pipeline patterns](#deep-dive--lasr-pipeline-patterns)
9. [Tutorial 5 — The round trip: lasR writes COPC too](#tutorial-5--the-round-trip)
10. [Tutorial 6 — Terrain metrics with adaptive chunking](#tutorial-6--terrain-metrics-with-adaptive-chunking)
11. [When to use which tool](#when-to-use-which-tool)
12. [VPC: the interoperability glue](#vpc-the-interoperability-glue)
13. [Reference cheat sheet](#reference-cheat-sheet)

---

## The big picture

R has a mature LiDAR ecosystem anchored by two packages from
Jean-Romain Roussel's [r-lidar](https://github.com/r-lidar) project:

| Package | Role |
|---------|------|
| **lidR** | Interactive analysis — loads LAS into R `data.frame`/`LAS` objects, rich API for DTM, CHM, segmentation, metrics, plotting |
| **lasR** | Production processing — C++ pipeline engine, multi-file, streaming, single-pass, minimal RAM |

Both packages assume **local files**.  They read `.las` / `.laz` tiles
that already live on disk.

**copc4R** fills a gap: it is a *cloud-native COPC I/O specialist*.
It speaks HTTP range requests to the
[COPC](https://copc.io/) octree, fetches **only the bytes you need**
from remote servers (S3, Azure Blob, Planetary Computer), and hands
clean data to lidR or lasR for processing.

```
  ┌────────────────────────────────────────┐
  │           Cloud / HTTP                 │
  │  Planetary Computer, AWS, Azure, ...   │
  └────────────┬───────────────────────────┘
               │  COPC range reads
               ▼
  ┌────────────────────────────────────────┐
  │             copc4R                     │
  │  discover · query · fetch · VPC write  │
  └──────────┬───────────────┬─────────────┘
             │               │
             ▼               ▼
        ┌──────────┐    ┌──────────┐
        │   lidR   │    │   lasR   │
        │ (as_las) │    │ (exec on │
        │ R-native │    │  .vpc)   │
        └──────────┘    └──────────┘
```

The **VPC** (Virtual Point Cloud) is the interoperability format —
a lightweight JSON catalog that copc4R writes and that lasR, QGIS,
and PDAL wrench all understand natively.

---

## Installation

```r
# copc4R (GitHub)
remotes::install_github("bi0m3trics/copc4R")

# lidR (CRAN)
install.packages("lidR")

# lasR (r-universe -- not on CRAN)
install.packages("lasR", repos = "https://r-lidar.r-universe.dev")

# supporting packages
install.packages(c("sf", "terra", "jsonlite", "httr", "future"))
```

---

## What each package does

### copc4R — cloud-native COPC I/O

- **HTTP range reads** — fetches only the octree nodes overlapping your
  AOI from remote COPC files; no full-tile downloads.
- **STAC integration** — pass a Planetary Computer (or any STAC) search
  URL to `read_copc()` and it discovers tiles, signs URLs, merges
  results automatically.
- **Spatial filters** — bbox, polygon, corridor, Z-range, attribute
  filters, LOD depth/resolution.
- **Octree estimation** — `copc_estimate_points()` counts hierarchy
  nodes without decompressing a single point.
- **Adaptive splitting** — `copc_adaptive_split()` subdivides large
  areas into balanced work units, useful with `copc_apply()`.
- **VPC read/write** — spec-conformant STAC ItemCollection files
  readable by lasR, QGIS, and PDAL wrench.
- **`as_las()`** — converts copc4R output to `lidR::LAS` in one call.
- **Cache** — in-memory + on-disk caching of hierarchy pages and
  headers for repeated queries.

### lidR — interactive R analysis

- Loads full point clouds into R memory as `LAS` objects.
- Rich API: `classify_ground()`, `normalize_height()`,
  `rasterize_canopy()`, `segment_trees()`, `pixel_metrics()`, …
- Interactive 3-D plotting (`plot(las)`).
- `LAScatalog` for large multi-file projects with `catalog_apply()`.
- Best for: exploratory work, small-to-medium datasets,
  custom R-based algorithms.

### lasR — production C++ pipelines

- Tiny R API wrapping a full C++ processing engine.
- *Pipelines*: chain stages with `+`; every stage sees each point
  once in a single pass — no intermediate files, minimal RAM.
- **Reading stages:** `reader()` (auto-detect), `reader_las()`,
  `reader_copc()`.  All three accept VPC files, multi-file vectors,
  or single files.
- **Processing stages:** `classify_with_csf()`, `classify_with_ivf()`,
  `triangulate()`, `normalize()`, `rasterize()`, `local_maximum()`,
  `hulls()`, `sampling_*()`, …
- **Output stages:** `dtm()`, `chm()`, `write_las()`,
  `write_copc()`, `write_vpc()`.  Rasters are returned as
  `terra::SpatRaster`; written LAS/COPC tiles go to disk.
- **Custom R logic:** `callback()` lets you inject an arbitrary R
  function into the pipeline to compute per-chunk metrics, run
  classification models, or fuse with external data.
- **Filters:** `keep_ground()`, `keep_first()`, `keep_class()`,
  `drop_noise()`, etc. — applied at the reader level or per-stage.
- Reads VPC files natively via `exec(pipeline, on = "file.vpc")`.
- Multi-core parallelism at the C++ level (`ncores`).
- **Can also write COPC and VPC** — so the data can round-trip
  between copc4R and lasR in both directions.
- Best for: batch production on local or shared storage.

---

## Tutorial 1 — Quick start: point + buffer → lidR {#tutorial-1--quick-start}

**Goal:** Grab a 100 m radius point cloud from USGS 3DEP and visualize
it in lidR.  One `read_copc()` call handles STAC search, tile
discovery, range reads, clipping, and merging.

```r
library(sf)
library(copc4R)

# 1. Define a point anywhere in the US and buffer to a circle
sf_use_s2(TRUE)
aoi <- st_buffer(
  st_sf(geometry = st_sfc(st_point(c(-111.6550, 35.1700)), crs = 4326)),
  dist = 100        # metres (s2 engine)
)

# 2. Stream points from Planetary Computer via STAC + COPC range reads
result <- read_copc(
  "https://planetarycomputer.microsoft.com/api/stac/v1/search",
  aoi    = aoi,
  select = "xyzic",          # X Y Z Intensity Classification
  filter = "-drop_withheld -drop_noise -keep_first",
  progress = TRUE
)

cat("Points:", nrow(result$data), "\n")

# 3. Convert to lidR and plot
las <- as_las(result)
lidR::plot(las, color = "Z")
```

**What happened under the hood:**

1. `read_copc()` detected the Planetary Computer STAC endpoint.
2. It searched for `3dep-lidar-copc` tiles intersecting the AOI bbox.
3. Each matching tile URL was signed with a SAS token automatically.
4. Only the COPC octree nodes overlapping the (reprojected) AOI were
   fetched via HTTP range requests — not the entire multi-GB tile.
5. Points were clipped to the circular polygon boundary.
6. Results from all tiles were merged into one `data.table`.

---

## Tutorial 2 — Canopy Height Model from Planetary Computer {#tutorial-2--chm-from-planetary-computer}

**Goal:** End-to-end pipeline — fetch a ~300 × 300 m forested area,
classify ground, normalize heights, and produce a 0.5 m CHM raster.

```r
library(copc4R)
library(lidR)
library(sf)

# 1. Rectangular AOI near NAU School of Forestry, Flagstaff AZ
aoi <- st_sf(
  geometry = st_sfc(st_polygon(list(matrix(c(
    -111.654, 35.1682,
    -111.649, 35.1682,
    -111.649, 35.1707,
    -111.654, 35.1707,
    -111.654, 35.1682
  ), ncol = 2, byrow = TRUE))),
  crs = 4326)
)

# 2. Fetch all points inside AOI from Planetary Computer
result <- read_copc(
  "https://planetarycomputer.microsoft.com/api/stac/v1/search",
  aoi    = aoi,
  select = "xyzicrnap",
  filter = "-drop_withheld",
  progress = TRUE
)

las <- as_las(result)
cat(sprintf("Points: %s | CRS: %s\n",
            format(npoints(las), big.mark = ","),
            st_crs(las)$input))

# 3. Classify ground (Cloth Simulation Filter)
las <- classify_ground(las, algorithm = csf(
  sloop_smooth     = FALSE,
  class_threshold  = 0.5,
  cloth_resolution = 0.5,
  rigidness        = 1L
))

# 4. Normalize heights (TIN interpolation)
nlas <- normalize_height(las, algorithm = tin())
cat(sprintf("Height range: %.1f – %.1f m\n", min(nlas$Z), max(nlas$Z)))

# 5. Build a pitfree CHM at 0.5 m resolution
chm <- rasterize_canopy(nlas, res = 0.5, algorithm = pitfree(
  thresholds = c(0, 2, 5, 10, 15),
  max_edge   = c(0, 1.5)
))

terra::plot(chm, col = height.colors(50),
            main = "CHM — copc4R + 3DEP + lidR")
```

---

## Tutorial 3 — Download and catalog with VPC {#tutorial-3--download-and-catalog-with-vpc}

**Goal:** Use `download_3dep_copc()` to grab multiple tiles, write a
VPC index, and then process the catalog in parallel chunks.

```r
library(copc4R)
library(sf)

# 1. Define a larger AOI (~1 km²)
aoi <- st_sf(
  geometry = st_sfc(st_polygon(list(matrix(c(
    -111.660, 35.165,
    -111.648, 35.165,
    -111.648, 35.175,
    -111.660, 35.175,
    -111.660, 35.165
  ), ncol = 2, byrow = TRUE))),
  crs = 4326)
)

# 2. Download -- only AOI points are fetched per tile (range reads)
dest <- file.path(tempdir(), "3dep_tiles")
download_3dep_copc(
  aoi      = aoi,
  dest_dir = dest,
  filter   = "-drop_withheld -drop_noise",
  merge    = FALSE,          # keep individual tile clips
  progress = TRUE
)

# 3. List the COPC files we just wrote
files <- list.files(dest, "\\.copc\\.laz$|\\.laz$", full.names = TRUE)
cat(length(files), "tile file(s) downloaded\n")

# 4. Write a VPC index (interoperable with lasR, QGIS, PDAL wrench)
vpc_path <- file.path(dest, "flagstaff.vpc")
copc4R::write_copc_vpc(files, vpc_path, statistics = TRUE, progress = TRUE) # Careful, there is a lasR::write_vpc
cat("VPC written:", vpc_path, "\n")

# 5. Read VPC as a catalog — fast path trusts VPC metadata
ctg <- read_copc_vpc(vpc_path, chunk_size = 200, chunk_buffer = 10)
print(ctg)

# 6. Process in parallel chunks with copc_apply()
library(future)
plan(multisession, workers = 4)

results <- copc_apply(ctg, function(las) {
  data.frame(
    n     = nrow(las@data),
    z_avg = mean(las$Z),
    z_max = max(las$Z)
  )
}, progress = TRUE)

plan(sequential)
# copc_apply() auto-merges data.frame results via rbind
print(results)
```

**What the VPC file looks like** (abridged):

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "stac_version": "1.0.0",
      "id": "USGS_LPC_AZ_VerdeKa686...",
      "bbox": [-111.66, 35.165, 1950.3, -111.648, 35.175, 2105.7],
      "properties": {
        "pc:count": 1847230,
        "pc:type": "lidar",
        "pc:indexed": true,
        "proj:wkt2": "PROJCRS[\"NAD83 / UTM zone 12N\", ...]",
        "proj:bbox": [421350.0, 3891200.0, 1950.3, 422400.0, 3892300.0, 2105.7]
      },
      "assets": {
        "data": { "href": "./USGS_LPC_AZ_...copc.laz", "roles": ["data"] }
      }
    }
  ]
}
```

Key VPC properties written by copc4R:
- **`pc:indexed = true`** — COPC always has a spatial index (octree).
  lasR reads this to skip redundant lax indexing.
- **`proj:wkt2`** and **`proj:projjson`** — robust CRS round-tripping,
  especially for compound CRS (horizontal + vertical datums).
- **3-D `proj:bbox`** — 6-element array with Z-range, as the VPC spec
  recommends.

---

## Tutorial 4 — Large-area processing with lasR {#tutorial-4--large-area-processing-with-lasr}

**Goal:** Use copc4R for data acquisition, then hand the VPC to lasR
for fast C++ pipeline-based processing (DTM + CHM in a single pass).

```r
library(copc4R)

# 1. We already have a VPC from Tutorial 3
vpc_path <- file.path(tempdir(), "3dep_tiles", "flagstaff.vpc")

# --- Option A: lasR processes the VPC directly ---
if (requireNamespace("lasR", quietly = TRUE)) {
  library(lasR)

  # Build a pipeline: ground classification → DTM → CHM
  pipeline <- classify_with_csf() +
    triangulate(filter = keep_ground()) +
    dtm(1) +                             # 1 m DTM
    normalize() +
    chm(0.5)                             # 0.5 m CHM

  # exec() reads the VPC natively — no copc4R needed at this stage
  result <- exec(pipeline, on = vpc_path, ncores = 4, progress = TRUE)

  # result is a named list of SpatRasters — dtm() and chm() are
  # convenience wrappers around rasterize(), so results are named
  # "rasterize" and "rasterize.1" (in pipeline order).
  cat("Result names:", paste(names(result), collapse = ", "), "\n")
  dtm_rast <- result[[1]]               # first rasterize stage = DTM
  chm_rast <- result[[2]]               # second rasterize stage = CHM

  # lasR temp GeoTIFFs lack min/max stats; terra::plot() needs them
  dtm_rast <- terra::setMinMax(dtm_rast)
  chm_rast <- terra::setMinMax(chm_rast)

  terra::plot(dtm_rast,  main = "DTM (lasR)", col = gray.colors(25))
  terra::plot(chm_rast, main = "CHM (lasR)",
              col = grDevices::colorRampPalette(
                c("blue", "cyan", "yellow", "red"))(25))
}
```

**Why this works well:**

- copc4R fetches *just the AOI* from huge remote COPC tiles via range
  reads — no multi-GB downloads.
- `write_copc_vpc()` creates a spec-conformant VPC with `pc:indexed = true`,
  3-D bboxes, and proj:wkt2 — everything lasR expects.
- lasR's `exec()` reads the VPC metadata in one shot (no header
  re-reads), processes each tile through the full pipeline in C++,
  writes tiled rasters, and merges seamlessly.

```r
# --- Option B: copc4R read → as_las() → lidR ---
# Better for small AOIs or when you need full R-level control

result <- read_copc(files[1], progress = FALSE)
las <- as_las(result)

las <- lidR::classify_ground(las, lidR::csf())
nlas <- lidR::normalize_height(las, lidR::tin())
chm <- lidR::rasterize_canopy(nlas, 0.5, lidR::pitfree())
```

---

## Deep dive — lasR pipeline patterns {#deep-dive--lasr-pipeline-patterns}

lasR's design is different from lidR's.  Where lidR gives you
imperative R functions you call one at a time, lasR lets you declare
a **pipeline** — a chain of processing stages connected with `+`.
The engine reads the data *once*, passes every point through every
stage in sequence, and writes all outputs at the end.  This
single-pass design is what makes lasR dramatically faster and more
memory-efficient than equivalent lidR code on large datasets.

### How pipelines work

```r
library(lasR)

# Each stage is a function call that returns a "stage" object.
# Chain them with + to build a pipeline.
pipeline <-
  reader() +                            # stage 1: read tiles
  classify_with_csf() +                 # stage 2: ground classification
  triangulate(filter = keep_ground()) + # stage 3: Delaunay TIN of ground
  normalize() +                         # stage 4: subtract TIN from Z
  dtm(1) +                              # stage 5: rasterize a 1 m DTM
  chm(0.5)                              # stage 6: rasterize a 0.5 m CHM

# exec() drives the pipeline across one or more input files.
# "on" can be: a single file, a vector of files, or a VPC path.
ans <- exec(pipeline, on = vpc_path, ncores = 4, progress = TRUE)
```

`exec()` returns a named list only for stages that produce R-side
output (rasters, data frames).  Stages like `write_las()` write
directly to disk and return nothing.

### Filtering

lasR filters follow the LAStools/LASlib convention and can be
applied globally (at the reader) or per-stage:

```r
# Global filter — only first returns enter the pipeline
pipeline <- reader(filter = keep_first()) +
  rasterize(1, "zmean")

# Per-stage filter — only ground points go into the TIN
pipeline <- reader() +
  triangulate(filter = keep_ground()) +
  dtm(1)
```

### Writing output files

lasR can write LAS, LAZ, COPC, and VPC:

```r
# Write normalized, noise-free tiles as compressed LAZ
pipeline <- classify_with_csf() +
  classify_with_ivf() +
  normalize() +
  write_las(paste0(tempdir(), "/*_norm.laz"),
            filter = drop_noise())

# Write COPC (cloud-optimized) output
pipeline <- classify_with_csf() +
  normalize() +
  write_copc(paste0(tempdir(), "/*_norm.copc.laz"))

exec(pipeline, on = vpc_path, ncores = 4)
```

Output filenames use `*` as a placeholder — lasR replaces it with
the source tile name, so each input tile produces its own output.

### Custom R callbacks

The `callback()` stage lets you inject any R function into the
pipeline.  The function receives a `data.frame` of the current
chunk's points and can return metrics, modify points, or trigger
side effects.

```r
# Per-chunk canopy metrics via callback
my_metrics <- function(data) {
  above2 <- data$Z[data$Z > 2]
  data.frame(
    n        = nrow(data),
    z_mean   = mean(data$Z),
    z_p95    = quantile(data$Z, 0.95, names = FALSE),
    cover    = length(above2) / nrow(data),
    gap_frac = 1 - length(above2) / nrow(data)
  )
}

pipeline <- classify_with_csf() +
  normalize() +
  callback(my_metrics, expose = "xyzc")

ans <- exec(pipeline, on = vpc_path, ncores = 4)
# ans$callback is a data.frame with one row per chunk
```

### Rasterize with expressions

`rasterize()` accepts R expressions as operators, giving you
fine-grained control over raster outputs:

```r
# Canopy cover raster at 10 m resolution
pipeline <- normalize() +
  rasterize(10, mean(Z > 2))

# Multi-band: mean height + point count
pipeline <- normalize() +
  rasterize(10, c("zmean", "count"))
```

### copc4R + lasR: the full picture

The key insight: **copc4R handles the cloud-to-local bridge;
lasR handles the local processing at speed.**

```
  ┌─── Cloud ──────────────────────────────────────────────┐
  │  Planetary Computer / S3 / Azure Blob / HTTP           │
  └────────────┬───────────────────────────────────────────┘
               │ COPC range reads (copc4R)
               ▼
  ┌────────────────────────────────────────────────────────┐
  │  copc4R                                                │
  │    read_copc()          → as_las() → lidR (interactive)│
  │    download_3dep_copc() → write_copc_vpc() → .vpc file ─────┤
  └────────────────────────────────────────────────────────┘
                                                    │
                                                    ▼
  ┌────────────────────────────────────────────────────────┐
  │  lasR                                                  │
  │    exec(pipeline, on = "file.vpc", ncores = 8)         │
  │    → DTM, CHM, normalized tiles, metrics, ...          │
  │    → write_copc() / write_vpc() back out               │
  └────────────────────────────────────────────────────────┘
```

Without copc4R, getting data *into* lasR from the cloud requires
manually downloading multi-GB tiles, unzipping, building indexes.
copc4R fetches just the AOI via range reads, writes clean local
COPC tiles, and produces a VPC index that lasR can `exec()` on
directly.

---

## Tutorial 5 — The round trip: lasR writes COPC too {#tutorial-5--the-round-trip}

**Goal:** Show the full round trip — copc4R fetches from the cloud,
lasR processes and writes new COPC output, and the cycle continues.

```r
library(copc4R)
library(sf)

# ── Step 1: copc4R fetches data from Planetary Computer ──
aoi <- st_sf(
  geometry = st_sfc(st_polygon(list(matrix(c(
    -111.654, 35.168,
    -111.649, 35.168,
    -111.649, 35.172,
    -111.654, 35.172,
    -111.654, 35.168
  ), ncol = 2, byrow = TRUE))),
  crs = 4326)
)

dest <- file.path(tempdir(), "round_trip")
dir.create(dest, showWarnings = FALSE)

download_3dep_copc(
  aoi      = aoi,
  dest_dir = dest,
  filter   = "-drop_withheld -drop_noise",
  progress = TRUE
)

# Write VPC for the downloaded tiles
raw_files <- list.files(dest, "\\.copc\\.laz$", full.names = TRUE)
vpc_raw   <- file.path(dest, "raw.vpc")
write_copc_vpc(raw_files, vpc_raw, statistics = TRUE)

cat("Raw VPC:", vpc_raw, "\n")
cat("Tiles:", length(raw_files), "\n\n")

# ── Step 2: lasR processes everything in one pass ──
if (requireNamespace("lasR", quietly = TRUE)) {
  library(lasR)

  out_dir <- file.path(dest, "processed")
  dir.create(out_dir, showWarnings = FALSE)

  pipeline <-
    classify_with_csf() +
    classify_with_ivf() +                      # isolated-voxel noise filter
    normalize() +
    dtm(1) +                                   # 1 m DTM raster
    chm(0.5) +                                 # 0.5 m CHM raster
    write_copc(file.path(out_dir, "*.copc.laz"),
               filter = drop_noise())          # clean COPC output

  result <- exec(pipeline, on = vpc_raw, ncores = 2, progress = TRUE)

  # dtm() and chm() are wrappers around rasterize(), so access by index
  dtm_rast <- result[[1]]               # first rasterize stage = DTM
  chm_rast <- result[[2]]               # second rasterize stage = CHM

  # lasR temp GeoTIFFs lack min/max stats; terra::plot() needs them
  dtm_rast <- terra::setMinMax(dtm_rast)
  chm_rast <- terra::setMinMax(chm_rast)

  terra::plot(dtm_rast, main = "DTM — lasR", col = gray.colors(25))
  terra::plot(chm_rast, main = "CHM — lasR",
              col = grDevices::colorRampPalette(
                c("blue", "green", "yellow", "red"))(50))

  # ── Step 3: lasR writes a new VPC for its output ──
  # Option A: use lasR's write_vpc() stage (inline in pipeline)
  # Option B: use copc4R's write_copc_vpc() on the output files
  processed_files <- list.files(out_dir, "\\.copc\\.laz$",
                                full.names = TRUE)
  vpc_processed <- file.path(out_dir, "processed.vpc")
  write_copc_vpc(processed_files, vpc_processed, statistics = TRUE)

  cat("\nProcessed VPC:", vpc_processed, "\n")
  cat("Processed tiles:", length(processed_files), "\n")

  # ── Step 4: further analysis with lidR on processed data ──
  # The processed COPC files are now ground-classified and normalized.
  # copc4R can read them back for targeted spatial queries.
  result2 <- read_copc(processed_files[1], progress = FALSE)
  las     <- as_las(result2)

  # Already normalized — go straight to analysis
  ttops <- lidR::locate_trees(las, lidR::lmf(ws = 5))
  cat(sprintf("Found %d trees\n", nrow(ttops)))
}
```

**The round-trip pattern:**

```
Cloud (MPC/S3)
  │  copc4R::download_3dep_copc()
  ▼
Local COPC tiles + raw.vpc
  │  lasR::exec(pipeline, on = "raw.vpc")
  ▼
Processed COPC tiles + processed.vpc   ← lasR wrote these
  │  copc4R::read_copc()  →  as_las()
  ▼
lidR analysis (tree detection, segmentation, custom metrics)
```

Each tool does what it does best:
- **copc4R** handles cloud access and COPC I/O
- **lasR** handles heavy batch processing
- **lidR** handles interactive analysis and custom R algorithms
- **VPC** is the contract that ties them together

---

## Tutorial 6 — Terrain metrics with adaptive chunking {#tutorial-6--terrain-metrics-with-adaptive-chunking}

**Goal:** Process a dense drone point cloud using copc4R's
octree-aware adaptive splitting to balance chunk sizes, then compute
per-chunk terrain metrics with lidR.

```r
library(copc4R)

# A dense UAS/drone COPC file hosted on CyVerse (public, no auth needed)
url <- "https://data.cyverse.org/dav-anon/iplant/home/jgillan/USGA/imagery_products/hole17_point_cloud.copc.laz"

# 1. Pre-flight: estimate point density across the full extent
est <- copc_estimate_points(url)
cat(sprintf("Full extent: ~%s points, density %.1f pts/m²\n",
            format(est$estimated_points, big.mark = ","),
            est$density))

# 2. Adaptive split — subdivide until each chunk has ≤ 500k points
tiles <- copc_adaptive_split(url, max_points = 500000,
                             min_size = 50, collar = 10)
cat(sprintf("Split into %d tiles (depth range: %d – %d)\n",
            nrow(tiles), min(tiles$depth), max(tiles$depth)))

# 3. Process each tile
library(lidR)

results <- lapply(seq_len(nrow(tiles)), function(i) {
  tile <- tiles[i, ]

  # Read with collar (buffered bbox)
  result <- read_copc(url,
    bbox = c(tile$read_xmin, tile$read_ymin,
             tile$read_xmax, tile$read_ymax),
    filter = "-drop_withheld -drop_noise",
    progress = FALSE
  )

  if (nrow(result$data) == 0) return(NULL)
  las <- as_las(result)

  # Classify ground + normalize
  las <- classify_ground(las, csf())
  nlas <- normalize_height(las, tin())

  # Clip buffer (use core bbox, not read bbox)
  nlas <- clip_rectangle(nlas,
    tile$xmin, tile$ymin, tile$xmax, tile$ymax)

  # Compute terrain metrics
  data.frame(
    tile_id = i,
    n_pts   = npoints(nlas),
    z_mean  = mean(nlas$Z),
    z_p95   = quantile(nlas$Z, 0.95),
    z_max   = max(nlas$Z),
    canopy_cover = mean(nlas$Z > 2)
  )
})

metrics <- do.call(rbind, Filter(Negate(is.null), results))
print(metrics)
```

The key advantage: `copc_adaptive_split()` uses the COPC octree to
*estimate* point counts without decompressing data. Dense areas
get more subdivisions; sparse areas get fewer. Each chunk stays
near the target size, balancing memory and compute time.

---

## When to use which tool

| Scenario | Tool | Why |
|----------|------|-----|
| Grab a small AOI from a remote dataset | **copc4R** `read_copc()` | HTTP range reads fetch only what you need |
| Browse Planetary Computer 3DEP | **copc4R** `download_3dep_copc()` | STAC search + signed URLs + polygon clipping in one call |
| Quick preview (thousands of points) | **copc4R** `read_copc(max_depth = 2)` | LOD sampling from octree — fast even on huge files |
| Estimate point density before downloading | **copc4R** `copc_estimate_points()` | Reads hierarchy only — zero point decompression |
| Build a shareable tile index | **copc4R** `write_copc_vpc()` | STAC-conformant JSON understood by lasR, QGIS, PDAL |
| Interactive analysis / plotting | **lidR** | Rich R API, `LAS` objects, ggplot, 3-D viewer |
| Ground classification / normalization (small area) | **lidR** | `classify_ground()`, `normalize_height()`, `tin()`, `csf()` |
| CHM, DTM, tree segmentation (small-medium) | **lidR** | `rasterize_canopy()`, `segment_trees()`, `pixel_metrics()` |
| Batch production over many tiles | **lasR** | Single-pass C++ pipeline, minimal RAM, multi-core |
| DTM + CHM + classification in one pass | **lasR** | `dtm(1) + chm(0.5) + classify_with_csf()` via `exec()` |
| Write normalized/cleaned COPC tiles | **lasR** | `write_copc()` stage in pipeline — noise-free, classified output |
| Custom per-chunk metrics at scale | **lasR** | `callback()` stage runs your R function inside the C++ pipeline |
| Point-cloud noise filtering | **lasR** | `classify_with_ivf()` isolated-voxel filter, fast and built-in |
| Process a VPC catalog | **lasR** `exec(pipeline, on = "file.vpc")` | Reads VPC natively, parallel tile processing |
| Round-trip: cloud → process → re-index | **copc4R** + **lasR** | copc4R fetches + writes VPC; lasR processes + writes COPC; copc4R re-indexes |

**The handoff pattern:**

```
copc4R  →  VPC  →  lasR
  │                  │
  └── as_las() ──→ lidR
```

---

## VPC: the interoperability glue

Virtual Point Cloud (`.vpc`) files are a STAC ItemCollection
(GeoJSON FeatureCollection) that index multiple point cloud tiles
with their spatial extents, CRS, and optional metadata.

**Who reads VPC?**

| Consumer | How |
|----------|-----|
| **copc4R** | `read_copc_vpc("file.vpc")` → `copc_catalog` |
| **lasR** | `exec(pipeline, on = "file.vpc")` |
| **QGIS** | Open as a native point cloud layer (3.32+) |
| **PDAL wrench** | `pdal wrench info file.vpc` |

**copc4R's `write_copc_vpc()` is spec-conformant:**

- `pc:indexed = true` — signals that tiles have spatial indexes
  (always true for COPC; lasR skips lax generation).
- `proj:wkt2` + `proj:projjson` — robust CRS encoding, handles
  compound CRS (horizontal + vertical).
- 3-D bounding boxes — `proj:bbox` is always 6 elements (with Z).
- Hierarchical overviews — optional pyramid of decimated tiles for
  seamless QGIS zoom.

**`read_copc_vpc()` fast path:**

```r
# Default: trust_vpc = TRUE
# Builds catalog from VPC metadata — no COPC header re-reads
ctg <- read_copc_vpc(vpc_path)  # near-instant for 10,000 tiles

# Override: re-read every header for maximum accuracy
ctg <- read_copc_vpc(vpc_path, trust_vpc = FALSE)
```

---

## Reference cheat sheet

### Data access (copc4R)

```r
# Read points from a bundled file with a bounding box
f <- system.file("extdata", "NoAZCampus_Zuni_Block.copc.laz", package = "copc4R")
result <- read_copc(f, bbox = c(440881, 3892042, 440900, 3892100))

# Read from a remote URL with an sf polygon AOI and LOD limit
result <- read_copc(url, aoi = aoi, max_depth = 3)

# STAC search — discovers + signs + fetches in one call
result <- read_copc(
  "https://planetarycomputer.microsoft.com/api/stac/v1/search",
  aoi = aoi, filter = "-keep_first")

# Sample / thin
result <- read_copc_sample(url, n = 10000)
result <- read_copc_sample(url, voxel_size = 2, mode = "random")

# Metadata only (no point decompression)
hdr     <- read_copc_header(f)
info    <- copc_info(f)
bounds  <- copc_bounds(f, as_sf = TRUE)
density <- copc_density(f)
est     <- copc_estimate_points(f)
```

### Conversion (copc4R → lidR)

```r
las <- as_las(result)      # copc4R result → lidR::LAS
```

### Catalog & VPC (copc4R)

```r
# files and vpc_path are from Tutorial 3
write_copc_vpc(files, vpc_path, statistics = TRUE)
ctg <- read_copc_vpc(vpc_path)                   # fast path (default)
ctg <- read_copc_catalog(files, chunk_size = 200)

# Parallel chunk processing (see Tutorial 3 for full example)
library(future)
plan(multisession, workers = 4)
copc_apply(ctg, function(las) {
  data.frame(n = nrow(las@data), z_avg = mean(las$Z))
}, progress = TRUE)
plan(sequential)
```

### Adaptive splitting (copc4R)

```r
tiles <- copc_adaptive_split(url, max_points = 500000, collar = 10)
```

### lidR processing

```r
las <- classify_ground(las, csf())
nlas <- normalize_height(las, tin())
chm  <- rasterize_canopy(nlas, 0.5, pitfree())
dtm  <- rasterize_terrain(las, 1, tin())
ttops <- locate_trees(nlas, lmf(ws = 5))
```

### lasR pipeline processing

```r
library(lasR)

# Basic DTM + CHM pipeline
pipeline <- classify_with_csf() +
  triangulate(filter = keep_ground()) +
  dtm(1) +
  normalize() +
  chm(0.5)

exec(pipeline, on = vpc_path, ncores = 8, progress = TRUE)

# Write clean COPC output (drop noise, keep classified + normalized)
pipeline <- classify_with_csf() +
  classify_with_ivf() +
  normalize() +
  write_copc(paste0(tempdir(), "/*_clean.copc.laz"),
             filter = drop_noise())

exec(pipeline, on = vpc_path, ncores = 4)

# Custom per-chunk metrics via callback
pipeline <- normalize() +
  callback(function(data) {
    data.frame(z_p95 = quantile(data$Z, 0.95),
               cover = mean(data$Z > 2))
  }, expose = "xyzc")

ans <- exec(pipeline, on = vpc_path)
# ans$callback => data.frame with one row per chunk

# Rasterize custom expressions
pipeline <- normalize() +
  rasterize(10, mean(Z > 2))               # canopy cover raster

# Global filter: only first returns
pipeline <- reader(filter = keep_first()) +
  rasterize(1, "zmax")

# Per-stage filter: ground-only TIN
pipeline <- reader() +
  triangulate(filter = keep_ground()) +
  dtm(1)
```

### Cache management (copc4R)

```r
copc_cache_config(mem_max_mb = 512)
copc_cache_stats()
copc_cache_clear()
```

---

*copc4R is developed by Andrew Sánchez Meador at Northern Arizona
University.*
*lidR and lasR are developed by Jean-Romain Roussel at
[r-lidar](https://www.r-lidar.com/).*
