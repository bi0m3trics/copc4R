# copc4R <img src="https://github.com/bi0m3trics/copc4R/blob/main/img/copc4R_logo.png" width="150" align="right"/>

![licence](https://img.shields.io/badge/Licence-MIT-blue.svg)

R package to read [Cloud Optimized Point Cloud (COPC)](https://copc.io/) files from local disk or HTTP endpoints using range reads.

`copc4R` relies on a modified version of `LASzip` that was adapted to be compatible with `R`. The library can therefore be compiled into `R` without any complaints from `R CMD check`. It enables R users to read `.copc.laz` binary files — a single-file, spatially indexed variant of LAZ commonly used for streaming LiDAR data from the cloud.

## Installation

```r
# Install the development version from GitHub:
# install.packages("remotes")
remotes::install_github("bi0m3trics/copc4R")
```

## Quick example

```r
library(sf)
library(copc4R)

# ── Define a point location and buffer to a 100 m circular AOI ───────────────
sf_use_s2(TRUE)   # st_buffer unit is metres for geographic CRS

aoi <- st_sf(
  id       = 1,
  geometry = st_sfc(st_point(c(-111.72671, 35.10700)), crs = 4326)
)
aoi <- st_buffer(aoi, dist = 100)

# ── Stream points from USGS 3DEP via Planetary Computer STAC ─────────────────
# read_copc() detects the STAC endpoint, discovers intersecting tiles,
# fetches only the octree nodes overlapping the AOI via HTTP range reads,
# clips to the circle boundary, and merges results across tile boundaries.
# The "3dep-lidar-copc" collection is selected automatically for this host.

result <- read_copc(
  "https://planetarycomputer.microsoft.com/api/stac/v1/search",
  aoi    = aoi,
  select = "xyzicrnap",    # X Y Z Intensity Classification ReturnNumber
                           # NumberOfReturns ScanAngle PointSourceID
  filter = paste(
    "-drop_withheld",      # remove withheld flag points
    "-keep_voxel 1.0",     # 3-D voxel thinning at ~1 m spacing
    "-keep_first",         # first returns only
    "-drop_noise"          # drop ASPRS noise classes (7, 18)
  ),
  progress = TRUE
)

# result$data   -- data.table of point attributes
# result$header -- named list of LAS header fields

# ── Convert to a lidR LAS object and visualise ───────────────────────────────
las <- as_las(result)
lidR::plot(las)
```

## Features

- Reads `.copc.laz` files from **local paths** or **HTTP(S) URLs** (via range reads).
- Spatial bounding box and Z-range filtering through the COPC octree hierarchy.
- Returns `data.table` output with [rlas](https://github.com/r-lidar/rlas)-compatible header conventions.
- Column selection via a lidR-style `select` string (`"xyz"`, `"xyzirc"`, `"*"`, etc.).
- Optional conversion to `lidR::LAS` objects via `as_las()`.
- Bundled `LASzip` for LAZ decompression — no external system dependencies required.

## Copyright Information

`copc4R` contains code written by Andrew Sánchez Meador as well as third-party code included for technical reasons. Details below.

- For `LASzip`:
  - (c) 2007-2021 martin.isenburg@rapidlasso.com - http://rapidlasso.com
  - Provided under GPL-3 license.
- For `rlas` code conventions enabling Martin Isenburg's code to be wrapped into R:
  - (c) 2016-2021 Jean-Romain Roussel
  - Provided under GPL-3 license.
- For the COPC specification:
  - (c) 2021 Andrew Bell, Howard Butler, and Connor Manning of [Hobu, Inc.](https://hobu.co/)
  - Provided under MIT license.
  - https://github.com/copcio/copcio.github.io/tree/main
- For `copc4R` R code:
  - (c) 2025-2026 Andrew Sánchez Meador
  - Provided under MIT license.
