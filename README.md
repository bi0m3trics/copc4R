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
library(copc4R)

# ── Read a remote COPC file ───────────────────────────────────────────────────
copc_file <- "https://s3.amazonaws.com/hobu-lidar/autzen-classified.copc.laz"
result    <- read_copc(copc_file)

# result$header is a named list of LAS header fields
# result$data   is a data.table of point attributes
head(result$data)

# ── Spatial subset with a bounding box ────────────────────────────────────────
result <- read_copc(copc_file, bbox = c(637000, 851000, 638000, 852000))

# ── Inspect header metadata ───────────────────────────────────────────────────
hdr <- read_copc_header(copc_file)
cat("LAS Version :", hdr[["Version Major"]], ".", hdr[["Version Minor"]], "\n")
cat("Point Format:", hdr[["Point Data Format ID"]], "\n")
cat("Total Points:", format(hdr[["Number of point records"]], big.mark = ","), "\n")

# ── Convert to a lidR LAS object (requires lidR) ─────────────────────────────
las <- as_las(result)
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
