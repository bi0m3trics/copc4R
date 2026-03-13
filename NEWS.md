# copc4R 0.2.0

## Build system

* **Replaced `configure` / `cleanup` scripts with a self-contained
  `src/Makevars`**.  The new Makevars detects libcurl at build time
  using Make's `$(shell ...)` — no executable shell scripts required.
  This fixes "not found" errors when building on Linux from filesystems
  that strip execute permissions (OneDrive, iRODS/FUSE, etc.).
  `Makevars.in` is retained for reference but no longer processed.

## New examples

* **`lidR_vs_copc4R_density.R`**: Apples-to-apples benchmark comparing
  a traditional lidR catalog workflow (`readLAScatalog()` → `clip_roi()`
  → `rasterize_density()`) against copc4R's polygon read
  (`read_copc_polygon()` → `as_las()` → `rasterize_density()`).  Times
  each stage, tracks memory usage via `gc()`, and prints a head-to-head
  summary table.

## New features

### Planetary Computer integration

* **`download_3dep_copc()`**: New convenience function to search and download
  USGS 3DEP COPC point-cloud data from Microsoft's Planetary Computer.
  Accepts any `sf` / `sfc` polygon as the area of interest, converts it to
  valid GeoJSON in WGS 84, queries the STAC API (`3dep-lidar-copc`
  collection), and uses COPC range reads to fetch **only the points inside
  the AOI** from each remote tile — no full-tile downloads.  Clipped points
  are written to local `.laz` files via `lidR::writeLAS()`.  Supports
  pagination, per-tile or merged output, overwrite control, `select` /
  `filter` pass-through, and progress reporting.

### Query & filtering

* **`-keep_voxel <size>` filter predicate**: The `filter` string in
  `read_copc()` (and all query wrappers) now accepts `-keep_voxel <size>`.
  This overlays a 3-D voxel grid of the specified cell size (in CRS
  units) and keeps the single point nearest to each voxel's center —
  equivalent to PDAL `filters.voxelcenternearestneighbor`.  Parsed in R
  (not part of the C++ per-point filter) and applied as a post-processing
  step after the AOI clip.

* **Expanded attribute filters**: The `filter` parameter on `read_copc()`
  now supports a much richer set of lidR/lastools-style predicates
  (inspired by PDAL `filters.*`, laslib, and lidR):
    - Classification: `-keep_class`, `-drop_class`, `-keep_ground`,
      `-drop_noise`
    - Return number: `-keep_first`, `-keep_last`, `-keep_single`,
      `-keep_return`, `-drop_return`
    - Flags: `-drop_withheld`, `-drop_overlap`
    - Intensity: `-keep_intensity_above`, `-keep_intensity_below`
    - Z-value: `-keep_z_above`, `-keep_z_below`
    - Scan angle: `-keep_scan_angle_above`, `-keep_scan_angle_below`
    - Thinning: `-keep_random_fraction`, `-keep_every_nth`

* **Resolution-based LOD**: New `resolution` parameter on `read_copc()`
  (and all query wrappers).  Specify a target point spacing in CRS units
  and the appropriate `max_depth` is computed automatically from the COPC
  octree spacing — like PDAL's `readers.copc` `resolution` option.

* **Voxel-based subsampling**: `read_copc_sample()` gains a `voxel_size`
  parameter and `mode` argument (`"random"`, `"first"`, `"center"`).
  Overlays a 3-D voxel grid and keeps one point per occupied cell —
  equivalent to PDAL's `filters.voxeldownsize`.  The optional `n`
  parameter applies a secondary random cap after voxel thinning.

### PDAL-inspired additions

  Features modelled after PDAL's COPC reader and filter stages:
    - `resolution` parameter (cf. PDAL `readers.copc` `resolution`)
    - Voxel thinning modes (cf. PDAL `filters.voxeldownsize`)
    - Decimation via `-keep_every_nth` (cf. PDAL `filters.decimation`)
    - Random fraction via `-keep_random_fraction`

---

# copc4R 0.1.1

## New features

### Data access & performance

* **Smart caching**: Persistent on-disk cache + in-memory LRU for HTTP
  range-read chunks, now fully wired into the read pipeline.  When reading
  remote COPC files, `read_copc()` automatically checks the chunk cache
  before fetching, and newly fetched chunks are stored for future use.
  Header reads are also cached (memoized) across calls — repeated
  `read_copc_header()`, `copc_info()`, `copc_bounds()`, and `copc_density()`
  calls hit the cache instead of issuing new HTTP requests.
  Configure via `copc_cache_config()`, inspect with `copc_cache_stats()`,
  and clear with `copc_cache_clear()`.

* **Parallel HTTP chunk fetch**: `read_copc()` gains a `threads` parameter
  (default 4).  For remote files, compressed chunks are fetched concurrently
  using C++-level `std::async` with persistent per-thread HTTP connections.
  In the cache-aware path, only uncached chunks are parallel-fetched via the
  new `cpp_fetch_raw_chunks_parallel()` C++ function; already-cached chunks
  are served instantly from memory/disk.

* **Streaming iterator**: `read_copc_iter()` returns a closure-based iterator
  that yields point batches on demand, suitable for out-of-core workflows.

* **Header-first metadata**: New lightweight functions that read only the
  file header (no point data):
    - `copc_info()` — COPC octree metadata (center, halfsize, spacing, GPS-time
      range, CRS, point count).
    - `copc_bounds()` — 3-D bounding box (or `sf::st_bbox()` with CRS).
    - `copc_density()` — estimate point density for an AOI without
      decompressing any points (uses COPC hierarchy node counts).

### Query capability

* **Polygon AOI**: `read_copc()` gains an `aoi` parameter accepting any `sf`
  polygon / multipolygon.  The octree is queried via the bbox, and points are
  clipped to the actual geometry boundary.

* **Corridor queries**: `read_copc_corridor()` buffers a LINESTRING by a
  specified width and reads the resulting corridor polygon.

* **Attribute filters**: New `filter` parameter on `read_copc()` with
  lidR-style predicates: `-keep_class`, `-drop_class`, `-keep_first`,
  `-keep_last`, `-keep_single`, `-drop_withheld`, `-drop_overlap`,
  `-keep_intensity_above`, `-keep_intensity_below`.

* **Multi-resolution / LOD sampling**: New `max_depth` parameter on
  `read_copc()`.  Setting `max_depth = 0` reads only the root node
  (coarsest resolution), `max_depth = 3` reads the first 4 octree
  levels, etc.  Combine with `read_copc_sample()` for fast previews.

* **Random subsampling**: `read_copc_sample()` reads from COPC and then
  randomly samples down to at most *n* points.

* **Deterministic tiling**: `read_copc_tiles()` splits an extent into a
  regular grid and reads each tile independently.  Supports configurable
  tile size and overlap buffer.

### lidR / rlas interop

* **Catalog ergonomics**: `read_copc_catalog()` builds a `lidR::LAScatalog`
  (or lightweight `copc_catalog` S3 object) from one or more COPC
  paths/URLs.  `catalog_apply()` processes chunks with a user function.

* **Extra Bytes / dimensions**: Files with Extra Bytes VLR
  (`LASF_Spec` record 4) now have their additional dimensions decoded
  and returned as extra numeric columns in the output `data.table`.
  Scale and offset from the VLR are applied automatically.

### Convenience wrappers

* `read_copc_polygon()` — read points within an `sf` polygon.
* `read_copc_corridor()` — read points along a buffered line.
* `read_copc_sample()` — read and thin to *n* points.

## Breaking changes

* `read_copc()` signature expanded: new parameters `aoi`, `filter`,
  `max_depth`, `threads` inserted before `max_points`/`progress`.
  Existing code using positional arguments may need updating;
  keyword arguments are unaffected.

* `cpp_read_copc()` now takes 10 arguments (was 6).  The R-side
  `RcppExports.R` and registration table are updated accordingly.

## Internal

* C++ `select_nodes()` now accepts `max_depth` for LOD pruning.
* C++ `count_nodes()` added for fast hierarchy-level point counting.
* New C++ exports: `cpp_select_nodes()`, `cpp_fetch_raw_chunk()`,
  `cpp_fetch_raw_chunks_parallel()` — expose node selection and chunk
  fetching to R for cache orchestration.
* Parallel HTTP fetching uses `std::async` with per-thread persistent
  `HttpRangeReader` instances (connection reuse via HTTP keep-alive).
* Three-mode fetch pipeline in `cpp_read_copc()`: prefetched → parallel
  HTTP → sequential, chosen automatically based on parameters.
* Header cache (`R/cache.R`) memoizes `read_copc_header()` results;
  chunk cache stores raw compressed bytes keyed by (URL, offset, size).
* Extra Bytes VLR parsing added to `las14_header.{h,cpp}`.
* `PointFilter` + `parse_filter()` in `rcpp_exports.cpp` for
  attribute-level point filtering during decompression.
* Registration table updated from 3 to 7 C++ entry points.

---

# copc4R 0.1.0

* Initial release.
* Read COPC (.copc.laz) files from local disk or HTTP(S) URLs.
* HTTP range-read support via libcurl (streams only relevant octree
  nodes, no full file download required).
* Spatial bounding box and Z-range queries via the COPC hierarchy.
* Column selection with lidR-style `select` string.
* `max_points` cap.
* Conversion to `lidR::LAS` objects via `as_las()`.
* Bundled LASzip library for in-process LAZ decompression.
