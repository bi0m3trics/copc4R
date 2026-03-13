# =============================================================================
# copc4R -- Example: Pick a point + buffer AOI + STAC read_copc()
# =============================================================================
#
# Demonstrates the most concise copc4R workflow:
#
#   1. Define a point location and buffer it to a circular AOI (100 m radius)
#   2. Pass a STAC search URL directly to read_copc()
#      - Collections are auto-detected for Planetary Computer
#      - Tiles are discovered, signed, and range-read automatically
#      - Points are clipped to the buffered circle across tile boundaries
#   3. Apply a combined filter string in one call:
#        -drop_withheld  remove withheld points
#        -keep_voxel 1.0 subsample to ~1 m point spacing (3-D voxel grid)
#        -keep_first     only first returns
#        -drop_noise     remove class 7 & 18 (ASPRS noise)
#   4. Convert to a lidR::LAS object and visualise
#
# Key point: read_copc() handles multi-tile STAC searches transparently --
# the AOI (in WGS 84) is automatically reprojected into each tile's native
# CRS, and results from all intersecting tiles are merged into one table.
#
# Required packages:
#   install.packages(c("sf", "lidR", "httr", "jsonlite"))
#
# Data: https://planetarycomputer.microsoft.com/dataset/3dep-lidar-copc
# STAC: https://planetarycomputer.microsoft.com/api/stac/v1
# =============================================================================

library(sf)
library(copc4R)

# =============================================================================
# 1. Define a point AOI and buffer to a 100 m circle
# =============================================================================
# Point near the southwest edge of Flagstaff, AZ (WGS 84 / EPSG:4326).
# Change coords to any WGS 84 longitude/latitude you like.

coords <- c(-111.72671, 35.10700)   # lon, lat

sf_use_s2(TRUE)   # spherical geometry on -- st_buffer unit is metres

aoi <- st_sf(
  id       = 1,
  geometry = st_sfc(st_point(coords), crs = 4326)
)

# Buffer in metres (s2 geometry engine interprets distance as metres for
# geographic CRS when sf_use_s2(TRUE))
aoi <- st_buffer(aoi, dist = 100)

# =============================================================================
# 2. Stream points from Planetary Computer 3DEP via STAC + COPC range reads
# =============================================================================
# read_copc() detects the STAC endpoint, searches for intersecting tiles,
# fetches only the octree nodes that overlap the AOI bounding box via HTTP
# range reads, then clips to the polygon boundary and merges across tiles.
# The "3dep-lidar-copc" collection is selected automatically for this host.

result <- read_copc(
  "https://planetarycomputer.microsoft.com/api/stac/v1/search",
  aoi    = aoi,
  select = "xyzicrnap",    # X Y Z Intensity Classification ReturnNumber
                           # NumberOfReturns ScanAngle PointSourceID
  filter = paste(
    "-drop_withheld",      # remove withheld flag points
    "-keep_voxel 1.0",     # 3-D voxel thinning at 1 m grid
    "-keep_first",         # first returns only
    "-drop_noise"          # drop ASPRS noise classes (7, 18)
  ),
  progress = TRUE
)

# =============================================================================
# 3. Convert to a lidR::LAS object
# =============================================================================
# as_las() builds a valid LAS 1.4 object from the data.table + header.
# The native tile CRS (NAD83 / UTM zone 12N, EPSG:26912) is preserved.

las <- as_las(result)

# =============================================================================
# 4. Visualise
# =============================================================================
# 3-D interactive point cloud coloured by height (default).

lidR::plot(las)
