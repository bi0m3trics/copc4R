# -- Virtual Point Cloud (VPC) support ----------------------------------------
#
# Read and write VPC files -- the STAC ItemCollection / GeoJSON
# FeatureCollection index used by the QGIS/PDAL ecosystem.
#
# Spec: https://github.com/PDAL/wrench/blob/main/vpc-spec.md
# ---------------------------------------------------------------------------

#' Read a Virtual Point Cloud (VPC) file
#'
#' Reads a `.vpc` file (a STAC ItemCollection / GeoJSON FeatureCollection)
#' and returns a `copc_catalog` or `lidR::LAScatalog` ready for processing
#' with [copc_catalog_apply()] or [copc_apply()].
#'
#' A VPC file is a lightweight JSON index referencing multiple point cloud
#' files (typically `.copc.laz`) with their spatial extents, CRS, and
#' optional statistics.
#' This is the format used by QGIS and PDAL wrench for virtual point
#' clouds.
#'
#' @note Currently only `.copc.laz` files are supported as referenced
#'   tiles. Non-COPC formats (plain `.las`/`.laz`) are detected and
#'   skipped with a warning. A fallback to `lidR::readLAS()` for full
#'   QGIS parity is planned for a future release.
#'
#' @param path Character. Path or URL to a `.vpc` file.
#' @param chunk_size Numeric. Default chunk size for catalog processing
#'   (0 = file-level). Passed to [read_copc_catalog()].
#' @param chunk_buffer Numeric. Buffer around each chunk.
#' @param select Character. Default column selection.
#' @param filter Character. Default attribute filter string.
#' @param trust_vpc Logical. If `TRUE` (the default), build the catalog
#'   directly from VPC metadata (bbox, CRS, point count) without
#'   re-reading each COPC file header.  This is dramatically faster for
#'   large collections (thousands of files).  Set to `FALSE` to force a
#'   full header read for every tile (original behaviour).
#' @param progress Logical. Print progress messages?
#'
#' @return A `copc_catalog` or `lidR::LAScatalog` object. The VPC
#'   metadata (statistics, overview paths, etc.) is stored in an
#'   attribute `"vpc_metadata"`.
#'
#' @seealso [write_copc_vpc()], [read_copc_catalog()], [copc_apply()]
#'
#' @examples
#' \donttest{
#' ctg <- read_copc_vpc("tiles.vpc")
#' print(ctg)
#'
#' # Use with copc_apply()
#' results <- copc_apply(ctg, function(las) mean(las$Z),
#'                       chunk_size = 500)
#'
#' # Fast path (default) -- trusts VPC metadata, no header re-reads
#' ctg_fast <- read_copc_vpc("big_project.vpc", trust_vpc = TRUE)
#'
#' # Slow path -- re-reads every header for maximum accuracy
#' ctg_strict <- read_copc_vpc("big_project.vpc", trust_vpc = FALSE)
#' }
#'
#' @export
read_copc_vpc <- function(path,
                     chunk_size   = 0,
                     chunk_buffer = 0,
                     select       = "*",
                     filter       = "",
                     trust_vpc    = TRUE,
                     progress     = TRUE) {

  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for VPC support.", call. = FALSE)

  stopifnot(is.character(path), length(path) == 1L, nzchar(path))

  # Read the VPC JSON
  vpc_json <- if (grepl("^https?://", path, ignore.case = TRUE)) {
    jsonlite::fromJSON(url(path), simplifyVector = FALSE)
  } else {
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  }

  # Validate it's a FeatureCollection (STAC ItemCollection)
  fc_type <- vpc_json$type %||% ""
  if (!identical(fc_type, "FeatureCollection"))
    stop("VPC file must be a GeoJSON FeatureCollection (STAC ItemCollection). ",
         "Got type: '", fc_type, "'", call. = FALSE)

  features <- vpc_json$features
  if (length(features) == 0L)
    stop("VPC file contains no features (tiles).", call. = FALSE)

  # Base directory for resolving relative paths
  base_dir <- if (grepl("^https?://", path, ignore.case = TRUE)) {
    sub("/[^/]*$", "", path)
  } else {
    dirname(normalizePath(path, mustWork = FALSE))
  }

  # Extract file paths from STAC assets
  tile_paths <- character(length(features))
  vpc_meta   <- vector("list", length(features))

  for (i in seq_along(features)) {
    feat <- features[[i]]

    # Get the data asset href
    href <- NULL
    assets <- feat$assets
    if (!is.null(assets)) {
      # Prefer "data" role, fall back to first asset
      for (aname in names(assets)) {
        asset <- assets[[aname]]
        roles <- asset$roles %||% list()
        if ("data" %in% roles || aname == "data") {
          href <- asset$href
          break
        }
      }
      if (is.null(href) && length(assets) > 0L)
        href <- assets[[1L]]$href
    }

    if (is.null(href) || !nzchar(href))
      stop("Feature ", i, " has no data asset href.", call. = FALSE)

    # Resolve relative paths
    if (!grepl("^https?://|^/|^[A-Za-z]:", href)) {
      href <- file.path(base_dir, href)
      if (!grepl("^https?://", href))
        href <- normalizePath(href, mustWork = FALSE)
    }

    tile_paths[i] <- href

    # Warn about non-COPC files
    if (!grepl("\\.copc\\.la[sz]$", href, ignore.case = TRUE)) {
      warning("Feature ", i, " references a non-COPC file (",
              basename(href), "). ",
              "Only .copc.laz files are currently supported; ",
              "skipping. A lidR::readLAS() fallback is planned.",
              call. = FALSE)
      tile_paths[i] <- NA_character_
      vpc_meta[[i]] <- list(id = feat$id %||% NA_character_, skipped = TRUE)
      next
    }

    # Collect VPC-specific metadata (stats, overview, CRS, etc.)
    props <- feat$properties %||% list()
    meta <- list(
      id     = feat$id %||% NA_character_,
      bbox   = feat$bbox,
      datetime = props$datetime %||% NA_character_
    )

    # STAC pointcloud extension
    meta$pc_count  <- props$`pc:count` %||% NA_integer_
    meta$pc_type   <- props$`pc:type` %||% NA_character_
    meta$pc_schema <- props$`pc:schemas`

    # STAC projection extension
    meta$proj_epsg     <- props$`proj:epsg` %||% NA_integer_
    meta$proj_wkt2     <- props$`proj:wkt2` %||% NA_character_
    meta$proj_projjson <- props$`proj:projjson`
    meta$proj_bbox     <- props$`proj:bbox`
    meta$proj_geometry <- props$`proj:geometry`

    # Statistics (optional)
    meta$pc_statistics <- props$`pc:statistics`

    # Overview assets (may be multiple for hierarchical overviews)
    if (!is.null(assets)) {
      overview_hrefs <- character(0)
      for (aname in names(assets)) {
        asset <- assets[[aname]]
        roles <- asset$roles %||% list()
        if ("overview" %in% roles) {
          ovw_href <- asset$href
          if (!grepl("^https?://|^/|^[A-Za-z]:", ovw_href))
            ovw_href <- file.path(base_dir, ovw_href)
          overview_hrefs <- c(overview_hrefs, ovw_href)
        }
      }
      if (length(overview_hrefs) > 0L)
        meta$overviews <- overview_hrefs
    }

    vpc_meta[[i]] <- meta
  }

  # Drop non-COPC tiles that were skipped
  keep <- !is.na(tile_paths)
  if (sum(keep) == 0L)
    stop("No COPC tiles found in VPC file. Only .copc.laz files are ",
         "currently supported.", call. = FALSE)
  tile_paths <- tile_paths[keep]
  vpc_meta   <- vpc_meta[keep]

  if (progress)
    message(sprintf("VPC: %d COPC tile(s) referenced", length(tile_paths)))

  # -- Fast path: build catalog from VPC metadata without re-reading
  #    headers.  For large collections (1000s of files) this avoids

  #    HTTP HEAD + header-read round-trips for every tile.
  if (isTRUE(trust_vpc)) {
    ctg <- .vpc_build_catalog_from_meta(
      tile_paths, vpc_meta,
      chunk_size   = chunk_size,
      chunk_buffer = chunk_buffer,
      select       = select,
      filter       = filter,
      progress     = progress
    )
  } else {
    ctg <- read_copc_catalog(
      paths        = tile_paths,
      chunk_size   = chunk_size,
      chunk_buffer = chunk_buffer,
      select       = select,
      filter       = filter,
      progress     = progress
    )
  }

  # Attach VPC metadata
  attr(ctg, "vpc_metadata") <- list(
    source     = path,
    features   = vpc_meta,
    n_features = length(features)
  )

  ctg
}


#' Write a Virtual Point Cloud (VPC) file
#'
#' Creates a `.vpc` file (STAC ItemCollection / GeoJSON FeatureCollection)
#' from one or more COPC files.
#' The resulting VPC can be opened directly in QGIS or processed with
#' PDAL wrench.
#'
#' @param paths Character vector. Paths or URLs to `.copc.laz` files.
#' @param output Character. Output `.vpc` file path.
#' @param boundary Logical. If `TRUE`, compute the exact boundary polygon
#'   for each file (from the convex hull of points at LOD 0). Default
#'   `FALSE` uses the rectangular bounding box.
#' @param statistics Logical. If `TRUE`, compute per-attribute statistics
#'   (min, max, mean, count) for each file.
#'   Requires reading point data (slower). Default `FALSE`.
#' @param overview Character, logical, or `NULL`. A path to an existing
#'   overview `.copc.laz`, `TRUE` to auto-generate a single overview
#'   (1 in 1000 points), or `NULL` (default) for no overview.
#' @param overview_levels Integer. Number of hierarchical overview levels
#'   to generate when `overview = TRUE`. Level 0 covers the full extent;
#'   each subsequent level tiles into a 2x2 grid at higher density.
#'   Default `1L` (single merged overview). Values > 1 produce a pyramid.
#' @param progress Logical. Print progress messages?
#'
#' @return Invisibly returns the output file path.
#'
#' @seealso [read_copc_vpc()], [read_copc_catalog()]
#'
#' @examples
#' \donttest{
#' # Build a VPC from local COPC tiles
#' files <- list.files("tiles/", pattern = "\\.copc\\.laz$",
#'                     full.names = TRUE)
#' write_copc_vpc(files, "my_project.vpc")
#'
#' # With statistics and boundary polygons
#' write_copc_vpc(files, "my_project.vpc",
#'           boundary = TRUE, statistics = TRUE)
#' }
#'
#' @export
write_copc_vpc <- function(paths,
                      output,
                      boundary        = FALSE,
                      statistics      = FALSE,
                      overview        = NULL,
                      overview_levels = 1L,
                      progress        = TRUE) {

  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for VPC support.", call. = FALSE)

  stopifnot(is.character(paths), length(paths) >= 1L)
  stopifnot(is.character(output), length(output) == 1L, nzchar(output))

  # Ensure .vpc extension
  if (!grepl("\\.vpc$", output, ignore.case = TRUE))
    output <- paste0(output, ".vpc")

  out_dir <- dirname(normalizePath(output, mustWork = FALSE))

  features <- vector("list", length(paths))

  for (i in seq_along(paths)) {
    p <- paths[i]
    if (progress) message(sprintf("  [%d/%d] %s", i, length(paths),
                                  basename(sub("[?#].*$", "", p))))

    hdr <- read_copc_header(p)

    # -- Spatial extent -------------------------------------------------
    xmin <- hdr[["Min X"]]; xmax <- hdr[["Max X"]]
    ymin <- hdr[["Min Y"]]; ymax <- hdr[["Max Y"]]
    zmin <- hdr[["Min Z"]]; zmax <- hdr[["Max Z"]]

    # -- CRS -----------------------------------------------------------
    crs_wkt  <- NULL
    crs_epsg <- NA_integer_
    crs_projjson <- NULL
    is_compound  <- FALSE

    for (v in hdr[["Variable Length Records"]]) {
      if (!is.null(v[["WKT OGC COORDINATE SYSTEM"]])) {
        crs_wkt <- v[["WKT OGC COORDINATE SYSTEM"]]
        break
      }
    }

    if (!is.null(crs_wkt) && requireNamespace("sf", quietly = TRUE)) {
      crs_obj <- tryCatch(sf::st_crs(crs_wkt), error = function(e) NULL)
      if (!is.null(crs_obj)) {
        # Detect compound CRS (COMPD_CS / COMPOUNDCRS)
        is_compound <- grepl("COMPD_CS|COMPOUNDCRS", crs_wkt,
                             ignore.case = TRUE)
        if (!is_compound && !is.na(crs_obj$epsg))
          crs_epsg <- crs_obj$epsg

        # Generate proj:projjson via sf (PROJ >= 6)
        crs_projjson <- tryCatch({
          pj <- sf::st_crs(crs_obj)$proj4string  # ensure valid
          # sf exposes projjson via $wkt which is WKT2; build from input_obj
          raw_json <- sf::st_crs(crs_obj)$input
          if (requireNamespace("jsonlite", quietly = TRUE)) {
            # Use sf's GDAL/PROJ backend to get PROJJSON
            pj_json <- tryCatch(
              jsonlite::fromJSON(
                system2("projinfo", c("-o", "projjson", "-q",
                        shQuote(crs_wkt)), stdout = TRUE,
                        stderr = FALSE) |> paste(collapse = "\n"),
                simplifyVector = FALSE),
              error = function(e) NULL
            )
            pj_json
          }
        }, error = function(e) NULL)
      }
    }

    # -- WGS 84 bounding box (required by GeoJSON) --------------------
    bbox_wgs84 <- c(xmin, ymin, zmin, xmax, ymax, zmax)
    geom_wgs84 <- list(
      type = "Polygon",
      coordinates = list(list(
        c(xmin, ymin), c(xmax, ymin), c(xmax, ymax),
        c(xmin, ymax), c(xmin, ymin)
      ))
    )

    if (!is.null(crs_wkt) && requireNamespace("sf", quietly = TRUE)) {
      # Transform bbox corners to WGS 84
      tr <- tryCatch({
        pts <- sf::st_sfc(
          sf::st_point(c(xmin, ymin)),
          sf::st_point(c(xmax, ymin)),
          sf::st_point(c(xmax, ymax)),
          sf::st_point(c(xmin, ymax)),
          crs = sf::st_crs(crs_wkt)
        )
        pts84 <- sf::st_transform(pts, 4326)
        coords84 <- sf::st_coordinates(pts84)
        list(
          bbox = c(min(coords84[,1]), min(coords84[,2]), zmin,
                   max(coords84[,1]), max(coords84[,2]), zmax),
          geom = list(
            type = "Polygon",
            coordinates = list(list(
              c(coords84[1,1], coords84[1,2]),
              c(coords84[2,1], coords84[2,2]),
              c(coords84[3,1], coords84[3,2]),
              c(coords84[4,1], coords84[4,2]),
              c(coords84[1,1], coords84[1,2])
            ))
          )
        )
      }, error = function(e) NULL)

      if (!is.null(tr)) {
        bbox_wgs84 <- tr$bbox
        geom_wgs84 <- tr$geom
      }
    }

    # -- Boundary polygon (optional) -----------------------------------
    boundary_geom <- NULL
    if (boundary && requireNamespace("sf", quietly = TRUE)) {
      boundary_geom <- tryCatch({
        preview <- read_copc(p, max_depth = 0L, progress = FALSE)
        if (nrow(preview$data) > 2L) {
          pts_sf <- sf::st_as_sf(preview$data, coords = c("X", "Y"),
                                 remove = FALSE)
          hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts_sf)))
          # Convert to GeoJSON-like list
          coords <- sf::st_coordinates(hull)
          list(
            type = "Polygon",
            coordinates = list(lapply(seq_len(nrow(coords)), function(j)
              c(coords[j, 1], coords[j, 2])))
          )
        }
      }, error = function(e) NULL)
    }

    # -- Point cloud schema --------------------------------------------
    pc_schema <- .vpc_schema_from_header(hdr)

    # -- Statistics (optional) -----------------------------------------
    pc_statistics <- NULL
    if (statistics) {
      pc_statistics <- tryCatch({
        result <- read_copc(p, max_depth = 3L, progress = FALSE)
        if (nrow(result$data) > 0L)
          .vpc_compute_statistics(result$data)
      }, error = function(e) NULL)
    }

    # -- Asset href: relative to VPC location --------------------------
    asset_href <- p
    if (!grepl("^https?://", p)) {
      p_norm <- normalizePath(p, mustWork = FALSE)
      rel <- tryCatch({
        # Compute relative path from VPC dir to file
        .relative_path(p_norm, out_dir)
      }, error = function(e) p_norm)
      asset_href <- gsub("\\\\", "/", rel)
    }

    # -- Properties ----------------------------------------------------
    properties <- list(
      datetime = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
    )

    # STAC pointcloud extension
    properties[["pc:count"]]  <- as.integer(hdr[["Number of point records"]])
    properties[["pc:type"]]   <- "lidar"
    properties[["pc:schemas"]] <- pc_schema
    properties[["pc:indexed"]] <- TRUE   # COPC files always have a spatial index (octree)
    if (!is.null(pc_statistics))
      properties[["pc:statistics"]] <- pc_statistics

    # STAC projection extension
    # For compound CRS the single EPSG code is ambiguous (spec recommends
    # proj:wkt2 or proj:projjson instead).
    if (!is.na(crs_epsg) && !is_compound)
      properties[["proj:epsg"]] <- crs_epsg
    if (!is.null(crs_wkt))
      properties[["proj:wkt2"]] <- crs_wkt
    if (!is.null(crs_projjson))
      properties[["proj:projjson"]] <- crs_projjson

    # Native-CRS bbox (3D — strongly recommended by VPC spec)
    properties[["proj:bbox"]] <- c(xmin, ymin, zmin, xmax, ymax, zmax)

    # Native-CRS geometry (optional precise boundary)
    if (!is.null(boundary_geom))
      properties[["proj:geometry"]] <- boundary_geom

    # -- Assets --------------------------------------------------------
    assets <- list(
      data = list(
        href  = asset_href,
        roles = list("data")
      )
    )

    # -- Build feature -------------------------------------------------
    feature <- list(
      type       = "Feature",
      stac_version = "1.0.0",
      stac_extensions = list(
        "https://stac-extensions.github.io/pointcloud/v1.0.0/schema.json",
        "https://stac-extensions.github.io/projection/v1.1.0/schema.json"
      ),
      id         = tools::file_path_sans_ext(basename(sub("[?#].*$", "", p))),
      geometry   = if (!is.null(boundary_geom)) boundary_geom else geom_wgs84,
      bbox       = bbox_wgs84,
      properties = properties,
      assets     = assets
    )

    features[[i]] <- feature
  }

  # -- Build the ItemCollection ----------------------------------------
  vpc <- list(
    type     = "FeatureCollection",
    features = features
  )

  # -- Overviews --------------------------------------------------------
  if (!is.null(overview) && !identical(overview, FALSE)) {
    overview_levels <- max(1L, as.integer(overview_levels))

    if (isTRUE(overview)) {
      overview_assets <- .vpc_generate_overviews(
        paths, output, out_dir, overview_levels, progress
      )
    } else if (is.character(overview) && nzchar(overview)) {
      # Single user-supplied overview path
      ovw_rel <- if (!grepl("^https?://", overview)) {
        gsub("\\\\", "/", .relative_path(
          normalizePath(overview, mustWork = FALSE), out_dir))
      } else {
        overview
      }
      overview_assets <- list(
        list(key = "overview", href = ovw_rel, roles = list("overview"))
      )
    } else {
      overview_assets <- list()
    }

    # Attach overview assets to every feature
    for (oa in overview_assets) {
      for (i in seq_along(vpc$features)) {
        vpc$features[[i]]$assets[[oa$key]] <- list(
          href  = oa$href,
          roles = oa$roles
        )
      }
    }
  }

  # Write JSON
  json_text <- jsonlite::toJSON(vpc, auto_unbox = TRUE, pretty = TRUE,
                                digits = 10)
  writeLines(json_text, output)

  if (progress)
    message(sprintf("VPC written: %s (%d tiles)", output, length(paths)))

  invisible(output)
}


# -- VPC internal helpers --------------------------------------------------

#' Generate hierarchical overviews for a VPC
#'
#' Level 0 = single merged overview (every 1000th point).
#' Level k = 4^k tiles, each sampling every max(1, 1000 / 4^k)th point
#' from the files that overlap that tile.
#' @noRd
.vpc_generate_overviews <- function(paths, output, out_dir,
                                    levels = 1L, progress = TRUE) {
  base_name <- sub("\\.vpc$", "", basename(output), ignore.case = TRUE)
  ovw_dir   <- file.path(out_dir, paste0(base_name, "_overviews"))
  if (!dir.exists(ovw_dir)) dir.create(ovw_dir, recursive = TRUE)

  # Collect global extent
  hdrs <- lapply(paths, read_copc_header)
  global_xmin <- min(vapply(hdrs, `[[`, 0, "Min X"))
  global_ymin <- min(vapply(hdrs, `[[`, 0, "Min Y"))
  global_xmax <- max(vapply(hdrs, `[[`, 0, "Max X"))
  global_ymax <- max(vapply(hdrs, `[[`, 0, "Max Y"))

  overview_assets <- list()

  for (lvl in seq_len(levels) - 1L) {
    n_side <- 2L^lvl          # 1, 2, 4, 8, ...
    n_tiles <- n_side * n_side
    dx <- (global_xmax - global_xmin) / n_side
    dy <- (global_ymax - global_ymin) / n_side
    thin_nth <- max(1L, as.integer(1000 / (4L^lvl)))

    if (progress) message(sprintf("  Overview level %d: %d tile(s), every %dth point",
                                  lvl, n_tiles, thin_nth))

    for (tx in seq_len(n_side) - 1L) {
      for (ty in seq_len(n_side) - 1L) {
        tile_xmin <- global_xmin + tx * dx
        tile_ymin <- global_ymin + ty * dy
        tile_xmax <- tile_xmin + dx
        tile_ymax <- tile_ymin + dy
        tile_bbox <- c(tile_xmin, tile_ymin, tile_xmax, tile_ymax)

        ovw_file <- file.path(ovw_dir,
          sprintf("%s-overview-L%d_%d_%d.copc.laz", base_name, lvl, tx, ty))

        # Read thinned data from all files overlapping this tile
        tile_data <- list()
        for (p in paths) {
          samp <- tryCatch(
            read_copc(p, bbox = tile_bbox,
                      filter = paste0("-keep_every_nth ", thin_nth),
                      progress = FALSE),
            error = function(e) NULL)
          if (!is.null(samp) && nrow(samp$data) > 0L)
            tile_data[[length(tile_data) + 1L]] <- samp
        }

        if (length(tile_data) == 0L) next

        merged <- tile_data[[1L]]
        if (length(tile_data) > 1L)
          merged$data <- do.call(rbind, lapply(tile_data, `[[`, "data"))

        tryCatch(
          write_copc(merged, ovw_file, progress = FALSE),
          error = function(e) {
            warning("Overview tile failed: ", conditionMessage(e),
                    call. = FALSE)
          })

        if (file.exists(ovw_file)) {
          ovw_rel <- gsub("\\\\", "/",
            .relative_path(normalizePath(ovw_file, mustWork = FALSE), out_dir))
          key <- if (lvl == 0L && n_tiles == 1L) "overview"
                 else sprintf("overview-L%d_%d_%d", lvl, tx, ty)
          overview_assets[[length(overview_assets) + 1L]] <- list(
            key = key, href = ovw_rel, roles = list("overview"))
        }
      }
    }
  }

  overview_assets
}

#' Build a pc:schemas array from a COPC header
#' @noRd
.vpc_schema_from_header <- function(hdr) {
  fmt <- hdr[["Point Data Format ID"]]

  # Base LAS 1.4 format 6 attributes
  schema <- list(
    list(name = "X", size = 8, type = "floating"),
    list(name = "Y", size = 8, type = "floating"),
    list(name = "Z", size = 8, type = "floating"),
    list(name = "Intensity", size = 2, type = "unsigned"),
    list(name = "ReturnNumber", size = 1, type = "unsigned"),
    list(name = "NumberOfReturns", size = 1, type = "unsigned"),
    list(name = "ScanDirectionFlag", size = 1, type = "unsigned"),
    list(name = "EdgeOfFlightLine", size = 1, type = "unsigned"),
    list(name = "Classification", size = 1, type = "unsigned"),
    list(name = "ScanAngleRank", size = 4, type = "floating"),
    list(name = "UserData", size = 1, type = "unsigned"),
    list(name = "PointSourceID", size = 2, type = "unsigned"),
    list(name = "GpsTime", size = 8, type = "floating")
  )

  # Format 7: add RGB

  if (fmt >= 7L) {
    schema <- c(schema, list(
      list(name = "Red", size = 2, type = "unsigned"),
      list(name = "Green", size = 2, type = "unsigned"),
      list(name = "Blue", size = 2, type = "unsigned")
    ))
  }

  # Format 8: add NIR
  if (fmt >= 8L) {
    schema <- c(schema, list(
      list(name = "NIR", size = 2, type = "unsigned")
    ))
  }

  schema
}


#' Compute per-attribute statistics for VPC pc:statistics
#' @noRd
.vpc_compute_statistics <- function(dt) {
  stats <- list()
  numeric_cols <- names(dt)[vapply(dt, is.numeric, FALSE)]
  for (col in numeric_cols) {
    vals <- dt[[col]]
    stats[[length(stats) + 1L]] <- list(
      name     = col,
      position = length(stats),
      average  = mean(vals, na.rm = TRUE),
      count    = sum(!is.na(vals)),
      maximum  = max(vals, na.rm = TRUE),
      minimum  = min(vals, na.rm = TRUE),
      stddev   = stats::sd(vals, na.rm = TRUE),
      variance = stats::var(vals, na.rm = TRUE)
    )
  }
  stats
}


#' Compute a relative path from target to base
#' @noRd
.relative_path <- function(target, base) {
  # Normalise separators
  target <- gsub("\\\\", "/", target)
  base   <- gsub("\\\\", "/", base)

  # Ensure base ends with /
  if (!grepl("/$", base)) base <- paste0(base, "/")

  # If the target starts with the base, strip it
  if (startsWith(target, base)) {
    return(substring(target, nchar(base) + 1L))
  }

  # Otherwise use a simple ./relative approach
  # Split paths and find common prefix
  t_parts <- strsplit(target, "/")[[1L]]
  b_parts <- strsplit(base, "/")[[1L]]
  # Remove empty trailing from base
  b_parts <- b_parts[nzchar(b_parts)]

  common <- 0L
  for (j in seq_len(min(length(t_parts), length(b_parts)))) {
    if (tolower(t_parts[j]) == tolower(b_parts[j])) {
      common <- j
    } else {
      break
    }
  }

  ups <- length(b_parts) - common
  rest <- t_parts[(common + 1L):length(t_parts)]
  paste(c(rep("..", ups), rest), collapse = "/")
}


#' Null-coalescing operator (safe re-definition)
#' @noRd
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a


#' Build a copc_catalog from VPC metadata (no header re-reads)
#'
#' The "fast path" used by `read_copc_vpc(trust_vpc = TRUE)`.  Constructs the
#' same catalog object that `read_copc_catalog()` would produce, but
#' extracts spatial extent and CRS from the already-parsed STAC
#' properties rather than performing a header read for every tile.
#' @noRd
.vpc_build_catalog_from_meta <- function(tile_paths, vpc_meta,
                                         chunk_size, chunk_buffer,
                                         select, filter, progress) {
  n <- length(tile_paths)

  # Extract bbox components from VPC metadata (proj:bbox is native CRS,
  # 6-element: xmin ymin zmin xmax ymax zmax).  Fall back to WGS 84
  # bbox (which may be geographic) if proj:bbox is missing.
  xmin <- ymin <- zmin <- xmax <- ymax <- zmax <- numeric(n)
  npts <- numeric(n)
  pdrf <- integer(n)

  for (i in seq_len(n)) {
    m <- vpc_meta[[i]]
    pb <- unlist(m$proj_bbox %||% m$bbox)
    if (!is.null(pb) && length(pb) >= 6L) {
      xmin[i] <- pb[[1]]; ymin[i] <- pb[[2]]; zmin[i] <- pb[[3]]
      xmax[i] <- pb[[4]]; ymax[i] <- pb[[5]]; zmax[i] <- pb[[6]]
    } else if (!is.null(pb) && length(pb) >= 4L) {
      xmin[i] <- pb[[1]]; ymin[i] <- pb[[2]]
      xmax[i] <- pb[[3]]; ymax[i] <- pb[[4]]
      zmin[i] <- NA_real_; zmax[i] <- NA_real_
    }
    npts[i] <- as.numeric(m$pc_count %||% NA_real_)

    # Attempt to recover PDRF from pc:schemas (presence of Red → fmt 7+)
    schema_names <- if (!is.null(m$pc_schema))
      vapply(m$pc_schema, function(s) s$name %||% "", "") else character(0)
    pdrf[i] <- if ("NIR" %in% schema_names) 8L
               else if ("Red" %in% schema_names) 7L
               else 6L
  }

  file_info <- data.frame(
    filename                = tile_paths,
    Min.X                   = xmin,
    Max.X                   = xmax,
    Min.Y                   = ymin,
    Max.Y                   = ymax,
    Min.Z                   = zmin,
    Max.Z                   = zmax,
    Number.of.point.records = npts,
    Point.Data.Format.ID    = pdrf,
    stringsAsFactors        = FALSE
  )

  # -- Attempt to build a lidR LAScatalog (same logic as read_copc_catalog)
  if (requireNamespace("lidR", quietly = TRUE) &&
      requireNamespace("sf", quietly = TRUE)) {
    ctg <- tryCatch({
      # Recover CRS from VPC metadata (use first tile's WKT or projjson)
      crs <- sf::NA_crs_
      for (m in vpc_meta) {
        wkt <- m$proj_wkt2
        if (!is.null(wkt) && nzchar(wkt)) {
          crs <- tryCatch(sf::st_crs(wkt), error = function(e) sf::NA_crs_)
          break
        }
      }

      polys <- lapply(seq_len(n), function(i) {
        fi <- file_info[i, ]
        sf::st_polygon(list(matrix(c(
          fi$Min.X, fi$Min.Y,
          fi$Max.X, fi$Min.Y,
          fi$Max.X, fi$Max.Y,
          fi$Min.X, fi$Max.Y,
          fi$Min.X, fi$Min.Y
        ), ncol = 2, byrow = TRUE)))
      })
      geom <- sf::st_sfc(polys, crs = crs)
      file_sf <- sf::st_sf(file_info, geometry = geom)

      ctg_obj <- methods::new("LAScatalog")
      ctg_obj@data <- file_sf
      lidR::opt_chunk_size(ctg_obj)   <- chunk_size
      lidR::opt_chunk_buffer(ctg_obj) <- chunk_buffer
      lidR::opt_select(ctg_obj)       <- select
      lidR::opt_filter(ctg_obj)       <- filter
      ctg_obj
    }, error = function(e) NULL)
    if (!is.null(ctg)) return(ctg)
  }

  # -- Fallback: lightweight S3 catalog
  structure(
    list(
      files   = file_info,
      headers = NULL,
      options = list(
        chunk_size   = chunk_size,
        chunk_buffer = chunk_buffer,
        select       = select,
        filter       = filter
      )
    ),
    class = "copc_catalog"
  )
}
