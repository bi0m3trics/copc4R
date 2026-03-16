# -- Planetary Computer: USGS 3DEP COPC tile download -------------------------
#
# Search and download Cloud Optimized Point Clouds (COPC) from
# Microsoft's Planetary Computer STAC API for the USGS 3DEP collection.
# Only the points within the user's AOI are fetched (via COPC range reads)
# and written to local .copc.laz files -- no full-tile downloads.
#
# Data info : https://planetarycomputer.microsoft.com/dataset/3dep-lidar-copc
# STAC API  : https://planetarycomputer.microsoft.com/api/stac/v1
# ---------------------------------------------------------------------------

# Null-coalescing helper (guarded against multiple definition)
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Download 3DEP COPC point clouds from Microsoft Planetary Computer
#'
#' Queries the
#' \href{https://planetarycomputer.microsoft.com/dataset/3dep-lidar-copc}{Microsoft Planetary Computer STAC API}
#' for USGS 3DEP LiDAR COPC tiles that intersect a user-supplied \code{sf}
#' geometry, reads \emph{only the points inside the AOI} from each tile
#' via COPC range reads, and writes the clipped results as \code{.copc.laz}
#' files in a user-specified directory.
#'
#' @details
#' The workflow is:
#' \enumerate{
#'   \item The \code{aoi} geometry is reprojected to WGS 84 (EPSG:4326)
#'     and converted to valid GeoJSON for the STAC \code{intersects}
#'     query parameter.
#'   \item All matching STAC items are retrieved (with automatic
#'     pagination).
#'   \item For each tile the AOI is reprojected into the tile's native
#'     CRS (read from the COPC header) and passed to
#'     \code{\link{read_copc}()} as the \code{aoi} argument.  Only the
#'     octree nodes that overlap the AOI bounding box are fetched via
#'     HTTP range reads, and points are then clipped to the polygon
#'     boundary -- so only the data you need is transferred.
#'   \item The clipped points are written to a local \code{.copc.laz}
#'     file via \code{\link{write_copc}()}.
#' }
#'
#' @param aoi An \code{sf} or \code{sfc} object (polygon / multipolygon).
#'   The geometry is unioned, converted to a single polygon/multipolygon,
#'   and sent as GeoJSON to the STAC \code{intersects} parameter.
#' @param dest_dir Character. Directory where \code{.copc.laz} files
#'   will be saved.  Created recursively if it does not exist.
#' @param select Character. Column selection string passed to
#'   \code{\link{read_copc}()} (default \code{"*"} = all attributes).
#'   See \code{\link{read_copc}} for the full syntax.
#' @param filter Character. Attribute filter string passed to
#'   \code{\link{read_copc}()} (default \code{NULL} = no filter).
#' @param max_points Numeric. Maximum number of points to read per tile
#'   (default \code{Inf}).
#' @param limit Integer. Maximum number of STAC items returned per
#'   request page (default 100, max 1000).
#' @param overwrite Logical. If \code{TRUE}, re-download tiles that
#'   already exist locally. Default \code{FALSE}.
#' @param merge Logical.  If \code{TRUE} (the default), all tiles are
#'   merged into a single output file named
#'   \code{3dep_copc_merge.copc.laz}.  If \code{FALSE}, one file per
#'   tile is written (named by STAC item ID).
#' @param progress Logical. Show progress during range reads
#'   (default \code{TRUE}).
#' @param verbose Logical. Print informational messages during the
#'   search and download steps (default \code{TRUE}).
#'
#' @return A \code{data.frame} (invisibly) with one row per output file
#'   and the following columns:
#'   \describe{
#'     \item{item_id}{STAC item identifier(s) -- comma-separated when
#'       \code{merge = TRUE}.}
#'     \item{file}{Full local path to the written \code{.copc.laz} file.}
#'     \item{n_points}{Number of points in the file.}
#'     \item{written}{Logical: \code{TRUE} if the file was written in
#'       this call.}
#'   }
#'
#' @examples
#' \donttest{
#' library(sf)
#' library(copc4R)
#'
#' # -- Small AOI near Flagstaff, AZ ------------------------------------
#' coords <- matrix(c(
#'   -111.672, 35.170,
#'   -111.664, 35.170,
#'   -111.664, 35.172,
#'   -111.672, 35.172,
#'   -111.672, 35.170
#' ), ncol = 2, byrow = TRUE)
#'
#' aoi <- st_sf(
#'   id       = 1,
#'   geometry = st_sfc(st_polygon(list(coords)), crs = 4326)
#' )
#'
#' # Fetch only the AOI portion and save as a single merged .laz
#' result <- download_3dep_copc(aoi, dest_dir = "~/3dep_tiles")
#' print(result)
#'
#' # Read the local file back in
#' las <- read_copc(result$file[1])
#'
#' # Or keep one file per source tile
#' result2 <- download_3dep_copc(aoi, dest_dir = "~/3dep_tiles",
#'                               merge = FALSE, overwrite = TRUE)
#' }
#'
#' @export
download_3dep_copc <- function(aoi,
                               dest_dir,
                               select     = "*",
                               filter     = NULL,
                               max_points = Inf,
                               limit      = 100L,
                               overwrite  = FALSE,
                               merge      = TRUE,
                               progress   = TRUE,
                               verbose    = TRUE) {
  # -- Check dependencies --------------------------------------------------
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required. Install with: install.packages('sf')",
         call. = FALSE)
  if (!requireNamespace("httr", quietly = TRUE))
    stop("Package 'httr' is required. Install with: install.packages('httr')",
         call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')",
         call. = FALSE)
  if (!requireNamespace("lidR", quietly = TRUE))
    stop("Package 'lidR' is required (for as_las() conversion). ",
         "Install with: install.packages('lidR')",
         call. = FALSE)

  # -- Validate inputs ----------------------------------------------------
  stopifnot(inherits(aoi, c("sf", "sfc")))
  stopifnot(is.character(dest_dir), length(dest_dir) == 1L, nzchar(dest_dir))
  limit <- as.integer(min(max(limit, 1L), 1000L))

  # -- Prepare output directory --------------------------------------------
  dest_dir <- normalizePath(dest_dir, mustWork = FALSE)
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
    if (verbose) message("Created directory: ", dest_dir)
  }

  # -- Convert AOI to WGS 84 GeoJSON geometry -----------------------------
  aoi_geom <- sf::st_geometry(aoi)
  if (is.na(sf::st_crs(aoi_geom)))
    stop("The 'aoi' object has no CRS. Assign one with sf::st_set_crs().",
         call. = FALSE)

  aoi_4326 <- sf::st_transform(aoi_geom, 4326)
  aoi_union <- sf::st_union(aoi_4326)
  aoi_union <- sf::st_cast(aoi_union, "MULTIPOLYGON")

  # Build GeoJSON geometry list from sf coordinates
  intersects_geom <- .sf_to_geojson_geometry(aoi_union)

  # -- Search STAC API (with pagination) ----------------------------------
  stac_url <- "https://planetarycomputer.microsoft.com/api/stac/v1/search"
  all_items <- list()

  body <- list(
    collections = list("3dep-lidar-copc"),
    intersects  = intersects_geom,
    limit       = limit
  )

  repeat {
    resp <- httr::POST(
      stac_url,
      httr::content_type_json(),
      body = jsonlite::toJSON(body, auto_unbox = TRUE, digits = NA)
    )
    httr::stop_for_status(resp, task = "search Planetary Computer STAC API")

    result <- jsonlite::fromJSON(
      httr::content(resp, as = "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )

    items <- result$features
    all_items <- c(all_items, items)
    if (verbose)
      message("STAC search: fetched ", length(items), " item(s) ",
              "(total so far: ", length(all_items), ")")

    # Check for a "next" link (pagination)
    next_link <- NULL
    for (lnk in result$links %||% list()) {
      if (identical(lnk$rel, "next")) {
        next_link <- lnk
        break
      }
    }
    if (is.null(next_link)) break

    if (!is.null(next_link$body)) {
      body <- modifyList(body, next_link$body)
    } else if (!is.null(next_link$href)) {
      resp <- httr::GET(next_link$href)
      httr::stop_for_status(resp)
      result <- jsonlite::fromJSON(
        httr::content(resp, as = "text", encoding = "UTF-8"),
        simplifyVector = FALSE
      )
      items <- result$features
      all_items <- c(all_items, items)
      if (verbose)
        message("STAC search: fetched ", length(items), " item(s) ",
                "(total so far: ", length(all_items), ")")
      next_link <- NULL
      for (lnk in result$links %||% list()) {
        if (identical(lnk$rel, "next")) {
          next_link <- lnk
          break
        }
      }
      if (is.null(next_link)) break
    } else {
      break
    }
  }

  n_items <- length(all_items)
  if (n_items == 0L) {
    message("No 3DEP COPC tiles found for the supplied AOI.")
    return(invisible(data.frame(
      item_id  = character(0), file = character(0),
      n_points = integer(0),   written = logical(0),
      stringsAsFactors = FALSE
    )))
  }

  if (verbose) message("Found ", n_items, " COPC tile(s) intersecting the AOI.")

  # -- Obtain a SAS token -------------------------------------------------
  sas_token <- .pc_get_sas_token()

  # -- Read only the AOI portion from each tile via COPC range reads ------
  all_las   <- list()
  item_ids  <- character(0)

  for (i in seq_len(n_items)) {
    item    <- all_items[[i]]
    item_id <- item$id

    href <- item$assets$data$href
    if (is.null(href)) {
      warning("Item '", item_id, "' has no 'data' asset \u2014 skipping.",
              call. = FALSE)
      next
    }

    sep <- if (grepl("\\?", href)) "&" else "?"
    signed_url <- paste0(href, sep, sas_token)

    if (verbose)
      message(sprintf("[%d/%d] Reading AOI from %s ...", i, n_items, item_id))

    # Reproject the AOI into the tile's native CRS so read_copc()
    # can use it directly for spatial querying + clipping.
    tile_hdr <- tryCatch(read_copc_header(signed_url),
                         error = function(e) NULL)
    if (is.null(tile_hdr)) {
      warning("Could not read header for '", item_id, "' \u2014 skipping.",
              call. = FALSE)
      next
    }

    tile_crs <- .crs_from_copc_header(tile_hdr)
    if (is.na(tile_crs)) {
      # Fallback: try STAC metadata
      tile_crs <- .crs_from_stac_item(item)
    }
    if (is.na(tile_crs)) {
      warning("Cannot determine CRS for '", item_id, "' \u2014 skipping.",
              call. = FALSE)
      next
    }

    aoi_proj <- sf::st_transform(aoi_union, tile_crs)

    res <- tryCatch(
      read_copc(
        path_or_url = signed_url,
        aoi         = aoi_proj,
        select      = select,
        filter      = filter,
        max_points  = max_points,
        progress    = progress
      ),
      error = function(e) {
        warning("read_copc() failed for '", item_id, "': ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )

    if (is.null(res) || nrow(res$data) == 0L) {
      if (verbose) message("  0 points in AOI \u2014 skipping.")
      next
    }

    las <- as_las(res)
    if (verbose)
      message(sprintf("  %s points read.",
                      format(lidR::npoints(las), big.mark = ",")))

    if (merge) {
      all_las[[length(all_las) + 1L]] <- las
      item_ids <- c(item_ids, item_id)
    } else {
      # Write one file per tile
      local_name <- paste0(item_id, ".copc.laz")
      local_path <- file.path(dest_dir, local_name)
      skip <- file.exists(local_path) && !overwrite
      if (skip) {
        if (verbose) message("  File exists \u2014 skipping write.")
      } else {
        write_copc(las, local_path, progress = FALSE)
        if (verbose) message("  Wrote ", local_path)
      }
      all_las[[length(all_las) + 1L]] <- data.frame(
        item_id  = item_id,
        file     = local_path,
        n_points = lidR::npoints(las),
        written  = !skip,
        stringsAsFactors = FALSE
      )
    }
  }

  # -- Write output -------------------------------------------------------
  if (merge) {
    if (length(all_las) == 0L) {
      message("No points found within the AOI across any tile.")
      return(invisible(data.frame(
        item_id = character(0), file = character(0),
        n_points = integer(0), written = logical(0),
        stringsAsFactors = FALSE
      )))
    }

    # Merge LAS objects
    merged <- do.call(rbind, all_las)
    local_path <- file.path(dest_dir, "3dep_copc_merge.copc.laz")
    skip <- file.exists(local_path) && !overwrite
    if (skip) {
      if (verbose) message("Merged file exists \u2014 skipping write.")
    } else {
      write_copc(merged, local_path, progress = FALSE)
      if (verbose)
        message(sprintf("Wrote %s points to %s",
                        format(lidR::npoints(merged), big.mark = ","),
                        local_path))
    }

    result_df <- data.frame(
      item_id  = paste(item_ids, collapse = ", "),
      file     = local_path,
      n_points = lidR::npoints(merged),
      written  = !skip,
      stringsAsFactors = FALSE
    )
  } else {
    result_df <- do.call(rbind, all_las)
  }

  if (is.null(result_df) || nrow(result_df) == 0L) {
    warning("No tiles produced output.", call. = FALSE)
    return(invisible(data.frame(
      item_id = character(0), file = character(0),
      n_points = integer(0), written = logical(0),
      stringsAsFactors = FALSE
    )))
  }

  rownames(result_df) <- NULL
  if (verbose) message("Done. Files saved to: ", dest_dir)
  invisible(result_df)
}


# -- Internal helpers ----------------------------------------------------------

#' @noRd
.pc_get_sas_token <- function() {
  sas_url <- "https://planetarycomputer.microsoft.com/api/sas/v1/token/3dep-lidar-copc"
  resp <- httr::GET(sas_url)
  httr::stop_for_status(resp, task = "obtain Planetary Computer SAS token")
  token_data <- jsonlite::fromJSON(
    httr::content(resp, as = "text", encoding = "UTF-8")
  )
  token_data$token
}

#' Extract CRS from a COPC header (WKT VLR)
#' @noRd
.crs_from_copc_header <- function(hdr) {
  vlrs <- hdr[["Variable Length Records"]]
  for (v in vlrs) {
    wkt <- v[["WKT OGC COORDINATE SYSTEM"]]
    if (!is.null(wkt) && nzchar(wkt)) {
      crs <- tryCatch(sf::st_crs(wkt), error = function(e) sf::st_crs(NA))
      if (!is.na(crs)) return(crs)
    }
  }
  sf::st_crs(NA)
}

#' Extract CRS from STAC item projection metadata
#' @noRd
.crs_from_stac_item <- function(item) {
  props <- item$properties

  # 1. proj:epsg

  epsg <- props[["proj:epsg"]]
  if (!is.null(epsg) && !is.na(epsg)) {
    crs <- tryCatch(sf::st_crs(as.integer(epsg)),
                    error = function(e) sf::st_crs(NA))
    if (!is.na(crs)) return(crs)
  }

  # 2. proj:projjson -- extract horizontal component EPSG
  projjson <- props[["proj:projjson"]]
  if (!is.null(projjson) && !is.null(projjson$components)) {
    horiz_epsg <- projjson$components[[1]]$id$code
    if (!is.null(horiz_epsg)) {
      crs <- tryCatch(sf::st_crs(as.integer(horiz_epsg)),
                      error = function(e) sf::st_crs(NA))
      if (!is.na(crs)) return(crs)
    }
  }

  # 3. proj:wkt2
  wkt2 <- props[["proj:wkt2"]]
  if (!is.null(wkt2) && nzchar(wkt2)) {
    crs <- tryCatch(sf::st_crs(wkt2), error = function(e) sf::st_crs(NA))
    if (!is.na(crs)) return(crs)
  }

  sf::st_crs(NA)
}

#' Convert an sfc geometry to a GeoJSON geometry list
#' @noRd
.sf_to_geojson_geometry <- function(sfc) {
  geom <- sfc[[1]]
  geom_type <- sf::st_geometry_type(geom)

  if (geom_type == "POLYGON") {
    coords <- unclass(geom)
    list(
      type        = "Polygon",
      coordinates = lapply(coords, function(ring) {
        lapply(seq_len(nrow(ring)), function(j) as.numeric(ring[j, ]))
      })
    )
  } else if (geom_type == "MULTIPOLYGON") {
    polys <- unclass(geom)
    list(
      type        = "MultiPolygon",
      coordinates = lapply(polys, function(poly) {
        lapply(poly, function(ring) {
          lapply(seq_len(nrow(ring)), function(j) as.numeric(ring[j, ]))
        })
      })
    )
  } else {
    stop("Unsupported geometry type '", geom_type,
         "' for GeoJSON conversion. Use a POLYGON or MULTIPOLYGON AOI.",
         call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# .normalise_aoi() -- robustly coerce any aoi input to a valid sfc polygon.
# Returns a list:
#   $type      : "bbox" | "polygon"
#   $sfc       : sfc (polygon/multipolygon), or NULL for bbox-only
#   $bbox      : numeric(4) c(xmin, ymin, xmax, ymax)
#   $needs_clip: logical
# @noRd
.normalise_aoi <- function(aoi, buffer = NULL, caller = "read_copc") {
  if (is.null(aoi)) return(NULL)

  # --- numeric(4) bbox shorthand ---
  if (is.numeric(aoi)) {
    if (length(aoi) != 4L || anyNA(aoi) || any(!is.finite(aoi)))
      stop(caller, "(): when 'aoi' is numeric it must be a finite vector of ",
           "length 4: c(xmin, ymin, xmax, ymax).", call. = FALSE)
    if (aoi[1L] >= aoi[3L] || aoi[2L] >= aoi[4L])
      stop(caller, "(): bbox xmin >= xmax or ymin >= ymax.", call. = FALSE)
    return(list(type = "bbox", sfc = NULL, bbox = aoi, needs_clip = FALSE))
  }

  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for spatial 'aoi' arguments.", call. = FALSE)

  # --- coerce to sfc ---
  sfc <- tryCatch({
    if (inherits(aoi, "sfg")) {
      sf::st_sfc(aoi)
    } else if (inherits(aoi, c("sf", "sfc"))) {
      sf::st_geometry(aoi)
    } else {
      stop("Cannot coerce 'aoi' to an sf geometry. ",
           "Provide an sf, sfc, sfg, or numeric(4) object.", call. = FALSE)
    }
  }, error = function(e)
    stop(caller, "(): invalid 'aoi' -- ", conditionMessage(e), call. = FALSE)
  )

  if (length(sfc) == 0L)
    stop(caller, "(): 'aoi' contains zero features.", call. = FALSE)

  # Drop empty geometries
  is_empty <- sf::st_is_empty(sfc)
  if (all(is_empty))
    stop(caller, "(): all geometries in 'aoi' are empty.", call. = FALSE)
  if (any(is_empty)) {
    warning(caller, "(): ", sum(is_empty),
            " empty geometry(ies) in 'aoi' dropped.", call. = FALSE)
    sfc <- sfc[!is_empty]
  }

  # CRS required
  crs <- sf::st_crs(sfc)
  if (is.na(crs))
    stop(caller, "(): 'aoi' has no CRS. ",
         "Set one with sf::st_set_crs() before passing to ", caller, "().",
         call. = FALSE)

  # Union multiple features
  if (length(sfc) > 1L)
    sfc <- sf::st_union(sfc)

  # Repair invalid geometries silently
  sfc <- tryCatch(sf::st_make_valid(sfc), error = function(e) sfc)

  geom_type <- as.character(sf::st_geometry_type(sfc, by_geometry = FALSE))

  # Extract polygons from GEOMETRYCOLLECTION
  if (grepl("GEOMETRYCOLLECTION", geom_type, ignore.case = TRUE)) {
    sfc <- tryCatch({
      p <- sf::st_collection_extract(sfc, "POLYGON")
      if (length(p) == 0L) sfc else p
    }, error = function(e) sfc)
    geom_type <- as.character(sf::st_geometry_type(sfc, by_geometry = FALSE))
  }

  # Buffer points / lines
  if (grepl("POINT|LINE", geom_type, ignore.case = TRUE)) {
    if (is.null(buffer) || !is.numeric(buffer) || length(buffer) != 1L ||
        !is.finite(buffer) || buffer <= 0)
      stop(caller, "(): 'buffer' must be a single positive number when ",
           "'aoi' is a ", geom_type, " geometry.", call. = FALSE)
    sfc       <- sf::st_buffer(sfc, dist = buffer)
    geom_type <- as.character(sf::st_geometry_type(sfc, by_geometry = FALSE))
  }

  if (!grepl("POLYGON", geom_type, ignore.case = TRUE))
    stop(caller, "(): 'aoi' resolved to unsupported geometry type '",
         geom_type, "'. Use POLYGON/MULTIPOLYGON, or supply 'buffer' ",
         "for point/line inputs.", call. = FALSE)

  list(type       = "polygon",
       sfc        = sfc,
       bbox       = as.numeric(sf::st_bbox(sfc)),
       needs_clip = TRUE)
}

# ---------------------------------------------------------------------------
# .stac_copc_href() -- find the COPC asset href from a STAC item.
# Different catalogs use different asset keys; we try common ones in order.
# @noRd
.stac_copc_href <- function(item) {
  assets <- item$assets
  if (is.null(assets)) return(NULL)

  # Common asset keys used across known STAC catalogs:
  #   "data"       -- Planetary Computer 3DEP
  #   "copc"       -- MAAP, some Element84 collections
  #   "pointcloud" -- some experimental catalogs
  #   "lidar"      -- NOAA, some NZ catalogs
  #   "asset"      -- generic fallback
  candidate_keys <- c("data", "copc", "pointcloud", "lidar", "asset")
  laz_pat <- "\\.copc\\.laz|\\.laz$"

  for (key in candidate_keys) {
    href <- assets[[key]]$href
    if (!is.null(href) && nzchar(href) &&
        grepl(laz_pat, href, ignore.case = TRUE))
      return(href)
  }

  # Last resort: scan all assets for any .copc.laz / .laz href
  for (key in names(assets)) {
    href <- assets[[key]]$href
    if (!is.null(href) && nzchar(href) &&
        grepl(laz_pat, href, ignore.case = TRUE))
      return(href)
  }

  NULL
}

# ---------------------------------------------------------------------------
# .stac_search_items() -- query a STAC search endpoint and return all items
# (with automatic pagination).
# @noRd
.stac_search_items <- function(stac_url,
                               aoi_sfc,
                               collections = NULL,
                               limit       = 100L,
                               progress    = TRUE) {
  # Guard against static cloud-storage URLs which are not STAC search APIs.
  # A POST to storage.googleapis.com or *.s3*.amazonaws.com returns 400/403/405
  # because they serve objects, not applications.
  static_patterns <- c(
    "storage\\.googleapis\\.com",
    "\\.s3[.-][a-z0-9-]*\\.amazonaws\\.com",
    "^https://s3\\.amazonaws\\.com/"
  )
  if (any(vapply(static_patterns, grepl,
                 logical(1L), x = stac_url, ignore.case = TRUE, perl = TRUE))) {
    stop(
      "The URL looks like a static cloud-storage bucket, which cannot serve ",
      "STAC search queries.\n",
      "  URL: ", stac_url, "\n",
      "  Static STAC catalogs (catalog.json on S3/GCS) do not support the ",
      "STAC Item Search API.\n",
      "  Look for a dedicated STAC search service that indexes the same data.",
      call. = FALSE
    )
  }

  # Guard: CRS must be set before reprojecting
  if (is.na(sf::st_crs(aoi_sfc)))
    stop("AOI has no CRS -- cannot reproject to WGS84 for STAC search.",
         call. = FALSE)

  aoi_4326 <- tryCatch(
    sf::st_transform(aoi_sfc, 4326),
    error = function(e)
      stop("Failed to reproject AOI to WGS84: ", conditionMessage(e),
           call. = FALSE)
  )
  aoi_union       <- sf::st_cast(sf::st_union(aoi_4326), "MULTIPOLYGON")
  intersects_geom <- .sf_to_geojson_geometry(aoi_union)

  body <- list(
    intersects = intersects_geom,
    limit      = as.integer(min(max(limit, 1L), 1000L))
  )
  if (!is.null(collections))
    body$collections <- as.list(collections)

  .stac_post <- function(url, bdy) {
    resp <- tryCatch(
      httr::POST(url, httr::content_type_json(),
                 body = jsonlite::toJSON(bdy, auto_unbox = TRUE, digits = NA)),
      error = function(e)
        stop("STAC HTTP request failed: ", conditionMessage(e),
             "\nCheck that the URL is a valid STAC search endpoint: ", url,
             call. = FALSE)
    )

    status <- httr::status_code(resp)
    if (status >= 400L) {
      # Extract server error body for richer diagnostics
      server_hint <- tryCatch({
        raw_err <- httr::content(resp, as = "text", encoding = "UTF-8")
        if (!nzchar(trimws(raw_err))) return("")
        parsed_err <- tryCatch(
          jsonlite::fromJSON(raw_err, simplifyVector = FALSE),
          error = function(e) NULL
        )
        msg <- parsed_err[["description"]] %||%
               parsed_err[["message"]]     %||%
               parsed_err[["detail"]]      %||% NULL
        if (!is.null(msg) && nzchar(msg))
          paste0(" -- ", msg)
        else if (nchar(raw_err) <= 400L)
          paste0(" -- ", trimws(gsub("\\s+", " ", raw_err)))
        else ""
      }, error = function(e) "")
      stop(httr::http_status(resp)$message,
           ". Failed to search STAC API: ", url, server_hint,
           call. = FALSE)
    }

    raw <- tryCatch(
      httr::content(resp, as = "text", encoding = "UTF-8"),
      error = function(e)
        stop("Could not read STAC response body: ", conditionMessage(e),
             call. = FALSE)
    )
    tryCatch(
      jsonlite::fromJSON(raw, simplifyVector = FALSE),
      error = function(e)
        stop("STAC response is not valid JSON: ", conditionMessage(e),
             "\nEndpoint: ", url, call. = FALSE)
    )
  }

  all_items <- list()
  repeat {
    parsed <- .stac_post(stac_url, body)
    items  <- parsed$features %||% list()
    all_items <- c(all_items, items)
    if (progress)
      message("STAC: fetched ", length(items),
              " item(s) (total: ", length(all_items), ")")

    next_link <- NULL
    for (lnk in parsed$links %||% list())
      if (identical(lnk$rel, "next")) { next_link <- lnk; break }
    if (is.null(next_link)) break

    if (!is.null(next_link$body)) {
      body <- modifyList(body, next_link$body)
    } else if (!is.null(next_link$href)) {
      resp2 <- tryCatch(httr::GET(next_link$href),
                        error = function(e) NULL)
      if (is.null(resp2)) break
      httr::stop_for_status(resp2)
      parsed2 <- tryCatch(
        jsonlite::fromJSON(
          httr::content(resp2, as = "text", encoding = "UTF-8"),
          simplifyVector = FALSE),
        error = function(e) list(features = list())
      )
      items2 <- parsed2$features %||% list()
      all_items <- c(all_items, items2)
      if (progress)
        message("STAC: fetched ", length(items2),
                " item(s) (total: ", length(all_items), ")")
      break
    } else {
      break
    }
  }
  all_items
}

# ---------------------------------------------------------------------------
# .read_copc_from_stac() -- search a STAC endpoint, stream AOI points from
# every matching tile, and return a merged read_copc()-style list
# (list(data = data.table, header = list)).
# @noRd
.read_copc_from_stac <- function(stac_url,
                                 aoi,
                                 buffer      = NULL,
                                 collections = NULL,
                                 zrange      = NULL,
                                 select      = NULL,
                                 filter      = NULL,
                                 max_depth   = -1L,
                                 resolution  = NULL,
                                 max_points  = Inf,
                                 threads     = 4L,
                                 limit       = 100L,
                                 progress    = TRUE) {
  for (pkg in c("sf", "httr", "jsonlite"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for STAC queries.", call. = FALSE)

  # AOI is required for STAC queries
  if (is.null(aoi))
    stop("'aoi' is required when reading from a STAC endpoint.\n",
         "Supply a polygon, point/line (with 'buffer'), or numeric(4) bbox.",
         call. = FALSE)

  # Normalise and validate AOI through the robust helper
  aoi_norm <- .normalise_aoi(aoi, buffer = buffer,
                              caller = ".read_copc_from_stac")

  if (aoi_norm$type == "bbox") {
    # Wrap bare bbox into a WGS84 polygon sfc for the STAC search
    bb <- aoi_norm$bbox
    aoi_sfc <- sf::st_sfc(
      sf::st_polygon(list(matrix(
        c(bb[1], bb[2],  bb[3], bb[2],
          bb[3], bb[4],  bb[1], bb[4],
          bb[1], bb[2]),
        ncol = 2L, byrow = TRUE
      ))),
      crs = 4326
    )
  } else {
    aoi_sfc <- aoi_norm$sfc
  }

  # Auto-default collections for known STAC endpoints
  if (is.null(collections)) {
    if (grepl("planetarycomputer\\.microsoft\\.com", stac_url, ignore.case = TRUE))
      collections <- "3dep-lidar-copc"
  }

  # Search -- wrapped for better error messages on bad/unsupported URLs
  all_items <- tryCatch(
    .stac_search_items(stac_url, aoi_sfc,
                       collections = collections,
                       limit = limit, progress = progress),
    error = function(e)
      stop("STAC search failed.\n",
           "  Endpoint : ", stac_url, "\n",
           "  Reason   : ", conditionMessage(e), "\n",
           "  Hint     : Ensure this is a STAC *search* endpoint ",
           "(e.g. .../stac/v1/search), not a static catalog.json or browser URL.",
           call. = FALSE)
  )

  n_items <- length(all_items)
  if (n_items == 0L) {
    msg <- "No COPC tiles found for the supplied AOI."
    if (grepl("planetarycomputer\\.microsoft\\.com", stac_url, ignore.case = TRUE))
      msg <- paste0(msg, "\n",
        "  Note: 'collection = 3dep-lidar-copc' covers only the contiguous US,\n",
        "        Alaska, and Hawaii. If your AOI is outside the US consider a\n",
        "        different STAC endpoint or dataset.")
    message(msg)
    return(list(data = data.table::data.table(), header = list()))
  }
  if (progress) message("Found ", n_items, " COPC tile(s) intersecting the AOI.")

  # SAS token (Planetary Computer); skip signing for other endpoints
  is_pc <- grepl("planetarycomputer\\.microsoft\\.com", stac_url,
                 ignore.case = TRUE)
  sas_token <- if (is_pc) {
    tryCatch(
      .pc_get_sas_token(),
      error = function(e) {
        warning("Could not obtain Planetary Computer SAS token: ",
                conditionMessage(e), "\nTile URLs will be unsigned.",
                call. = FALSE)
        NULL
      }
    )
  } else NULL

  all_data   <- list()
  last_hdr   <- list()
  n_total    <- 0

  for (i in seq_len(n_items)) {
    item    <- all_items[[i]]
    item_id <- item$id %||% paste0("item_", i)

    href <- .stac_copc_href(item)
    if (is.null(href)) {
      warning("Item '", item_id, "' has no recognisable COPC/LAZ asset -- skipping.",
              call. = FALSE)
      next
    }

    signed_url <- if (!is.null(sas_token)) {
      sep <- if (grepl("\\?", href)) "&" else "?"
      paste0(href, sep, sas_token)
    } else {
      href
    }

    if (progress)
      message(sprintf("[%d/%d] Fetching AOI from %s ...", i, n_items, item_id))

    # Pass aoi_sfc -- read_copc() auto-reprojects into the tile CRS
    res <- tryCatch(
      read_copc(
        path_or_url = signed_url,
        aoi         = sf::st_sf(geometry = aoi_sfc),
        zrange      = zrange,
        select      = select,
        filter      = filter,
        max_depth   = max_depth,
        resolution  = resolution,
        max_points  = Inf,
        threads     = threads,
        progress    = FALSE   # per-tile progress is noisy; we print per-tile above
      ),
      error = function(e) {
        warning("read_copc() failed for '", item_id, "': ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )

    if (is.null(res) || nrow(res$data) == 0L) {
      if (progress) message("  0 points in AOI -- skipping.")
      next
    }

    n_pts <- nrow(res$data)
    if (progress)
      message(sprintf("  %s point(s) read.", format(n_pts, big.mark = ",")))

    all_data[[length(all_data) + 1L]] <- res$data
    last_hdr <- res$header
    n_total  <- n_total + n_pts

    if (is.finite(max_points) && n_total >= max_points) {
      if (progress) message("max_points reached -- stopping early.")
      break
    }
  }

  if (length(all_data) == 0L) {
    message("No points found within the AOI across any tile.")
    return(list(data = data.table::data.table(), header = list()))
  }

  merged <- data.table::rbindlist(all_data, fill = TRUE)
  if (is.finite(max_points) && nrow(merged) > max_points)
    merged <- merged[seq_len(as.integer(max_points)), ]

  list(data = merged, header = last_hdr)
}
