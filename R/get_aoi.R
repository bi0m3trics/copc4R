#' Interactively draw an Area of Interest on a leaflet map
#'
#' Opens an interactive leaflet map with satellite imagery and lets the
#' user draw a rectangle AOI.
#' The map attempts browser geolocation on load; if denied it falls back
#' to the supplied centre and zoom.
#'
#' @param fallback_center Numeric length-2 vector `c(lng, lat)` used when
#'   geolocation is unavailable.
#' @param fallback_zoom Integer. Initial zoom level for the fallback view.
#'
#' @return An `sf` data frame containing a single polygon feature, or
#'   `NULL` if no AOI was drawn.
#'
#' @details
#' Requires the \pkg{sf}, \pkg{leaflet}, \pkg{mapedit}, and
#' \pkg{htmlwidgets} packages.
#' Only the first drawn feature is retained.
#'
#' @examples
#' \donttest{
#' aoi <- get_aoi()
#' tiles <- download_3dep_copc(aoi, dest_dir = "tiles/")
#' }
#'
#' @export
get_aoi <- function(
    fallback_center = c(-111.6513, 35.1983),
    fallback_zoom   = 12
) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required.", call. = FALSE)
  if (!requireNamespace("leaflet", quietly = TRUE))
    stop("Package 'leaflet' is required.", call. = FALSE)
  if (!requireNamespace("mapedit", quietly = TRUE))
    stop("Package 'mapedit' is required.", call. = FALSE)
  if (!requireNamespace("htmlwidgets", quietly = TRUE))
    stop("Package 'htmlwidgets' is required.", call. = FALSE)

  # build leaflet map
  m <- leaflet::leaflet() |>
    leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery) |>
    leaflet::setView(
      lng  = fallback_center[1],
      lat  = fallback_center[2],
      zoom = fallback_zoom
    )

  # use browser geolocation on load, if allowed
  m <- htmlwidgets::onRender(
    m,
    "
    function(el, x) {
      var map = this;
      function onLocationFound(e) { map.setView(e.latlng, 15); }
      function onLocationError(e) {
        console.log('Location not available, using fallback view.');
      }
      map.on('locationfound', onLocationFound);
      map.on('locationerror', onLocationError);
      map.locate({ setView: true, maxZoom: 15, enableHighAccuracy: true });
    }
    "
  )

  # launch drawing tool (rectangle only)
  aoi <- mapedit::drawFeatures(
    map   = m,
    sf    = TRUE,
    title = "Draw AOI box",
    editorOptions = list(
      drawRectangle = TRUE,
      drawPolygon   = FALSE,
      drawMarker    = FALSE,
      drawCircle    = FALSE,
      drawPolyline  = FALSE,
      cutPolygon    = FALSE,
      removalMode   = FALSE,
      editMode      = FALSE
    )
  )

  if (is.null(aoi) || nrow(aoi) == 0) {
    message("No AOI drawn.")
    return(NULL)
  }

  # keep only the first feature
  aoi <- aoi[1, , drop = FALSE]
  aoi
}
