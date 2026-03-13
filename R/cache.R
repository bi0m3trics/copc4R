# -- Range-read cache ---------------------------------------------------------
#
# Persistent on-disk cache (and in-memory LRU) keyed by URL + byte-range so
# repeated AOI queries don't re-download the same COPC chunks.                 #
# The on-disk cache lives in `tools::R_user_dir("copc4R", "cache")` by default
# and can be overridden via `options(copc4R.cache_dir = "/some/path")`.
#
# Internal helpers are prefixed with `.cache_` and are not exported.
# ---------------------------------------------------------------------------

# -- In-memory LRU cache (environment-based) ---------------------------------
.mem_cache <- new.env(parent = emptyenv())
.mem_cache$store   <- list()
.mem_cache$keys    <- character(0)
.mem_cache$max_mb  <- 256L            # default 256 MB

# -- Header cache (lightweight memoization) ----------------------------------
.header_cache <- new.env(parent = emptyenv())
.header_cache$store       <- list()
.header_cache$max_entries <- 50L

# -- Public: configure cache -------------------------------------------------

#' Configure the COPC range-read cache
#'
#' Controls the in-memory LRU size and on-disk cache directory.
#'
#' @param mem_max_mb Numeric scalar. Maximum in-memory cache size in MB.
#'   Default 256.
#' @param disk_dir Character.  Path to the on-disk cache directory.
#'   \code{NULL} (default) uses `tools::R_user_dir("copc4R", "cache")`.
#'   Set to \code{FALSE} to disable on-disk caching.
#' @param enabled Logical.  Set to \code{FALSE} to disable all caching.
#'
#' @return Invisibly returns a list with the previous settings.
#'
#' @examples
#' # Increase in-memory cache to 512 MB
#' prev <- copc_cache_config(mem_max_mb = 512)
#'
#' # Restore previous settings
#' copc_cache_config(mem_max_mb = prev$mem_max_mb)
#'
#' @export
copc_cache_config <- function(mem_max_mb = NULL,
                              disk_dir   = NULL,
                              enabled    = NULL) {
  prev <- list(
    mem_max_mb = .mem_cache$max_mb,
    disk_dir   = getOption("copc4R.cache_dir"),
    enabled    = getOption("copc4R.cache_enabled", TRUE)
  )

  if (!is.null(mem_max_mb)) {
    stopifnot(is.numeric(mem_max_mb), length(mem_max_mb) == 1L)
    .mem_cache$max_mb <- as.integer(mem_max_mb)
  }
  if (!is.null(disk_dir)) {
    if (isFALSE(disk_dir)) {
      options(copc4R.cache_dir = FALSE)
    } else {
      stopifnot(is.character(disk_dir), length(disk_dir) == 1L)
      options(copc4R.cache_dir = disk_dir)
    }
  }
  if (!is.null(enabled)) {
    stopifnot(is.logical(enabled), length(enabled) == 1L)
    options(copc4R.cache_enabled = enabled)
  }
  invisible(prev)
}


#' Report cache statistics
#'
#' @return A named list with \code{mem_entries}, \code{mem_mb},
#'   \code{disk_entries}, \code{disk_mb}, \code{disk_dir}.
#'
#' @examples
#' copc_cache_stats()
#'
#' @export
copc_cache_stats <- function() {
  mem_entries <- length(.mem_cache$keys)
  mem_bytes   <- sum(vapply(.mem_cache$store,
                            function(x) as.double(object.size(x)), 0.0))
  ddir <- .cache_disk_dir()
  disk_entries <- 0L
  disk_bytes   <- 0
  if (!is.null(ddir) && dir.exists(ddir)) {
    files <- list.files(ddir, full.names = TRUE, recursive = TRUE)
    disk_entries <- length(files)
    disk_bytes   <- sum(file.size(files), na.rm = TRUE)
  }

  list(
    mem_entries     = mem_entries,
    mem_mb          = round(mem_bytes / 1e6, 2),
    header_entries  = length(.header_cache$store),
    disk_entries    = disk_entries,
    disk_mb         = round(disk_bytes / 1e6, 2),
    disk_dir        = ddir %||% "(disabled)"
  )
}


#' Clear the COPC cache
#'
#' @param mem  Logical.  Clear the in-memory LRU cache? Default \code{TRUE}.
#' @param disk Logical.  Clear the on-disk persistent cache? Default \code{FALSE}.
#'
#' @return Invisibly returns the cache stats \emph{before} clearing.
#'
#' @examples
#' copc_cache_clear(mem = TRUE, disk = FALSE)
#'
#' @export
copc_cache_clear <- function(mem = TRUE, disk = FALSE) {
  prev <- copc_cache_stats()
  if (mem) {
    .mem_cache$store <- list()
    .mem_cache$keys  <- character(0)
    .header_cache$store <- list()
  }
  if (disk) {
    ddir <- .cache_disk_dir()
    if (!is.null(ddir) && dir.exists(ddir)) {
      unlink(ddir, recursive = TRUE)
    }
  }
  invisible(prev)
}


# -- Internal helpers ---------------------------------------------------------

# Null-coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

.cache_enabled <- function() {
  isTRUE(getOption("copc4R.cache_enabled", TRUE))
}

.cache_disk_dir <- function() {
  dir <- getOption("copc4R.cache_dir")
  if (isFALSE(dir)) return(NULL)
  if (is.null(dir)) {
    dir <- tools::R_user_dir("copc4R", "cache")
  }
  dir
}

.cache_key <- function(url, offset, length) {
  digest_input <- paste0(url, "@", offset, "+", length)
  # Simple hash using rlang-free method (R's built-in)
  sprintf("%s_%s_%s",
          substr(gsub("[^a-zA-Z0-9]", "", basename(sub("[?#].*$", "", url))), 1, 20),
          as.character(offset),
          as.character(length))
}

.cache_get <- function(url, offset, length) {
  if (!.cache_enabled()) return(NULL)

  key <- .cache_key(url, offset, length)

  # 1. Try in-memory LRU
  if (key %in% .mem_cache$keys) {
    # Move to front (most recently used)
    .mem_cache$keys <- c(key, setdiff(.mem_cache$keys, key))
    return(.mem_cache$store[[key]])
  }

  # 2. Try on-disk
  ddir <- .cache_disk_dir()
  if (!is.null(ddir)) {
    fpath <- file.path(ddir, paste0(key, ".rds"))
    if (file.exists(fpath)) {
      val <- readRDS(fpath)
      .cache_mem_put(key, val)
      return(val)
    }
  }

  NULL
}

.cache_put <- function(url, offset, length, data) {
  if (!.cache_enabled()) return(invisible(NULL))

  key <- .cache_key(url, offset, length)

  # In-memory LRU
  .cache_mem_put(key, data)

  # On-disk
  ddir <- .cache_disk_dir()
  if (!is.null(ddir)) {
    if (!dir.exists(ddir)) dir.create(ddir, recursive = TRUE)
    fpath <- file.path(ddir, paste0(key, ".rds"))
    tryCatch(saveRDS(data, fpath), error = function(e) NULL)
  }

  invisible(NULL)
}

.cache_mem_put <- function(key, data) {
  .mem_cache$store[[key]] <- data
  .mem_cache$keys <- c(key, setdiff(.mem_cache$keys, key))

  # Evict LRU entries if over limit
  total_bytes <- sum(vapply(.mem_cache$store,
                            function(x) as.double(object.size(x)), 0.0))
  max_bytes <- .mem_cache$max_mb * 1e6
  while (total_bytes > max_bytes && length(.mem_cache$keys) > 1L) {
    evict_key <- .mem_cache$keys[length(.mem_cache$keys)]
    total_bytes <- total_bytes -
      as.double(object.size(.mem_cache$store[[evict_key]]))
    .mem_cache$store[[evict_key]] <- NULL
    .mem_cache$keys <- .mem_cache$keys[-length(.mem_cache$keys)]
  }
}


# -- Header cache helpers -----------------------------------------------------

.header_cache_get <- function(path_or_url) {
  .header_cache$store[[path_or_url]]
}

.header_cache_put <- function(path_or_url, header) {
  .header_cache$store[[path_or_url]] <- header
  # Evict oldest if over limit
  if (length(.header_cache$store) > .header_cache$max_entries) {
    .header_cache$store[[names(.header_cache$store)[1L]]] <- NULL
  }
}
