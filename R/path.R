#' Coerce spectral data to an ordered two-dimensional path
#'
#' `as_spectral_path()` converts a two-column matrix or data frame into the
#' path convention used by `spectropath`: an ordered coordinate `u` and a
#' measured signal or line flux `F`.
#'
#' The coordinate `u` may be wavelength, rest-frame wavelength, log-wavelength,
#' phase, time, or velocity relative to a transition. For spectral-line
#' morphology, `F` is usually continuum-subtracted or continuum-normalized
#' flux, but the path-signature calculations only require an ordered
#' two-dimensional curve.
#'
#' @param x A two-column numeric matrix, or a data frame containing the columns
#'   named by `u` and `flux`.
#' @param u,flux Column names to read when `x` is a data frame.
#' @param sort Logical; if `TRUE`, sort rows by the ordered coordinate.
#' @param remove_nonfinite Logical; if `TRUE`, remove rows containing `NA`,
#'   `NaN`, or infinite values.
#'
#' @return A two-column numeric matrix with columns `u` and `flux`.
#' @examples
#' d <- data.frame(lambda = c(5002, 5001, 5003), flux = c(0.1, 0.0, 0.2))
#' as_spectral_path(d, u = "lambda", flux = "flux")
#' @export
as_spectral_path <- function(x, u = NULL, flux = NULL,
                             sort = TRUE, remove_nonfinite = TRUE) {
  if (is.matrix(x) && is.numeric(x) && ncol(x) == 2) {
    path <- x
  } else if (is.data.frame(x)) {
    if (is.null(u) || is.null(flux)) {
      stop("For data frames, provide `u` and `flux` column names.", call. = FALSE)
    }
    path <- cbind(u = as.numeric(x[[u]]), flux = as.numeric(x[[flux]]))
  } else {
    stop("`x` must be a two-column numeric matrix or a data frame.", call. = FALSE)
  }

  if (remove_nonfinite) {
    keep <- is.finite(path[, 1]) & is.finite(path[, 2])
    path <- path[keep, , drop = FALSE]
  }

  if (nrow(path) < 2) {
    stop("A spectral path requires at least two finite points.", call. = FALSE)
  }

  if (sort) {
    path <- path[order(path[, 1]), , drop = FALSE]
  }

  colnames(path) <- c("u", "flux")
  validate_path_2d(path)
  path
}

validate_path_2d <- function(path) {
  if (!is.matrix(path) || !is.numeric(path)) {
    stop("`path` must be a numeric matrix.", call. = FALSE)
  }
  if (ncol(path) != 2) {
    stop("`path` must have exactly two columns: ordered coordinate and flux.", call. = FALSE)
  }
  if (nrow(path) < 2) {
    stop("`path` must have at least two rows.", call. = FALSE)
  }
  if (!all(is.finite(path))) {
    stop("`path` must contain only finite values.", call. = FALSE)
  }
  invisible(TRUE)
}

.path_ranges <- function(path) {
  list(
    du = max(path[, 1]) - min(path[, 1]),
    df = max(path[, 2]) - min(path[, 2])
  )
}

.trapezoid_integral <- function(x, y) {
  if (length(x) < 2) return(0)
  sum(0.5 * (y[-length(y)] + y[-1]) * diff(x))
}

.cumtrapz <- function(x, y) {
  if (length(x) < 2) return(rep(0, length(x)))
  c(0, cumsum(0.5 * (y[-length(y)] + y[-1]) * diff(x)))
}
