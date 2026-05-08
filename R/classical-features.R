#' Classical line-profile descriptors
#'
#' Convenience functions for benchmarking and machine-learning workflows.
#' These quantities are not path-signature features; they are provided so users
#' can compare classical summaries against the path-signature diagnostics in a
#' consistent table.
#'
#' `equivalent_width()` returns a signed pseudo-equivalent width for a
#' continuum-subtracted or continuum-normalized profile. If `continuum = 1`,
#' this is simply the integrated line flux in the coordinate units. For
#' conventional equivalent widths, provide the continuum level in the same units
#' as the flux.
#'
#' Moment-like quantities use non-negative profile weights. For emission lines,
#' the default `weight = "positive"` uses `pmax(F, 0)`. For mixed
#' emission--absorption profiles, `weight = "absolute"` can be useful.
#'
#' @param path Two-column numeric matrix or data frame coercible by
#'   [as_spectral_path()].
#' @param continuum Scalar or vector continuum level used by
#'   `equivalent_width()`.
#' @param weight Weighting convention for moments and non-parametric widths:
#'   `"positive"` uses positive flux only, `"absolute"` uses absolute flux,
#'   and `"negative"` uses absorption depth `pmax(-F, 0)`.
#' @return Numeric scalar, except `line_moments()` and `classical_features()`,
#'   which return data frames.
#' @name classical_line_features
NULL

.profile_weight <- function(path, weight = c("positive", "absolute", "negative")) {
  weight <- match.arg(weight)
  flux <- path[, 2]
  w <- switch(
    weight,
    positive = pmax(flux, 0),
    absolute = abs(flux),
    negative = pmax(-flux, 0)
  )
  area <- .trapezoid_integral(path[, 1], w)
  if (!is.finite(area) || area <= 0) {
    return(rep(NA_real_, length(w)))
  }
  w / area
}

#' @rdname classical_line_features
#' @export
line_flux <- function(path) {
  path <- as_spectral_path(path)
  .scalar(.trapezoid_integral(path[, 1], path[, 2]))
}

#' @rdname classical_line_features
#' @export
peak_flux <- function(path) {
  path <- as_spectral_path(path)
  .scalar(max(path[, 2], na.rm = TRUE))
}

#' @rdname classical_line_features
#' @export
equivalent_width <- function(path, continuum = 1) {
  path <- as_spectral_path(path)
  if (length(continuum) == 1) {
    continuum <- rep(continuum, nrow(path))
  }
  if (length(continuum) != nrow(path)) {
    stop("`continuum` must be a scalar or have one value per path row.", call. = FALSE)
  }
  if (any(!is.finite(continuum)) || any(continuum == 0)) {
    stop("`continuum` must contain finite non-zero values.", call. = FALSE)
  }
  .scalar(.trapezoid_integral(path[, 1], path[, 2] / continuum))
}

#' @rdname classical_line_features
#' @export
line_moments <- function(path, weight = c("positive", "absolute", "negative")) {
  path <- as_spectral_path(path)
  u <- path[, 1]
  w <- .profile_weight(path, weight)

  if (all(is.na(w))) {
    return(data.frame(
      centroid = NA_real_,
      line_dispersion = NA_real_,
      skewness = NA_real_,
      kurtosis_excess = NA_real_
    ))
  }

  centroid <- .trapezoid_integral(u, u * w)
  variance <- .trapezoid_integral(u, (u - centroid)^2 * w)
  sigma <- sqrt(max(variance, 0))

  if (!is.finite(sigma) || sigma <= 0) {
    skew <- NA_real_
    kurt <- NA_real_
  } else {
    z <- (u - centroid) / sigma
    skew <- .trapezoid_integral(u, z^3 * w)
    kurt <- .trapezoid_integral(u, z^4 * w) - 3
  }

  data.frame(
    centroid = .scalar(centroid),
    line_dispersion = .scalar(sigma),
    skewness = .scalar(skew),
    kurtosis_excess = .scalar(kurt),
    row.names = NULL
  )
}

#' @rdname classical_line_features
#' @export
fwhm <- function(path, weight = c("positive", "absolute", "negative")) {
  path <- as_spectral_path(path)
  weight <- match.arg(weight)
  u <- path[, 1]
  flux <- switch(
    weight,
    positive = pmax(path[, 2], 0),
    absolute = abs(path[, 2]),
    negative = pmax(-path[, 2], 0)
  )

  mx <- max(flux, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) return(NA_real_)

  idx <- which(flux >= 0.5 * mx)
  if (length(idx) < 2) return(NA_real_)

  .scalar(max(u[idx]) - min(u[idx]))
}

#' @rdname classical_line_features
#' @export
w80 <- function(path, weight = c("positive", "absolute", "negative")) {
  path <- as_spectral_path(path)
  u <- path[, 1]
  w <- .profile_weight(path, weight)
  if (all(is.na(w))) return(NA_real_)

  cdf <- .cumtrapz(u, w)
  total <- max(cdf, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  cdf <- cdf / total

  q10 <- stats::approx(cdf, u, xout = 0.10, ties = "ordered")$y
  q90 <- stats::approx(cdf, u, xout = 0.90, ties = "ordered")$y

  .scalar(q90 - q10)
}

#' @rdname classical_line_features
#' @export
classical_features <- function(path, continuum = 1,
                               weight = c("positive", "absolute", "negative")) {
  path <- as_spectral_path(path)
  weight <- match.arg(weight)
  moments <- line_moments(path, weight = weight)
  data.frame(
    line_flux = line_flux(path),
    equivalent_width = equivalent_width(path, continuum = continuum),
    peak_flux = peak_flux(path),
    fwhm = fwhm(path, weight = weight),
    w80 = w80(path, weight = weight),
    moments,
    check.names = FALSE,
    row.names = NULL
  )
}
