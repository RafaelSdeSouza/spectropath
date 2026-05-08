#' Path-signature diagnostics for spectral-line morphology
#'
#' These functions expose the low-order path-signature quantities used in the
#' accompanying paper. They are computed from the log-signature of the
#' piecewise-linear path `X = (u, F)`, where `u` is an ordered coordinate and
#' `F` is the measured signal or line flux.
#'
#' The names follow the paper notation:
#'
#' - `levy_area()` returns \eqn{A_{uF}}, the signed area.
#' - `curl_u()` returns \eqn{C_u}, an ordered-coordinate weighted area.
#' - `curl_f()` returns \eqn{C_F}, a flux-weighted area.
#' - `skew_signature()` returns \eqn{S_{\rm skew}}, an order-sensitive
#'   asymmetry contrast.
#' - `jerk_u()`, `twist_uf()`, and `jerk_f()` return fourth-order contrasts
#'   used to describe higher-order profile modulation.
#' - `emabs_area()` returns \eqn{A_{+-}}, the signed area of the cumulative
#'   emission--absorption path.
#'
#' @param path Two-column numeric matrix or data frame coercible by
#'   [as_spectral_path()].
#' @param normalize Logical; if `TRUE`, divide each contrast by its natural
#'   coordinate and flux ranges.
#' @return A numeric scalar.
#' @examples
#' u <- seq(-5, 5, length.out = 200)
#' f <- exp(-0.5 * (u / 0.8)^2) + 0.2 * exp(-0.5 * ((u - 2) / 1.2)^2)
#' path <- cbind(u, f)
#' levy_area(path)
#' curl_f(path)
#' @name path_signature_scalars
NULL

.logsig_for_depth <- function(path, depth) {
  if (depth == 2) {
    sig <- .stream2sig_d2_raw(path)
    list(L1 = sig$L1, L2 = sig$L2 - 0.5 * .kron2(sig$L1, sig$L1))
  } else if (depth == 3) {
    sig <- .stream2sig_d3_raw(path)
    .logsig_tensor_d4(list(L1 = sig$L1, L2 = sig$L2, L3 = sig$L3, L4 = numeric(16)))
  } else if (depth == 4) {
    .logsig_tensor_d4(.stream2sig_d4_raw(path))
  } else {
    stop("`depth` must be 2, 3, or 4.", call. = FALSE)
  }
}

.maybe_norm <- function(path, value, pow_u, pow_f, normalize) {
  if (!normalize) return(.scalar(value))
  ranges <- .path_ranges(path)
  .norm_scales(value, ranges$du, ranges$df, pow_u, pow_f)
}

#' @rdname path_signature_scalars
#' @export
levy_area <- function(path, normalize = TRUE) {
  path <- as_spectral_path(path)
  L <- .logsig_for_depth(path, 2)
  .maybe_norm(path, 0.5 * (L$L2[1, 2] - L$L2[2, 1]), 1, 1, normalize)
}

#' @rdname path_signature_scalars
#' @export
curl_u <- function(path, normalize = TRUE) {
  path <- as_spectral_path(path)
  L <- .logsig_for_depth(path, 3)
  .maybe_norm(path, 0.5 * (L$L3[2] - L$L3[3]), 2, 1, normalize)
}

#' @rdname path_signature_scalars
#' @export
curl_f <- function(path, normalize = TRUE) {
  path <- as_spectral_path(path)
  L <- .logsig_for_depth(path, 3)
  .maybe_norm(path, 0.5 * (L$L3[6] - L$L3[7]), 1, 2, normalize)
}

#' @rdname path_signature_scalars
#' @export
skew_signature <- function(path, normalize = TRUE) {
  path <- as_spectral_path(path)
  L <- .logsig_for_depth(path, 3)
  value <- 0.5 * (L$L3[2] + L$L3[3] - L$L3[6] - L$L3[7])
  .maybe_norm(path, value, 2, 1, normalize)
}

#' @rdname path_signature_scalars
#' @export
jerk_u <- function(path, normalize = TRUE) {
  path <- as_spectral_path(path)
  L <- .logsig_for_depth(path, 4)
  value <- 0.25 * (L$L4[2] - L$L4[3] + L$L4[5] - L$L4[9])
  .maybe_norm(path, value, 3, 1, normalize)
}

#' @rdname path_signature_scalars
#' @export
twist_uf <- function(path, normalize = TRUE) {
  path <- as_spectral_path(path)
  L <- .logsig_for_depth(path, 4)
  value <- 0.25 * (L$L4[6] - L$L4[7] + L$L4[10] - L$L4[13])
  .maybe_norm(path, value, 2, 2, normalize)
}

#' @rdname path_signature_scalars
#' @export
jerk_f <- function(path, normalize = TRUE) {
  path <- as_spectral_path(path)
  L <- .logsig_for_depth(path, 4)
  value <- 0.25 * (L$L4[15] - L$L4[14] + L$L4[11] - L$L4[8])
  .maybe_norm(path, value, 1, 3, normalize)
}

#' Cumulative emission--absorption path
#'
#' Builds the cumulative path \eqn{(C_+, C_-)}, where positive and negative
#' flux are accumulated separately along the ordered coordinate. This is useful
#' for P-Cygni-like or mixed emission--absorption morphology.
#'
#' @param path Two-column numeric matrix or data frame coercible by
#'   [as_spectral_path()].
#' @param orient If `"increasing"`, force the path to run from smaller to
#'   larger coordinate. If `"as_is"`, preserve input order.
#' @return A two-column matrix with columns `C_plus` and `C_minus`.
#' @export
emabs_path <- function(path, orient = c("increasing", "as_is")) {
  orient <- match.arg(orient)
  path <- as_spectral_path(path, sort = orient == "increasing")

  u <- path[, 1]
  flux <- path[, 2]

  if (orient == "increasing" && any(diff(u) <= 0)) {
    stop("For orient = 'increasing', the ordered coordinate must be strictly increasing.", call. = FALSE)
  }
  if (orient == "as_is" && any(diff(u) == 0)) {
    stop("The ordered coordinate contains repeated values.", call. = FALSE)
  }

  out <- cbind(
    C_plus = .cumtrapz(u, pmax(flux, 0)),
    C_minus = .cumtrapz(u, pmax(-flux, 0))
  )
  rownames(out) <- NULL
  out
}

#' @rdname path_signature_scalars
#' @param orient If `"increasing"`, force the path to run from smaller to
#'   larger coordinate before computing the cumulative emission--absorption
#'   area. If `"as_is"`, preserve input order.
#' @param eps Numerical tolerance used to detect absent emission or absorption.
#' @export
emabs_area <- function(path, normalize = TRUE,
                       orient = c("increasing", "as_is"),
                       eps = sqrt(.Machine$double.eps)) {
  orient <- match.arg(orient)
  Xpm <- emabs_path(path, orient = orient)
  d_plus <- diff(range(Xpm[, 1]))
  d_minus <- diff(range(Xpm[, 2]))

  if (!is.finite(d_plus) || !is.finite(d_minus) ||
      d_plus <= eps || d_minus <= eps) {
    return(0)
  }

  levy_area(Xpm, normalize = normalize)
}

#' Compute the paper-facing path-signature feature set
#'
#' `path_features()` returns a one-row data frame containing the low-order
#' path-signature diagnostics discussed in the paper. At depth 2 it returns the
#' signed area and emission--absorption area. At depth 3 it adds the weighted
#' area and skew-signature contrasts. At depth 4 it adds the higher-order
#' modulation terms.
#'
#' @param path Two-column numeric matrix or data frame coercible by
#'   [as_spectral_path()].
#' @param depth Integer depth: 2, 3, or 4.
#' @param normalize Logical; if `TRUE`, divide each contrast by its natural
#'   coordinate and flux ranges.
#' @return A one-row data frame.
#' @examples
#' u <- seq(-5, 5, length.out = 200)
#' f <- exp(-0.5 * (u / 0.8)^2) + 0.2 * exp(-0.5 * ((u - 2) / 1.2)^2)
#' path_features(cbind(u, f), depth = 4)
#' @export
path_features <- function(path, depth = 4, normalize = TRUE) {
  path <- as_spectral_path(path)
  if (!depth %in% 2:4) stop("`depth` must be 2, 3, or 4.", call. = FALSE)

  out <- list(
    levy_area = levy_area(path, normalize = normalize),
    emabs_area = emabs_area(path, normalize = normalize)
  )

  if (depth >= 3) {
    out <- c(out, list(
      curl_u = curl_u(path, normalize = normalize),
      curl_f = curl_f(path, normalize = normalize),
      skew_signature = skew_signature(path, normalize = normalize)
    ))
  }

  if (depth >= 4) {
    out <- c(out, list(
      jerk_u = jerk_u(path, normalize = normalize),
      twist_uf = twist_uf(path, normalize = normalize),
      jerk_f = jerk_f(path, normalize = normalize)
    ))
  }

  data.frame(out, check.names = FALSE, row.names = NULL)
}
