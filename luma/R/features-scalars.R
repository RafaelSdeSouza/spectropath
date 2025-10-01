#' @include algebra-helpers.R chen.R streamsig.R logsig.R
NULL
# =========================
# Public, one-scalar APIs
# =========================

#' Scalar log-signature features
#'
#' Standalone estimators for drift, Levy area, curls, jerks, and twist.
#' All functions take a 2-column numeric matrix `path` with columns (t, f).
#'
#' @param path Numeric matrix with columns (t, f).
#' @param normalize Logical; if TRUE, return the scale-invariant version.
#' @return Numeric scalar.
#' @examples
#' path <- cbind(t = c(0, 1, 2), f = c(0, 0.5, 0.2))
#' drift_f(path)
#' levy_tf(path)
#' curl_t(path)
#' curl_f(path)
#' jerk_t(path)
#' twist_tf(path)
#' jerk_f(path)
#' @name scalar_features
NULL

# ---- Drift --------------------------------------------------------------

#' @rdname scalar_features
#' @title Drift in f (endpoint displacement)
#' @export
drift_f <- function(path, normalize = FALSE) {
  validate_path_2d(path)
  df <- path[nrow(path), 2] - path[1, 2]
  if (normalize) {
    if (abs(df) < .Machine$double.eps^0.5) return(0)
    return(df / df)
  }
  df
}

# ---- Level 2: Levy area -------------------------------------------------

#' @rdname scalar_features
#' @title Levy area [t,f]
#' @export
levy_tf <- function(path, normalize = TRUE) {
  validate_path_2d(path)
  sig2 <- .stream2sig_d2_raw(path)
  L2 <- sig2$L2 - 0.5 * .kron2(sig2$L1, sig2$L1)
  val <- 0.5 * (L2[1, 2] - L2[2, 1])  # (tf - ft)/2
  if (!normalize) return(val)
  dt <- max(path[, 1]) - min(path[, 1])
  df <- max(path[, 2]) - min(path[, 2])
  .norm_scales(val, dt, df, pow_t = 1, pow_f = 1)
}

# ---- Level 3: curls -----------------------------------------------------

#' @rdname scalar_features
#' @title Curl in t: [t,[t,f]]
#' @export
curl_t <- function(path, normalize = TRUE) {
  validate_path_2d(path)
  sig3 <- .stream2sig_d3_raw(path)
  Lt <- .logsig_tensor_d4(list(L1 = sig3$L1, L2 = sig3$L2, L3 = sig3$L3, L4 = numeric(16)))
  val <- 0.5 * (Lt$L3[2] - Lt$L3[3])  # (ttf - tft)/2
  if (!normalize) return(val)
  dt <- max(path[, 1]) - min(path[, 1])
  df <- max(path[, 2]) - min(path[, 2])
  .norm_scales(val, dt, df, pow_t = 2, pow_f = 1)
}

#' @rdname scalar_features
#' @title Curl in f: [f,[t,f]]
#' @export
curl_f <- function(path, normalize = TRUE) {
  validate_path_2d(path)
  sig3 <- .stream2sig_d3_raw(path)
  Lt <- .logsig_tensor_d4(list(L1 = sig3$L1, L2 = sig3$L2, L3 = sig3$L3, L4 = numeric(16)))
  val <- 0.5 * (Lt$L3[6] - Lt$L3[7])  # (ftf - fft)/2
  if (!normalize) return(val)
  dt <- max(path[, 1]) - min(path[, 1])
  df <- max(path[, 2]) - min(path[, 2])
  .norm_scales(val, dt, df, pow_t = 1, pow_f = 2)
}

# ---- Level 4: jerks and twist ------------------------------------------

#' @rdname scalar_features
#' @title Jerk in t: [t,[t,[t,f]]]
#' @export
jerk_t <- function(path, normalize = TRUE) {
  validate_path_2d(path)
  sig4 <- .stream2sig_d4_raw(path)
  Lt <- .logsig_tensor_d4(sig4)
  # indices in degree-4 layout: 2,3,5,9 => tttf - ttft + tftt - fttt
  val <- 0.25 * (Lt$L4[2] - Lt$L4[3] + Lt$L4[5] - Lt$L4[9])
  if (!normalize) return(val)
  dt <- max(path[, 1]) - min(path[, 1])
  df <- max(path[, 2]) - min(path[, 2])
  .norm_scales(val, dt, df, pow_t = 3, pow_f = 1)
}

#' @rdname scalar_features
#' @title Twist tf: [t,[f,[t,f]]]
#' @export
twist_tf <- function(path, normalize = TRUE) {
  validate_path_2d(path)
  sig4 <- .stream2sig_d4_raw(path)
  Lt <- .logsig_tensor_d4(sig4)
  # indices: 6,7,10,13 => tftf - tfft + fttf - fftt
  val <- 0.25 * (Lt$L4[6] - Lt$L4[7] + Lt$L4[10] - Lt$L4[13])
  if (!normalize) return(val)
  dt <- max(path[, 1]) - min(path[, 1])
  df <- max(path[, 2]) - min(path[, 2])
  .norm_scales(val, dt, df, pow_t = 2, pow_f = 2)
}

#' @rdname scalar_features
#' @title Jerk in f: [f,[f,[t,f]]]
#' @export
jerk_f <- function(path, normalize = TRUE) {
  validate_path_2d(path)
  sig4 <- .stream2sig_d4_raw(path)
  Lt <- .logsig_tensor_d4(sig4)
  # indices: 15,14,11,8 => ffft - fftf + ftft - tfff
  val <- 0.25 * (Lt$L4[15] - Lt$L4[14] + Lt$L4[11] - Lt$L4[8])
  if (!normalize) return(val)
  dt <- max(path[, 1]) - min(path[, 1])
  df <- max(path[, 2]) - min(path[, 2])
  .norm_scales(val, dt, df, pow_t = 1, pow_f = 3)
}
