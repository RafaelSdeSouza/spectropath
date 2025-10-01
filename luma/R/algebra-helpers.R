#' Small tensor helpers (internal)
#'
#' Consistent Kronecker-based layouts used across the package.
#' - .kron2(a,b): 2x2 with [1,1]=tt, [1,2]=tf, [2,1]=ft, [2,2]=ff
#' - .kron3(a,b,c): 8-vector (ttt, ttf, tft, tff, ftt, ftf, fft, fff)
#' - .kron4(a,b,c,d): 16-vector (tttt, tttf, ttft, ttff, tftt, tftf, tfft, tfff, fttt, fttf, ftft, ftff, fftt, fftf, ffft, ffff)
#'
#' @keywords internal
#' @noRd
# ---------- Small tensor helpers (consistent layouts) ----------
.kron2 <- function(a, b) {
  outer(a, b)
}
#' @keywords internal
#' @noRd
.kron3 <- function(a, b, c) {
  # 8-vector layout: (ttt, ttf, tft, tff, ftt, ftf, fft, fff)
  as.vector(kronecker(kronecker(a, b), c))
}
#' @keywords internal
#' @noRd
.kron4 <- function(a, b, c, d) {
  # 16-vector layout:
  # (tttt, tttf, ttft, ttff, tftt, tftf, tfft, tfff, fttt, fttf, ftft, ftff, fftt, fftf, ffft, ffff)
  as.vector(kronecker(kronecker(kronecker(a, b), c), d))
}
# ---------- Scale normalization helper ----------
#' @keywords internal
.norm_scales <- function(val, dt, df, pow_t, pow_f) {
  eps <- .Machine$double.eps^0.5
  dt <- if (dt > 0) dt else eps
  df <- if (df > 0) df else eps
  val / (dt^pow_t * df^pow_f)
}
