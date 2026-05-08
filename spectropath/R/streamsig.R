# ---------- Raw signatures for piecewise-linear paths ----------
#' @keywords internal
.stream2sig_d2_raw <- function(path) {
  validate_path_2d(path)
  inc <- diff(path)
  acc <- list(L1 = c(0, 0), L2 = matrix(0, 2, 2))
  if (nrow(inc) == 0) return(acc)
  for (r in seq_len(nrow(inc))) acc <- .chen_prod_d2(acc, .sig_inc_d2(inc[r, ]))
  acc
}

#' @keywords internal
.stream2sig_d3_raw <- function(path) {
  validate_path_2d(path)
  inc <- diff(path)
  acc <- list(L1 = c(0, 0), L2 = matrix(0, 2, 2), L3 = numeric(8))
  if (nrow(inc) == 0) return(acc)
  for (r in seq_len(nrow(inc))) acc <- .chen_prod_d3(acc, .sig_inc_d3(inc[r, ]))
  acc
}

#' @keywords internal
.stream2sig_d4_raw <- function(path) {
  validate_path_2d(path)
  inc <- diff(path)
  acc <- list(L1 = c(0,0), L2 = matrix(0,2,2), L3 = numeric(8), L4 = numeric(16))
  if (nrow(inc) == 0) return(acc)
  for (r in seq_len(nrow(inc))) acc <- .chen_prod_d4(acc, .sig_inc_d4(inc[r, ]))
  acc
}
