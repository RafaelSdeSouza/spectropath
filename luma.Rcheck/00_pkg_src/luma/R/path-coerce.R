#' @title Ensure (t,f) matrix
#' @export
as_path2d <- function(x, t = NULL, f = NULL, sort = TRUE) {
  if (is.matrix(x) && is.numeric(x) && ncol(x) == 2) {
    m <- x
  } else if (is.data.frame(x)) {
    stopifnot(!is.null(t), !is.null(f))
    m <- cbind(t = as.numeric(x[[t]]), f = as.numeric(x[[f]]))
  } else {
    stop("as_path2d: supply a 2-col numeric matrix or a data.frame with t and f columns.")
  }
  keep <- is.finite(m[,1]) & is.finite(m[,2])
  m <- m[keep, , drop = FALSE]
  if (nrow(m) < 2L) stop("as_path2d: need at least 2 finite rows.")
  if (sort) {
    o <- order(m[,1])
    m <- m[o, , drop = FALSE]
  }
  colnames(m) <- c("t","f")
  m
}

#' @export
is_path2d <- function(x) is.matrix(x) && is.numeric(x) && identical(colnames(x), c("t","f")) && ncol(x)==2L
