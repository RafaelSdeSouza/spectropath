validate_path_2d <- function(path) {
  if (!is.matrix(path) || !is.numeric(path))
    stop("`path` must be a numeric matrix.", call. = FALSE)
  if (ncol(path) != 2)
    stop("`path` must have exactly 2 columns: t and f.", call. = FALSE)
  if (nrow(path) < 2)
    stop("`path` must have at least 2 rows.", call. = FALSE)
  invisible(TRUE)
}
