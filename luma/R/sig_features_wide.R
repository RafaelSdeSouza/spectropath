#' Batch signature features for a "wide" matrix (rows = objects, cols = samples)
#'
#' @param mat numeric matrix, rows = objects, cols = samples along a common axis
#' @param axis numeric vector for the horizontal axis (length = ncol(mat)).
#'   If NULL, tries as.numeric(colnames(mat)); if that fails, uses seq_len(ncol(mat)).
#' @param depth integer in {2,3,4} for feature depth
#' @param normalize logical; pass through to scalar feature fns
#' @param id_col character; name of the id column in the output
#' @return data.frame with one row per object and signature features as columns
#' @export
sig_features_wide <- function(mat, axis = NULL, depth = 3, normalize = TRUE, id_col = "id") {
  if (!is.matrix(mat) || !is.numeric(mat)) stop("mat must be a numeric matrix.")
  if (is.null(axis)) {
    cn <- colnames(mat)
    if (!is.null(cn)) {
      suppressWarnings(ax_try <- as.numeric(cn))
    } else ax_try <- NA_real_
    if (all(is.finite(ax_try))) axis <- ax_try else axis <- seq_len(ncol(mat))
  }
  if (length(axis) != ncol(mat)) stop("length(axis) must equal ncol(mat).")
  if (!depth %in% c(2,3,4)) stop("depth must be 2, 3, or 4.")

  # build an id from rownames or 1..n
  ids <- rownames(mat)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(mat)))

  out_list <- vector("list", nrow(mat))
  for (i in seq_len(nrow(mat))) {
    f <- as.numeric(mat[i, ])
    keep <- is.finite(axis) & is.finite(f)
    t <- axis[keep]; f <- f[keep]
    if (length(t) < 2L) {
      base <- c(drift_f = NA_real_, levy_tf = NA_real_)
      if (depth >= 3) base <- c(base, curl_t = NA_real_, curl_f = NA_real_)
      if (depth >= 4) base <- c(base, jerk_t = NA_real_, twist_tf = NA_real_, jerk_f = NA_real_)
      out_list[[i]] <- base
      next
    }
    o <- order(t); t <- t[o]; f <- f[o]
    path <- cbind(t = t, f = f)

    feats <- c(
      drift_f = drift_f(path, normalize = normalize),
      levy_tf = levy_tf(path, normalize = normalize)
    )
    if (depth >= 3) {
      feats <- c(feats,
                 curl_t = curl_t(path, normalize = normalize),
                 curl_f = curl_f(path, normalize = normalize)
      )
    }
    if (depth >= 4) {
      feats <- c(feats,
                 jerk_t   = jerk_t(path,  normalize = normalize),
                 twist_tf = twist_tf(path, normalize = normalize),
                 jerk_f   = jerk_f(path,  normalize = normalize)
      )
    }
    out_list[[i]] <- feats
  }

  feat_mat <- do.call(rbind, out_list)
  feat_df  <- data.frame(id = ids, feat_mat, check.names = FALSE, row.names = NULL)
  names(feat_df)[1] <- id_col
  feat_df
}

