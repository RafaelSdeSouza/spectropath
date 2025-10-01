#' @export
sig_features <- function(path, depth = 3,
                         normalize =  TRUE) {
  validate_path_2d(path)
  if (depth == 2) sig <- ._stream2sig_d2_raw(path)
  else if (depth == 3) sig <- .stream2sig_d3_raw(path)
  else if (depth == 4) sig <- .stream2sig_d4_raw(path)
  else stop("depth in {2,3,4}")

  Lt <- if (depth == 2) {
    list(L1 = sig$L1, L2 = sig$L2 - 0.5 * .kron2(sig$L1, sig$L1))
  } else {
    .logsig_tensor_d4(sig)
  }

  dt <- max(path[,1]) - min(path[,1]); df <- max(path[,2]) - min(path[,2])
  feats <- list(
    drift_f = if (normalize) (path[nrow(path),2]-path[1,2])/(ifelse(df>0,df,Inf)) else (path[nrow(path),2]-path[1,2]),
    levy_tf = if (depth>=2) .norm_scales(0.5*(Lt$L2[1,2]-Lt$L2[2,1]), dt, df, 1, 1) else NA_real_,
    curl_t  = if (depth>=3) .norm_scales(0.5*(Lt$L3[2]-Lt$L3[3]), dt, df, 2, 1) else NA_real_,
    curl_f  = if (depth>=3) .norm_scales(0.5*(Lt$L3[6]-Lt$L3[7]), dt, df, 1, 2) else NA_real_,
    jerk_t  = if (depth>=4) .norm_scales(0.25*(Lt$L4[2]-Lt$L4[3]+Lt$L4[5]-Lt$L4[9]), dt, df, 3, 1) else NA_real_,
    twist_tf= if (depth>=4) .norm_scales(0.25*(Lt$L4[6]-Lt$L4[7]+Lt$L4[10]-Lt$L4[13]), dt, df, 2, 2) else NA_real_,
    jerk_f  = if (depth>=4) .norm_scales(0.25*(Lt$L4[15]-Lt$L4[14]+Lt$L4[11]-Lt$L4[8]), dt, df, 1, 3) else NA_real_
  )
  as.data.frame(feats)
}
