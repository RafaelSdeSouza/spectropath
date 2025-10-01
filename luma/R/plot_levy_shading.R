# -----------------------------
# Internal: build Levy polygons
# -----------------------------

#' Build Levy-area polygons and diagnostics (internal)
#'
#' Constructs per-segment polygons between the curve and its end-to-end chord,
#' along with tidy per-segment data, resampled lines, and Levy totals.
#'
#' @param t Numeric time vector (strictly increasing after cleaning).
#' @param f Numeric flux vector, same length as `t`.
#' @param mode Resampling mode: "native", "uniform", or "max_dt".
#' @param grid_n Integer grid size for "uniform" mode.
#' @param max_dt Maximum segment length for "max_dt" mode.
#' @param normalize Logical; also compute normalized Levy area (levy_raw / (Delta_t * Delta_f)).
#' @importFrom stats approx
#' @return A list with components `polys`, `segments`, `lines`, `totals`.
#' @keywords internal
#' @noRd
shade_levy_polys <- function(t, f,
                             mode      = c("native", "uniform", "max_dt"),
                             grid_n    = 256,
                             max_dt    = NULL,
                             normalize = TRUE) {
  mode <- match.arg(mode)
  stopifnot(is.numeric(t), is.numeric(f), length(t) == length(f), length(t) >= 2)

  # sort / clean
  o <- order(t); t <- t[o]; f <- f[o]
  keep <- is.finite(t) & is.finite(f)
  t <- t[keep]; f <- f[keep]
  if ((t[length(t)] - t[1]) <= 0) stop("t must increase.")

  # optional resampling
  if (mode == "uniform") {
    t_new <- seq(t[1], t[length(t)], length.out = grid_n)
    f     <- approx(t, f, xout = t_new, rule = 2)$y
    t     <- t_new
  } else if (mode == "max_dt") {
    if (is.null(max_dt) || max_dt <= 0) stop("Provide positive max_dt for mode='max_dt'.")
    tt <- t[1]; ff <- f[1]
    for (k in seq_len(length(t) - 1)) {
      nsub <- max(1L, ceiling((t[k + 1] - t[k]) / max_dt))
      tk   <- seq(t[k], t[k + 1], length.out = nsub + 1L)
      fk   <- approx(c(t[k], t[k + 1]), c(f[k], f[k + 1]), xout = tk, rule = 2)$y
      tt <- c(tt, tk[-1]); ff <- c(ff, fk[-1])
    }
    t <- tt; f <- ff
  }

  # chord through endpoints
  dt_range <- t[length(t)] - t[1]
  slope <- (f[length(f)] - f[1]) / dt_range
  chord <- f[1] + slope * (t - t[1])

  # per-segment quantities
  i0 <- seq_len(length(t) - 1L); i1 <- i0 + 1L
  t0 <- t[i0]; t1 <- t[i1]
  f0 <- f[i0]; f1 <- f[i1]
  c0 <- chord[i0]; c1 <- chord[i1]
  dt <- t1 - t0

  # drop zero-length segments (duplicate t)
  ok <- dt > 0
  if (!all(ok)) {
    t0 <- t0[ok]; t1 <- t1[ok]
    f0 <- f0[ok]; f1 <- f1[ok]
    c0 <- c0[ok]; c1 <- c1[ok]
    dt <- dt[ok]
    i0 <- i0[ok]
  }

  # exact trapezoid area between curve and chord on each linear segment
  area_seg <- (((f0 + f1) - (c0 + c1)) / 2) * dt
  sign_seg <- ifelse(area_seg >= 0, "above", "below")

  # polygons for shading: (t0,f0)-(t1,f1)-(t1,c1)-(t0,c0)
  polys <- lapply(seq_along(i0), function(k) {
    data.frame(
      seg_id   = k,
      t        = c(t0[k], t1[k], t1[k], t0[k]),
      y        = c(f0[k], f1[k], c1[k], c0[k]),
      area_seg = area_seg[k],
      sign     = sign_seg[k],
      stringsAsFactors = FALSE
    )
  })
  polys <- do.call(rbind, polys)
  polys$sign <- factor(polys$sign, levels = c("above", "below"), ordered = TRUE)

  # tidy segments table
  segments <- data.frame(
    seg_id = seq_along(i0),
    t0 = t0, t1 = t1,
    f0 = f0, f1 = f1,
    c0 = c0, c1 = c1,
    dt = dt,
    area_seg = area_seg,
    sign = factor(sign_seg, levels = c("above", "below"), ordered = TRUE),
    stringsAsFactors = FALSE
  )

  # lines for plotting
  lines <- data.frame(t = t, f = f, chord = chord)

  # totals & scales
  levy_raw <- sum(area_seg)
  Delta_t  <- max(t) - min(t)
  Delta_f  <- max(f) - min(f)
  levy_norm <- if (normalize) levy_raw / (Delta_t * Delta_f) else NA_real_

  totals <- list(
    levy_raw  = levy_raw,
    levy_norm = levy_norm,
    Delta_t   = Delta_t,
    Delta_f   = Delta_f,
    area_pos  = sum(area_seg[area_seg > 0]),
    area_neg  = sum(area_seg[area_seg < 0])
  )

  list(polys = polys, segments = segments, lines = lines, totals = totals)
}

# -----------------------------
# Public: plotting wrapper
# -----------------------------

#' Plot Levy-area shading between curve and chord
#'
#' Draws the curve, its end-to-end chord, and shades signed areas between them.
#' The subtitle reports the raw and normalized Levy area and, when the
#' \pkg{luma} package is available, cross-checks against `luma::levy_tf()`.
#'
#' @param path Optional 2-column numeric matrix with columns (t, f).
#' @param t,f  Optional numeric vectors (ignored if `path` is provided).
#' @param mode Resampling mode for shading: "native", "uniform", or "max_dt".
#' @param grid_n Integer grid size for "uniform" mode.
#' @param max_dt Maximum segment length for "max_dt" mode.
#' @param normalize Logical; also compute normalized Levy area in subtitle.
#' @param fill_above,fill_below Fill colors for positive/negative strips.
#' @param alpha_fill Polygon alpha.
#' @param show_legend Logical; show fill legend.
#' @param subtitle Logical; show a subtitle with Levy values.
#' @return A list with `plot` (ggplot object) and `data` (the list returned by `shade_levy_polys()`).
#' @export
plot_levy_shading <- function(path = NULL, t = NULL, f = NULL,
                              mode = c("native","uniform","max_dt"),
                              grid_n = 256, max_dt = NULL,
                              normalize = TRUE,
                              fill_above = "orange", fill_below = "steelblue",
                              alpha_fill = 0.8, show_legend = FALSE,
                              subtitle = TRUE) {
  mode <- match.arg(mode)
  if (!is.null(path)) {
    stopifnot(is.matrix(path), ncol(path) == 2)
    t <- path[,1]; f <- path[,2]
  } else {
    stopifnot(is.numeric(t), is.numeric(f), length(t) == length(f))
  }

  # build shading data
  out <- shade_levy_polys(t, f, mode = mode, grid_n = grid_n,
                          max_dt = max_dt, normalize = normalize)

  # optional comparison with luma::levy_tf if available
  have_luma <- requireNamespace("luma", quietly = TRUE)
  if (have_luma) {
    path_mat <- cbind(t = out$lines$t, f = out$lines$f)
    levy_raw_pkg  <- luma::levy_tf(path_mat, normalize = FALSE)
    levy_norm_pkg <- luma::levy_tf(path_mat, normalize = TRUE)
  } else {
    levy_raw_pkg  <- NA_real_
    levy_norm_pkg <- NA_real_
  }

  # subtitle text
  sub_txt <- NULL
  if (subtitle) {
    if (have_luma) {
      sub_txt <- sprintf("Levy raw: %.6f (pkg: %.6f)   Levy norm: %.6g (pkg: %.6g)",
                         out$totals$levy_raw, levy_raw_pkg,
                         out$totals$levy_norm, levy_norm_pkg)
    } else {
      sub_txt <- sprintf("Levy raw: %.6f   Levy norm: %.6g",
                         out$totals$levy_raw, out$totals$levy_norm)
    }
  }

  # plot
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required. Install with install.packages('ggplot2').")
  }
  gg <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = out$polys,
                          ggplot2::aes(x = t, y = y, group = seg_id, fill = sign),
                          alpha = alpha_fill, color = NA) +
    ggplot2::geom_path(data = out$lines,
                       ggplot2::aes(x = t, y = f),
                       linewidth = 0.8, linetype = "dashed", color = "gray30") +
    ggplot2::geom_path(data = out$lines,
                       ggplot2::aes(x = t, y = chord),
                       linewidth = 0.7, linetype = 3, color = "black") +
    ggplot2::scale_fill_manual(values = c(above = fill_above, below = fill_below)) +
    ggplot2::labs(x = "time", y = "flux",
                  fill = "signed area",
                  subtitle = sub_txt) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = if (isTRUE(show_legend)) "right" else "none")

  list(plot = gg, data = out)
}
