#' Plot Local Projection Impulse Responses
#'
#' Plots impulse responses and confidence bands from an object returned by
#' \code{fLP()}. Handles both single-band and multi-band fits transparently:
#' if \code{fLP()} was called with \code{conf = c(0.68, 0.95)} the plot
#' automatically draws one ribbon per level (widest ribbon is drawn first with
#' lowest opacity so the tighter band is visually on top).
#'
#' @param x An object of class \code{"fLP"} returned by \code{fLP()}.
#' @param variables Optional character vector of response variables to plot.
#'   Defaults to all responses.
#' @param scale Numeric scalar applied to IRFs and confidence bands before
#'   plotting. Default is 1.
#' @param return_data Logical. If \code{TRUE}, return the plotting data
#'   (long format when multi-band) instead of a ggplot object.
#' @param ribbon_fill Color for the confidence band fill.
#' @param ribbon_alpha Transparency for the confidence band. For multi-band
#'   fits, either a single number (interpreted as the alpha of the widest band;
#'   inner bands get proportionally darker) or a numeric vector of the same
#'   length as the confidence levels, matched widest-to-tightest.
#' @param line_color Color for the IRF line.
#' @param zero_line_color Color for the horizontal zero line.
#' @param facet_scales Facet scales option: \code{"free_y"}, \code{"free"},
#'   \code{"free_x"}, or \code{"fixed"}.
#' @param facet_ncol Number of facet columns. Defaults to an automatic layout.
#'
#' @return A ggplot object, or a tibble if \code{return_data = TRUE}.
#'
#' @export
fPlotLP <- function(x,
                    variables       = NULL,
                    scale           = 1,
                    return_data     = FALSE,
                    ribbon_fill     = "#407EC9",
                    ribbon_alpha    = 0.2,
                    line_color      = "#910048",
                    zero_line_color = "#707372",
                    facet_scales    = "free_y",
                    facet_ncol      = NULL) {
  if (!inherits(x, "fLP")) {
    stop("fPlotLP: 'x' must be an object returned by fLP().", call. = FALSE)
  }
  if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) ||
      !is.finite(scale)) {
    stop("fPlotLP: 'scale' must be a finite numeric scalar.", call. = FALSE)
  }

  varnames <- colnames(x$irfs)
  if (is.null(varnames)) varnames <- x$lhs_vars
  if (is.null(varnames)) varnames <- paste0("response_", seq_len(ncol(x$irfs)))

  if (!is.null(variables)) {
    missing_vars <- setdiff(variables, varnames)
    if (length(missing_vars) > 0L) {
      stop(sprintf(
        "fPlotLP: variable(s) not found in LP result: %s",
        paste(missing_vars, collapse = ", ")
      ), call. = FALSE)
    }
    keep <- match(variables, varnames)
    varnames <- varnames[keep]
  } else {
    keep <- seq_along(varnames)
  }

  if (is.null(facet_ncol)) {
    facet_ncol <- ceiling(sqrt(length(varnames)))
  }

  horizon <- x$horizons
  if (is.null(horizon)) horizon <- seq_len(nrow(x$irfs)) - 1L

  multi <- is.list(x$irfs_upper)

  # Confidence keys ordered widest to tightest (fLP sorts conf descending).
  conf_keys <- if (multi) names(x$irfs_upper) else format(x$conf, trim = TRUE)

  # Resolve alphas: one per conf level, widest first.
  n_bands <- length(conf_keys)
  if (length(ribbon_alpha) == 1L) {
    if (n_bands == 1L) {
      alphas <- ribbon_alpha
    } else {
      # widest = base alpha; each tighter band gets +alpha_step
      alphas <- ribbon_alpha + (seq_len(n_bands) - 1L) *
        min(ribbon_alpha, (1 - ribbon_alpha) / n_bands)
    }
  } else {
    if (length(ribbon_alpha) != n_bands) {
      stop("fPlotLP: 'ribbon_alpha' must be scalar or match the number of confidence bands.",
           call. = FALSE)
    }
    alphas <- ribbon_alpha
  }
  names(alphas) <- conf_keys

  # ---- Build long-format plot data --------------------------------------
  get_band <- function(which, key, j) {
    m <- if (multi) {
      if (which == "upper") x$irfs_upper[[key]] else x$irfs_lower[[key]]
    } else {
      if (which == "upper") x$irfs_upper else x$irfs_lower
    }
    as.numeric(m[, j]) * scale
  }

  rows <- vector("list", length(varnames) * n_bands)
  idx  <- 1L
  for (i in seq_along(varnames)) {
    j <- keep[i]
    point_i <- as.numeric(x$irfs[, j]) * scale
    for (k in conf_keys) {
      rows[[idx]] <- tibble::tibble(
        variable = varnames[i],
        shock    = x$shock,
        horizon  = horizon,
        conf     = as.numeric(k),
        point    = point_i,
        upper    = get_band("upper", k, j),
        lower    = get_band("lower", k, j)
      )
      idx <- idx + 1L
    }
  }

  plot_data <- do.call(rbind, rows)
  plot_data$variable <- factor(plot_data$variable, levels = varnames)
  plot_data$conf     <- factor(plot_data$conf,
                               levels = as.numeric(conf_keys))  # widest first

  if (return_data) {
    # Single-band: keep the historical shape (no 'conf' column).
    if (!multi && n_bands == 1L) {
      plot_data$conf <- NULL
    }
    return(tibble::as_tibble(plot_data))
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = horizon))
  for (k in conf_keys) {
    band_k <- plot_data[plot_data$conf == as.numeric(k), , drop = FALSE]
    p <- p + ggplot2::geom_ribbon(
      data  = band_k,
      ggplot2::aes(ymin = lower, ymax = upper),
      fill  = ribbon_fill,
      alpha = alphas[[k]]
    )
  }

  p +
    ggplot2::geom_line(
      data = plot_data[plot_data$conf == as.numeric(conf_keys[1]), , drop = FALSE],
      ggplot2::aes(y = point),
      color = line_color,
      linewidth = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = zero_line_color,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::facet_wrap(~ variable, scales = facet_scales, ncol = facet_ncol) +
    ggplot2::labs(x = "Horizon", y = NULL)
}
