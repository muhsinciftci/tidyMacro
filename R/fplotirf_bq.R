#' Plot Blanchard-Quah Long-Run Identified Impulse Response Functions
#'
#' Creates a grid of IRF plots with shaded confidence bands using ggplot2.
#' Designed for BQ long-run identified VARs. Supports one or multiple shocks;
#' when a single shock is supplied the facet shows only response variables.
#'
#' @param point Array of point estimates (N x N x horizon+1), typically from
#'   \code{fbqIRF()}. Variables listed in \code{cumulate} will be cumulated
#'   inside this function to recover level responses.
#' @param boot_result List returned by \code{fbootstrapBQ()}. Must contain
#'   \code{upper} and \code{lower}. If \code{upper2} and \code{lower2} are
#'   present, a second (inner) confidence band is drawn automatically.
#' @param varnames Character vector of variable names (length N).
#' @param shocknames Character vector of shock names (length equal to
#'   \code{length(shocks)}).
#' @param cumulate Integer vector (1-based) of variable indices whose point
#'   IRFs should be cumulated. Must match the \code{cumulate} argument passed
#'   to \code{fbootstrapBQ()}.
#' @param shocks Integer vector (1-based) of shock indices to plot.
#' @param return_data Logical, if TRUE returns a tibble instead of a plot.
#'   Default FALSE.
#' @param ribbon_fill Color for confidence band fill (default: "#407EC9")
#' @param ribbon_alpha Transparency for outer confidence band (default: 0.2)
#' @param ribbon_alpha2 Transparency for inner confidence band (default: 0.35)
#' @param line_color Color for point estimate line (default: "#910048")
#' @param zero_line_color Color for horizontal zero line (default: "#707372")
#' @param facet_scales Scale option for \code{facet_wrap}: one of
#'   \code{"free"} (default), \code{"free_y"}, \code{"free_x"}, or
#'   \code{"fixed"}.
#' @param facet_ncol Number of columns in the facet layout. Default NULL.
#'
#' @return A ggplot object (if \code{return_data = FALSE}) or a tibble.
#'
#' @export
fplotirf_bq <- function(point, boot_result, varnames, shocknames,
                        cumulate, shocks,
                        return_data     = FALSE,
                        ribbon_fill     = "#407EC9",
                        ribbon_alpha    = 0.2,
                        ribbon_alpha2   = 0.35,
                        line_color      = "#910048",
                        zero_line_color = "#707372",
                        facet_scales    = c("free", "free_y", "free_x", "fixed"),
                        facet_ncol      = NULL) {

    upper <- boot_result$upper
    lower <- boot_result$lower

    has_second_band <- !is.null(boot_result$upper2) && !is.null(boot_result$lower2)
    if (has_second_band) {
        upper2 <- boot_result$upper2
        lower2 <- boot_result$lower2
    }

    facet_scales <- match.arg(facet_scales)

    dims     <- dim(point)
    N        <- dims[1]
    horizon  <- dims[3]
    n_shocks <- length(shocks)

    if (is.null(facet_ncol))
        facet_ncol <- if (n_shocks == 1L) ceiling(sqrt(N)) else n_shocks

    if (length(varnames) != N)
        stop(sprintf("varnames must have length %d (number of variables)", N))
    if (length(shocknames) != n_shocks)
        stop("shocknames must have the same length as shocks")
    if (any(shocks < 1) || any(shocks > N))
        stop(sprintf("All shocks must be between 1 and %d", N))
    if (any(cumulate < 1) || any(cumulate > N))
        stop(sprintf("All cumulate indices must be between 1 and %d", N))

    # Cumulate point IRF (upper/lower already cumulated inside bootstrap loop)
    point_plot <- point
    for (s in shocks) {
        for (ci in cumulate) {
            point_plot[ci, s, ] <- cumsum(point_plot[ci, s, ])
        }
    }

    irf_list <- vector("list", N * n_shocks)
    idx <- 1L
    for (j in seq_along(shocks)) {
        s <- shocks[j]
        for (i in seq_len(N)) {
            row <- tibble::tibble(
                variable = varnames[i],
                shock    = shocknames[j],
                horizon  = 0L:(horizon - 1L),
                point    = as.numeric(point_plot[i, s, ]),
                upper    = as.numeric(upper[i, s, ]),
                lower    = as.numeric(lower[i, s, ])
            )
            if (has_second_band) {
                row$upper2 <- as.numeric(upper2[i, s, ])
                row$lower2 <- as.numeric(lower2[i, s, ])
            }
            irf_list[[idx]] <- row
            idx <- idx + 1L
        }
    }

    irf_data <- do.call(rbind, irf_list)
    irf_data$variable <- factor(irf_data$variable, levels = varnames)
    irf_data$shock    <- factor(irf_data$shock,    levels = shocknames)

    if (return_data)
        return(tibble::as_tibble(irf_data))

    p <- ggplot2::ggplot(irf_data, ggplot2::aes(x = horizon)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                             fill = ribbon_fill, alpha = ribbon_alpha)

    if (has_second_band) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower2, ymax = upper2),
                                      fill = ribbon_fill, alpha = ribbon_alpha2)
    }

    p <- p +
        ggplot2::geom_line(ggplot2::aes(y = point),
                           color = line_color, linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = 0,
                            color      = zero_line_color,
                            linetype   = "dashed",
                            linewidth  = 0.6) +
        ggplot2::labs(x = NULL, y = NULL, caption = NULL)

    if (n_shocks == 1L) {
        p <- p + ggplot2::facet_wrap(~ variable,
                                     scales = facet_scales,
                                     ncol   = facet_ncol)
    } else {
        p <- p + ggplot2::facet_wrap(~ shock + variable,
                                     scales = facet_scales,
                                     ncol   = facet_ncol)
    }

    return(p)
}
