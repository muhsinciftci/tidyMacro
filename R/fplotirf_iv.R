#' Plot IV/Proxy SVAR Impulse Response Functions with Confidence Bands
#'
#' Creates a grid of IRF plots with shaded confidence bands using ggplot2.
#' Accepts the full list returned by \code{fbootstrapIV_mbb()} directly.
#'
#' @param result_mbb List returned by \code{fbootstrapIV_mbb()}. Must contain
#'   \code{meanirf}, \code{upper}, and \code{lower}. If \code{upper2} and
#'   \code{lower2} are present, a second (inner) confidence band is drawn.
#' @param varnames Character vector of variable names (length N)
#' @param shockname Character string for the shock label
#' @param scale Numeric scalar applied to all IRF values before plotting.
#'   Useful for unit rescaling (e.g. \code{scale = 10} to express as percent
#'   of a 10-unit shock). Default 1.
#' @param return_data Logical, if TRUE returns tibble instead of plot. Default FALSE.
#' @param ribbon_fill Color for confidence band fill (default: "#407EC9")
#' @param ribbon_alpha Transparency for outer confidence band (default: 0.2)
#' @param ribbon_alpha2 Transparency for inner confidence band (default: 0.35)
#' @param line_color Color for point estimate line (default: "#910048")
#' @param zero_line_color Color for horizontal zero line (default: "#707372")
#' @param facet_scales Facet scales option: "free", "free_y", "free_x", or
#'   "fixed" (default: "free")
#' @param facet_ncol Number of columns in facet grid (default: 2)
#'
#' @return A ggplot object (if \code{return_data = FALSE}) or a tibble with
#'   columns \code{variable}, \code{horizon}, \code{point}, \code{upper},
#'   \code{lower} (and optionally \code{upper2}, \code{lower2}).
#'
#' @export
fplotirf_iv <- function(result_mbb, varnames, shockname,
                        scale           = 1,
                        return_data     = FALSE,
                        ribbon_fill     = "#407EC9",
                        ribbon_alpha    = 0.2,
                        ribbon_alpha2   = 0.35,
                        line_color      = "#910048",
                        zero_line_color = "#707372",
                        facet_scales    = "free",
                        facet_ncol      = 2) {

    point  <- result_mbb$meanirf * scale
    upper  <- result_mbb$upper   * scale
    lower  <- result_mbb$lower   * scale

    has_second_band <- !is.null(result_mbb$upper2) && !is.null(result_mbb$lower2)
    if (has_second_band) {
        upper2 <- result_mbb$upper2 * scale
        lower2 <- result_mbb$lower2 * scale
    }

    N       <- nrow(point)
    horizon <- ncol(point)

    irf_list <- vector("list", N)
    for (i in seq_len(N)) {
        row <- tibble::tibble(
            variable = varnames[i],
            horizon  = 0L:(horizon - 1L),
            point    = as.numeric(point[i, ]),
            upper    = as.numeric(upper[i, ]),
            lower    = as.numeric(lower[i, ]),
            shock    = shockname
        )
        if (has_second_band) {
            row$upper2 <- as.numeric(upper2[i, ])
            row$lower2 <- as.numeric(lower2[i, ])
        }
        irf_list[[i]] <- row
    }

    irf_data <- do.call(rbind, irf_list)
    irf_data$variable <- factor(irf_data$variable, levels = varnames)

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
                            color = zero_line_color, linetype = "dashed", linewidth = 0.6) +
        ggplot2::facet_wrap(~ variable, scales = facet_scales, ncol = facet_ncol)

    return(p)
}
