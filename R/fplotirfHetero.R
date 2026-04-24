#' Plot Heteroskedasticity Bootstrap IRFs with Confidence Bands
#'
#' Creates a grid of IRF plots for results from \code{fbootstrapHetero()}.
#' Uses the plug-in point estimate (\code{result$point}) as the center line,
#' which is the correct center for Hall (1992) recentered confidence bands.
#'
#' @param result      List returned by \code{fbootstrapHetero()}. Must contain
#'   \code{point}, \code{upper}, \code{lower}. If \code{upper2} and
#'   \code{lower2} are present, a second (inner) band is drawn.
#' @param varnames    Character vector of variable names (length N).
#' @param shockname   Character string for the shock label.
#' @param scale       Numeric scalar applied to all IRF values. Default 1.
#' @param return_data Logical; if \code{TRUE} returns a tibble. Default \code{FALSE}.
#' @param ribbon_fill    Fill colour for confidence bands. Default \code{"#407EC9"}.
#' @param ribbon_alpha   Transparency for outer band. Default 0.2.
#' @param ribbon_alpha2  Transparency for inner band. Default 0.35.
#' @param line_color     Colour for point estimate line. Default \code{"#910048"}.
#' @param zero_line_color Colour for zero reference line. Default \code{"#707372"}.
#' @param facet_scales   Facet scales: \code{"free"}, \code{"free_y"},
#'   \code{"free_x"}, or \code{"fixed"}. Default \code{"free"}.
#' @param facet_ncol     Number of columns in facet grid. Default 2.
#'
#' @return A ggplot object (or a tibble if \code{return_data = TRUE}).
#'
#' @export
fplotirfHetero <- function(result, varnames, shockname,
                           scale           = 1,
                           return_data     = FALSE,
                           ribbon_fill     = "#407EC9",
                           ribbon_alpha    = 0.2,
                           ribbon_alpha2   = 0.35,
                           line_color      = "#910048",
                           zero_line_color = "#707372",
                           facet_scales    = "free",
                           facet_ncol      = 2) {

    point  <- result$point  * scale
    upper  <- result$upper  * scale
    lower  <- result$lower  * scale

    has_second_band <- !is.null(result$upper2) && !is.null(result$lower2)
    if (has_second_band) {
        upper2 <- result$upper2 * scale
        lower2 <- result$lower2 * scale
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
                            color = zero_line_color, linetype = "dashed",
                            linewidth = 0.6) +
        ggplot2::facet_wrap(~ variable, scales = facet_scales, ncol = facet_ncol)

    return(p)
}
