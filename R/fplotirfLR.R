#' Plot Long-Run Identified Impulse Response Functions with Confidence Bands
#'
#' Creates a grid of IRF plots with shaded confidence bands using ggplot2.
#' Designed for long-run identified VARs where the point estimate is an N x H
#' matrix — works for both LR-Max (\code{fmaxIRF()}) and Uhlig max-share
#' (\code{fuhligIRF()}) identifications.
#'
#' @param point N x (horizon+1) matrix of point estimates, from
#'   \code{fmaxIRF()} or \code{fuhligIRF()}.
#' @param boot_result List returned by \code{fbootstrapMax()},
#'   \code{fbootstrapMaxCorrected()}, \code{fbootstrapUhlig()}, or
#'   \code{fbootstrapUhligCorrected()}. Must contain \code{upper} and
#'   \code{lower}. If \code{upper2} and \code{lower2} are present, a second
#'   (inner) confidence band is drawn automatically.
#' @param varnames Character vector of variable names (length N).
#' @param return_data Logical, if TRUE returns tibble instead of plot. Default FALSE.
#' @param ribbon_fill Color for confidence band fill (default: "#407EC9")
#' @param ribbon_alpha Transparency for outer confidence band (default: 0.2)
#' @param ribbon_alpha2 Transparency for inner confidence band (default: 0.35)
#' @param line_color Color for point estimate line (default: "#910048")
#' @param zero_line_color Color for horizontal zero line (default: "#707372")
#' @param facet_scales Facet scales option: "free", "free_y", "free_x", or
#'   "fixed" (default: "free")
#' @param facet_ncol Number of columns in facet grid (default: NULL, auto-calculated)
#'
#' @return A ggplot object (if \code{return_data = FALSE}) or a tibble.
#'
#' @export
fplotirfLR <- function(point, boot_result = NULL, varnames,
                       return_data     = FALSE,
                       ribbon_fill     = "#407EC9",
                       ribbon_alpha    = 0.2,
                       ribbon_alpha2   = 0.35,
                       line_color      = "#910048",
                       zero_line_color = "#707372",
                       facet_scales    = "free",
                       facet_ncol      = NULL) {

    has_bands <- !is.null(boot_result)
    if (has_bands) {
        upper <- boot_result$upper
        lower <- boot_result$lower
    }

    has_second_band <- has_bands &&
        !is.null(boot_result$upper2) && !is.null(boot_result$lower2)
    if (has_second_band) {
        upper2 <- boot_result$upper2
        lower2 <- boot_result$lower2
    }

    N       <- nrow(point)
    horizon <- ncol(point)

    if (length(varnames) != N)
        stop(sprintf("varnames must have length %d (number of variables)", N))
    if (is.null(facet_ncol))
        facet_ncol <- ceiling(sqrt(N))

    irf_list <- vector("list", N)
    for (i in seq_len(N)) {
        row <- tibble::tibble(
            variable = varnames[i],
            horizon  = 0L:(horizon - 1L),
            point    = as.numeric(point[i, ])
        )
        if (has_bands) {
            row$upper <- as.numeric(upper[i, ])
            row$lower <- as.numeric(lower[i, ])
        }
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

    p <- ggplot2::ggplot(irf_data, ggplot2::aes(x = horizon))

    if (has_bands) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                                      fill = ribbon_fill, alpha = ribbon_alpha)
    }

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
