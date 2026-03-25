#' Plot Blanchard-Quah Long-Run Identified Impulse Response Functions
#'
#' Creates a grid of IRF plots with shaded confidence bands using ggplot2.
#' Designed for BQ long-run identified VARs. Supports one or multiple shocks;
#' when a single shock is supplied the facet shows only response variables.
#'
#' @param point Array of point estimates (N x N x horizon+1), typically from
#'   \code{fbqIRF()}. Variables listed in \code{cumulate} will be cumulated
#'   inside this function to recover level responses.
#' @param upper Array of upper confidence bounds (N x N x horizon+1), from
#'   \code{fbootstrapBQ()$upper}. Already cumulated inside the bootstrap loop.
#' @param lower Array of lower confidence bounds (N x N x horizon+1), from
#'   \code{fbootstrapBQ()$lower}. Already cumulated inside the bootstrap loop.
#' @param varnames Character vector of variable names (length N)
#' @param shocknames Character vector of shock names (length equal to
#'   \code{length(shocks)})
#' @param cumulate Integer vector (1-based) of variable indices whose point
#'   IRFs should be cumulated. Must match the \code{cumulate} argument passed
#'   to \code{fbootstrapBQ()}.
#' @param shocks Integer vector (1-based) of shock indices to plot.
#' @param prc Numeric, confidence level percentage (e.g. 68 or 95)
#' @param return_data Logical, if TRUE returns a tibble instead of a plot.
#'   Default FALSE.
#' @param ribbon_fill Color for confidence band fill (default: "#407EC9")
#' @param ribbon_alpha Transparency for confidence band (default: 0.3)
#' @param line_color Color for point estimate line (default: "#910048")
#' @param zero_line_color Color for horizontal zero line (default: "#707372")
#' @param facet_scales Scale option for \code{facet_wrap}: one of
#'   \code{"free"} (default), \code{"free_y"}, \code{"free_x"},
#'   or \code{"fixed"}.
#' @param facet_ncol Number of columns in the facet layout. Default NULL.
#'   For a single shock defaults to \code{ceiling(sqrt(N))}; for multiple
#'   shocks defaults to \code{length(shocks)} (one column per shock).
#'
#' @return A ggplot object (if \code{return_data = FALSE}) or a tibble with
#'   columns \code{variable}, \code{shock}, \code{horizon}, \code{point},
#'   \code{upper}, \code{lower} (if \code{return_data = TRUE}).
#'
#' @details
#' The point IRF from \code{fbqIRF()} is in first-difference units. This
#' function cumulates the rows specified in \code{cumulate} to recover level
#' responses, matching what \code{fbootstrapBQ()} does internally so that
#' the point estimate and confidence bands are on the same scale.
#'
#' When a single shock is supplied the facet formula is \code{~ variable},
#' keeping the layout clean. When multiple shocks are supplied the formula
#' becomes \code{~ shock + variable} so each shock gets its own column.
#'
#' The plot uses the currently active ggplot2 theme. For a consistent look
#' with other tidyMacro plots, use \code{ftheme_tidyMacro()} or set it
#' globally with \code{set_tidyMacro_theme()}.
#'
#' @examples
#' \dontrun{
#' var_result <- fVAR(y, p = p, c = 1)
#' wold       <- fwoldIRF(var_result, horizon = 40)
#' Sigma      <- var_result$sigma_full
#' C1         <- apply(wold, c(1, 2), sum)
#' D1         <- t(chol(C1 %*% Sigma %*% t(C1)))
#' K          <- solve(C1, D1)
#' point_irf  <- fbqIRF(wold, K)
#'
#' boot <- fbootstrapBQ(y, var_result, nboot = 1000, horizon = 40,
#'                      prc = 68, bootscheme = "residual",
#'                      cumulate = c(1, 2))
#'
#' # Single shock — facets by variable only
#' fplotirf_bq(
#'   point      = point_irf,
#'   upper      = boot$upper,
#'   lower      = boot$lower,
#'   varnames   = c("LABPROD", "HOURS"),
#'   shocknames = "Technology",
#'   cumulate   = c(1, 2),
#'   shocks     = 1,
#'   prc        = 68
#' )
#'
#' # Multiple shocks — facets by shock x variable
#' fplotirf_bq(
#'   point      = point_irf,
#'   upper      = boot$upper,
#'   lower      = boot$lower,
#'   varnames   = c("LABPROD", "HOURS"),
#'   shocknames = c("Technology", "Non-Technology"),
#'   cumulate   = c(1, 2),
#'   shocks     = c(1, 2),
#'   prc        = 68
#' )
#' }
#'
#' @export
fplotirf_bq <- function(point, upper, lower, varnames, shocknames,
                        cumulate, shocks, prc,
                        return_data     = FALSE,
                        ribbon_fill     = "#407EC9",
                        ribbon_alpha    = 0.3,
                        line_color      = "#910048",
                        zero_line_color = "#707372",
                        facet_scales    = c("free", "free_y", "free_x", "fixed"),
                        facet_ncol      = NULL) {

    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("Package 'ggplot2' is required. Please install it.", call. = FALSE)
    if (!requireNamespace("tibble", quietly = TRUE))
        stop("Package 'tibble' is required. Please install it.", call. = FALSE)

    facet_scales <- match.arg(facet_scales)

    dims     <- dim(point)
    N        <- dims[1]
    horizon  <- dims[3]
    n_shocks <- length(shocks)

    # Default ncol: sqrt layout for single shock, one col per shock otherwise
    if (is.null(facet_ncol))
        facet_ncol <- if (n_shocks == 1L) ceiling(sqrt(N)) else n_shocks

    # Validate inputs
    if (length(varnames) != N)
        stop(sprintf("varnames must have length %d (number of variables)", N))
    if (length(shocknames) != n_shocks)
        stop("shocknames must have the same length as shocks")
    if (any(shocks < 1) || any(shocks > N))
        stop(sprintf("All shocks must be between 1 and %d", N))
    if (any(cumulate < 1) || any(cumulate > N))
        stop(sprintf("All cumulate indices must be between 1 and %d", N))

    # Cumulate point IRF for selected rows (upper/lower already cumulated
    # inside fbootstrapBQ loop)
    point_plot <- point
    for (s in shocks) {
        for (ci in cumulate) {
            point_plot[ci, s, ] <- cumsum(point_plot[ci, s, ])
        }
    }

    # Build long-format tibble over all variable x shock combinations
    irf_list <- vector("list", N * n_shocks)
    idx <- 1L
    for (j in seq_along(shocks)) {
        s <- shocks[j]
        for (i in seq_len(N)) {
            irf_list[[idx]] <- tibble::tibble(
                variable = varnames[i],
                shock    = shocknames[j],
                horizon  = 0L:(horizon - 1L),
                point    = as.numeric(point_plot[i, s, ]),
                upper    = as.numeric(upper[i, s, ]),
                lower    = as.numeric(lower[i, s, ])
            )
            idx <- idx + 1L
        }
    }

    irf_data <- do.call(rbind, irf_list)

    # Preserve ordering in facets
    irf_data$variable <- factor(irf_data$variable, levels = varnames)
    irf_data$shock    <- factor(irf_data$shock,    levels = shocknames)

    if (return_data)
        return(tibble::as_tibble(irf_data))

    p <- ggplot2::ggplot(irf_data, ggplot2::aes(x = horizon)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                             fill = ribbon_fill, alpha = ribbon_alpha) +
        ggplot2::geom_line(ggplot2::aes(y = point),
                           color = line_color, linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = 0,
                            color      = zero_line_color,
                            linetype   = "dashed",
                            linewidth  = 0.6) +
        ggplot2::labs(x = NULL, y = NULL, caption = NULL)

    # Single shock: facet by variable only (cleaner, no redundant shock label)
    # Multiple shocks: facet by shock x variable (one column per shock)
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
