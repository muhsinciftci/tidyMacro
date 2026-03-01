#' Plot Historical Decomposition
#'
#' Plots the historical decomposition for multiple variables, showing the
#' contribution of the chosen structural shock, the residual (all other shocks),
#' and the total explained variation as a stacked bar chart with a total line
#' overlay. Mirrors the MATLAB Bloom replication figure with facets per variable.
#'
#' @param histdec_list Named list of (T-p) x N matrices, one per variable to
#'   plot. Each element is the \code{$histdec} output of \code{fhistdec()} called
#'   for that variable. Names are used as facet labels.
#' @param shock Integer (1-indexed) selecting the shock of interest.
#' @param shockname Character string label for the selected shock (legend entry).
#' @param dates Optional vector of dates (Date, numeric, or integer). Can be the
#'   full T-length vector (trimmed internally using \code{p}) or the already-
#'   trimmed (T-p)-length vector. Default \code{NULL} uses an integer index.
#' @param p Integer lag order of the VAR. Used only when \code{dates} has length
#'   T (full sample) to trim to \code{dates[(p+1):T]}. Default \code{0L}.
#' @param return_data Logical. If \code{TRUE} returns the wide-format tibble
#'   (columns: \code{date}, \code{variable}, shockname, \code{Residual},
#'   \code{Total}) instead of the plot. Default \code{FALSE}.
#' @param facet_ncol Number of columns in the facet grid. Default \code{2}.
#' @param shock_fill Fill color for the shock contribution area
#'   (default: \code{"#407EC9"}).
#' @param resid_fill Fill color for the residual (other shocks) area
#'   (default: \code{"#E87722"}).
#' @param total_color Color for the total line overlay
#'   (default: \code{"#1a1a1a"}).
#' @param total_linewidth Line width for the total line (default: \code{0.8}).
#' @param area_alpha Transparency for the stacked areas (default: \code{1}).
#'
#' @return A ggplot object (if \code{return_data = FALSE}) or a tibble with
#'   columns \code{date}, \code{variable}, \code{series}, and \code{value}.
#'
#' @details
#' For each variable and each time period t the plot shows:
#' \itemize{
#'   \item \strong{Shock bar} — \code{histdec[t, shock]}: contribution of the
#'     selected structural shock.
#'   \item \strong{Residual bar} — sum of \code{histdec[t, -shock]}: combined
#'     contribution of all other shocks.
#'   \item \strong{Total line} — sum of both bars: total explained variation.
#' }
#'
#' @examples
#' \dontrun{
#' var_result <- fVAR(y, p = 12, c = 1)
#' K      <- t(chol(var_result$sigma_full))
#' shock  <- match("UNCERT", colnames(y))
#' series <- setdiff(seq_len(ncol(y)), shock)
#'
#' # Call fhistdec for each variable and collect into a named list
#' histdec_list <- setNames(
#'     lapply(series, function(i) fhistdec(y, var_result, K, i)$histdec),
#'     colnames(y)[series]
#' )
#'
#' fplot_histdec(
#'     histdec_list = histdec_list,
#'     shock        = shock,
#'     shockname    = "Uncertainty",
#'     dates        = dates_vec,
#'     p            = 12
#' ) + ftheme_tidyMacro()
#' }
#'
#' @seealso \code{\link{fhistdec}}
#'
#' @importFrom ggplot2 ggplot aes geom_area geom_line geom_hline labs
#'   scale_fill_manual scale_color_manual facet_wrap
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#'
#' @export
fplot_histdec <- function(histdec_list,
                          shock,
                          shockname,
                          dates           = NULL,
                          p               = 0L,
                          return_data     = FALSE,
                          facet_ncol      = 2,
                          shock_fill      = "#407EC9",
                          resid_fill      = "#E87722",
                          total_color     = "#1a1a1a",
                          total_linewidth = 0.8,
                          area_alpha      = 1) {

    # ------------------------------------------------------------------ #
    # Validate list input
    # ------------------------------------------------------------------ #
    if (!is.list(histdec_list) || length(histdec_list) == 0) {
        stop("histdec_list must be a non-empty named list of (T-p) x N matrices")
    }
    if (is.null(names(histdec_list))) {
        names(histdec_list) <- paste0("Var", seq_along(histdec_list))
    }

    # All matrices must have the same dimensions
    T_eff <- nrow(as.matrix(histdec_list[[1]]))
    N     <- ncol(as.matrix(histdec_list[[1]]))

    for (nm in names(histdec_list)) {
        m <- as.matrix(histdec_list[[nm]])
        if (nrow(m) != T_eff || ncol(m) != N) {
            stop(sprintf(
                "All matrices in histdec_list must be %d x %d; '%s' is %d x %d",
                T_eff, N, nm, nrow(m), ncol(m)
            ))
        }
    }

    if (shock < 1L || shock > N) {
        stop(sprintf("shock must be between 1 and %d", N))
    }

    # ------------------------------------------------------------------ #
    # Date axis — accept full T-length or pre-trimmed (T-p)-length vector
    # ------------------------------------------------------------------ #
    if (is.null(dates)) {
        date_vec <- seq_len(T_eff)
    } else if (length(dates) == T_eff) {
        date_vec <- dates
    } else if (length(dates) == T_eff + p) {
        date_vec <- dates[(p + 1L):length(dates)]
    } else {
        stop(sprintf(
            "dates length (%d) must equal nrow(histdec) (%d) or nrow(histdec) + p (%d)",
            length(dates), T_eff, T_eff + p
        ))
    }

    # ------------------------------------------------------------------ #
    # Build combined long-format tibble (wide wrt series, then pivot)
    # ------------------------------------------------------------------ #
    resid_label <- "Residual"
    var_names   <- names(histdec_list)

    df_list <- vector("list", length(histdec_list))

    for (i in seq_along(histdec_list)) {
        hd <- as.matrix(histdec_list[[i]])

        shock_contrib <- as.numeric(hd[, shock])
        resid_contrib <- rowSums(hd[, -shock, drop = FALSE])
        total         <- shock_contrib + resid_contrib

        df_list[[i]] <- tibble::tibble(
            date            = date_vec,
            variable        = var_names[i],
            !!shockname     := shock_contrib,
            Residual        = resid_contrib,
            Total           = total
        )
    }

    # Wide-format tibble: one row per (date, variable)
    df_wide <- dplyr::bind_rows(df_list)
    df_wide$variable <- factor(df_wide$variable, levels = var_names)

    if (return_data) {
        return(df_wide)
    }

    # ------------------------------------------------------------------ #
    # Plot — geom_area per column, geom_line for Total
    # ------------------------------------------------------------------ #
    fill_map  <- stats::setNames(c(resid_fill,   shock_fill),
                                 c("Residual",   shockname))
    color_map <- c(Total = total_color)

    p <- ggplot2::ggplot(df_wide, ggplot2::aes(x = .data$date)) +
        ggplot2::geom_area(
            ggplot2::aes(y = .data$Residual, fill = "Residual"),
            alpha = area_alpha
        ) +
        ggplot2::geom_area(
            ggplot2::aes(y = .data[[shockname]], fill = shockname),
            alpha = area_alpha
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = .data$Total, color = "Total"),
            linewidth = total_linewidth
        ) +
        ggplot2::facet_wrap(
            ~ .data$variable,
            scales = "free",
            ncol   = facet_ncol
        ) +
        ggplot2::scale_fill_manual(values  = fill_map) +
        ggplot2::scale_color_manual(values = color_map) +
        ggplot2::labs(
            x     = NULL,
            y     = NULL,
            fill  = NULL,
            color = NULL
        )

    return(p)
}
