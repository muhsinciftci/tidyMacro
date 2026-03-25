#' Plot Historical Decomposition
#'
#' Plots the historical decomposition for multiple variables, showing the
#' contribution of each structural shock as individual stacked areas, plus
#' a total line overlay. Mirrors the MATLAB Bloom replication figure with
#' facets per variable.
#'
#' @param histdec_list Named list of (T-p) x N matrices, one per variable to
#'   plot. Each element is the \code{$histdec} output of \code{fhistdec()}
#'   called for that variable. Names are used as facet labels.
#' @param shock Integer (1-indexed) selecting the primary shock of interest.
#' @param shockname Character string label for the selected shock (legend entry).
#' @param shocknames Optional character vector of length N with names for
#'   \emph{all} shocks. When supplied the non-selected shocks are shown as
#'   individually named stacked areas instead of a single "Residual" bar.
#'   Default \code{NULL} (falls back to "Residual").
#' @param dates Optional vector of dates (Date, numeric, or integer). Can be
#'   the full T-length vector (trimmed internally using \code{p}) or the
#'   already-trimmed (T-p)-length vector. Default \code{NULL} uses an integer
#'   index.
#' @param p Integer lag order of the VAR. Used only when \code{dates} has
#'   length T to trim to \code{dates[(p+1):T]}. Default \code{0}.
#' @param return_data Logical. If \code{TRUE} returns the long-format tibble
#'   instead of the plot. Default \code{FALSE}.
#' @param facet_ncol Number of columns in the facet grid. Default \code{2}.
#' @param facet_scales Scale option for \code{facet_wrap}: one of
#'   \code{"free"} (default), \code{"free_y"}, \code{"free_x"},
#'   or \code{"fixed"}.
#' @param shock_fill Fill color for the primary shock area
#'   (default: \code{"#407EC9"}).
#' @param other_fills Character vector of fill colors for the remaining shocks.
#'   Recycled if shorter than \code{N - 1}. Default \code{"#E87722"} (single
#'   color used for the "Residual" bar or cycled across named other shocks).
#' @param total_color Color for the total line overlay
#'   (default: \code{"#1a1a1a"}).
#' @param total_linewidth Line width for the total line (default: \code{0.8}).
#' @param area_alpha Transparency for the stacked areas (default: \code{1}).
#'
#' @return A ggplot object or a tibble (if \code{return_data = TRUE}) with
#'   columns \code{date}, \code{variable}, \code{series}, \code{value}.
#'
#' @details
#' When \code{shocknames} is \code{NULL} all non-selected shocks are summed
#' into a single "Residual" bar (original behaviour, suitable for the Bloom
#' replication). When \code{shocknames} is supplied each shock gets its own
#' named stacked area — useful for the BQ case where the second shock has a
#' meaningful label (e.g. "Non-Technology").
#'
#' @examples
#' \dontrun{
#' # BQ case — name both shocks
#' fplot_histdec(
#'     histdec_list = histdec_all,
#'     shock        = 1,
#'     shockname    = "Technology",
#'     shocknames   = c("Technology", "Non-Technology"),
#'     dates        = dates,
#'     p            = p,
#'     facet_ncol   = 1
#' )
#'
#' # Cholesky case — collapse others into "Residual" (default)
#' fplot_histdec(
#'     histdec_list = histdec_list,
#'     shock        = shock,
#'     shockname    = "Uncertainty",
#'     dates        = dates_vec,
#'     p            = 12
#' )
#' }
#'
#' @seealso \code{\link{fhistdec}}
#'
#' @importFrom ggplot2 ggplot aes geom_area geom_line labs
#'   scale_fill_manual scale_color_manual facet_wrap
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#'
#' @export
fplot_histdec <- function(histdec_list,
                          shock,
                          shockname,
                          shocknames      = NULL,
                          dates           = NULL,
                          p               = 0,
                          return_data     = FALSE,
                          facet_ncol      = 2,
                          facet_scales    = c("free", "free_y", "free_x", "fixed"),
                          shock_fill      = "#407EC9",
                          other_fills     = "#E87722",
                          total_color     = "#1a1a1a",
                          total_linewidth = 0.8,
                          area_alpha      = 1) {

    facet_scales <- match.arg(facet_scales)

    # ------------------------------------------------------------------ #
    # Validate list input
    # ------------------------------------------------------------------ #
    if (!is.list(histdec_list) || length(histdec_list) == 0)
        stop("histdec_list must be a non-empty named list of (T-p) x N matrices")
    if (is.null(names(histdec_list)))
        names(histdec_list) <- paste0("Var", seq_along(histdec_list))

    T_eff <- nrow(as.matrix(histdec_list[[1]]))
    N     <- ncol(as.matrix(histdec_list[[1]]))

    for (nm in names(histdec_list)) {
        m <- as.matrix(histdec_list[[nm]])
        if (nrow(m) != T_eff || ncol(m) != N)
            stop(sprintf(
                "All matrices must be %d x %d; '%s' is %d x %d",
                T_eff, N, nm, nrow(m), ncol(m)
            ))
    }

    if (shock < 1 || shock > N)
        stop(sprintf("shock must be between 1 and %d", N))

    if (!is.null(shocknames) && length(shocknames) != N)
        stop(sprintf("shocknames must have length %d (one per shock)", N))

    # ------------------------------------------------------------------ #
    # Date axis
    # ------------------------------------------------------------------ #
    if (is.null(dates)) {
        date_vec <- seq_len(T_eff)
    } else if (length(dates) == T_eff) {
        date_vec <- dates
    } else if (length(dates) == T_eff + p) {
        date_vec <- dates[(p + 1L):length(dates)]
    } else {
        stop(sprintf(
            "dates length (%d) must equal nrow(histdec) (%d) or nrow(histdec)+p (%d)",
            length(dates), T_eff, T_eff + p
        ))
    }

    var_names    <- names(histdec_list)
    other_shocks <- setdiff(seq_len(N), shock)   # indices of non-selected shocks

    # Labels for non-selected shocks
    if (!is.null(shocknames)) {
        other_labels <- shocknames[other_shocks]
    } else {
        other_labels <- "Residual"               # collapse all others
    }

    # ------------------------------------------------------------------ #
    # Build long-format tibble
    # ------------------------------------------------------------------ #
    df_list <- vector("list", length(histdec_list))

    for (i in seq_along(histdec_list)) {
        hd    <- as.matrix(histdec_list[[i]])
        total <- rowSums(hd)

        if (!is.null(shocknames)) {
            # One column per shock
            row_df <- tibble::tibble(date = date_vec, variable = var_names[i])
            for (s in seq_len(N)) {
                row_df[[shocknames[s]]] <- as.numeric(hd[, s])
            }
            row_df[["Total"]] <- as.numeric(total)
        } else {
            # Selected shock + collapsed Residual
            row_df <- tibble::tibble(
                date             = date_vec,
                variable         = var_names[i],
                !!shockname      := as.numeric(hd[, shock]),
                Residual         = as.numeric(rowSums(hd[, other_shocks, drop = FALSE])),
                Total            = as.numeric(total)
            )
        }

        df_list[[i]] <- row_df
    }

    df_wide          <- dplyr::bind_rows(df_list)
    df_wide$variable <- factor(df_wide$variable, levels = var_names)

    if (return_data)
        return(tibble::as_tibble(df_wide))

    # ------------------------------------------------------------------ #
    # Build fill map
    # ------------------------------------------------------------------ #
    if (!is.null(shocknames)) {
        # Cycle other_fills across non-selected shocks
        n_other     <- length(other_labels)
        other_cols  <- rep_len(other_fills, n_other)
        fill_map    <- stats::setNames(
            c(other_cols,          shock_fill),
            c(other_labels,        shockname)
        )
        # All shock columns in stacking order: others first, selected on top
        stack_cols <- c(other_labels, shockname)
    } else {
        fill_map   <- stats::setNames(
            c(other_fills[1], shock_fill),
            c("Residual",     shockname)
        )
        stack_cols <- c("Residual", shockname)
    }

    color_map <- c(Total = total_color)

    # ------------------------------------------------------------------ #
    # Plot — pivot to long for stacked geom_area
    # ------------------------------------------------------------------ #
    df_long <- tidyr::pivot_longer(
        df_wide,
        cols      = tidyselect::all_of(stack_cols),
        names_to  = "series",
        values_to = "value"
    )
    df_long$series <- factor(df_long$series, levels = stack_cols)

    p_out <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$date)) +
        ggplot2::geom_area(
            ggplot2::aes(y = .data$value, fill = .data$series),
            position = "stack",
            alpha    = area_alpha
        ) +
        ggplot2::geom_line(
            data      = df_wide,
            ggplot2::aes(y = .data$Total, color = "Total"),
            linewidth = total_linewidth
        ) +
        ggplot2::facet_wrap(
            ~ .data$variable,
            scales = facet_scales,
            ncol   = facet_ncol
        ) +
        ggplot2::scale_fill_manual(values  = fill_map) +
        ggplot2::scale_color_manual(values = color_map) +
        ggplot2::labs(x = NULL, y = NULL, fill = NULL, color = NULL)

    return(p_out)
}
