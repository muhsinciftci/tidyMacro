#' Plot MSW Weak-IV Robust IRFs with Two Confidence Bands
#'
#' Creates a faceted IRF plot showing both the plug-in (delta method) and
#' Anderson-Rubin weak-IV robust confidence sets from \code{fMSW()}.
#' The outer (lighter) band is the Anderson-Rubin set; the inner (darker) band
#' is the delta-method interval.
#'
#' @param msw_result List returned by \code{fMSW()}. Must contain
#'   \code{IRF}, \code{Dmlbound}, \code{Dmubound}, \code{MSWlbound},
#'   \code{MSWubound}.
#' @param varnames Character vector of variable names (length N).
#' @param shockname Character string for the shock label.
#' @param scale Numeric scalar applied to all values before plotting. Default 1.
#' @param return_data Logical; if \code{TRUE} returns a tibble instead of a
#'   ggplot object. Default \code{FALSE}.
#' @param ribbon_fill_ar  Fill colour for Anderson-Rubin band. Default "#407EC9".
#' @param ribbon_fill_dm  Fill colour for delta-method band. Default "#407EC9".
#' @param ribbon_alpha_ar Transparency for AR outer band. Default 0.18.
#' @param ribbon_alpha_dm Transparency for delta-method inner band. Default 0.35.
#' @param line_color      Colour for point estimate line. Default "#910048".
#' @param zero_line_color Colour for zero reference line. Default "#707372".
#' @param facet_scales    Facet scales: "free", "free_y", "free_x", or "fixed".
#'   Default "free".
#' @param facet_ncol Number of columns in the facet grid. Default 2.
#'
#' @return A ggplot object (or a tibble if \code{return_data = TRUE}).
#'
#' @export
fplotirfMSW <- function(msw_result, varnames, shockname,
                         scale           = 1,
                         return_data     = FALSE,
                         ribbon_fill_ar  = "#407EC9",
                         ribbon_fill_dm  = "#407EC9",
                         ribbon_alpha_ar = 0.18,
                         ribbon_alpha_dm = 0.35,
                         line_color      = "#910048",
                         zero_line_color = "#707372",
                         facet_scales    = "free",
                         facet_ncol      = 2) {

    irf     <- msw_result$IRF       * scale
    dm_lo   <- msw_result$Dmlbound  * scale
    dm_hi   <- msw_result$Dmubound  * scale
    ar_lo   <- msw_result$MSWlbound * scale
    ar_hi   <- msw_result$MSWubound * scale

    N       <- nrow(irf)
    horizon <- ncol(irf)

    irf_list <- vector("list", N)
    for (i in seq_len(N)) {
        irf_list[[i]] <- tibble::tibble(
            variable = varnames[i],
            horizon  = 0L:(horizon - 1L),
            point    = as.numeric(irf[i, ]),
            dm_lo    = as.numeric(dm_lo[i, ]),
            dm_hi    = as.numeric(dm_hi[i, ]),
            ar_lo    = as.numeric(ar_lo[i, ]),
            ar_hi    = as.numeric(ar_hi[i, ]),
            shock    = shockname
        )
    }

    irf_data <- do.call(rbind, irf_list)
    irf_data$variable <- factor(irf_data$variable, levels = varnames)

    if (return_data)
        return(tibble::as_tibble(irf_data))

    ggplot2::ggplot(irf_data, ggplot2::aes(x = horizon)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = ar_lo, ymax = ar_hi, fill = "Anderson-Rubin"),
            alpha = ribbon_alpha_ar) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = dm_lo, ymax = dm_hi, fill = "Delta Method"),
            alpha = ribbon_alpha_dm) +
        ggplot2::geom_line(ggplot2::aes(y = point, color = "Point estimate"),
                           linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = 0,
                            color = zero_line_color, linetype = "dashed",
                            linewidth = 0.6) +
        ggplot2::scale_fill_manual(
            name   = NULL,
            values = c("Anderson-Rubin" = ribbon_fill_ar,
                       "Delta Method"   = ribbon_fill_dm),
            guide  = ggplot2::guide_legend(
                override.aes = list(alpha = c(ribbon_alpha_ar, ribbon_alpha_dm)))
        ) +
        ggplot2::scale_color_manual(
            name   = NULL,
            values = c("Point estimate" = line_color)
        ) +
        ggplot2::facet_wrap(~ variable, scales = facet_scales, ncol = facet_ncol)
}
