#' Plot Impulse Responses from a Sign-Restricted SVAR
#'
#' Draws the posterior median and credible band of the responses to one
#' structural shock, for every variable in the system. Two models can be
#' overlaid, which is the standard way of showing what narrative restrictions
#' add to plain sign restrictions.
#'
#' @param sr An object of class \code{"fSignRestr"} from
#'   \code{\link{fSignRestr}}.
#' @param shock Integer (1-indexed) structural shock to plot (default 1).
#' @param compare Optional second \code{"fSignRestr"} object drawn on the same
#'   panels for comparison (default NULL).
#' @param labels Character vector of length 1 or 2 naming the models in the
#'   legend (default \code{c("Sign", "Sign + narrative")}).
#' @param varnames Optional character vector overriding the variable names.
#' @param scale Numeric vector of length 1 or N multiplying the responses, e.g.
#'   to express a variable in basis points (default 1).
#' @param scale_compare Numeric scaling applied to \code{compare}; defaults to
#'   \code{scale}.
#' @param return_data Logical; return the plotting tibble instead of the plot
#'   (default FALSE).
#' @param colors Length-2 character vector of colours for the two models.
#' @param ribbon_alpha Transparency of the credible bands (default 0.18).
#' @param conf Optional numeric vector of band levels to draw, in percent. Each
#'   must have been requested in the \code{fSignRestr} call. Defaults to the
#'   first level of the fitted object. Several levels are drawn as nested
#'   ribbons, the widest most transparent.
#' @param facet_scales Facet scale option passed to
#'   \code{ggplot2::facet_wrap} (default "free_y").
#' @param facet_ncol Number of facet columns (default NULL, auto).
#'
#' @return A \code{ggplot} object, or a tibble when \code{return_data = TRUE}.
#'
#' @seealso \code{\link{fSignRestr}}, \code{\link{fPlotVarDec}}
#'
#' @export
fPlotIRFSign <- function(sr, shock = 1, compare = NULL,
                         labels        = c("Sign", "Sign + narrative"),
                         varnames      = NULL,
                         scale         = 1,
                         scale_compare = NULL,
                         return_data   = FALSE,
                         colors        = c("#407EC9", "#E4002B"),
                         ribbon_alpha  = 0.18,
                         conf          = NULL,
                         facet_scales  = "free_y",
                         facet_ncol    = NULL) {

    if (!inherits(sr, "fSignRestr")) stop("`sr` must be an 'fSignRestr' object.")

    if (is.null(conf)) conf <- sr$conf[1]
    conf <- sort(as.numeric(conf), decreasing = TRUE)

    ## Bands come from the `bands` list so that a level the model was not asked
    ## for fails loudly instead of being silently drawn at the wrong coverage.
    band_of <- function(obj, level) {
        key <- format(level)
        b   <- obj$bands[[key]]
        if (is.null(b))
            stop(sprintf(
                "Band %s%% is not available; refit with conf = c(%s).",
                key, paste(unique(c(obj$conf, level)), collapse = ", ")))
        b
    }

    grab <- function(obj, lab, sc) {
        N  <- dim(obj$IRmed)[1]
        H  <- dim(obj$IRmed)[3]
        if (shock < 1 || shock > N)
            stop(sprintf("`shock` must be between 1 and %d.", N))
        vn <- if (!is.null(varnames)) varnames else obj$varnames
        sc <- rep_len(sc, N)
        do.call(rbind, lapply(conf, function(level) {
            b <- band_of(obj, level)
            do.call(rbind, lapply(seq_len(N), function(i) tibble::tibble(
                model    = lab,
                level    = level,
                variable = vn[i],
                horizon  = 0L:(H - 1L),
                median   = as.numeric(obj$IRmed[i, shock, ]) * sc[i],
                lower    = as.numeric(b$IRinf[i, shock, ]) * sc[i],
                upper    = as.numeric(b$IRsup[i, shock, ]) * sc[i])))
        }))
    }

    dat <- grab(sr, labels[1], scale)
    if (!is.null(compare)) {
        if (!inherits(compare, "fSignRestr"))
            stop("`compare` must be an 'fSignRestr' object.")
        if (is.null(scale_compare)) scale_compare <- scale
        lab2 <- if (length(labels) >= 2) labels[2] else "Comparison"
        dat  <- rbind(dat, grab(compare, lab2, scale_compare))
    }

    vn <- if (!is.null(varnames)) varnames else sr$varnames
    dat$variable <- factor(dat$variable, levels = vn)
    dat$model    <- factor(dat$model, levels = unique(dat$model))

    if (return_data) return(tibble::as_tibble(dat))

    if (is.null(facet_ncol)) facet_ncol <- ceiling(sqrt(length(vn)))
    pal <- stats::setNames(colors[seq_along(levels(dat$model))], levels(dat$model))

    ## One ribbon per (model, level); the median is drawn once per model.
    dat$band <- interaction(dat$model, dat$level, drop = TRUE)
    med <- dat[dat$level == conf[1], , drop = FALSE]

    ## Widest band faintest, so nested levels stay readable.
    alphas <- if (length(conf) == 1L) ribbon_alpha
              else ribbon_alpha * seq(0.75, 1.45, length.out = length(conf))
    names(alphas) <- format(conf)

    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$horizon)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                         fill = .data$model, group = .data$band,
                         alpha = factor(format(.data$level), levels = format(conf))),
            colour = NA) +
        ggplot2::scale_alpha_manual(values = alphas, name = NULL,
                                    labels = paste0(format(conf), "%"),
                                    guide = if (length(conf) == 1L) "none"
                                            else ggplot2::guide_legend()) +
        ggplot2::geom_line(
            data = med,
            ggplot2::aes(y = .data$median, colour = .data$model),
            linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = 0, colour = "#707372",
                            linetype = "dashed", linewidth = 0.5) +
        ggplot2::scale_colour_manual(values = pal, name = NULL) +
        ggplot2::scale_fill_manual(values = pal, name = NULL) +
        ggplot2::facet_wrap(~ .data$variable, scales = facet_scales,
                            ncol = facet_ncol) +
        ggplot2::labs(x = "Horizon", y = NULL)

    if (is.null(compare)) p <- p + ggplot2::guides(colour = "none", fill = "none")
    p
}
