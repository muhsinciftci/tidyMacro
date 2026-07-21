# -----------------------------------------------------------------------
# fLPDID.R — LP-DiD (Dube-Girardi-Jorda-Taylor) with fixest-like syntax
#
# Provides:
#   fLPDID()      — main estimator (calls the C++ engine fLPDID_cpp)
#   fPlotLPDID()  — event-study plot; pass one result or several named
#                   results directly: fPlotLPDID(`Spec A` = a, `Spec B` = b)
#
# Formula syntax (mirrors fLPpanel):
#   fLPDID(y ~ treat + controls, data, panel_id = c("id", "year"),
#          treat = "treat", post = 9, pre = 9)
#
#   - LHS: single outcome variable.
#   - RHS: the treatment LEVEL column plus any controls. The term named in
#     `treat` is the treatment; everything else is a control.
#   - Time fixed effects are ALWAYS absorbed (they are part of the
#     estimator); unit effects are removed by the long difference. An
#     optional "| time" part is accepted but redundant.
#
# Macros (like fixest ..vars):
#   ctrl <- c("grgsp", "corptax", "unionmem")
#   fLPDID(lshare ~ bank + l(..ctrl, 1:4), ...)
#
# Panel-aware operators (respect unit boundaries, preserve time gaps):
#   l(x, k), l(x, k1:k2)  — lags        → columns "x_lk"
#   f(x, k)               — leads       → columns "x_fk"
#   d(x)                  — first diff  → column  "x_d1"
#   l(d(x), 1:4)          — lags of the first difference
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------


# =======================================================================
# Internal helpers
# =======================================================================

# Panel-aware shift: lag (k > 0) or lead (k < 0) within units, robust to
# time gaps. Returns a vector aligned with the original row order.
.flpdid_shift <- function(x, id, time, k) {
  if (k == 0L) return(x)
  ord <- order(id, time)
  xo  <- x[ord]; io <- id[ord]; to <- time[ord]
  n   <- length(x)
  res <- rep(NA_real_, n)
  a   <- abs(k)
  if (a < n) {
    if (k > 0L) {                       # lag
      idx <- (a + 1L):n
      ok  <- io[idx] == io[idx - a] & to[idx] == to[idx - a] + a
      res[idx[ok]] <- xo[(idx - a)[ok]]
    } else {                            # lead
      idx <- 1L:(n - a)
      ok  <- io[idx] == io[idx + a] & to[idx] == to[idx + a] - a
      res[idx[ok]] <- xo[(idx + a)[ok]]
    }
  }
  out <- rep(NA_real_, n)
  out[ord] <- res
  out
}

# Expand d(x) tokens: creates column "x_d1" = x - l(x, 1) and rewrites the
# formula string. Runs BEFORE the l()/f() expansion, so l(d(x), 1:4)
# becomes l(x_d1, 1:4) and then expands to x_d1_l1 ... x_d1_l4.
.fLPDID_expand_diff_terms <- function(formula_chr, data, id_var, time_var) {
  pattern <- "d\\(([A-Za-z_.][A-Za-z0-9_.]*)\\)"
  repeat {
    m <- regexpr(pattern, formula_chr, perl = TRUE)
    if (m == -1L) break
    full_match <- regmatches(formula_chr, m)
    var <- sub(pattern, "\\1", full_match, perl = TRUE)
    if (!var %in% names(data))
      stop(sprintf(
        "fLPDID: variable '%s' used in difference term '%s' not found in data.",
        var, full_match), call. = FALSE)
    col_name <- sprintf("%s_d1", var)
    if (!col_name %in% names(data)) {
      x <- as.numeric(data[[var]])
      data[[col_name]] <- x - .flpdid_shift(x, data[[id_var]],
                                            data[[time_var]], 1L)
    }
    formula_chr <- sub(pattern, col_name, formula_chr, perl = TRUE)
  }
  list(formula_chr = formula_chr, data = data)
}


# =======================================================================
# Main function
# =======================================================================

#' Local Projections Difference-in-Differences (fixest-like syntax)
#'
#' Estimates event-study treatment effects with the LP-DiD estimator of
#' Dube, Girardi, Jordà and Taylor (2025, \emph{Journal of Applied
#' Econometrics}, \doi{10.1002/jae.70000}): long-difference local
#' projections on clean-control samples with absorbed time fixed effects
#' and cluster-robust (CR1, Stata `reghdfe` convention) standard errors.
#' The horizon loop runs in parallel C++ (OpenMP). Formula syntax mirrors
#' \code{\link{fLPpanel}}: supports \code{..macro} expansion and the
#' panel-aware operators \code{l()}, \code{f()} and \code{d()}.
#'
#' @param formula A formula \code{y ~ treat + controls}, e.g.
#'   \code{lshare ~ bank + l(lshare, 1:4) + l(d(lshare), 1:4)}.
#'   \itemize{
#'     \item \strong{LHS}: single outcome variable.
#'     \item \strong{RHS}: the treatment \emph{level} column (0/1;
#'       \code{NA} allowed) plus any controls. The term named in
#'       \code{treat} is the treatment; all other terms are controls.
#'     \item \strong{Macros}: \code{..name} where \code{name} is a
#'       character vector in the calling environment.
#'     \item \strong{Operators}: \code{l(x, k)}, \code{l(x, k1:k2)},
#'       \code{f(x, k)}, \code{d(x)}, \code{l(d(x), 1:4)},
#'       \code{l(..vars, 1:4)}. Never cross unit boundaries; time gaps
#'       are preserved.
#'     \item Time fixed effects are always absorbed; an optional
#'       \code{| time} part is accepted but redundant. Unit effects are
#'       removed by the long difference.
#'   }
#' @param data A \code{data.frame} or \code{tibble} in long panel form.
#' @param panel_id Length-2 character vector \code{c(id_col, time_col)}.
#' @param treat Character. The RHS term label of the treatment level.
#' @param post,pre Number of post- and pre-treatment horizons
#'   (\code{event_time = -1} is the normalization unless \code{pmd}).
#' @param cluster Column name for clustering. Default: the unit id.
#' @param nonabsorbing Set \code{TRUE} when treatment can switch on and
#'   off (entries and reversals). Requires \code{L}.
#' @param L Stabilization window for the clean-control sets
#'   (non-absorbing case).
#' @param ccc Clean-control condition flavor for the non-absorbing case:
#'   \code{0} = none (ANRR-style), \code{1} = clean past for treated and
#'   controls, \code{2} = controls additionally clean through \code{t+h}.
#' @param pmd Use the pre-mean-differenced baseline (average of all past
#'   outcomes) instead of \code{y[t-1]}.
#' @param reweight Recover the equally-weighted ATT via inverse implicit
#'   weights \code{1/(1-p_t)} (exogenous/absorbing case; with covariates
#'   this is an approximation — the regression-adjustment version is not
#'   yet implemented in C++).
#' @param conf Confidence level for the reported bands (default 0.95).
#' @param n_threads OpenMP threads (\code{-1} = all available).
#'
#' @return A tibble of class \code{"fLPDID"} with columns
#'   \code{event_time}, \code{estimate}, \code{se}, \code{conf_low},
#'   \code{conf_high}, \code{nobs}, \code{nclust}.
#'
#' @examples
#' \dontrun{
#' data("BankingDeregulation")
#' m <- fLPDID(lshare ~ bank, data = BankingDeregulation,
#'             panel_id = c("id", "year"), treat = "bank",
#'             post = 9, pre = 9)
#' fPlotLPDID(m)
#'
#' # with controls via panel-aware operators and macros
#' ctrl <- c("grgsp", "corptax", "unionmem")
#' m2 <- fLPDID(lshare ~ bank + l(lshare, 1:4) + l(d(lshare), 1:4) +
#'                l(..ctrl, 1:4),
#'              data = BankingDeregulation, panel_id = c("id", "year"),
#'              treat = "bank", post = 9, pre = 9)
#' }
#' @export
fLPDID <- function(formula, data,
                   panel_id,
                   treat,
                   post = 10L, pre = 5L,
                   cluster = NULL,
                   nonabsorbing = FALSE, L = NULL, ccc = 1L,
                   pmd = FALSE, reweight = FALSE,
                   conf = 0.95, n_threads = -1L) {

  env <- parent.frame()

  # ---- required arguments ----------------------------------------------
  if (missing(treat))
    stop("fLPDID: 'treat' is required (the RHS term of the treatment level).",
         call. = FALSE)
  if (missing(panel_id))
    stop("fLPDID: 'panel_id' is required (c(id_col, time_col)).",
         call. = FALSE)
  if (!is.character(panel_id) || length(panel_id) != 2L || anyNA(panel_id))
    stop("fLPDID: 'panel_id' must be a length-2 character vector.",
         call. = FALSE)
  id_var   <- panel_id[1L]
  time_var <- panel_id[2L]
  if (!is.character(treat) || length(treat) != 1L)
    stop("fLPDID: 'treat' must be a single character string.", call. = FALSE)
  if (!is.data.frame(data)) data <- as.data.frame(data)
  for (v in c(id_var, time_var))
    if (!v %in% names(data))
      stop(sprintf("fLPDID: column '%s' not found in data.", v),
           call. = FALSE)
  if (nonabsorbing && ccc > 0 && is.null(L))
    stop("fLPDID: `nonabsorbing = TRUE` requires the stabilization window `L`.",
         call. = FALSE)
  if (!pmd && pre < 2L)
    stop("fLPDID: `pre` must be >= 2 (event_time -1 is the normalization).",
         call. = FALSE)
  if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1)
    stop("fLPDID: 'conf' must be a single number in (0, 1).", call. = FALSE)

  # ---- normalized integer time index (preserves gaps, as in fLPpanel) --
  t_raw <- data[[time_var]]
  if (inherits(t_raw, "Date") || inherits(t_raw, "POSIXct"))
    t_raw <- as.numeric(t_raw)
  if (!is.numeric(t_raw))
    stop("fLPDID: 'time' column must be numeric, Date, or POSIXct.",
         call. = FALSE)
  steps <- diff(sort(unique(t_raw)))
  step  <- if (length(steps)) min(steps) else 1
  if (!is.finite(step) || step <= 0)
    stop("fLPDID: could not determine the time step.", call. = FALSE)
  t_index <- as.integer(round((t_raw - min(t_raw, na.rm = TRUE)) / step))

  # ---- formula -> string; macro, d(), l()/f() expansion ----------------
  formula_chr <- paste(deparse(formula), collapse = " ")
  formula_chr <- .fLP_expand_macro_lag(formula_chr, env)
  formula_chr <- .fLP_expand_macro(formula_chr, env)

  # optional "| fe" part: accepted only if it is the time variable
  parts    <- strsplit(formula_chr, "\\|")[[1]]
  main_str <- trimws(parts[1L])
  if (length(parts) > 1L) {
    fe_vars <- attr(stats::terms(
      stats::as.formula(paste("~", trimws(parts[2L])), env = env)),
      "term.labels")
    if (!all(fe_vars %in% time_var))
      stop("fLPDID: only time fixed effects are meaningful in LP-DiD ",
           "(they are always absorbed); unit effects are removed by the ",
           "long difference. Drop '| ", paste(fe_vars, collapse = " + "),
           "' from the formula.", call. = FALSE)
  }

  diff_ex  <- .fLPDID_expand_diff_terms(main_str, data, id_var, time_var)
  main_ex  <- .fLPpanel_expand_lag_terms(diff_ex$formula_chr, diff_ex$data,
                                         id_var, time_var)
  main_str <- main_ex$formula_chr
  data     <- main_ex$data

  main_formula <- tryCatch(
    stats::as.formula(main_str, env = env),
    error = function(e) stop("fLPDID: could not parse the formula:\n  ",
                             main_str, "\n  ", conditionMessage(e),
                             call. = FALSE))
  lhs_vars <- all.vars(main_formula[[2L]])
  if (length(lhs_vars) != 1L)
    stop("fLPDID: exactly one LHS variable is required.", call. = FALSE)

  rhs_str    <- paste(deparse(main_formula[[3L]]), collapse = " ")
  rhs_terms  <- stats::terms(stats::as.formula(paste("~", rhs_str), env = env),
                             keep.order = TRUE)
  rhs_labels <- attr(rhs_terms, "term.labels")

  if (!treat %in% rhs_labels)
    stop(sprintf(
      "fLPDID: treatment term '%s' not found in RHS. Available: %s",
      treat, paste(rhs_labels, collapse = ", ")), call. = FALSE)
  if (!treat %in% names(data))
    stop("fLPDID: 'treat' must be a bare column (the treatment level), ",
         "not an interaction or transformed term.", call. = FALSE)
  control_labels <- setdiff(rhs_labels, treat)

  # ---- numeric checks --------------------------------------------------
  base_vars <- unique(c(lhs_vars, all.vars(rhs_terms)))
  for (nm in base_vars)
    if (!is.numeric(data[[nm]]))
      stop(sprintf("fLPDID: column '%s' must be numeric.", nm),
           call. = FALSE)

  # ---- control design matrix (NA rows preserved; C++ handles them) -----
  if (length(control_labels)) {
    design_formula <- stats::as.formula(
      paste("~", paste(control_labels, collapse = " + "), "- 1"), env = env)
    mf <- stats::model.frame(design_formula, data = data,
                             na.action = stats::na.pass)
    X  <- stats::model.matrix(design_formula, mf)
    if (ncol(X) != length(control_labels))
      stop("fLPDID: internal error — design matrix has ", ncol(X),
           " columns but ", length(control_labels), " control terms.",
           call. = FALSE)
    colnames(X) <- control_labels
    storage.mode(X) <- "double"
  } else {
    X <- matrix(0.0, nrow(data), 0L)
  }

  cl_var <- if (is.null(cluster)) id_var else cluster
  if (!cl_var %in% names(data))
    stop(sprintf("fLPDID: cluster column '%s' not found in data.", cl_var),
         call. = FALSE)

  # ---- C++ engine -------------------------------------------------------
  out <- fLPDID_cpp(
    y            = as.numeric(data[[lhs_vars]]),
    treat        = as.numeric(data[[treat]]),
    X            = X,
    i_index      = as.integer(as.factor(data[[id_var]])),
    t_index      = t_index,
    cl_index     = as.integer(as.factor(data[[cl_var]])),
    pre_window   = as.integer(pre),
    post_window  = as.integer(post),
    nonabsorbing = isTRUE(nonabsorbing),
    Lwin         = as.integer(if (is.null(L)) 0L else L),
    ccc          = as.integer(ccc),
    pmd          = isTRUE(pmd),
    reweight     = isTRUE(reweight),
    n_threads    = as.integer(n_threads)
  )

  res <- tibble::tibble(
    event_time = as.integer(out$event_time),
    estimate   = as.numeric(out$estimate),
    se         = as.numeric(out$se),
    nobs       = as.integer(out$nobs),
    nclust     = as.integer(out$nclust)
  )
  if (!pmd)
    res <- dplyr::bind_rows(
      res,
      tibble::tibble(event_time = -1L, estimate = 0, se = 0,
                     nobs = NA_integer_, nclust = NA_integer_))

  z   <- stats::qnorm(1 - (1 - conf) / 2)
  res <- dplyr::arrange(res, .data$event_time)
  res <- dplyr::mutate(res,
                       conf_low  = .data$estimate - z * .data$se,
                       conf_high = .data$estimate + z * .data$se)
  res <- dplyr::select(res, "event_time", "estimate", "se",
                       "conf_low", "conf_high", "nobs", "nclust")

  attr(res, "lpdid_spec") <- list(
    formula = formula_chr, y = lhs_vars, treat = treat,
    controls = control_labels, post = post, pre = pre,
    nonabsorbing = nonabsorbing, L = L, ccc = ccc,
    pmd = pmd, reweight = reweight, conf = conf, cluster = cl_var
  )
  class(res) <- c("fLPDID", class(res))
  res
}


# =======================================================================
# Plotting
# =======================================================================

#' Plot LP-DiD event-study estimates
#'
#' Pass a single \code{fLPDID} result, or several named results directly —
#' one facet per result:
#' \code{fPlotLPDID(`Only FEs` = a, `With controls` = b)}.
#'
#' Horizons that could not be estimated (non-finite estimates) are dropped
#' from the plot with a message, so no ggplot warnings are emitted.
#'
#' @param ... One or more \code{fLPDID} results (ideally named).
#' @param title,xlab,ylab Plot labels.
#' @param ncol Number of facet columns when several results are given.
#' @param scales Facet scales: \code{"fixed"}, \code{"free"},
#'   \code{"free_x"} or \code{"free_y"} (default).
#' @param ci How to draw the confidence intervals: \code{"ribbon"}
#'   (shaded band, DGJT Figure-4 style; default) or \code{"bars"}
#'   (point estimates with capped error bars, coefplot / Figure-3 style).
#' @return A ggplot object.
#' @export
fPlotLPDID <- function(..., title = NULL,
                       xlab = "Event time", ylab = "Estimate",
                       ncol = 2L, scales = "free_y",
                       ci = c("ribbon", "bars")) {
  ci <- match.arg(ci)
  if (!is.character(scales) || length(scales) != 1L ||
      !scales %in% c("fixed", "free", "free_x", "free_y"))
    stop("fPlotLPDID: 'scales' must be one of \"fixed\", \"free\", ",
         "\"free_x\", \"free_y\".", call. = FALSE)
  specs <- list(...)
  # backward compatibility: a single plain list of results is unwrapped
  if (length(specs) == 1L && !inherits(specs[[1L]], "fLPDID") &&
      is.list(specs[[1L]]))
    specs <- specs[[1L]]
  if (!length(specs))
    stop("Provide at least one fLPDID result.", call. = FALSE)
  ok <- vapply(specs, inherits, logical(1L), what = "fLPDID")
  if (!all(ok))
    stop("All inputs must be fLPDID results.", call. = FALSE)
  nm <- names(specs)
  if (is.null(nm)) nm <- rep("", length(specs))
  nm[nm == ""] <- paste("Spec", which(nm == ""))
  if (length(specs) == 1L && nm[1L] == "Spec 1") nm[1L] <- "LP-DiD"

  df <- dplyr::bind_rows(
    lapply(seq_along(specs), function(i) {
      d <- as.data.frame(specs[[i]])
      d$spec <- nm[i]
      d
    })
  )
  df$spec <- factor(df$spec, levels = nm)

  bad <- !is.finite(df$estimate)
  if (any(bad)) {
    message("fPlotLPDID: dropping ", sum(bad),
            " horizon(s) that could not be estimated.")
    df <- df[!bad, ]
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$event_time,
                                        y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55",
                        linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = -0.5, colour = "grey80",
                        linetype = 2, linewidth = 0.3)

  if (ci == "ribbon") {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf_low,
                                        ymax = .data$conf_high),
                           alpha = 0.15, fill = "#2C6E91") +
      ggplot2::geom_line(colour = "#2C6E91", linewidth = 0.6) +
      ggplot2::geom_point(colour = "#2C6E91", size = 1)
  } else {
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf_low,
                                          ymax = .data$conf_high),
                             width = 0.35, colour = "#2C6E91",
                             linewidth = 0.4) +
      ggplot2::geom_point(colour = "#2C6E91", size = 1.4)
  }

  p <- p + ggplot2::labs(title = title, x = xlab, y = ylab)
  if (length(specs) > 1L)
    p <- p + ggplot2::facet_wrap(~spec, ncol = ncol, scales = scales)
  p
}
