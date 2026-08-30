# -----------------------------------------------------------------------
# fLPDID.R — LP-DiD (Dube-Girardi-Jorda-Taylor) with fixest-like syntax
#
# Provides:
#   fLPDID()      — main estimator (calls the C++ engine fLPDID_cpp)
#   fPlotLPDID()  — event-study plot; pass one result or several named
#                   results directly: fPlotLPDID(`Spec A` = a, `Spec B` = b)
#
# Formula syntax (mirrors fLPPanel):
#   fLPDID(y ~ treat + controls, data, panel_id = c("id", "year"),
#          treat = "treat", post = 9, pre = 9)
#
#   - LHS: the outcome; transformations such as log(y) are evaluated.
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

# Panel-aware shift on the normalized integer time grid: lag (k > 0) or
# lead (k < 0) within units. `id` and `t_norm` are aligned with `x` and
# already share the same integer step, so this correctly handles Date /
# POSIXct / non-unit-step numeric time indexes.
.flpdid_shift <- function(x, id, t_norm, k) {
  if (k == 0L) return(x)
  key      <- paste(id, t_norm, sep = "|")
  key_look <- paste(id, t_norm - k, sep = "|")
  idx <- match(key_look, key)
  out <- rep(NA_real_, length(x))
  ok  <- !is.na(idx)
  out[ok] <- x[idx[ok]]
  out
}

# Expand d(x) tokens: creates column "x_d1" = x - l(x, 1) and rewrites the
# formula string. Runs BEFORE the l()/f() expansion so `l(d(x), 1:4)`
# becomes `l(x_d1, 1:4)` and then expands to x_d1_l1 ... x_d1_l4.
# `t_norm` is the normalized integer time; `id_vec` the id column.
.fLPDID_expand_diff_terms <- function(formula_chr, data, id_vec, t_norm) {
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
      data[[col_name]] <- x - .flpdid_shift(x, id_vec, t_norm, 1L)
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
#' \code{\link{fLPPanel}}: supports \code{..macro} expansion and the
#' panel-aware operators \code{l()}, \code{f()} and \code{d()}.
#'
#' @param formula A formula \code{y ~ treat + controls}, e.g.
#'   \code{lshare ~ bank + l(lshare, 1:4) + l(d(lshare), 1:4)}.
#'   \itemize{
#'     \item \strong{LHS}: the outcome. Transformations are evaluated,
#'       so \code{log(y) ~ treat} estimates the model for \code{log(y)}.
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
#'   \code{conf_high}, \code{nobs}, \code{nclust}, \code{ndrop}.
#'   Confidence bands use a cluster-t critical value with
#'   \code{nclust - 1} degrees of freedom (the \code{reghdfe}
#'   convention); \code{ndrop} counts control columns omitted for
#'   collinearity at that horizon.
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
  # Missing identifiers cannot be honoured: `as.integer(as.factor(NA))` is
  # NA_integer_, which the C++ hash maps would treat as one ordinary shared
  # unit — silently merging unrelated rows. Reject instead.
  if (anyNA(data[[id_var]]))
    stop(sprintf("fLPDID: %d missing value(s) in unit id '%s'; ",
                 sum(is.na(data[[id_var]])), id_var),
         "remove those rows before calling fLPDID.", call. = FALSE)
  # ---- scalar / range validation for numeric + logical args -----------
  for (nm_flag in c("nonabsorbing", "pmd", "reweight")) {
    v <- get(nm_flag)
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop(sprintf("fLPDID: '%s' must be a single TRUE/FALSE value.", nm_flag),
           call. = FALSE)
  }
  if (!is.numeric(pre)  || length(pre)  != 1L || is.na(pre)  ||
      pre  < 0 || pre  != floor(pre))
    stop("fLPDID: 'pre' must be a non-negative integer.", call. = FALSE)
  if (!is.numeric(post) || length(post) != 1L || is.na(post) ||
      post < 0 || post != floor(post))
    stop("fLPDID: 'post' must be a non-negative integer.", call. = FALSE)
  if (!is.numeric(ccc)  || length(ccc)  != 1L || is.na(ccc)  ||
      !ccc %in% c(0L, 1L, 2L))
    stop("fLPDID: 'ccc' must be one of 0, 1, or 2.", call. = FALSE)
  if (!is.null(L) && (!is.numeric(L) || length(L) != 1L || is.na(L) ||
                      L < 0 || L != floor(L)))
    stop("fLPDID: 'L' must be a non-negative integer or NULL.", call. = FALSE)
  if (!is.numeric(n_threads) || length(n_threads) != 1L || is.na(n_threads) ||
      n_threads != floor(n_threads))
    stop("fLPDID: 'n_threads' must be an integer (<=0 = all cores).",
         call. = FALSE)
  if (nonabsorbing && ccc > 0 && is.null(L))
    stop("fLPDID: `nonabsorbing = TRUE` requires the stabilization window `L`.",
         call. = FALSE)
  if (!pmd && pre < 2L)
    stop("fLPDID: `pre` must be >= 2 (event_time -1 is the normalization).",
         call. = FALSE)
  if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1)
    stop("fLPDID: 'conf' must be a single number in (0, 1).", call. = FALSE)

  # Reweight + controls is NOT the DDCG `teffects ra` regression-adjustment
  # estimator — surface that at runtime, not only in docs.
  if (isTRUE(reweight)) {
    # controls existence is only known after formula parsing; emit a
    # tentative warning here that becomes accurate once controls are
    # parsed (we recheck below and downgrade the message if unneeded).
    reweight_needs_note <- TRUE
  } else {
    reweight_needs_note <- FALSE
  }

  # ---- normalized integer time index (preserves gaps, as in fLPPanel) --
  t_raw <- data[[time_var]]
  if (inherits(t_raw, "Date") || inherits(t_raw, "POSIXct"))
    t_raw <- as.numeric(t_raw)
  if (!is.numeric(t_raw))
    stop("fLPDID: 'time' column must be numeric, Date, or POSIXct.",
         call. = FALSE)
  # A non-finite time would become NA_integer_ in `t_index` and reach the C++
  # engine as an ordinary (garbage) period, silently corrupting every shift.
  if (!all(is.finite(t_raw)))
    stop(sprintf("fLPDID: %d missing or non-finite value(s) in time column '%s'; ",
                 sum(!is.finite(t_raw)), time_var),
         "remove those rows before calling fLPDID.", call. = FALSE)
  steps <- diff(sort(unique(t_raw)))
  step  <- if (length(steps)) min(steps) else 1
  if (!is.finite(step) || step <= 0)
    stop("fLPDID: could not determine the time step.", call. = FALSE)
  t_index <- as.integer(round((t_raw - min(t_raw, na.rm = TRUE)) / step))

  # Duplicate (id, time) rows silently corrupt the panel time->row map
  # used by the C++ engine — reject them up-front.
  dup <- duplicated(paste(data[[id_var]], t_index, sep = "|"))
  if (any(dup))
    stop(sprintf(
      "fLPDID: found %d duplicate (%s, %s) row(s); each unit-time pair must be unique.",
      sum(dup), id_var, time_var), call. = FALSE)

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

  diff_ex  <- .fLPDID_expand_diff_terms(main_str, data,
                                        data[[id_var]], t_index)
  main_ex  <- .fLPPanel_expand_lag_terms(diff_ex$formula_chr, diff_ex$data,
                                         diff_ex$data[[id_var]], t_index)
  main_str <- main_ex$formula_chr
  data     <- main_ex$data

  main_formula <- tryCatch(
    stats::as.formula(main_str, env = env),
    error = function(e) stop("fLPDID: could not parse the formula:\n  ",
                             main_str, "\n  ", conditionMessage(e),
                             call. = FALSE))
  # The LHS is an expression, not just a column name: RHS controls already
  # go through model.matrix() and honour transformations, so `log(y) ~ treat`
  # must be evaluated too rather than silently collapsing to `y`.
  lhs_expr  <- main_formula[[2L]]
  lhs_vars  <- all.vars(lhs_expr)
  lhs_label <- paste(deparse(lhs_expr), collapse = "")
  if (!length(lhs_vars))
    stop("fLPDID: the LHS must reference at least one column of `data`.",
         call. = FALSE)
  for (nm in lhs_vars)
    if (!nm %in% names(data))
      stop(sprintf("fLPDID: LHS variable '%s' not found in data.", nm),
           call. = FALSE)

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

  # The engine defines an entry as dD == 1 and a candidate control as
  # dD == 0 with treatment level 0. Any other coding (e.g. 0/2) selects no
  # rows at all and would return NaN for every horizon without explanation.
  tr_obs <- data[[treat]][!is.na(data[[treat]])]
  if (length(tr_obs) && !all(tr_obs %in% c(0, 1))) {
    offending <- unique(tr_obs[!tr_obs %in% c(0, 1)])
    stop(sprintf(
      "fLPDID: treatment '%s' must be coded 0/1 (NA allowed); found %s.",
      treat, paste(utils::head(sort(offending), 5L), collapse = ", ")),
      call. = FALSE)
  }
  if (!any(tr_obs == 1))
    warning(sprintf(
      "fLPDID: treatment '%s' is never 1 — no treatment entries exist, so ",
      treat), "every horizon will be non-estimable.", call. = FALSE)

  # LHS evaluated on the (operator-expanded) data.
  y_vec <- tryCatch(eval(lhs_expr, data, env),
                    error = function(e)
                      stop("fLPDID: could not evaluate the LHS '", lhs_label,
                           "':\n  ", conditionMessage(e), call. = FALSE))
  if (!is.numeric(y_vec) || length(y_vec) != nrow(data))
    stop("fLPDID: the LHS '", lhs_label, "' must evaluate to a numeric ",
         "vector of length nrow(data).", call. = FALSE)

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
  # As for the unit id: NA cluster codes would be pooled into one cluster,
  # changing G and the CR1 factor without any diagnostic.
  if (anyNA(data[[cl_var]]))
    stop(sprintf("fLPDID: %d missing value(s) in cluster column '%s'; ",
                 sum(is.na(data[[cl_var]])), cl_var),
         "remove those rows before calling fLPDID.", call. = FALSE)

  # Finalize the reweight+controls note now that controls are parsed.
  if (reweight_needs_note && length(control_labels) > 0L) {
    warning(
      "fLPDID: reweight = TRUE with controls uses simple inverse-implicit ",
      "reweighting, NOT the DDCG `teffects ra` regression-adjustment ",
      "estimator (see Dube-Girardi-Jorda-Taylor Figure 4). Interpret the ",
      "reweighted ATT accordingly.",
      call. = FALSE)
  }

  # ---- C++ engine -------------------------------------------------------
  out <- fLPDID_cpp(
    y            = as.numeric(y_vec),
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
    n_threads    = as.integer(n_threads),
    verbose      = FALSE
  )

  res <- tibble::tibble(
    event_time = as.integer(out$event_time),
    estimate   = as.numeric(out$estimate),
    se         = as.numeric(out$se),
    nobs       = as.integer(out$nobs),
    nclust     = as.integer(out$nclust),
    ndrop      = as.integer(out$ndrop),
    .norm      = FALSE
  )
  if (!pmd)
    res <- dplyr::bind_rows(
      res,
      tibble::tibble(event_time = -1L, estimate = 0, se = 0,
                     nobs = NA_integer_, nclust = NA_integer_,
                     ndrop = NA_integer_, .norm = TRUE))

  res <- dplyr::arrange(res, .data$event_time)

  # Cluster-robust inference uses a t critical value with G - 1 degrees of
  # freedom (the `reghdfe` / `fixest` convention that matches the CR1 factor
  # applied in C++), not a normal quantile.
  crit <- rep(NA_real_, nrow(res))
  ok   <- !is.na(res$nclust) & res$nclust > 1L
  crit[ok] <- stats::qt(1 - (1 - conf) / 2, res$nclust[ok] - 1L)
  res$conf_low  <- res$estimate - crit * res$se
  res$conf_high <- res$estimate + crit * res$se
  # The normalized event_time == -1 point is exactly zero by construction.
  res$conf_low[res$.norm]  <- 0
  res$conf_high[res$.norm] <- 0

  bad <- !res$.norm & !is.finite(res$estimate)
  if (any(bad))
    warning(sprintf(
      "fLPDID: %d horizon(s) not estimable (event_time %s) — too few clean-control rows, no treatment entries, or a degenerate treatment column after within-year demeaning.",
      sum(bad), paste(res$event_time[bad], collapse = ", ")), call. = FALSE)
  if (any(res$ndrop > 0L, na.rm = TRUE))
    warning(sprintf(
      "fLPDID: collinear control column(s) omitted at event_time %s (as `reghdfe` does); see the `ndrop` column.",
      paste(res$event_time[!is.na(res$ndrop) & res$ndrop > 0L],
            collapse = ", ")), call. = FALSE)

  res <- dplyr::select(res, "event_time", "estimate", "se",
                       "conf_low", "conf_high", "nobs", "nclust", "ndrop")

  attr(res, "lpdid_spec") <- list(
    formula = formula_chr, y = lhs_label, treat = treat,
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
