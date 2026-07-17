# -----------------------------------------------------------------------
# fLPIV.R — Local Projections with an External Instrument (LP-IV)
#
# Syntax mirrors fLP():
#   fLPIV(
#     formula     = y ~ D + l(..controls, 1:6),
#     instruments = ~ l(z, 0:6),
#     data        = df,
#     endog       = "D",
#     horizons    = 48,
#     conf        = c(0.68, 0.95),
#     nw_lags_iv  = 6,          # fixed HAC bandwidth (matches Stata vce(hac nw 6))
#     cumulative  = TRUE
#   )
#
# 'formula' is the main outcome equation: LHS = outcome(s), RHS =
# endogenous treatment + exogenous controls (with l()/f()/..macro support).
# 'instruments' is a one-sided formula whose RHS is the full external
# instrument set (again with l()/f()/..macro support).
#
# 'endog' names the single RHS term that is the endogenous treatment. The
# remaining RHS terms become the exogenous control matrix C.
# -----------------------------------------------------------------------

#' Local Projections with an External Instrument (LP-IV)
#'
#' Estimates 2SLS local projections horizon-by-horizon using the delta-method
#' Newey-West HAC sandwich for the endogenous treatment effect. Matches
#' Cesa-Bianchi's VAR Toolbox \code{LPmodel.m} in LP-IV mode when
#' \code{nw_lags_iv} is set to the toolbox's \code{nlag_iv}.
#'
#' @param formula A formula \code{y ~ D + controls}. LHS may be
#'   \code{c(y1, y2, ...)} for multi-outcome LP-IV. RHS lists the endogenous
#'   treatment plus all exogenous controls. \code{l()} / \code{f()} lag/lead
#'   operators and \code{..macro} expansion behave exactly as in \code{fLP()}.
#' @param instruments A one-sided formula \code{~ z + ...} whose RHS is the
#'   full external instrument matrix. Use \code{l(z, 0:6)} to include the
#'   contemporaneous instrument plus six lags (over-identified 2SLS).
#' @param data A \code{data.frame} or \code{tibble}.
#' @param endog Character. The single RHS term (of \code{formula}) that is
#'   the endogenous treatment. The remaining RHS terms are the exogenous
#'   control matrix.
#' @param horizons Integer or the sequence \code{0:H}. Default 12.
#' @param conf Numeric in (0, 1). Scalar or vector of confidence levels.
#'   Vectors trigger multi-band output; see \code{\link{fLP}}.
#' @param nw_lags_iv Integer. Bartlett HAC bandwidth for the IV score
#'   variance. \code{> 0}: fixed bandwidth across horizons — this matches
#'   Stata's \code{vce(hac nw <k>)} and the toolbox's LP-IV bands.
#'   \code{= 0}: horizon-varying bandwidth \code{= h} (just-identified rule).
#'   Default \code{NULL} → use \code{ncol(Z) - 1}, i.e. the lag order implied
#'   by the instrument formula, which reproduces the toolbox default.
#' @param cumulative Logical. \code{TRUE} regresses
#'   \eqn{y_{t+h} - y_{t-1}} on the RHS. Default \code{TRUE} — the standard
#'   long-difference LP-IV convention.
#' @param n_threads Integer. OpenMP threads for the horizon loop. Default 0
#'   (= all cores minus one).
#'
#' @return An object of class \code{c("fLPIV", "fLP")} so all existing
#'   \code{fLP}-family methods (\code{print}, \code{tidy.fLP}, \code{fPlotLP})
#'   work unchanged. Additional fields: \code{Fstat_fs} and \code{rsqr_fs}
#'   (both length \code{H + 1}: per-horizon first-stage F-statistic and R²),
#'   plus \code{endog} and \code{instrument_vars}.
#'
#' @examples
#' \dontrun{
#' controls <- c("Unemployment", "Inflation", "FFRates")
#' fLPIV(
#'   Unemployment ~ FFRates + l(..controls, 1:6),
#'   instruments = ~ l(RRShock, 0:6),
#'   data        = jt_iv,
#'   endog       = "FFRates",
#'   horizons    = 48,
#'   conf        = c(0.68, 0.95),
#'   nw_lags_iv  = 6,
#'   cumulative  = TRUE
#' )
#' }
#'
#' @export
fLPIV <- function(formula, instruments, data,
                  endog,
                  horizons   = 12L,
                  conf       = 0.90,
                  nw_lags_iv = NULL,
                  cumulative = TRUE,
                  n_threads  = 0L) {

  env <- parent.frame()

  # -------------------------------------------------------------------
  # 0. Basic checks
  # -------------------------------------------------------------------
  if (missing(endog)) {
    stop("fLPIV: 'endog' is required (name of the endogenous treatment RHS term).")
  }
  if (!is.character(endog) || length(endog) != 1L || is.na(endog)) {
    stop("fLPIV: 'endog' must be a single character string.")
  }
  if (missing(instruments) || !inherits(instruments, "formula")) {
    stop("fLPIV: 'instruments' must be a one-sided formula, e.g. ~ l(z, 0:6).")
  }

  H    <- .fLP_validate_horizon_max(horizons)
  conf <- .fLP_validate_conf_vector(conf)
  if (!is.data.frame(data)) data <- as.data.frame(data)

  # -------------------------------------------------------------------
  # 1. Parse 'formula' RHS/LHS with the fLP text-based expander
  # -------------------------------------------------------------------
  fmla_chr <- paste(deparse(formula), collapse = "")
  fmla_chr <- .fLP_expand_macro_lag(fmla_chr, env)
  fmla_chr <- .fLP_expand_macro(fmla_chr, env)
  exp_main <- .fLP_expand_lag_terms(fmla_chr, data)
  fmla_chr <- exp_main$formula_chr
  data     <- exp_main$data

  fmla_full <- tryCatch(as.formula(fmla_chr, env = env),
    error = function(e) stop("fLPIV: cannot parse expanded formula:\n  ", fmla_chr,
                             "\n  ", conditionMessage(e)))

  lhs_vars <- all.vars(fmla_full[[2L]])
  rhs_str  <- paste(deparse(fmla_full[[3L]]), collapse = " ")
  rhs_tms  <- terms(as.formula(paste("~", rhs_str), env = env), keep.order = TRUE)
  rhs_vars <- attr(rhs_tms, "term.labels")

  if (length(lhs_vars) == 0L) stop("fLPIV: formula must have at least one LHS variable.")
  if (length(rhs_vars) == 0L) stop("fLPIV: formula RHS must have at least the endogenous treatment.")

  # -------------------------------------------------------------------
  # 2. Parse 'instruments' (one-sided) with the same expander
  # -------------------------------------------------------------------
  iv_chr <- paste(deparse(instruments), collapse = "")
  iv_chr <- .fLP_expand_macro_lag(iv_chr, env)
  iv_chr <- .fLP_expand_macro(iv_chr, env)
  exp_iv <- .fLP_expand_lag_terms(iv_chr, data)
  iv_chr <- exp_iv$formula_chr
  data   <- exp_iv$data

  iv_full <- tryCatch(as.formula(iv_chr, env = env),
    error = function(e) stop("fLPIV: cannot parse expanded instruments formula:\n  ", iv_chr,
                             "\n  ", conditionMessage(e)))
  iv_rhs  <- if (length(iv_full) == 2L) iv_full[[2L]] else iv_full[[3L]]
  iv_str  <- paste(deparse(iv_rhs), collapse = " ")
  iv_tms  <- terms(as.formula(paste("~", iv_str), env = env), keep.order = TRUE)
  iv_vars <- attr(iv_tms, "term.labels")
  if (length(iv_vars) == 0L)
    stop("fLPIV: 'instruments' must contribute at least one instrument.")

  # -------------------------------------------------------------------
  # 3. Validate & split endog / exog controls
  # -------------------------------------------------------------------
  if (!endog %in% rhs_vars) {
    stop(sprintf("fLPIV: endogenous treatment '%s' not on the RHS. RHS terms: %s",
                 endog, paste(rhs_vars, collapse = ", ")))
  }
  ctrl_vars <- setdiff(rhs_vars, endog)  # can be empty

  # -------------------------------------------------------------------
  # 4. Check all needed vars exist; drop NA rows
  # -------------------------------------------------------------------
  all_needed <- unique(c(lhs_vars, endog, ctrl_vars, iv_vars))
  miss_v <- setdiff(all_needed, names(data))
  if (length(miss_v) > 0L)
    stop(sprintf("fLPIV: variable(s) not in data: %s", paste(miss_v, collapse = ", ")))

  raw <- data[, all_needed, drop = FALSE]
  na_rows <- which(!complete.cases(raw))
  if (length(na_rows) > 0L) {
    message(sprintf("fLPIV: removing %d row(s) with NA.", length(na_rows)))
    raw <- raw[-na_rows, , drop = FALSE]
  }
  non_num <- names(raw)[!vapply(raw, is.numeric, logical(1))]
  if (length(non_num) > 0L)
    stop(sprintf("fLPIV: all variables must be numeric. Non-numeric: %s",
                 paste(non_num, collapse = ", ")))

  Y <- as.matrix(raw[, lhs_vars, drop = FALSE]); storage.mode(Y) <- "double"
  D <- as.numeric(raw[[endog]])
  Z <- as.matrix(raw[, iv_vars, drop = FALSE]); storage.mode(Z) <- "double"
  C <- if (length(ctrl_vars) > 0L)
         as.matrix(raw[, ctrl_vars, drop = FALSE])
       else matrix(0.0, nrow = nrow(raw), ncol = 0)
  storage.mode(C) <- "double"

  T_eff <- nrow(Y)

  # -------------------------------------------------------------------
  # 5. Sanity / defaults
  # -------------------------------------------------------------------
  if (is.null(nw_lags_iv)) {
    nw_lags_iv <- max(0L, ncol(Z) - 1L)   # toolbox default
  } else {
    nw_lags_iv <- .fLP_validate_scalar_int(nw_lags_iv, "nw_lags_iv")
  }
  if (!is.logical(cumulative) || length(cumulative) != 1L || is.na(cumulative))
    stop("fLPIV: 'cumulative' must be TRUE or FALSE.")
  n_threads <- .fLP_validate_scalar_int(n_threads, "n_threads")

  min_obs <- ncol(C) + ncol(Z) + 3L + H + as.integer(cumulative)
  if (T_eff < min_obs) {
    stop(sprintf("fLPIV: not enough observations. Need >= %d rows after NA removal; have %d.",
                 min_obs, T_eff))
  }

  # -------------------------------------------------------------------
  # 6. Call C++ engine (with conf[1]; bands rebuilt in R from irfs_se)
  # -------------------------------------------------------------------
  res_cpp <- fLPIV_cpp(
    Y          = Y,
    D          = D,
    Z          = Z,
    C          = C,
    H          = H,
    conf_level = as.double(conf[1L]),
    nw_lags_iv = as.integer(nw_lags_iv),
    cumulative = isTRUE(cumulative),
    n_threads  = as.integer(n_threads)
  )

  # -------------------------------------------------------------------
  # 7. Rebuild multi-band bands (matching fLP behavior)
  # -------------------------------------------------------------------
  .build_band <- function(cl) {
    z <- stats::qnorm(0.5 * (1.0 + cl))
    list(upper = res_cpp$irfs + z * res_cpp$irfs_se,
         lower = res_cpp$irfs - z * res_cpp$irfs_se)
  }
  if (length(conf) == 1L) {
    b <- .build_band(conf)
    res_cpp$irfs_upper <- b$upper
    res_cpp$irfs_lower <- b$lower
  } else {
    keys <- format(conf, trim = TRUE)
    upper_list <- vector("list", length(conf)); names(upper_list) <- keys
    lower_list <- vector("list", length(conf)); names(lower_list) <- keys
    for (i in seq_along(conf)) {
      b <- .build_band(conf[i])
      upper_list[[i]] <- b$upper
      lower_list[[i]] <- b$lower
    }
    res_cpp$irfs_upper <- upper_list
    res_cpp$irfs_lower <- lower_list
  }

  # -------------------------------------------------------------------
  # 8. Labeling
  # -------------------------------------------------------------------
  h_seq <- 0L:H
  rownames(res_cpp$irfs)    <- as.character(h_seq); colnames(res_cpp$irfs)    <- lhs_vars
  rownames(res_cpp$irfs_se) <- as.character(h_seq); colnames(res_cpp$irfs_se) <- lhs_vars
  names(res_cpp$Fstat_fs)   <- as.character(h_seq)
  names(res_cpp$rsqr_fs)    <- as.character(h_seq)

  .label_band <- function(m) {
    rownames(m) <- as.character(h_seq); colnames(m) <- lhs_vars; m
  }
  if (is.list(res_cpp$irfs_upper)) {
    res_cpp$irfs_upper <- lapply(res_cpp$irfs_upper, .label_band)
    res_cpp$irfs_lower <- lapply(res_cpp$irfs_lower, .label_band)
  } else {
    res_cpp$irfs_upper <- .label_band(res_cpp$irfs_upper)
    res_cpp$irfs_lower <- .label_band(res_cpp$irfs_lower)
  }

  res_cpp$lhs_vars        <- lhs_vars
  res_cpp$rhs_vars        <- rhs_vars
  res_cpp$endog           <- endog
  res_cpp$ctrl_vars       <- ctrl_vars
  res_cpp$instrument_vars <- iv_vars
  res_cpp$shock           <- endog                # for fPlotLP labelling
  res_cpp$horizons        <- h_seq
  res_cpp$conf            <- conf
  res_cpp$nw_lags_iv      <- nw_lags_iv
  res_cpp$cumulative      <- isTRUE(cumulative)
  res_cpp$n_threads       <- as.integer(n_threads)
  res_cpp$nobs            <- T_eff
  res_cpp$formula         <- fmla_full
  res_cpp$instruments     <- iv_full
  res_cpp$formula_orig    <- formula
  res_cpp$call            <- match.call()

  class(res_cpp) <- c("fLPIV", "fLP")
  res_cpp
}


# =======================================================================
# S3 methods
# =======================================================================

#' @export
print.fLPIV <- function(x, digits = 4L, ...) {
  cat("\nLocal Projections — IV (fLPIV)\n")
  cat(strrep("-", 45L), "\n")
  cat("Original formula : ", deparse(x$formula_orig), "\n")
  cat("Expanded formula : ", deparse(x$formula),      "\n")
  cat("Instruments      : ", deparse(x$instruments),  "\n")
  cat("Endogenous       : ", x$endog,                 "\n")
  cat("Controls         : ", if (length(x$ctrl_vars)) paste(x$ctrl_vars, collapse = ", ") else "(none)", "\n")
  cat("LHS variables    : ", paste(x$lhs_vars, collapse = ", "), "\n")
  cat("Horizons         : 0 to", max(x$horizons),     "\n")
  cat("Cumulative       : ", isTRUE(x$cumulative),     "\n")
  cat("Confidence       : ",
      paste0(format(x$conf * 100, trim = TRUE), "%", collapse = ", "), "\n")
  cat("HAC bandwidth    :  ")
  if (isTRUE(x$nw_lags_iv > 0)) {
    cat(sprintf("fixed = %d (vce(hac nw %d))\n", x$nw_lags_iv, x$nw_lags_iv))
  } else {
    cat("horizon-varying = h (just-identified rule)\n")
  }
  cat("Observations     : ", x$nobs, "\n")
  cat("\nIRF (endog = '", x$endog, "'):\n", sep = "")
  print(round(x$irfs, digits))
  if (!is.null(x$Fstat_fs)) {
    cat("\nFirst-stage F (h = 0..H): ",
        paste(round(x$Fstat_fs, 2), collapse = ", "), "\n")
  }
  invisible(x)
}
