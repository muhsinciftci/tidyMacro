# -----------------------------------------------------------------------
# fLP.R — Classical Local Projections with fixest-like formula syntax
#
# Provides:
#   fLP()        — main user-facing function
#   print.fLP()  — print method
#   tidy.fLP()   — broom-compatible tidy method
#   coef.fLP()   — extract IRF coefficients
#
# Syntax mirrors fixest::feols():
#   fLP(y ~ shock + control1 + control2, data, horizons = 12, shock = "shock")
#
# Macro expansion (like fixest's ..vars):
#   controls <- c("ip_l1", "cpi_l1")
#   fLP(ip ~ mp_shock + ..controls, data, horizons = 12, shock = "mp_shock")
#
# Lag / lead operators:
#   l(var, k)     — lag var by k periods  → column "var_lk"
#   l(var, k1:k2) — lags k1 through k2    → columns "var_lk1", ..., "var_lk2"
#   f(var, k)     — lead var by k periods  → column "var_fk"
#
#   Works inside macros too:
#   l(..controls, 1:3) — lags all variables in 'controls' by 1, 2, 3
#
# Multi-equation:
#   fLP(c(ip, cpi, rate) ~ mp_shock + ..controls, ..., shock = "mp_shock")
#   (runs one LP per LHS variable; LHS vars are NOT added to the RHS)
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------


# =======================================================================
# Internal helpers
# =======================================================================

.fLP_validate_horizon_max <- function(horizons) {
  if (!is.numeric(horizons) || length(horizons) == 0L ||
      anyNA(horizons) || any(!is.finite(horizons)) ||
      any(horizons < 0) || any(horizons != floor(horizons))) {
    stop("fLP: 'horizons' must be a non-negative integer or the sequence 0:H.")
  }

  h <- as.integer(horizons)
  if (length(h) == 1L) {
    return(h)
  }

  max_h <- max(h)
  if (!identical(sort(unique(h)), 0L:max_h)) {
    stop("fLP: vector 'horizons' must be the complete sequence 0:H.")
  }
  max_h
}

.fLP_validate_scalar_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0 || x != floor(x)) {
    stop(sprintf("fLP: '%s' must be a non-negative integer.", name))
  }
  as.integer(x)
}

# Like .fLP_validate_scalar_int but allows negative values (e.g. nw_offset,
# which is legitimately -1 for the classic Jorda h-1 NW bandwidth rule).
.fLP_validate_scalar_int_signed <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x != floor(x)) {
    stop(sprintf("fLP: '%s' must be an integer.", name))
  }
  as.integer(x)
}

.fLP_validate_scalar_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= lower || x >= upper) {
    stop(sprintf("fLP: '%s' must be a scalar in (%s, %s).",
                 name, lower, upper))
  }
  as.double(x)
}

# Validates 'conf': one or more confidence levels in (0, 1), e.g. 0.90 or
# c(0.68, 0.95) for a double-band plot. Returned sorted descending (widest
# first) so downstream band-building and plotting order is deterministic.
.fLP_validate_conf_vector <- function(conf) {
  if (!is.numeric(conf) || length(conf) == 0L || anyNA(conf) ||
      any(!is.finite(conf)) || any(conf <= 0) || any(conf >= 1)) {
    stop("fLP: 'conf' must be one or more confidence levels in (0, 1), e.g. 0.90 or c(0.68, 0.95).")
  }
  conf <- sort(unique(as.double(conf)), decreasing = TRUE)
  conf
}

# -----------------------------------------------------------------------
# Step 1-helper: expand l(..macro, lags) / f(..macro, lags)
#   Called BEFORE general ..macro expansion so the inner macro token
#   is still intact and can be found in the parent environment.
#
#   l(..controls, 1:3)  →  l(ip, 1:3) + l(cpi, 1:3) + l(rate, 1:3)
# -----------------------------------------------------------------------
.fLP_expand_macro_lag <- function(formula_chr, env) {

  # Pattern: l or f, then (..macroname, anything)
  pattern <- "(l|f)\\(\\.\\.[A-Za-z_.][A-Za-z0-9_.]*,([^)]+)\\)"

  repeat {
    m <- regexpr(pattern, formula_chr, perl = TRUE)
    if (m == -1L) break

    full_match <- regmatches(formula_chr, m)
    parsed     <- regmatches(full_match,
                             regexec("(l|f)\\(\\.\\.([A-Za-z_.][A-Za-z0-9_.]*),([^)]+)\\)",
                                     full_match, perl = TRUE))[[1]]
    op       <- parsed[2L]            # "l" or "f"
    mac_name <- parsed[3L]            # macro name (without ..)
    lags_str <- trimws(parsed[4L])   # e.g. "1:3"

    if (!exists(mac_name, envir = env, inherits = TRUE)) {
      stop(sprintf("fLP: macro '..%s' not found in the calling environment.", mac_name))
    }
    vars <- get(mac_name, envir = env, inherits = TRUE)
    if (!is.character(vars) || length(vars) == 0L) {
      stop(sprintf("fLP: macro '..%s' must be a non-empty character vector.", mac_name))
    }

    # Build e.g. l(ip, 1:3) + l(cpi, 1:3) + l(rate, 1:3)
    expansion <- paste(
      sprintf("%s(%s, %s)", op, vars, lags_str),
      collapse = " + "
    )
    formula_chr <- sub(pattern, expansion, formula_chr, perl = TRUE)
  }

  formula_chr
}


# -----------------------------------------------------------------------
# Step 2-helper: expand plain ..macro tokens (not inside l/f calls)
# -----------------------------------------------------------------------
.fLP_expand_macro <- function(formula_chr, env) {

  pattern <- "\\.\\.[A-Za-z_.][A-Za-z0-9_.]*"
  matches <- unique(regmatches(formula_chr,
                               gregexpr(pattern, formula_chr))[[1]])

  for (macro in matches) {
    mac_name <- sub("^\\.\\.", "", macro)
    if (!exists(mac_name, envir = env, inherits = TRUE)) {
      stop(sprintf("fLP: macro '..%s' not found in the calling environment.", mac_name))
    }
    vars <- get(mac_name, envir = env, inherits = TRUE)
    if (!is.character(vars) || length(vars) == 0L) {
      stop(sprintf("fLP: macro '..%s' must be a non-empty character vector.", mac_name))
    }
    formula_chr <- gsub(macro,
                        paste(vars, collapse = " + "),
                        formula_chr, fixed = TRUE)
  }

  formula_chr
}


# -----------------------------------------------------------------------
# Step 3-helper: expand l(var, lags) / f(var, leads) term-by-term.
#   Creates new columns in (a copy of) data, returns updated formula
#   string and data.
#
#   Supported lag specs:
#     l(ip, 1)      → ip_l1
#     l(ip, 1:3)    → ip_l1 + ip_l2 + ip_l3
#     f(ip, 1:2)    → ip_f1 + ip_f2
# -----------------------------------------------------------------------
.fLP_expand_lag_terms <- function(formula_chr, data) {

  # Match l(...) or f(...) anywhere in the string.
  # We iterate to handle multiple occurrences.
  pattern <- "(l|f)\\(([A-Za-z_.][A-Za-z0-9_.]*),([^)]+)\\)"

  repeat {
    m <- regexpr(pattern, formula_chr, perl = TRUE)
    if (m == -1L) break

    full_match <- regmatches(formula_chr, m)
    parsed     <- regmatches(full_match,
                             regexec("(l|f)\\(([A-Za-z_.][A-Za-z0-9_.]*),([^)]+)\\)",
                                     full_match, perl = TRUE))[[1]]
    op       <- parsed[2L]            # "l" or "f"
    var      <- trimws(parsed[3L])   # variable name
    lags_str <- trimws(parsed[4L])   # lag spec, e.g. "1:3"

    # Evaluate lag indices safely
    lag_idx <- tryCatch(
      eval(parse(text = lags_str)),
      error = function(e) stop(sprintf(
        "fLP: cannot parse lag/lead spec '%s' in term '%s': %s",
        lags_str, full_match, conditionMessage(e)
      ))
    )
    if (!is.numeric(lag_idx) || anyNA(lag_idx) || any(!is.finite(lag_idx)) ||
        any(lag_idx < 0L) || any(lag_idx != floor(lag_idx))) {
      stop(sprintf(
        "fLP: lag/lead indices must be non-negative integers; got '%s'.", lags_str
      ))
    }
    lag_idx <- as.integer(lag_idx)

    if (!var %in% names(data)) {
      stop(sprintf(
        "fLP: variable '%s' used in lag/lead term '%s' not found in data.",
        var, full_match
      ))
    }

    x_raw <- data[[var]]
    T_dat <- length(x_raw)

    new_cols <- character(length(lag_idx))
    for (i in seq_along(lag_idx)) {
      k        <- lag_idx[i]
      col_name <- if (op == "l") sprintf("%s_l%d", var, k) else
                                 sprintf("%s_f%d", var, k)

      # Create the column only if it doesn't already exist
      if (!col_name %in% names(data)) {
        if (k == 0L) {
          data[[col_name]] <- x_raw
        } else if (k >= T_dat) {
          data[[col_name]] <- rep(NA_real_, T_dat)
        } else if (op == "l") {
          # lag: shift forward k rows, fill head with NA
          data[[col_name]] <- c(rep(NA_real_, k), x_raw[seq_len(T_dat - k)])
        } else {
          # lead: shift backward k rows, fill tail with NA
          data[[col_name]] <- c(x_raw[seq.int(k + 1L, T_dat)], rep(NA_real_, k))
        }
      }
      new_cols[i] <- col_name
    }

    # Replace the lag/lead term with expanded column names
    replacement <- paste(new_cols, collapse = " + ")
    formula_chr <- sub(pattern, replacement, formula_chr, perl = TRUE)
  }

  list(formula_chr = formula_chr, data = data)
}


# =======================================================================
# Main function
# =======================================================================

#' Classical Local Projections (fixest-like syntax)
#'
#' Estimates impulse responses via the Jordà (2005) local projection method.
#' The syntax mirrors \pkg{fixest}: use a formula for the model,
#' \code{..macroname} tokens to expand character-vector variable lists, and
#' \code{l()} / \code{f()} operators to create lags and leads on-the-fly.
#'
#' @param formula A formula of the form \code{y ~ shock + controls}.
#'   \itemize{
#'     \item \strong{LHS}: one variable, or \code{c(y1, y2)} for multi-equation.
#'       LHS variables are \emph{never} automatically placed on the RHS.
#'     \item \strong{RHS}: exactly what you want in the regression. No automatic
#'       lag augmentation. Use \code{l()} / \code{f()} for that.
#'     \item \strong{Macros}: write \code{..name} where \code{name} is a
#'       character vector in the calling environment.
#'     \item \strong{Lag operator}: \code{l(var, k)} or \code{l(var, k1:k2)}.
#'     \item \strong{Lead operator}: \code{f(var, k)} or \code{f(var, k1:k2)}.
#'     \item Macros inside lag/lead: \code{l(..controls, 1:3)}.
#'   }
#' @param data A \code{data.frame} or \code{tibble}. Must contain all variables
#'   referenced in the formula (before lag/lead expansion).
#' @param horizons Integer. Maximum horizon H, or the complete sequence
#'   \code{0:H}; projections run for \eqn{h = 0, 1, \ldots, H}. Default: 12.
#' @param shock Character. \strong{Required.} Name of the shock/impulse variable
#'   (must appear on the RHS). This is the variable whose coefficient at each
#'   horizon is returned as the IRF.
#' @param conf Numeric in (0, 1). One or more confidence levels for bands, e.g.
#'   \code{0.90} or \code{c(0.68, 0.95)}. When more than one level is supplied,
#'   \code{irfs_upper} and \code{irfs_lower} become named lists of matrices
#'   (one per level, keyed by the level as a string) and \code{fPlotLP()} draws
#'   one ribbon per level. Bands are always rebuilt in R from \code{irfs_se},
#'   so multi-level fits cost the same as single-level fits. Default: 0.90.
#' @param nw_lags Integer or \code{NULL}. Base Newey-West bandwidth.
#'   At horizon \eqn{h} the effective bandwidth is
#'   \code{max(nw_lags + h + nw_offset, 0)} (Bartlett kernel).
#'   Defaults to \code{floor(T^(1/3))}.
#' @param nw_offset Integer. Shift applied to the NW bandwidth rule. Default
#'   \code{1L}, the Miranda-Agrippino & Ricco convention (bandwidth grows as
#'   \eqn{h+1}). Set to \code{0L} to reproduce the classic Jordà (2005) rule,
#'   bandwidth = \eqn{h} (LP residuals at horizon \eqn{h} are MA(h) by
#'   construction). \code{fLP}'s \code{h} is 0-indexed (h=0 is the impact
#'   horizon), which makes this the exact equivalent of Cesa-Bianchi's VAR
#'   Toolbox (\code{LPmodel.m}: \code{OLSmodel(Y,X,0,hh-1)} with 1-indexed
#'   \code{hh}, so \code{hh-1 == h}) — needed to match that toolbox's SEs
#'   exactly, e.g. in the Jordà-Taylor (2025) replication.
#' @param store_full Logical. If \code{TRUE}, return full coefficient and
#'   standard-error matrices for every horizon. Defaults to \code{FALSE} to
#'   keep the result small and fast.
#' @param cumulative Logical. If \code{TRUE}, the LHS at horizon \eqn{h} is
#'   \eqn{y_{t+h} - y_{t-1}} (cumulative response from \eqn{t-1} to \eqn{t+h}).
#'   Default \code{FALSE}.
#' @param n_threads Integer. Number of OpenMP threads for the horizon loop.
#'   \code{0} = all cores minus one. Default \code{0}.
#'
#' @return An object of class \code{"fLP"}.
#'
#' @examples
#' \dontrun{
#' ## 1. Basic usage — explicit shock, explicit controls
#' fLP(ip ~ mp_shock + ip_l1 + ip_l2 + cpi_l1,
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 2. Macro expansion
#' controls <- c("ip_l1", "ip_l2", "cpi_l1", "rate_l1")
#' fLP(ip ~ mp_shock + ..controls,
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 3. Lag operator — auto-creates ip_l1, ip_l2, ip_l3 from data
#' fLP(ip ~ mp_shock + l(ip, 1:3) + l(cpi, 1:2),
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 4. Macro + lag operator
#' vars <- c("ip", "cpi", "rate")
#' fLP(ip ~ mp_shock + l(..vars, 1:3),
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 5. Multi-equation (ip, cpi, rate as LHS — separate regressions)
#' fLP(c(ip, cpi, rate) ~ mp_shock + l(..vars, 1:3),
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 6. Lead operator
#' fLP(ip ~ mp_shock + l(ip, 1:12) + f(sentiment, 1),
#'     data = df, horizons = 24, shock = "mp_shock")
#' }
#'
#' @export
fLP <- function(formula, data, horizons = 12L,
                shock,
                conf         = 0.90,
                nw_lags      = NULL,
                nw_offset    = 1L,
                store_full   = FALSE,
                cumulative   = FALSE,
                n_threads    = 0L) {

  # Capture caller environment for macro lookup
  env <- parent.frame()

  # =====================================================================
  # 0. Require 'shock' argument
  # =====================================================================
  if (missing(shock)) {
    stop(
      "fLP: 'shock' is required. ",
      "Pass the name of the impulse variable, e.g. shock = \"mp_shock\"."
    )
  }
  if (!is.character(shock) || length(shock) != 1L || is.na(shock)) {
    stop("fLP: 'shock' must be a single character string.")
  }

  H <- .fLP_validate_horizon_max(horizons)
  conf <- .fLP_validate_conf_vector(conf)   # sorted descending (widest first)
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # =====================================================================
  # 1. Stringify the formula for text-based processing
  # =====================================================================
  formula_chr <- paste(deparse(formula), collapse = "")

  # =====================================================================
  # 2a. Expand l(..macro, lags) / f(..macro, lags)
  #     Must happen BEFORE general ..macro expansion.
  # =====================================================================
  formula_chr <- .fLP_expand_macro_lag(formula_chr, env)

  # =====================================================================
  # 2b. Expand plain ..macro tokens
  # =====================================================================
  formula_chr <- .fLP_expand_macro(formula_chr, env)

  # =====================================================================
  # 3. Expand l(var, lags) / f(var, leads) — creates columns in data
  # =====================================================================
  expanded      <- .fLP_expand_lag_terms(formula_chr, data)
  formula_chr   <- expanded$formula_chr
  data          <- expanded$data   # may have new lag/lead columns

  # =====================================================================
  # 4. Rebuild formula and parse LHS / RHS variable names
  # =====================================================================
  formula_expanded <- tryCatch(
    as.formula(formula_chr, env = env),
    error = function(e) stop(
      "fLP: could not parse the expanded formula:\n  ", formula_chr,
      "\n  ", conditionMessage(e)
    )
  )

  lhs_expr <- formula_expanded[[2L]]
  rhs_expr <- formula_expanded[[3L]]

  # LHS: handles plain 'y' and 'c(y1, y2)'
  lhs_vars <- all.vars(lhs_expr)

  # RHS: one element per term (no intercept; we add it internally).
  # deparse() can split long expressions across multiple strings; collapse
  # to a single string so as.formula() doesn't emit its length>1 warning.
  rhs_str   <- paste(deparse(rhs_expr), collapse = " ")
  rhs_terms <- terms(as.formula(paste("~", rhs_str), env = env),
                     keep.order = TRUE)
  rhs_vars  <- attr(rhs_terms, "term.labels")

  if (length(lhs_vars) == 0L)
    stop("fLP: formula must have at least one LHS variable.")
  if (length(rhs_vars) == 0L)
    stop("fLP: formula must have at least one RHS variable.")

  # =====================================================================
  # 5. Validate shock variable
  # =====================================================================
  if (!shock %in% rhs_vars) {
    stop(sprintf(
      "fLP: shock variable '%s' not found on the RHS.\n  RHS terms: %s",
      shock, paste(rhs_vars, collapse = ", ")
    ))
  }

  # 0-indexed column position in X (C++ convention)
  shock_col <- which(rhs_vars == shock) - 1L

  # =====================================================================
  # 6. Validate all variables exist in data
  # =====================================================================
  all_needed <- unique(c(lhs_vars, rhs_vars))
  missing_v  <- setdiff(all_needed, names(data))
  if (length(missing_v) > 0L) {
    stop(sprintf(
      "fLP: the following variables are not in data: %s",
      paste(missing_v, collapse = ", ")
    ))
  }

  # =====================================================================
  # 7. Build numeric matrices Y (T x ny) and X (T x k)
  #    Remove rows with NA in any used variable.
  # =====================================================================
  raw     <- data[, all_needed, drop = FALSE]
  na_rows <- which(!complete.cases(raw))
  if (length(na_rows) > 0L) {
    message(sprintf(
      "fLP: removing %d row(s) with NA (from lag/lead or missing data).",
      length(na_rows)
    ))
    raw <- raw[-na_rows, , drop = FALSE]
  }

  non_numeric <- names(raw)[!vapply(raw, is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop(sprintf(
      "fLP: all model variables must be numeric. Non-numeric variable(s): %s",
      paste(non_numeric, collapse = ", ")
    ))
  }

  Y <- as.matrix(raw[, lhs_vars, drop = FALSE])
  X <- as.matrix(raw[, rhs_vars, drop = FALSE])

  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"

  T_eff <- nrow(Y)

  # =====================================================================
  # 8. Default Newey-West base bandwidth
  # =====================================================================
  if (is.null(nw_lags)) {
    nw_lags <- as.integer(floor(T_eff^(1.0 / 3.0)))
  } else {
    nw_lags <- .fLP_validate_scalar_int(nw_lags, "nw_lags")
  }
  nw_offset <- .fLP_validate_scalar_int_signed(nw_offset, "nw_offset")

  # =====================================================================
  # 9. Sanity checks
  # =====================================================================
  if (!is.logical(store_full) || length(store_full) != 1L || is.na(store_full))
    stop("fLP: 'store_full' must be TRUE or FALSE.")
  if (!is.logical(cumulative) || length(cumulative) != 1L || is.na(cumulative))
    stop("fLP: 'cumulative' must be TRUE or FALSE.")
  n_threads <- .fLP_validate_scalar_int(n_threads, "n_threads")
  min_obs <- ncol(X) + 2L + H + as.integer(cumulative)
  if (T_eff < min_obs) {
    stop(sprintf(
      "fLP: not enough observations. Need at least %d rows after NA removal; have %d.",
      min_obs, T_eff
    ))
  }

  # =====================================================================
  # 10. Call C++ engine (all horizon loops in native C++)
  # =====================================================================
  res_cpp <- fLP_cpp(
    Y            = Y,
    X            = X,
    H            = H,
    shock_col    = as.integer(shock_col),
    conf_level   = as.double(conf[1L]),   # any level works; bands rebuilt below
    nw_lags_base = as.integer(nw_lags),
    store_full   = isTRUE(store_full),
    cumulative   = isTRUE(cumulative),
    n_threads    = as.integer(n_threads),
    nw_offset    = as.integer(nw_offset)
  )

  # ---------------------------------------------------------------------
  # Rebuild bands from irfs_se so multiple confidence levels come from a
  # single fit. Single conf: irfs_upper/lower stay matrices (unchanged
  # from the C++ output shape). Multi conf: they become named lists of
  # matrices keyed by the level printed as a string ("0.68", "0.95", ...).
  # ---------------------------------------------------------------------
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
    upper_list <- vector("list", length(conf))
    lower_list <- vector("list", length(conf))
    names(upper_list) <- keys
    names(lower_list) <- keys
    for (i in seq_along(conf)) {
      b <- .build_band(conf[i])
      upper_list[[i]] <- b$upper
      lower_list[[i]] <- b$lower
    }
    res_cpp$irfs_upper <- upper_list
    res_cpp$irfs_lower <- lower_list
  }

  # =====================================================================
  # 11. Annotate and return
  # =====================================================================
  h_seq <- 0L:H

  rownames(res_cpp$irfs)    <- as.character(h_seq)
  colnames(res_cpp$irfs)    <- lhs_vars
  rownames(res_cpp$irfs_se) <- as.character(h_seq)
  colnames(res_cpp$irfs_se) <- lhs_vars

  .label_band <- function(m) {
    rownames(m) <- as.character(h_seq)
    colnames(m) <- lhs_vars
    m
  }
  if (is.list(res_cpp$irfs_upper)) {
    res_cpp$irfs_upper <- lapply(res_cpp$irfs_upper, .label_band)
    res_cpp$irfs_lower <- lapply(res_cpp$irfs_lower, .label_band)
  } else {
    res_cpp$irfs_upper <- .label_band(res_cpp$irfs_upper)
    res_cpp$irfs_lower <- .label_band(res_cpp$irfs_lower)
  }

  if (isTRUE(store_full)) {
    rhs_with_const <- c("(Intercept)", rhs_vars)
    for (i in seq_along(res_cpp$betas)) {
      rownames(res_cpp$betas[[i]]) <- rhs_with_const
      colnames(res_cpp$betas[[i]]) <- lhs_vars
      rownames(res_cpp$ses[[i]])   <- rhs_with_const
      colnames(res_cpp$ses[[i]])   <- lhs_vars
    }
    names(res_cpp$betas) <- paste0("h", h_seq)
    names(res_cpp$ses)   <- paste0("h", h_seq)
  }

  res_cpp$lhs_vars        <- lhs_vars
  res_cpp$rhs_vars        <- rhs_vars
  res_cpp$shock           <- shock
  res_cpp$horizons        <- h_seq
  res_cpp$conf            <- conf
  res_cpp$nw_lags         <- nw_lags
  res_cpp$nw_offset       <- nw_offset
  res_cpp$store_full      <- isTRUE(store_full)
  res_cpp$cumulative      <- isTRUE(cumulative)
  res_cpp$n_threads       <- as.integer(n_threads)
  res_cpp$nobs            <- T_eff
  res_cpp$formula         <- formula_expanded
  res_cpp$formula_orig    <- formula
  res_cpp$call            <- match.call()

  class(res_cpp) <- "fLP"
  return(res_cpp)
}


# =======================================================================
# S3 methods
# =======================================================================

#' @export
print.fLP <- function(x, digits = 4L, ...) {
  cat("\nLocal Projections (fLP)\n")
  cat(strrep("-", 45L), "\n")
  cat("Original formula : ", deparse(x$formula_orig), "\n")
  cat("Expanded formula : ", deparse(x$formula),      "\n")
  cat("Shock            : ", x$shock,                 "\n")
  cat("LHS variables    : ", paste(x$lhs_vars, collapse = ", "), "\n")
  cat("Horizons         : 0 to", max(x$horizons),     "\n")
  cat("Cumulative       : ", isTRUE(x$cumulative),     "\n")
  cat("Confidence       : ",
      paste0(format(x$conf * 100, trim = TRUE), "%", collapse = ", "), "\n")
  nw_off <- if (!is.null(x$nw_offset)) x$nw_offset else 1L
  cat("NW lags          : base =", x$nw_lags,
      sprintf(", offset = %d -> effective at horizon h: max(base + h + %d, 0)\n",
              nw_off, nw_off))
  cat("HAC variance     :  full Newey-West sandwich (X'X)^-1 X'\u03a9X (X'X)^-1\n")
  cat("Observations     : ", x$nobs, "\n")
  cat("\nIRF (shock = '", x$shock, "'):\n", sep = "")
  print(round(x$irfs, digits))
  invisible(x)
}


#' @export
coef.fLP <- function(object, ...) {
  object$irfs
}


#' @export
tidy.fLP <- function(x, ...) {
  H  <- max(x$horizons)
  ny <- length(x$lhs_vars)

  multi <- is.list(x$irfs_upper)
  h_seq <- 0L:H

  base <- data.frame(
    horizon  = rep(h_seq, times = ny),
    lhs      = rep(x$lhs_vars, each = length(h_seq)),
    shock    = x$shock,
    estimate = as.vector(x$irfs),
    se       = as.vector(x$irfs_se),
    stringsAsFactors = FALSE
  )

  if (!multi) {
    base$lower <- as.vector(x$irfs_lower)
    base$upper <- as.vector(x$irfs_upper)
    return(base)
  }

  keys <- names(x$irfs_upper)
  for (k in keys) {
    base[[paste0("lower_", k)]] <- as.vector(x$irfs_lower[[k]])
    base[[paste0("upper_", k)]] <- as.vector(x$irfs_upper[[k]])
  }
  base
}
