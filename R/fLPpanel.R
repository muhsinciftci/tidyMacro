# -----------------------------------------------------------------------
# fLPpanel.R — Panel Local Projections with fixest-like formula syntax
#
# Provides:
#   fLPpanel()       — main user-facing function
#   print.fLPpanel() — print method
#   coef.fLPpanel()  — extract IRF matrix
#   tidy.fLPpanel()  — broom-compatible tidy method
#
# Formula syntax (fixest-style):
#   fLPpanel(y ~ shock + control1 + control2 | id + time, data,
#            id = "id", time = "time", shock = "shock", horizons = 12)
#
#   - LHS: single response variable.
#   - RHS-main (before "|"): shock + additional controls.
#   - RHS-FE (after "|"): categorical fixed effects (any number).
#     Omit "|" for no fixed effects.
#
# Macros (like fixest ..vars):
#   controls <- c("gdp_l1", "inf_l1", "rate_l1")
#   fLPpanel(y ~ mp_shock + ..controls | id + time, ..., shock = "mp_shock")
#
# Panel-aware lag / lead operators (respect unit boundaries):
#   l(x, k)     lag x by k periods within each unit  → column "x_lk"
#   l(x, k1:k2) → columns "x_lk1", ..., "x_lk2"
#   f(x, k)     lead x by k periods within each unit → column "x_fk"
#   l(..vars, 1:3) works too.
#
# Heterogeneous responses (Almuzara-Sancibrián s ⊗ X):
#   fLPpanel(y ~ mp_shock | id + time, ..., shock = "mp_shock",
#            interact = c("size", "leverage"))
#   returns one IRF column per interaction variable.
#
# Compatible with fPlotLP() via class = c("fLPpanel", "fLP").
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------


# =======================================================================
# Internal helpers
# =======================================================================

# Panel-aware lag/lead expander — like .fLP_expand_lag_terms but respects
# unit boundaries. Creates one column per (var, k) pair on a *copy* of
# data, replacing the l() / f() token in the formula string.
.fLPpanel_expand_lag_terms <- function(formula_chr, data, id_var, time_var) {

  pattern <- "(l|f)\\(([A-Za-z_.][A-Za-z0-9_.]*),([^)]+)\\)"

  # Row order for the shift: sort by (id, time)
  ord     <- order(data[[id_var]], data[[time_var]])
  inv_ord <- order(ord)
  id_s    <- data[[id_var]][ord]
  # Boundary: TRUE when the row starts a new unit
  new_unit <- c(TRUE, id_s[-1L] != id_s[-length(id_s)])
  # Cumulative unit index within the sorted vector: 1, 1, 1, 2, 2, ...
  unit_id  <- cumsum(new_unit)
  T_dat    <- nrow(data)

  panel_lag <- function(x_raw, k, op) {
    # x_raw is in original order; sort by (id, time), shift within unit,
    # then invert the ordering to restore original row indexing.
    x_s <- x_raw[ord]
    y_s <- rep(NA_real_, T_dat)
    if (op == "l") {
      # y_s[i] = x_s[i - k] iff row i-k belongs to the same unit
      if (k >= T_dat) return(rep(NA_real_, T_dat))
      idx <- seq_len(T_dat - k) + k
      valid <- unit_id[idx] == unit_id[idx - k]
      y_s[idx[valid]] <- x_s[(idx - k)[valid]]
    } else {
      if (k >= T_dat) return(rep(NA_real_, T_dat))
      idx <- seq_len(T_dat - k)
      valid <- unit_id[idx] == unit_id[idx + k]
      y_s[idx[valid]] <- x_s[(idx + k)[valid]]
    }
    y_s[inv_ord]
  }

  repeat {
    m <- regexpr(pattern, formula_chr, perl = TRUE)
    if (m == -1L) break

    full_match <- regmatches(formula_chr, m)
    parsed     <- regmatches(full_match,
                             regexec("(l|f)\\(([A-Za-z_.][A-Za-z0-9_.]*),([^)]+)\\)",
                                     full_match, perl = TRUE))[[1]]
    op       <- parsed[2L]
    var      <- trimws(parsed[3L])
    lags_str <- trimws(parsed[4L])

    lag_idx <- tryCatch(
      eval(parse(text = lags_str)),
      error = function(e) stop(sprintf(
        "fLPpanel: cannot parse lag/lead spec '%s' in term '%s': %s",
        lags_str, full_match, conditionMessage(e)
      ))
    )
    if (!is.numeric(lag_idx) || anyNA(lag_idx) || any(!is.finite(lag_idx)) ||
        any(lag_idx < 0L) || any(lag_idx != floor(lag_idx))) {
      stop(sprintf(
        "fLPpanel: lag/lead indices must be non-negative integers; got '%s'.",
        lags_str
      ))
    }
    lag_idx <- as.integer(lag_idx)

    if (!var %in% names(data)) {
      stop(sprintf(
        "fLPpanel: variable '%s' used in lag/lead term '%s' not found in data.",
        var, full_match
      ))
    }
    x_raw <- data[[var]]

    new_cols <- character(length(lag_idx))
    for (i in seq_along(lag_idx)) {
      k <- lag_idx[i]
      col_name <- if (op == "l") sprintf("%s_l%d", var, k) else
                                 sprintf("%s_f%d", var, k)
      if (!col_name %in% names(data)) {
        if (k == 0L) {
          data[[col_name]] <- x_raw
        } else {
          data[[col_name]] <- panel_lag(x_raw, k, op)
        }
      }
      new_cols[i] <- col_name
    }

    # Wrap in parens so downstream ":"-interactions distribute correctly:
    # l(x, 1:2):z  →  (x_l1 + x_l2):z  →  x_l1:z + x_l2:z.
    replacement <- if (length(new_cols) > 1L)
      paste0("(", paste(new_cols, collapse = " + "), ")")
    else
      new_cols[1L]
    formula_chr <- sub(pattern, replacement, formula_chr, perl = TRUE)
  }

  list(formula_chr = formula_chr, data = data)
}


# =======================================================================
# Main function
# =======================================================================

#' Panel Local Projections (fixest-like syntax)
#'
#' Estimates panel impulse responses via the local-projection estimator of
#' Almuzara & Sancibrián (2024, NY Fed SR 1090). Formula syntax mirrors
#' \pkg{fixest}: RHS variables before \code{|} are the shock and controls,
#' variables after \code{|} are absorbed fixed effects. Interactions
#' (e.g. \code{shock:size}) are computed explicitly from the formula — the
#' estimator uses \emph{exactly} the specification the user writes.
#' Supports \code{..macro} expansion and panel-aware \code{l()} / \code{f()}
#' operators.
#'
#' @param formula A formula of the form
#'   \code{y ~ shock + controls | fe1 + fe2}, e.g.
#'   \code{y ~ shock:size + l(y, 1:2) | unit + tt}.
#'   \itemize{
#'     \item \strong{LHS}: single response variable.
#'     \item \strong{RHS-main} (before \code{|}): every term you want in the
#'       regression. Interactions with \code{:} are computed as elementwise
#'       products (this is where the Almuzara-Sancibrián \eqn{s_{it} \cdot X_t}
#'       appears — write \code{shock:size} in the formula).
#'     \item \strong{RHS-FE} (after \code{|}): categorical fixed-effect
#'       variables (omit \code{|} for no fixed effects).
#'     \item \strong{Macros}: \code{..name} where \code{name} is a
#'       character vector in the calling environment.
#'     \item \strong{Panel-aware lag/lead}: \code{l(x, k)}, \code{l(x, k1:k2)},
#'       \code{f(x, k)}, \code{l(..vars, 1:3)}. Never crosses unit boundaries.
#'   }
#' @param data A \code{data.frame} or \code{tibble} in long panel form.
#' @param panel_id Length-2 character vector \code{c(id_col, time_col)}
#'   naming the unit-id column and the time column, respectively.
#' @param horizons Integer. Maximum horizon H (projections run for
#'   \eqn{h = 0, 1, \ldots, H}) or the full sequence \code{0:H}. Default 12.
#' @param shock Character (single or vector). The RHS-main term label(s)
#'   whose coefficient is the impulse response. Everything else in the
#'   RHS-main is treated as a control. Must match the term label exactly
#'   (after macro / lag expansion); interactions are labelled with
#'   \code{":"}, e.g. \code{"shock:size"}.
#' @param conf Numeric in (0, 1). One or more confidence levels, e.g.
#'   \code{0.90} or \code{c(0.68, 0.95)}. Bands are rebuilt from the raw
#'   SEs, so multi-level fits cost the same as single-level fits.
#' @param small_sample Logical. If \code{TRUE}, apply the Imbens-Kolesar
#'   (2016, REStat) small-sample refinement. Default \code{FALSE}.
#' @param cumulative Logical. If \code{TRUE}, project
#'   \eqn{\sum_{k=0}^h y_{i,t+k}} instead of \eqn{y_{i,t+h}}. Default
#'   \code{FALSE}.
#' @param p_max Integer. Maximum number of lags of the regressand and of
#'   the shock terms added as controls (at horizon h the effective number
#'   of lags is \eqn{\min(h, p_{\max})}). Default 0.
#' @param n_threads Integer. OpenMP threads for the horizon loop. 0 = all
#'   cores minus one. Default 0.
#'
#' @return An object of class \code{c("fLPpanel", "fLP")} that plugs into
#'   \code{\link{fPlotLP}} and \code{\link{tidy.fLP}}.
#'
#' @examples
#' \dontrun{
#' # Homogeneous shock
#' fLPpanel(y ~ shock + l(y, 1:2) | unit + tt,
#'          data = df, panel_id = c("unit", "tt"),
#'          shock = "shock", horizons = 12)
#'
#' # Heterogeneous shock via explicit interaction (Almuzara-Sancibrián)
#' fLPpanel(y ~ shock:size | unit + tt,
#'          data = df, panel_id = c("unit", "tt"),
#'          shock = "shock:size", horizons = 12,
#'          small_sample = TRUE, cumulative = TRUE)
#'
#' # Multiple heterogeneity margins → multiple IRFs
#' fLPpanel(y ~ shock:size + shock:leverage | unit + tt,
#'          data = df, panel_id = c("unit", "tt"),
#'          shock = c("shock:size", "shock:leverage"), horizons = 12)
#' }
#'
#' @export
fLPpanel <- function(formula, data,
                     panel_id,
                     horizons     = 12L,
                     shock,
                     conf         = 0.90,
                     small_sample = FALSE,
                     cumulative   = FALSE,
                     p_max        = 0L,
                     n_threads    = 0L) {

  env <- parent.frame()

  # ---- required arguments ----------------------------------------------
  if (missing(shock))    stop("fLPpanel: 'shock' is required.")
  if (missing(panel_id)) stop("fLPpanel: 'panel_id' is required (c(id_col, time_col)).")

  if (!is.character(panel_id) || length(panel_id) != 2L || anyNA(panel_id))
    stop("fLPpanel: 'panel_id' must be a length-2 character vector, e.g. c('unit', 'tt').")
  id_var   <- panel_id[1L]
  time_var <- panel_id[2L]

  if (!is.character(shock) || length(shock) == 0L || anyNA(shock))
    stop("fLPpanel: 'shock' must be a non-empty character vector.")

  H    <- .fLP_validate_horizon_max(horizons)
  conf <- .fLP_validate_conf_vector(conf)

  if (!is.logical(small_sample) || length(small_sample) != 1L)
    stop("fLPpanel: 'small_sample' must be TRUE or FALSE.")
  if (!is.logical(cumulative)   || length(cumulative)   != 1L)
    stop("fLPpanel: 'cumulative' must be TRUE or FALSE.")

  p_max     <- .fLP_validate_scalar_int(p_max, "p_max")
  n_threads <- .fLP_validate_scalar_int(n_threads, "n_threads")

  if (!is.data.frame(data)) data <- as.data.frame(data)

  if (!id_var   %in% names(data))
    stop(sprintf("fLPpanel: id column '%s' not found in data.", id_var))
  if (!time_var %in% names(data))
    stop(sprintf("fLPpanel: time column '%s' not found in data.", time_var))

  # ---- formula → string -------------------------------------------------
  formula_chr <- paste(deparse(formula), collapse = " ")

  # Macro helpers from fLP (pure string operations)
  formula_chr <- .fLP_expand_macro_lag(formula_chr, env)
  formula_chr <- .fLP_expand_macro(formula_chr, env)

  # ---- split on "|" ----------------------------------------------------
  parts    <- strsplit(formula_chr, "\\|", fixed = FALSE)[[1]]
  main_str <- trimws(parts[1])
  fe_str   <- if (length(parts) > 1L) trimws(paste(parts[-1L], collapse = "|")) else ""

  # ---- panel-aware lag/lead expansion (main part only) ------------------
  main_ex  <- .fLPpanel_expand_lag_terms(main_str, data, id_var, time_var)
  main_str <- main_ex$formula_chr
  data     <- main_ex$data

  # ---- parse main formula ----------------------------------------------
  main_formula <- tryCatch(
    as.formula(main_str, env = env),
    error = function(e) stop(
      "fLPpanel: could not parse the main formula:\n  ", main_str,
      "\n  ", conditionMessage(e)
    )
  )
  lhs_expr <- main_formula[[2L]]
  rhs_expr <- main_formula[[3L]]

  lhs_vars <- all.vars(lhs_expr)
  if (length(lhs_vars) != 1L)
    stop("fLPpanel: exactly one LHS variable is required.")

  rhs_str   <- paste(deparse(rhs_expr), collapse = " ")
  rhs_terms <- terms(as.formula(paste("~", rhs_str), env = env),
                     keep.order = TRUE)
  rhs_labels <- attr(rhs_terms, "term.labels")
  # Raw base columns (e.g. 'shock', 'size' for term 'shock:size')
  rhs_base_vars <- all.vars(rhs_terms)

  # ---- validate shock terms -------------------------------------------
  bad_shocks <- setdiff(shock, rhs_labels)
  if (length(bad_shocks))
    stop(sprintf(
      "fLPpanel: shock term(s) not found in RHS-main: %s. Available: %s",
      paste(bad_shocks, collapse = ", "),
      paste(rhs_labels, collapse = ", ")))

  control_labels <- setdiff(rhs_labels, shock)

  # ---- parse FE variables ----------------------------------------------
  if (nzchar(fe_str)) {
    fe_expr <- tryCatch(
      as.formula(paste("~", fe_str), env = env),
      error = function(e) stop(
        "fLPpanel: could not parse the fixed-effects part '", fe_str,
        "': ", conditionMessage(e)
      )
    )
    fe_vars <- attr(terms(fe_expr, keep.order = TRUE), "term.labels")
  } else {
    fe_vars <- character(0L)
  }

  # ---- collect all needed base columns and drop NA rows ---------------
  all_needed <- unique(c(lhs_vars, rhs_base_vars, fe_vars, id_var, time_var))
  miss <- setdiff(all_needed, names(data))
  if (length(miss))
    stop(sprintf("fLPpanel: columns not found in data: %s",
                 paste(miss, collapse = ", ")))

  raw <- data[, all_needed, drop = FALSE]
  na_rows <- which(!complete.cases(raw))
  if (length(na_rows)) {
    message(sprintf("fLPpanel: dropping %d row(s) with NA.", length(na_rows)))
    raw <- raw[-na_rows, , drop = FALSE]
  }

  # ---- numeric checks on all vars used in the model matrix -------------
  for (nm in unique(c(lhs_vars, rhs_base_vars))) {
    if (!is.numeric(raw[[nm]]))
      stop(sprintf("fLPpanel: column '%s' must be numeric.", nm))
  }

  # ---- build the RHS design matrix from the formula --------------------
  # One column per term label (interactions computed as elementwise
  # products). No intercept — HDFE absorbs it.
  design_formula <- as.formula(paste("~", rhs_str, "- 1"), env = env)
  design <- stats::model.matrix(design_formula, data = raw)

  # model.matrix may rename term labels (e.g. reorder factors). Since all
  # our terms are numeric, colnames should match rhs_labels one-for-one,
  # in order.
  if (ncol(design) != length(rhs_labels)) {
    stop("fLPpanel: internal error — design matrix has ",
         ncol(design), " columns but ", length(rhs_labels), " term labels.")
  }
  colnames(design) <- rhs_labels

  shock_cols   <- match(shock, rhs_labels)
  control_cols <- setdiff(seq_along(rhs_labels), shock_cols)

  y_mat <- matrix(as.double(raw[[lhs_vars]]), ncol = 1L)
  W_mat <- if (length(control_cols)) design[, control_cols, drop = FALSE]
           else matrix(0, nrow(raw), 0L)
  storage.mode(W_mat) <- "double"

  # ---- shock / heterogeneity split -----------------------------------
  # Convention:
  #   * Bare shock term ("shock")  → X = data[["shock"]], s = ones
  #   * Single interaction "A:B"   → X = data[[A]], s = data[[B]]
  #     (engine forms s ⊗ X internally; mathematically identical to
  #      passing X = A * B and s = ones for point estimates and the
  #      time-clustered LAHR SE, and matches the paper's Imbens-Kolesar
  #      refinement — which uses the heterogeneity vector s as the
  #      per-period weighting variable)
  #   * Multiple shock terms       → X = design's shock cols, s = ones
  is_pairwise <- function(lbl) {
    p <- strsplit(lbl, ":", fixed = TRUE)[[1]]
    length(p) == 2L && all(p %in% names(raw)) &&
      all(vapply(p, \(nm) is.numeric(raw[[nm]]), logical(1L)))
  }

  if (length(shock) == 1L && is_pairwise(shock)) {
    parts  <- strsplit(shock, ":", fixed = TRUE)[[1]]
    X_mat  <- matrix(as.double(raw[[parts[1L]]]), ncol = 1L)
    s_mat  <- matrix(as.double(raw[[parts[2L]]]), ncol = 1L)
  } else {
    X_mat  <- design[, shock_cols, drop = FALSE]
    storage.mode(X_mat) <- "double"
    s_mat  <- matrix(1.0, nrow(raw), 1L)
  }

  fe_mat <- if (length(fe_vars)) {
    m <- matrix(0L, nrow(raw), length(fe_vars))
    for (j in seq_along(fe_vars))
      m[, j] <- as.integer(factor(raw[[fe_vars[j]]]))
    m
  } else {
    matrix(0L, nrow(raw), 0L)
  }

  i_index <- as.integer(factor(raw[[id_var]]))

  # Time index: normalize so that one unit-step = one period. This mirrors
  # panel_LP.m and handles unbalanced panels (including globally missing
  # periods) correctly: gaps in the observed time grid are PRESERVED, so
  # `l(x, 1)` at t = 5 correctly returns NaN if t = 4 is unobserved.
  # Only rank-based mapping would (wrongly) close the gap and skip over
  # the missing period.
  t_raw <- raw[[time_var]]
  if (inherits(t_raw, "Date") || inherits(t_raw, "POSIXct")) {
    t_raw <- as.numeric(t_raw)
  }
  if (!is.numeric(t_raw)) {
    stop("fLPpanel: 'time' column must be numeric, Date, or POSIXct.")
  }
  t_uniq <- sort(unique(t_raw))
  if (length(t_uniq) < 2L) {
    stop("fLPpanel: 'time' has fewer than 2 unique values.")
  }
  t_diff  <- min(diff(t_uniq))
  t_scaled <- (t_raw - t_uniq[1L]) / t_diff + 1
  on_grid  <- abs(t_scaled - round(t_scaled)) < 1e-2
  if (!all(on_grid)) {
    message(sprintf(
      "fLPpanel: dropping %d row(s) with off-grid time index.",
      sum(!on_grid)))
    raw     <- raw[on_grid, , drop = FALSE]
    y_mat   <- y_mat[on_grid, , drop = FALSE]
    X_mat   <- X_mat[on_grid, , drop = FALSE]
    W_mat   <- W_mat[on_grid, , drop = FALSE]
    s_mat   <- s_mat[on_grid, , drop = FALSE]
    fe_mat  <- fe_mat[on_grid, , drop = FALSE]
    i_index <- i_index[on_grid]
    t_scaled <- t_scaled[on_grid]
  }
  t_index <- as.integer(round(t_scaled))

  # ---- call C++ engine --------------------------------------------------
  res <- fLPpanel_cpp(
    y            = y_mat,
    s            = s_mat,
    X            = X_mat,
    W            = W_mat,
    FE           = fe_mat,
    i_index      = i_index,
    t_index      = t_index,
    H            = H,
    p_max        = p_max,
    small_sample = isTRUE(small_sample),
    cumulative   = isTRUE(cumulative),
    n_threads    = n_threads
  )

  # ---- build fLP-shaped result -----------------------------------------
  h_seq <- 0L:H

  # One IRF column per shock term (label = term label from the formula).
  irf_names <- shock

  colnames(res$estimate) <- irf_names
  colnames(res$SE)       <- irf_names
  colnames(res$df)       <- irf_names
  rownames(res$estimate) <- as.character(h_seq)
  rownames(res$SE)       <- as.character(h_seq)
  rownames(res$df)       <- as.character(h_seq)

  # Build bands at each requested conf level using per-cell t quantiles
  build_band <- function(cl) {
    q <- matrix(0, H + 1L, length(irf_names))
    for (i in seq_len(H + 1L))
      for (j in seq_along(irf_names))
        q[i, j] <- stats::qt(0.5 * (1 + cl), df = res$df[i, j])
    list(upper = res$estimate + q * res$SE,
         lower = res$estimate - q * res$SE)
  }

  if (length(conf) == 1L) {
    b <- build_band(conf)
    irfs_upper <- b$upper
    irfs_lower <- b$lower
    rownames(irfs_upper) <- rownames(irfs_lower) <- as.character(h_seq)
    colnames(irfs_upper) <- colnames(irfs_lower) <- irf_names
  } else {
    keys <- format(conf, trim = TRUE)
    irfs_upper <- vector("list", length(conf))
    irfs_lower <- vector("list", length(conf))
    names(irfs_upper) <- keys
    names(irfs_lower) <- keys
    for (i in seq_along(conf)) {
      b <- build_band(conf[i])
      rownames(b$upper) <- rownames(b$lower) <- as.character(h_seq)
      colnames(b$upper) <- colnames(b$lower) <- irf_names
      irfs_upper[[i]] <- b$upper
      irfs_lower[[i]] <- b$lower
    }
  }

  out <- list(
    irfs         = res$estimate,
    irfs_se      = res$SE,
    irfs_upper   = irfs_upper,
    irfs_lower   = irfs_lower,
    df           = res$df,
    pval         = res$pval,
    CI90         = res$CI90,   # engine-native fixed-level bands
    CI95         = res$CI95,
    CI99         = res$CI99,

    lhs_vars     = lhs_vars,
    rhs_vars     = rhs_labels,
    fe_vars      = fe_vars,
    shock        = shock,
    id           = id_var,
    time         = time_var,
    panel_id     = panel_id,
    horizons     = h_seq,
    conf         = conf,
    small_sample = isTRUE(small_sample),
    cumulative   = isTRUE(cumulative),
    p_max        = p_max,
    n_threads    = as.integer(n_threads),
    nobs         = nrow(raw),
    formula      = main_formula,
    formula_orig = formula,
    call         = match.call()
  )

  class(out) <- c("fLPpanel", "fLP")
  out
}


# =======================================================================
# S3 methods
# =======================================================================

#' @export
print.fLPpanel <- function(x, digits = 4L, ...) {
  cat("\nPanel Local Projections (fLPpanel)\n")
  cat(strrep("-", 45L), "\n")
  cat("Original formula : ", deparse(x$formula_orig), "\n")
  cat("Panel index      :  id = ", x$id, ", time = ", x$time, "\n", sep = "")
  cat("Shock term(s)    : ", paste(x$shock, collapse = ", "), "\n")
  cat("Response         : ", x$lhs_vars, "\n")
  if (length(x$fe_vars))
    cat("Fixed effects    : ", paste(x$fe_vars, collapse = " + "), "\n")
  else
    cat("Fixed effects    :  (none)\n")
  cat("Horizons         :  0 to ", max(x$horizons), "\n", sep = "")
  cat("Cumulative       : ", isTRUE(x$cumulative), "\n")
  cat("SE               : ",
      if (x$small_sample) "Imbens-Kolesar (2016) small-sample"
      else                "asymptotic time-clustered (LAHR)",
      "\n")
  cat("Confidence       : ",
      paste0(format(x$conf * 100, trim = TRUE), "%", collapse = ", "), "\n")
  cat("Observations     : ", x$nobs, "\n")
  cat("Lag order p_max  : ", x$p_max, "\n")
  cat("\nIRF (shock = '", x$shock, "'):\n", sep = "")
  print(round(x$irfs, digits))
  invisible(x)
}

#' @export
coef.fLPpanel <- function(object, ...) object$irfs

#' @export
tidy.fLPpanel <- function(x, ...) {
  H  <- max(x$horizons)
  ny <- ncol(x$irfs)
  irf_names <- colnames(x$irfs)
  multi <- is.list(x$irfs_upper)
  h_seq <- 0L:H

  base <- data.frame(
    horizon  = rep(h_seq, times = ny),
    shock    = rep(irf_names, each = length(h_seq)),
    estimate = as.vector(x$irfs),
    se       = as.vector(x$irfs_se),
    df       = as.vector(x$df),
    pval     = as.vector(x$pval),
    stringsAsFactors = FALSE
  )

  if (!multi) {
    base$lower <- as.vector(x$irfs_lower)
    base$upper <- as.vector(x$irfs_upper)
    return(base)
  }
  for (k in names(x$irfs_upper)) {
    base[[paste0("lower_", k)]] <- as.vector(x$irfs_lower[[k]])
    base[[paste0("upper_", k)]] <- as.vector(x$irfs_upper[[k]])
  }
  base
}
