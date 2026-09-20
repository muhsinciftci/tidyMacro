#' Sign, Narrative and Sign-Plus-IV Identification of a Structural VAR
#'
#' Front-end for \code{\link{fSR_cpp}}. Estimates a VAR, then identifies it by
#' sign restrictions, optionally combined with Antolin-Diaz and Rubio-Ramirez
#' narrative restrictions and with an external instrument that pins down the
#' first structural impact column. Inference is Bayesian throughout: reduced-form
#' parameters are drawn from their flat-prior Normal-inverse-Wishart posterior and
#' each draw is paired with a Haar-uniform rotation.
#'
#' @param y A T x N numeric matrix (or data frame) of endogenous variables.
#' @param p Integer lag order.
#' @param c Integer intercept indicator: 1 to include an intercept (default), 0
#'   to exclude it.
#' @param sign An N x ds matrix of sign restrictions: \code{+1} the response must
#'   be non-negative, \code{-1} non-positive, \code{0} unrestricted. \code{ds}
#'   must equal N when no instrument is used, and N - 1 with one.
#' @param nsteps Integer number of horizons to report, impact included
#'   (default 40).
#' @param ndraws Integer number of accepted draws (default 500).
#' @param sr_hor Integer. Restrictions are imposed at horizons
#'   \code{0, ..., sr_hor - 1} (default 1, impact only).
#' @param sr_rot Integer maximum rotations attempted per parameter draw
#'   (default 500).
#' @param max_post_draws Integer maximum parameter draws per accepted draw
#'   (default 1000).
#' @param conf Numeric coverage of the credible bands, in percent (default 90;
#'   Uhlig (2005) and Antolin-Diaz and Rubio-Ramirez (2018) report 68, so pass
#'   \code{conf = 68} to match them).
#'   A vector requests several levels at once, e.g. \code{c(68, 90)}: the first
#'   is reported in \code{IRinf}/\code{IRsup} and every level is returned in
#'   \code{bands}. Extra levels need \code{store_draws = TRUE}.
#' @param inference Integer: 1 (default) draws reduced-form parameters from the
#'   posterior; 0 holds them at OLS, leaving only rotation uncertainty.
#' @param narrative Optional list of narrative restrictions with components
#'   \code{sign} (a list of \code{shock}, \code{period}, \code{sign}) and
#'   \code{dom} (a list of \code{shock}, \code{period}, \code{var}). Periods may
#'   be integer rows of \code{y} or, when \code{dates} is supplied, date strings
#'   or \code{Date} values.
#' @param dates Optional vector of length \code{nrow(y)} — character labels or
#'   \code{Date} values — used to resolve narrative periods given as dates
#'   rather than as row numbers.
#' @param instrument Optional list with component \code{Z}, a T x k instrument
#'   matrix aligned with the rows of \code{y}, missing values allowed. Selects
#'   the \code{sign+iv} scheme. The identified column is held at its OLS value
#'   across draws, as the VAR Toolbox does, so the bands do not carry the
#'   instrument's own sampling uncertainty.
#' @param narr_weight_mc Integer. Positive values switch on the ADRR importance
#'   reweighting, using this many Monte Carlo replications per draw. \code{0}
#'   (default) uses plain rejection sampling, as the VAR Toolbox does.
#' @param resid_from_draw Logical. \code{FALSE} (default) evaluates narrative
#'   restrictions on OLS residuals, matching the VAR Toolbox; \code{TRUE} uses
#'   the residuals implied by each parameter draw.
#' @param store_draws Logical; keep the full IRF and FEVD draw distributions
#'   (default TRUE).
#' @param exog Optional T x M matrix of exogenous regressors.
#' @param varnames Optional character vector of length N.
#' @param n_threads Integer; 0 (default) uses all cores but one.
#' @param seed Integer base seed (default 42).
#' @param verbose Logical; print progress diagnostics (default FALSE).
#'
#' @return An object of class \code{"fSignRestr"}: the list returned by
#'   \code{\link{fSR_cpp}}, with \code{IRall} and \code{VDall} reshaped to
#'   \code{c(N, N, nsteps, ndraws)} and with \code{var}, \code{varnames},
#'   \code{ident}, \code{p}, \code{c}, \code{nsteps} and \code{iv} attached.
#'
#' @details
#' Impulse responses are indexed \code{[variable, shock, horizon]} and FEVD
#' shares are on the \code{[0, 1]} scale.
#'
#' With an instrument, the first column of \code{sign} corresponds to the second
#' structural shock: shock 1 is the instrumented one and is not sign-restricted.
#'
#' @references
#' Uhlig, H. (2005). \emph{Journal of Monetary Economics}, 52(2), 381--419.
#'
#' Antolin-Diaz, J., & Rubio-Ramirez, J. F. (2018). \emph{American Economic
#' Review}, 108(10), 2802--2829.
#'
#' @seealso \code{\link{fSR_cpp}}, \code{\link{fPlotIRFSign}},
#'   \code{\link{fRecoverBIV_cpp}}, \code{\link{fHDShock_cpp}}
#'
#' @export
fSignRestr <- function(y, p, c = 1, sign,
                       nsteps          = 40,
                       ndraws          = 500,
                       sr_hor          = 1,
                       sr_rot          = 500,
                       max_post_draws  = 1000,
                       conf            = 90,
                       inference       = 1,
                       narrative       = NULL,
                       dates           = NULL,
                       instrument      = NULL,
                       narr_weight_mc  = 0,
                       resid_from_draw = FALSE,
                       store_draws     = TRUE,
                       exog            = NULL,
                       varnames        = NULL,
                       n_threads       = 0,
                       seed            = 42,
                       verbose         = FALSE) {

    y <- as.matrix(y)
    if (!is.numeric(y)) stop("`y` must be numeric.")
    if (anyNA(y))       stop("`y` must not contain missing values.")

    N <- ncol(y)
    p <- as.integer(p)
    c <- as.integer(c)
    if (p < 1L)              stop("`p` must be at least 1.")
    if (!c %in% c(0L, 1L))   stop("`c` must be 0 or 1.")

    sign <- as.matrix(sign)
    if (nrow(sign) != N)
        stop(sprintf("`sign` must have %d rows, one per variable.", N))
    if (!all(sign %in% c(-1, 0, 1)))
        stop("`sign` entries must be -1, 0 or 1.")

    conf <- as.numeric(conf)
    if (length(conf) < 1L || anyNA(conf) || any(conf <= 0) || any(conf >= 100))
        stop("`conf` entries must lie strictly between 0 and 100.")
    ## A fraction would otherwise pass silently as a sub-1% band.
    if (any(conf < 1))
        stop("`conf` is in percent, not a fraction. Use ",
             paste(format(conf[conf < 1] * 100, trim = TRUE), collapse = ", "),
             " instead of ",
             paste(format(conf[conf < 1], trim = TRUE), collapse = ", "), ".")
    conf <- unique(conf)
    if (length(conf) > 1L && !isTRUE(store_draws))
        stop("Several `conf` levels require `store_draws = TRUE`.")

    if (is.null(varnames)) {
        varnames <- colnames(y)
        if (is.null(varnames)) varnames <- paste0("V", seq_len(N))
    }
    if (length(varnames) != N)
        stop(sprintf("`varnames` must have length %d.", N))

    if (!is.null(exog)) exog <- as.matrix(exog)

    var_ols <- fVAR(y, p = p, c = c, exog = exog)
    nobs    <- nrow(var_ols$residuals)

    ## ---- external instrument -------------------------------------------
    Bfix <- NULL; iv_info <- NULL
    ident <- "sign"

    if (!is.null(instrument)) {
        ident <- "sign+iv"
        ## `iv_refit` was removed: it re-estimated the first stage on each
        ## draw's residuals, which prices in reduced-form estimation error
        ## rather than instrument weakness, and so invited misreading as a
        ## weak-instrument correction. Fail loudly instead of ignoring it.
        if (!is.null(instrument$refit))
            stop("`instrument$refit` is no longer supported. The instrument column ",
                 "is always held at its OLS estimate, as the VAR Toolbox does; ",
                 "propagating first-stage uncertainty requires a proxy-augmented ",
                 "likelihood (Caldara and Herbst 2019), which is not implemented.")

        Z <- as.matrix(instrument$Z)
        if (nrow(Z) != nrow(y))
            stop("`instrument$Z` must have as many rows as `y`.")

        ## Align the instrument with the residual sample and keep the longest
        ## run of rows where every instrument column is observed.
        Zres  <- Z[(p + 1L):nrow(Z), , drop = FALSE]
        okrow <- stats::complete.cases(Zres)
        if (!any(okrow)) stop("`instrument$Z` has no observations inside the VAR sample.")
        idx   <- which(okrow)
        lo    <- min(idx); hi <- max(idx)
        if (!all(okrow[lo:hi]))
            stop("`instrument$Z` has gaps inside its sample; only a contiguous span is supported.")

        Z_sub <- Zres[lo:hi, , drop = FALSE]

        iv_info <- fRecoverBIV_cpp(
            resid_sub = var_ols$residuals[lo:hi, , drop = FALSE],
            Z_sub     = Z_sub,
            sigma     = var_ols$sigma,
            ntotcoeff = nrow(var_ols$beta))
        Bfix <- matrix(iv_info$b1, ncol = 1L)

        if (ncol(sign) != N - 1L)
            stop(sprintf("With an instrument, `sign` must have %d columns (shock 1 is the instrumented one).",
                         N - 1L))
    } else if (ncol(sign) != N) {
        stop(sprintf("`sign` must have %d columns.", N))
    }

    ## ---- narrative restrictions ----------------------------------------
    nr <- .fSR_narrative(narrative, dates, p, nobs, N)

    res <- fSR_cpp(y = y, p = p, c = c, SIGN = sign,
                   nsteps = as.integer(nsteps), ndraws = as.integer(ndraws),
                   sr_hor = as.integer(sr_hor), sr_rot = as.integer(sr_rot),
                   max_post_draws = as.integer(max_post_draws),
                   conf = conf[1], inference = as.integer(inference),
                   Bfix = Bfix,
                   narr_sign_shock  = nr$ns_shock,
                   narr_sign_period = nr$ns_period,
                   narr_sign_sign   = nr$ns_sign,
                   narr_dom_shock   = nr$nd_shock,
                   narr_dom_period  = nr$nd_period,
                   narr_dom_var     = nr$nd_var,
                   narr_weight_mc   = as.integer(narr_weight_mc),
                   resid_from_draw  = isTRUE(resid_from_draw),
                   store_draws      = isTRUE(store_draws),
                   exog = exog, n_threads = as.integer(n_threads),
                   seed = as.integer(seed), verbose = isTRUE(verbose))

    nkeep <- dim(res$Ball)[3]
    if (!is.null(res$IRall)) {
        dim(res$IRall) <- c(N, N, nsteps, nkeep)
        dim(res$VDall) <- c(N, N, nsteps, nkeep)
    }

    ## arma vectors arrive as one-column matrices; flatten for R ergonomics.
    res$weights <- as.numeric(res$weights)
    res$n_tried <- as.integer(res$n_tried)

    res$var      <- var_ols
    res$varnames <- varnames
    res$ident    <- ident
    res$p        <- p
    res$c        <- c
    res$nsteps   <- as.integer(nsteps)
    res$conf    <- conf
    res$bands    <- .fSR_bands(res, conf)
    res$iv       <- iv_info
    res$narrative_active <- nr$active

    class(res) <- c("fSignRestr", "list")
    res
}


#' @keywords internal
#' @noRd
.fSR_narrative <- function(narrative, dates, p, nobs, N) {

    empty <- list(ns_shock = NULL, ns_period = NULL, ns_sign = NULL,
                  nd_shock = NULL, nd_period = NULL, nd_var = NULL,
                  active = FALSE)
    if (is.null(narrative)) return(empty)
    if (!is.list(narrative)) stop("`narrative` must be a list.")

    ## Periods are expressed against the rows of `y`; the first p rows are
    ## absorbed as initial conditions, so the residual-sample row is period - p.
    ## A period may be a row number, a date string, or a Date; the latter two
    ## are matched against `dates`.
    resolve <- function(period, what) {
        if (is.character(period) || inherits(period, "Date")) {
            if (is.null(dates))
                stop(sprintf("`narrative$%s$period` is a date, so `dates` must be supplied.", what))
            key <- if (inherits(period, "Date")) format(period) else period
            ref <- if (inherits(dates, "Date")) format(dates) else as.character(dates)
            pos <- match(key, ref)
            if (anyNA(pos))
                stop(sprintf("Date(s) %s not found in `dates`.",
                             paste(sQuote(key[is.na(pos)]), collapse = ", ")))
        } else {
            pos <- as.integer(period)
        }
        idx <- pos - p
        if (any(idx < 1L))
            stop(sprintf("`narrative$%s$period` falls inside the %d-observation lag window.", what, p))
        if (any(idx > nobs))
            stop(sprintf("`narrative$%s$period` lies beyond the estimation sample.", what))
        as.integer(idx)
    }

    out <- empty
    if (!is.null(narrative$sign)) {
        s <- narrative$sign
        if (!all(c("shock", "period", "sign") %in% names(s)))
            stop("`narrative$sign` needs components `shock`, `period` and `sign`.")
        if (!all(s$sign %in% c(-1, 1)))
            stop("`narrative$sign$sign` entries must be -1 or 1.")
        if (any(s$shock < 1L | s$shock > N)) stop("`narrative$sign$shock` is out of range.")
        out$ns_shock  <- as.integer(s$shock)
        out$ns_period <- resolve(s$period, "sign")
        out$ns_sign   <- as.numeric(s$sign)
    }
    if (!is.null(narrative$dom)) {
        d <- narrative$dom
        if (!all(c("shock", "period", "var") %in% names(d)))
            stop("`narrative$dom` needs components `shock`, `period` and `var`.")
        if (any(d$shock < 1L | d$shock > N)) stop("`narrative$dom$shock` is out of range.")
        if (any(d$var   < 1L | d$var   > N)) stop("`narrative$dom$var` is out of range.")
        out$nd_shock  <- as.integer(d$shock)
        out$nd_period <- resolve(d$period, "dom")
        out$nd_var    <- as.integer(d$var)
    }
    out$active <- !is.null(out$ns_shock) || !is.null(out$nd_shock)
    out
}


#' @describeIn fSignRestr Print a compact summary of the identified model.
#' @param x An object of class \code{"fSignRestr"}.
#' @param ... Ignored.
#' @export
print.fSignRestr <- function(x, ...) {
    N <- length(x$varnames)
    cat("Sign-restricted SVAR (", x$ident, ")\n", sep = "")
    cat("  Variables      : ", paste(x$varnames, collapse = ", "), "\n", sep = "")
    cat("  Lags           : ", x$p, "   Intercept: ", x$c, "\n", sep = "")
    cat("  Accepted draws : ", dim(x$Ball)[3], " of ", dim(x$Ball)[3] + x$n_failed,
        " slots (", x$n_failed, " failed)\n", sep = "")
    cat("  Rotations tried: ", format(x$ndraws_tried, big.mark = ","),
        "   acceptance: ", sprintf("%.3f%%", 100 * x$accept_rate), "\n", sep = "")
    cat("  Horizons       : ", x$nsteps, "   bands: ",
        paste0(format(x$conf), collapse = ", "), "%\n", sep = "")
    if (x$narrative_active) cat("  Narrative restrictions: active\n")
    iv <- x[["iv"]]
    if (!is.null(iv)) {
        cat("  Instrument     : first-stage F = ", sprintf("%.2f", iv$fs_F),
            ", R2 = ", sprintf("%.3f", iv$fs_r2),
            ", n = ", iv$n_iv, "\n", sep = "")
        cat("  IV column      : fixed at OLS (VAR Toolbox convention)\n", sep = "")
    }
    invisible(x)
}


#' Percentile bands at one or more coverage levels
#'
#' Mirrors the percentile rules used inside \code{fSR_cpp} so that a level
#' recomputed here is identical to the one the C++ driver reports: plain
#' quantiles (type 7) when every draw carries the same weight, and the
#' mass-interpolating weighted percentile when the ADRR importance weights are
#' active.
#'
#' @keywords internal
#' @noRd
.fSR_bands <- function(res, conf) {

    if (is.null(res$IRall)) {
        if (length(conf) > 1L)
            stop("Several `conf` levels require `store_draws = TRUE`.")
        return(stats::setNames(
            list(list(IRinf = res$IRinf, IRsup = res$IRsup,
                      VDinf = res$VDinf, VDsup = res$VDsup)),
            format(conf)))
    }

    w        <- as.numeric(res$weights)
    weighted <- length(w) > 1L && stats::sd(w) > 0

    ## Weighted percentile of `nth_pct_w`: walk the sorted draws until the
    ## cumulative weight reaches the target mass, then interpolate inside the
    ## straddling cell.
    wq <- function(x, prob) {
        o     <- order(x)
        xs    <- x[o]
        ws    <- w[o]
        total <- sum(ws)
        if (!(total > 0)) return(xs[length(xs) %/% 2 + 1L])
        target <- prob / 100 * total
        cum    <- cumsum(ws)
        i      <- which(cum >= target)[1L]
        if (is.na(i)) return(xs[length(xs)])
        if (i == 1L || ws[i] <= 0) return(xs[i])
        frac <- (target - cum[i - 1L]) / ws[i]
        xs[i - 1L] + frac * (xs[i] - xs[i - 1L])
    }

    quant <- function(arr, prob) {
        d   <- dim(arr)
        mat <- matrix(arr, nrow = prod(d[1:3]), ncol = d[4])
        out <- if (weighted) apply(mat, 1L, wq, prob = prob)
               else          apply(mat, 1L, stats::quantile, probs = prob / 100,
                                   names = FALSE, type = 7)
        array(out, dim = d[1:3])
    }

    stats::setNames(lapply(conf, function(level) {
        lo <- (100 - level) / 2
        hi <- 100 - lo
        list(IRinf = quant(res$IRall, lo), IRsup = quant(res$IRall, hi),
             VDinf = quant(res$VDall, lo), VDsup = quant(res$VDall, hi))
    }), format(conf))
}
