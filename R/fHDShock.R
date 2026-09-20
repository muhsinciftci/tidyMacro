#' Historical Decomposition of a Structural VAR
#'
#' Splits each observed series into the contribution of every structural shock
#' plus the deterministic and initial-condition terms. The components sum back
#' to the data over the estimation sample.
#'
#' @param x Either an \code{"fSignRestr"} object from \code{\link{fSignRestr}},
#'   or a T x N numeric matrix of endogenous variables. When an object is
#'   supplied, \code{beta}, \code{p} and \code{c} are taken from it and only
#'   \code{y} or \code{B} need overriding.
#' @param B Structural impact matrix. Defaults to the Fry-Pagan draw
#'   \code{Bfp} when \code{x} is an \code{"fSignRestr"} object, since that is a
#'   genuine structural representation; the element-wise median \code{Bmed} is
#'   not, and its components would not sum back to the data.
#' @param y Optional T x N matrix of endogenous variables, required when
#'   \code{x} is not a matrix and the data are not recoverable from it.
#' @param beta Reduced-form coefficient matrix, \code{[const | lags | exog]}.
#' @param p Integer lag order.
#' @param c Integer intercept indicator, 1 or 0.
#' @param exog Optional T x M matrix of exogenous regressors.
#'
#' @return A list with \code{shock} (a T x N x N array: time, variable, shock),
#'   \code{init}, \code{const}, \code{exo} and \code{endo}, the last being the
#'   reconstructed data.
#'
#' @seealso \code{\link{fSignRestr}}, \code{\link{fPlotHistDec}}
#'
#' @examples
#' data("Uhlig2005")
#' y <- Uhlig2005 |>
#'   dplyr::mutate(dplyr::across(-c(Date, `Fed. Funds Rate`), \(v) 100 * v)) |>
#'   dplyr::select(-Date) |> as.matrix()
#'
#' SIGN <- matrix(0, 6, 6)
#' SIGN[, 1] <- c(0, -1, -1, 0, -1, 1)
#'
#' fit <- fSignRestr(y, p = 12, c = 1, sign = SIGN, nsteps = 24,
#'                   ndraws = 100, sr_hor = 6, seed = 42)
#' hd <- fHDShock(fit, y = y)
#'
#' # Components reproduce the data over the estimation sample.
#' max(abs(hd$endo[-(1:12), ] - y[-(1:12), ]))
#'
#' @export
fHDShock <- function(x, B = NULL, y = NULL, beta = NULL, p = NULL, c = NULL,
                     exog = NULL) {

    if (inherits(x, "fSignRestr")) {
        if (is.null(B))    B    <- x$Bfp
        if (is.null(beta)) beta <- x$var$beta
        if (is.null(p))    p    <- x$p
        if (is.null(c))    c    <- x$c
        if (is.null(y))    y    <- x$var$y
        if (is.null(y))
            stop("`y` must be supplied: the fitted object does not carry the data.")
    } else {
        y <- x
    }

    y <- as.matrix(y)
    if (!is.numeric(y)) stop("`y` must be numeric.")
    if (is.null(B) || is.null(beta) || is.null(p) || is.null(c))
        stop("`B`, `beta`, `p` and `c` are all required.")

    B <- as.matrix(B)
    if (nrow(B) != ncol(B) || nrow(B) != ncol(y))
        stop("`B` must be square with one row per variable.")

    fHDShock_cpp(y = y, beta = as.matrix(beta), B = B,
                 p = as.integer(p), c = as.integer(c), exog = exog)
}
