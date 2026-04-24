#' Invertibility Test for External Instrument SVAR
#'
#' Tests whether the structural shock is invertible (fundamental) by regressing
#' the cleaned instrument on current and \eqn{r} future VAR residuals and
#' applying an F-test for the joint significance of the future residual
#' coefficients (Null = invertible).
#'
#' @param z Numeric vector. Cleaned instrument of length \eqn{T - p}.
#' @param eps \eqn{(T-p) \times N} matrix of VAR residuals.
#' @param r Integer. Number of future residual lags to include in the test.
#'
#' @return A list with:
#' \describe{
#'   \item{F_stat}{F statistic for the joint significance of future residuals.}
#'   \item{pval}{p-value. Small values reject invertibility.}
#'   \item{instr_relev}{Correlation between fitted and actual instrument.}
#'   \item{delta}{Coefficient vector from regressing \eqn{z} on
#'     \eqn{[1, \varepsilon_t, \varepsilon_{t+1}, \ldots, \varepsilon_{t+r}]}.
#'     Length \eqn{N(r+1) + 1} (intercept first).}
#'   \item{fit}{Fitted instrument values \eqn{\hat{\delta}(F)\hat{\varepsilon}_t}.
#'     Length \eqn{T - p - r}.}
#'   \item{vt}{Noise residuals \eqn{z - \hat{z}}. Same length as \code{fit}.}
#' }
#'
#' @details
#' The test regresses the cleaned instrument \eqn{z_t} on current and
#' \eqn{r} future VAR residuals:
#' \deqn{z_t = c + \delta_0 \hat{\varepsilon}_t + \delta_1 \hat{\varepsilon}_{t+1}
#'   + \cdots + \delta_r \hat{\varepsilon}_{t+r} + u_t.}
#' The F-test checks \eqn{\delta_1 = \cdots = \delta_r = 0}. Rejection implies
#' the shock is non-invertible and the generalised SVAR-IV procedure of
#' Forni, Gambetti & Ricco (2022) should be applied. The returned \code{delta},
#' \code{fit}, and \code{vt} are needed as inputs to \code{fBootstrapIVRecover}.
#'
#' @seealso \code{\link{fTestRecoverability}}, \code{\link{fBootstrapIVRecover}}
#'
#' @references
#' Forni, M., Gambetti, L., & Ricco, G. (2022). Recovering the structural shock
#' from non-invertible structural VARs. \emph{Review of Economics and Statistics}.
#'
#' @examples
#' \dontrun{
#' var_result <- fVAR(X, p = 7, c = 1)
#' z   <- clean_instrument(instr, p = 7, X)
#' inv <- fTestInvertibility(z, var_result$residuals, r = 3)
#' cat(sprintf("F = %.3f  p = %.4f\n", inv$F_stat, inv$pval))
#' }
#'
#' @export
fTestInvertibility <- function(z, eps, r) {
  N     <- ncol(eps)
  T_eps <- nrow(eps)

  z_f         <- head(z, -r)
  eps_cur_fut <- do.call(cbind, lapply(0:r, function(j) eps[seq_len(T_eps - r) + j, ]))
  X_ftest     <- cbind(1, eps_cur_fut)

  # Unrestricted model
  n    <- nrow(X_ftest)
  N1   <- ncol(X_ftest) - 1L
  r2_u <- fOLS(as.matrix(z_f), X_ftest[, -1L, drop = FALSE], c = 1L)$r2

  # Restricted model: drop future residual columns
  constr <- (N + 2L):ncol(X_ftest)
  keep   <- setdiff(seq_len(ncol(X_ftest)), c(1L, constr))
  r2_r   <- fOLS(as.matrix(z_f), X_ftest[, keep, drop = FALSE], c = 1L)$r2
  N2     <- length(keep)

  F_stat <- ((r2_u - r2_r) / (N1 - N2)) / ((1 - r2_u) / (n - N1 - 1))
  pval   <- pf(F_stat, N1 - N2, n - N1 - 1, lower.tail = FALSE)

  delta <- as.vector(solve(crossprod(X_ftest), crossprod(X_ftest, z_f)))
  fit   <- as.vector(X_ftest %*% delta)
  vt    <- as.vector(z_f) - fit

  list(
    F_stat      = F_stat,
    pval        = pval,
    instr_relev = cor(fit, as.vector(z_f)),
    delta       = delta,
    fit         = fit,
    vt          = vt
  )
}
