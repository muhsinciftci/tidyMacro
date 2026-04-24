#' Recoverability Test for External Instrument SVAR
#'
#' Tests whether a non-invertible structural shock is recoverable by applying
#' a Ljung-Box autocorrelation test to the fitted instrument from the
#' invertibility regression (Null = recoverable, i.e. no autocorrelation).
#'
#' @param fit Numeric vector of fitted instrument values returned by
#'   \code{fTestInvertibility} (\code{$fit}).
#' @param lags Integer. Number of lags for the Ljung-Box test. Default \code{2L}.
#'
#' @return A list with:
#' \describe{
#'   \item{statistic}{Ljung-Box Q statistic.}
#'   \item{pval}{p-value. Large values fail to reject recoverability.}
#'   \item{lags}{Number of lags used.}
#' }
#'
#' @details
#' Implements the recoverability test of Forni, Gambetti & Ricco (2022). After
#' rejecting invertibility via \code{fTestInvertibility}, the fitted instrument
#' \eqn{\hat{z}_t = \hat{\delta}(F)\hat{\varepsilon}_t} is tested for
#' autocorrelation. Under recoverability the fitted series has no serial
#' correlation (it is a white noise projection); under non-recoverability it
#' is autocorrelated. Failure to reject the null means absolute IRFs, historical
#' decompositions, and spectral FEVDs can all be recovered via the generalised
#' SVAR-IV procedure.
#'
#' @seealso \code{\link{fTestInvertibility}}, \code{\link{fBootstrapIVRecover}}
#'
#' @references
#' Forni, M., Gambetti, L., & Ricco, G. (2022). Recovering the structural shock
#' from non-invertible structural VARs. \emph{Review of Economics and Statistics}.
#'
#' @examples
#' \dontrun{
#' inv <- fTestInvertibility(z, var_result$residuals, r = 3)
#' rec <- fTestRecoverability(inv$fit, lags = 2)
#' cat(sprintf("Ljung-Box p-value (lag %d): %.4f\n", rec$lags, rec$pval))
#' }
#'
#' @export
fTestRecoverability <- function(fit, lags = 2L) {
  box <- Box.test(fit, lag = lags, type = "Ljung-Box")
  list(
    statistic = unname(box$statistic),
    pval      = box$p.value,
    lags      = lags
  )
}
