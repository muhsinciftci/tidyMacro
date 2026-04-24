\name{fLP}
\alias{fLP}
\title{Simple Local Projections with Formula Syntax}
\usage{
fLP(
  formula,
  data,
  horizons = 0:20,
  lags = 0L,
  vcov = c("nw", "iid"),
  nw_lag = NULL,
  confidence = 0.95,
  cumulative = FALSE
)
}
\arguments{
\item{formula}{A two-sided formula like \code{y ~ shock + x1 + x2}.}

\item{data}{A data frame containing the variables in \code{formula}.}

\item{horizons}{Non-negative integer vector of forecast horizons.}

\item{lags}{Number of lags of all variables appearing in the formula to add
to the regression. Default is \code{0}.}

\item{vcov}{Character string specifying the covariance estimator:
\code{"nw"} for Newey-West HAC or \code{"iid"} for homoskedastic OLS.}

\item{nw_lag}{Newey-West lag length. If \code{NULL} and \code{vcov = "nw"},
it defaults to \code{max(horizons, 1)}.}

\item{confidence}{Confidence level used for the reported
normal-approximation interval.}

\item{cumulative}{Logical; if \code{TRUE}, the dependent variable at horizon
\code{h} is \eqn{\sum_{j=0}^h y_{t+j}} instead of \eqn{y_{t+h}}.}
}
\value{
A list with an \code{irf} tibble containing one row per horizon and columns
for the point estimate, standard error, test statistic, p-value, confidence
interval, and number of observations used.
}
\description{
Estimates a sequence of local projections using a formula interface such as
\code{y ~ shock + x1 + x2}. The first right-hand-side term is treated as the
shock variable of interest, while any remaining terms are included as
controls. Optional lags are appended in R before estimation, and the
horizon-by-horizon regressions are then fit in C++.
}
\details{
This first implementation is intentionally simple:
\itemize{
\item the data are assumed to be ordered in time already;
\item the first RHS term must expand to exactly one model-matrix column;
\item missing values are not supported in the variables used by the model.
}
}
\examples{
\dontrun{
set.seed(123)
n <- 200
shock <- rnorm(n)
x <- rnorm(n)
y <- numeric(n)

for (t in 2:n) {
  y[t] <- 0.4 * y[t - 1] + 0.8 * shock[t] + 0.2 * x[t] + rnorm(1)
}

dat <- data.frame(y = y, shock = shock, x = x)
lp <- fLP(y ~ shock + x, data = dat, horizons = 0:8, lags = 2)
lp$irf
}
}
