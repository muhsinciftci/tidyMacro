#' Galí (1999) Data
#'
#' Quarterly US data on labour productivity and hours used to replicate
#' Galí (1999). Both series enter the VAR in first differences (growth rates).
#'
#' @format A data frame with 183 rows and 3 columns:
#' \describe{
#'   \item{Date}{Quarter end date (\code{Date}), 1948 Q1 – 1994 Q4}
#'   \item{Productivity}{Growth rate of labour productivity (GDP per hour,
#'     all businesses), in percent}
#'   \item{Hours}{Growth rate of total employee hours in non-agricultural
#'     establishments, in percent}
#' }
#'
#' @details
#' Data sourced from Fernald's quarterly TFP series
#' (\url{https://www.johnfernald.net/TFP}). Productivity is defined as
#' output per hour for all businesses; hours are total employee hours in
#' non-agricultural establishments.
#'
#' The identifying assumption (Blanchard-Quah long-run restriction) is that
#' only technology shocks can have a permanent effect on the level of labour
#' productivity.
#'
#' @references
#' Galí, J. (1999). Technology, employment, and the business cycle: Do
#' technology shocks explain aggregate fluctuations?
#' \emph{American Economic Review}, 89(1), 249–271.
#'
#' @examples
#' data("Gali1999")
#' head(Gali1999)
#' y <- Gali1999 |> dplyr::select(Productivity, Hours) |> as.matrix()
"Gali1999"


#' Bloom (2009) Data
#'
#' Monthly US data used to replicate the uncertainty shock VAR of Bloom (2009).
#' Data sourced from \url{https://nbloom.people.stanford.edu/research},
#' covering July 1962 – June 2013.
#'
#' @format A data frame with 612 rows and 9 columns:
#' \describe{
#'   \item{Date}{Month start date (\code{Date}), 1962-07 – 2013-06}
#'   \item{STOCK}{Log S&P 500 stock price index, in percent}
#'   \item{FFR}{Federal funds rate, in percent}
#'   \item{WAGE}{Log average hourly earnings (manufacturing), in percent}
#'   \item{CPI}{Log Consumer Price Index, in percent}
#'   \item{HOURSM}{Average weekly hours in manufacturing}
#'   \item{EMPM}{Log employment in manufacturing, in percent}
#'   \item{IPM}{Log industrial production (manufacturing), in percent}
#'   \item{UNCERT}{Uncertainty indicator: equals 1 in months coinciding with
#'     one of 17 identified uncertainty events (e.g. Cuban Missile Crisis,
#'     Black Monday), 0 otherwise.}
#' }
#'
#' @details
#' The uncertainty shock is identified via Cholesky ordering with
#' \code{UNCERT} ordered first. The 17 event dates follow Bloom (2009)
#' Table 1.
#'
#' @references
#' Bloom, N. (2009). The impact of uncertainty shocks.
#' \emph{Econometrica}, 77(3), 623–685.
#'
#' @examples
#' data("Bloom2009")
#' head(Bloom2009)
#' y <- Bloom2009 |> dplyr::select(-Date) |> as.matrix()
"Bloom2009"


#' Kaenzig (2021) Data
#'
#' Monthly data on global oil markets and US macroeconomic aggregates used
#' to replicate the OPEC production cut shock identification of Kaenzig (2021).
#' Includes a high-frequency external instrument based on oil futures price
#' changes around OPEC announcements.
#'
#' @format A data frame with rows and 8 columns:
#' \describe{
#'   \item{Date}{Month start date (\code{Date})}
#'   \item{Oil_Price}{Log real oil price, in percent}
#'   \item{World_Oil_Prod}{Log world oil production, in percent}
#'   \item{World_Oil_Inven}{Log world oil inventories, in percent}
#'   \item{World_IP}{Log world industrial production, in percent}
#'   \item{US_IP}{Log US industrial production, in percent}
#'   \item{US_CPI}{Log US Consumer Price Index, in percent}
#'   \item{iv_kanzig_final}{High-frequency external instrument: oil futures
#'     price change in a narrow window around OPEC announcements. Contains
#'     \code{NA} outside announcement dates.}
#' }
#'
#' @details
#' The structural shock is identified using the external instrument
#' \code{iv_kanzig_final} via the proxy-SVAR / IV approach. The instrument
#' isolates exogenous variation in oil supply caused by OPEC production
#' decisions.
#'
#' @references
#' Kaenzig, D. A. (2021). The macroeconomic effects of oil supply news:
#' Evidence from OPEC announcements.
#' \emph{American Economic Review}, 111(4), 1092–1125.
#'
#' @examples
#' data("Kaenzig2021")
#' head(Kaenzig2021)
#' y <- Kaenzig2021 |>
#'   dplyr::select(Oil_Price, World_Oil_Prod, World_Oil_Inven,
#'                 World_IP, US_IP, US_CPI) |>
#'   as.matrix()
"Kaenzig2021"


#' Beaudry & Portier (2014) Data
#'
#' Quarterly US data used to identify news shocks following Beaudry and
#' Portier (2014) and Barsky and Sims (2011/2012). The news shock is
#' identified as the shock with no contemporaneous effect on TFP that
#' maximises the long-run forecast error variance of TFP.
#'
#' @format A data frame with rows and 4 columns:
#' \describe{
#'   \item{logTFP}{Log total factor productivity (Fernald index, base ≈ 100)}
#'   \item{logSP500}{Log real S&P 500 stock price}
#'   \item{logConsumption}{Log real private consumption}
#'   \item{logHours}{Log total hours worked in the business sector}
#' }
#'
#' @details
#' Data sourced from the replication package of Beaudry and Portier (2014),
#' available at
#' \url{https://www.openicpsr.org/openicpsr/project/113916/version/V1/view}.
#' The VAR is estimated with \eqn{p = 2} lags and a constant. The news shock
#' is ordered second (after TFP) and identified via the max-share criterion
#' over a 40-quarter horizon.
#'
#' @references
#' Beaudry, P., & Portier, F. (2014). News-driven business cycles: Insights
#' and challenges. \emph{Journal of Economic Literature}, 52(4), 993–1074.
#'
#' Barsky, R. B., & Sims, E. R. (2012). Information, animal spirits, and the
#' meaning of innovations in consumer confidence.
#' \emph{American Economic Review}, 102(4), 1343–1377.
#'
#' @examples
#' data("BeaudryPortier2014")
#' head(BeaudryPortier2014)
#' y <- BeaudryPortier2014 |>
#'   dplyr::select(logTFP, logSP500, logConsumption, logHours) |>
#'   as.matrix()
"BeaudryPortier2014"
