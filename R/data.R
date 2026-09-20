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


#' US Banking Deregulation and the Labor Share
#'
#' Annual state-level panel (46 US states, 1970–1996) on the private labor
#' share and staggered banking deregulation, from the replication package of
#' Dube, Girardi, Jordà and Taylor (LP-DiD). Treatment is absorbing.
#'
#' @format A data frame with 1,242 rows and 12 columns:
#' \describe{
#'   \item{id}{Integer state identifier}
#'   \item{state}{State abbreviation}
#'   \item{year}{Year, 1970–1996}
#'   \item{lshare}{Private labor share (compensation / private GSP)}
#'   \item{bank}{Inter-state banking deregulation in force (0/1)}
#'   \item{branch}{Intra-state branching deregulation in force (0/1)}
#'   \item{grgsp}{Growth rate of gross state product}
#'   \item{corptax}{State corporate tax rate}
#'   \item{unionmem}{Union membership share}
#'   \item{popgrowth}{Population growth}
#'   \item{lnunemp}{Log unemployment rate}
#'   \item{lnhpi}{Log house price index}
#' }
#'
#' @references
#' Dube, A., Girardi, D., Jordà, Ò., & Taylor, A. M. (2025). A Local
#' Projections Approach to Difference-in-Differences.
#' \emph{Journal of Applied Econometrics}, 40(7), 741–758.
#' \doi{10.1002/jae.70000}
#'
#' @examples
#' data("BankingDeregulation")
#' m <- fLPDID(lshare ~ bank, data = BankingDeregulation,
#'             panel_id = c("id", "year"), treat = "bank",
#'             post = 9, pre = 9)
"BankingDeregulation"


#' Democracy and Growth (ANRR / DDCG Panel)
#'
#' Annual country-level panel on democracy and GDP per capita from Acemoglu,
#' Naidu, Restrepo and Robinson (2019), as shipped in the Dube-Girardi-
#' Jordà-Taylor LP-DiD replication package. Treatment (democracy) is
#' non-absorbing: countries democratize and revert.
#'
#' @format A data frame with 5 columns:
#' \describe{
#'   \item{id}{Integer country identifier (World Bank code)}
#'   \item{country}{Country name}
#'   \item{year}{Year}
#'   \item{lgdp}{Log GDP per capita, multiplied by 100}
#'   \item{dem}{Democracy indicator (0/1; \code{NA} where undefined)}
#' }
#'
#' @references
#' Acemoglu, D., Naidu, S., Restrepo, P., & Robinson, J. A. (2019).
#' Democracy Does Cause Growth. \emph{Journal of Political Economy},
#' 127(1), 47–100.
#'
#' Dube, A., Girardi, D., Jordà, Ò., & Taylor, A. M. (2025). A Local
#' Projections Approach to Difference-in-Differences.
#' \emph{Journal of Applied Econometrics}, 40(7), 741–758.
#' \doi{10.1002/jae.70000}
#'
#' @examples
#' data("DemocracyGrowth")
#' m <- fLPDID(lgdp ~ dem + l(lgdp, 1:4), data = DemocracyGrowth,
#'             panel_id = c("id", "year"), treat = "dem",
#'             post = 30, pre = 20, nonabsorbing = TRUE, L = 20, ccc = 1)
"DemocracyGrowth"


#' Uhlig (2005) Data
#'
#' Monthly US data used to replicate the agnostic sign-restriction VAR of
#' Uhlig (2005), covering January 1965 - December 2003.
#'
#' @format A data frame with 468 rows and 7 columns:
#' \describe{
#'   \item{Date}{Month start date (\code{Date}), 1965-01 - 2003-12}
#'   \item{Real GDP}{Log real GDP, interpolated to monthly frequency}
#'   \item{GDP Deflator}{Log GDP deflator, interpolated to monthly frequency}
#'   \item{Commodity Price Idx.}{Log commodity price index}
#'   \item{Total Reserves}{Log total reserves}
#'   \item{Non-Borrowed Reserves}{Log non-borrowed reserves}
#'   \item{Fed. Funds Rate}{Federal funds rate, in percent}
#' }
#'
#' @details
#' The five log-level series are stored unscaled. Multiply them by 100 before
#' estimation so that impulse responses read in percent, as in the original
#' replication code; the funds rate is already in percent and must be left
#' alone.
#'
#' Identification is by sign restrictions only: a contractionary monetary
#' policy shock raises the funds rate and lowers the GDP deflator, commodity
#' prices and non-borrowed reserves for six months. Real GDP and total
#' reserves are deliberately left unrestricted - the response of output is
#' the object of the exercise, so it is not assumed.
#'
#' Data taken from the replication files of the VAR Toolbox of
#' Ambrogio Cesa-Bianchi (\url{https://github.com/ambropo/VAR-Toolbox}).
#'
#' @references
#' Uhlig, H. (2005). What are the effects of monetary policy on output?
#' Results from an agnostic identification procedure.
#' \emph{Journal of Monetary Economics}, 52(2), 381-419.
#' \doi{10.1016/j.jmoneco.2004.05.007}
#'
#' @seealso \code{\link{fSignRestr}} for the estimator, \code{\link{ADRR2018}}
#'   for the same variables over a longer sample.
#'
#' @examples
#' data("Uhlig2005")
#'
#' y <- Uhlig2005 |>
#'   dplyr::mutate(dplyr::across(-c(Date, `Fed. Funds Rate`), \(x) 100 * x)) |>
#'   dplyr::select(-Date) |>
#'   as.matrix()
#'
#' SIGN <- matrix(0, 6, 6)
#' SIGN[, 1] <- c(0, -1, -1, 0, -1, 1)
#'
#' fit <- fSignRestr(y, p = 12, c = 1, sign = SIGN, nsteps = 60,
#'                   ndraws = 100, sr_hor = 6, seed = 42)
"Uhlig2005"


#' Antolin-Diaz & Rubio-Ramirez (2018) Data
#'
#' Monthly US data used to replicate the narrative sign-restriction VAR of
#' Antolin-Diaz and Rubio-Ramirez (2018), covering January 1965 - November
#' 2007. The variables are those of \code{\link{Uhlig2005}}; the sample runs
#' four years longer.
#'
#' @format A data frame with 515 rows and 7 columns:
#' \describe{
#'   \item{Date}{Month start date (\code{Date}), 1965-01 - 2007-11}
#'   \item{Real GDP}{Log real GDP, interpolated to monthly frequency}
#'   \item{GDP Deflator}{Log GDP deflator, interpolated to monthly frequency}
#'   \item{Commodity Price Idx.}{Log commodity price index}
#'   \item{Total Reserves}{Log total reserves}
#'   \item{Non-Borrowed Reserves}{Log non-borrowed reserves}
#'   \item{Fed. Funds Rate}{Federal funds rate, in percent}
#' }
#'
#' @details
#' The paper keeps Uhlig's sign restrictions and adds two narrative
#' constraints tied to the Volcker announcement of October 1979: the monetary
#' policy shock was positive that month (a Type 1, sign restriction), and it
#' was the dominant contributor to the unexpected move in the funds rate (a
#' Type 2, dominance restriction). Both are imposed by
#' \code{\link{fSignRestr}} through its \code{narrative} argument.
#'
#' Estimation follows footnote 8 of the paper and uses no constant, the data
#' having been demeaned. Unlike \code{\link{Uhlig2005}}, the log-level series
#' are used as stored, without rescaling by 100.
#'
#' Data taken from the replication files of the VAR Toolbox of
#' Ambrogio Cesa-Bianchi (\url{https://github.com/ambropo/VAR-Toolbox}).
#'
#' @references
#' Antolin-Diaz, J., & Rubio-Ramirez, J. F. (2018). Narrative Sign
#' Restrictions for SVARs. \emph{American Economic Review}, 108(10),
#' 2802-2829. \doi{10.1257/aer.20161852}
#'
#' @seealso \code{\link{fSignRestr}}, \code{\link{Uhlig2005}}
#'
#' @examples
#' data("ADRR2018")
#'
#' y <- ADRR2018 |> dplyr::select(-Date) |> as.matrix()
#'
#' SIGN <- matrix(0, 6, 6)
#' SIGN[, 1] <- c(0, -1, -1, 0, -1, 1)
#'
#' # Volcker, October 1979: positive MP shock that dominates the funds-rate
#' # forecast error (variable 6).
#' narr <- list(
#'   sign = list(shock = 1, period = as.Date("1979-10-01"), sign = 1),
#'   dom  = list(shock = 1, period = as.Date("1979-10-01"), var = 6)
#' )
#'
#' fit <- fSignRestr(y, p = 12, c = 0, sign = SIGN, nsteps = 60, ndraws = 100,
#'                   sr_hor = 6, narrative = narr, dates = ADRR2018$Date,
#'                   seed = 42)
"ADRR2018"


#' Gertler & Karadi (2015) Data with External Instrument
#'
#' Monthly US data and the high-frequency monetary policy surprise used to
#' replicate the proxy-SVAR of Gertler and Karadi (2015), covering July 1979 -
#' June 2012.
#'
#' @format A data frame with 396 rows and 6 columns:
#' \describe{
#'   \item{Date}{Month start date (\code{Date}), 1979-07 - 2012-06}
#'   \item{Consumer Price Index}{Log CPI}
#'   \item{Industrial Production}{Log industrial production}
#'   \item{1-year T-Bill}{One-year Treasury rate, in percent}
#'   \item{EBP}{Excess bond premium of Gilchrist and Zakrajsek (2012), in
#'     percentage points}
#'   \item{FF4}{Three-month-ahead fed funds futures surprise in a 30-minute
#'     window around FOMC announcements, cumulated over the month, in
#'     percentage points. \code{NA} before January 1990, where the futures
#'     series does not exist.}
#' }
#'
#' @details
#' This is a different object from \code{\link{GK2015}}, which holds a
#' four-variable system on the shorter 1990m1 sample and carries no
#' instrument. \code{GK2015_2} adds the \code{FF4} proxy and starts in
#' 1979m7, so the VAR is estimated on the full sample while the first stage
#' uses only the 270 months where the instrument is observed.
#'
#' The columns are stored in source order. \code{\link{fSignRestr}}
#' instruments the residual of the \emph{first} endogenous variable, so the
#' policy rate must be reordered to the front before estimation - see the
#' example. Note that 49 of the 270 observed \code{FF4} values are genuine
#' zeros (months with no FOMC announcement surprise), not missing data.
#'
#' Data taken from the replication files of the VAR Toolbox of
#' Ambrogio Cesa-Bianchi (\url{https://github.com/ambropo/VAR-Toolbox}), where
#' missing instrument values are flagged with the sentinel \code{123456789};
#' they are stored here as \code{NA}.
#'
#' @references
#' Gertler, M., & Karadi, P. (2015). Monetary Policy Surprises, Credit Costs,
#' and Economic Activity. \emph{American Economic Journal: Macroeconomics},
#' 7(1), 44-76. \doi{10.1257/mac.20130329}
#'
#' Gilchrist, S., & Zakrajsek, E. (2012). Credit Spreads and Business Cycle
#' Fluctuations. \emph{American Economic Review}, 102(4), 1692-1720.
#' \doi{10.1257/aer.102.4.1692}
#'
#' @seealso \code{\link{fSignRestr}}, \code{\link{GK2015}}
#'
#' @examples
#' data("GK2015_2")
#'
#' # Policy rate first: the estimator instruments the first variable.
#' endo <- GK2015_2 |>
#'   dplyr::select(`1-year T-Bill`, `Consumer Price Index`,
#'                 `Industrial Production`, EBP) |>
#'   as.matrix()
#'
#' iv <- GK2015_2 |> dplyr::select(FF4) |> as.matrix()
#'
#' # Shock 1 is instrument-identified, so `sign` has k - 1 columns.
#' SIGN <- matrix(0, 4, 3)
#' SIGN[, 1] <- c(0, -1, -1, 0)
#'
#' fit <- fSignRestr(endo, p = 12, c = 1, sign = SIGN, nsteps = 48,
#'                   ndraws = 100, sr_hor = 3, instrument = list(Z = iv),
#'                   seed = 42)
"GK2015_2"


#' Gertler & Karadi (2015) Data, Forni-Gambetti-Ricco Sample
#'
#' Monthly US data and the Gertler-Karadi high-frequency monetary policy
#' surprise, restricted to the January 1990 - June 2012 window over which the
#' instrument is observed. This is the sample used by Forni, Gambetti and
#' Ricco (2024) to test invertibility and recoverability.
#'
#' @format A data frame with 270 rows and 5 columns:
#' \describe{
#'   \item{Date}{Month start date (\code{Date}), 1990-01 - 2012-06}
#'   \item{BY1}{One-year Treasury rate, in percent}
#'   \item{CPI_Inflation}{Monthly CPI inflation, in percent}
#'   \item{IP_Growth}{Monthly industrial production growth, in percent}
#'   \item{instr}{Three-month-ahead fed funds futures surprise (FF4) in a
#'     30-minute window around FOMC announcements, in percentage points}
#' }
#'
#' @details
#' Distinct from \code{\link{GK2015_2}}, which starts in 1979m7, keeps the
#' price and quantity series in log levels, adds the excess bond premium, and
#' leaves \code{FF4} as \code{NA} before 1990. Use this one for the
#' invertibility and recoverability tests, which need a balanced sample of
#' macro data and instrument; use \code{GK2015_2} for the proxy-SVAR, which
#' estimates the reduced form on the longer sample.
#'
#' The 49 zero values in \code{instr} are genuine - months without an FOMC
#' announcement surprise - not missing observations.
#'
#' @references
#' Gertler, M., & Karadi, P. (2015). Monetary Policy Surprises, Credit Costs,
#' and Economic Activity. \emph{American Economic Journal: Macroeconomics},
#' 7(1), 44-76. \doi{10.1257/mac.20130329}
#'
#' Forni, M., Gambetti, L., & Ricco, G. (2024). External Instrument SVAR
#' Analysis for Noninvertible Shocks. \emph{Journal of Applied Econometrics},
#' 39(7), 1173-1193. \doi{10.1002/jae.3072}
#'
#' @seealso \code{\link{fTestInvertibility}}, \code{\link{fTestRecoverability}},
#'   \code{\link{GK2015_2}}
#'
#' @examples
#' data("GK2015")
#'
#' X     <- GK2015 |> dplyr::select(BY1, CPI_Inflation, IP_Growth) |> as.matrix()
#' instr <- GK2015$instr
"GK2015"


#' Arias, Caldara & Rubio-Ramirez (2019) Data with the Romer-Romer Instrument
#'
#' Monthly US data for the six-variable monetary policy VAR of Arias, Caldara
#' and Rubio-Ramirez (2019), together with the Romer and Romer (2004) narrative
#' monetary policy shock used as an external instrument. Covers January 1965 -
#' December 2007.
#'
#' @format A data frame with 516 rows and 8 columns:
#' \describe{
#'   \item{Date}{Month start date (\code{Date}), 1965-01 - 2007-12}
#'   \item{Real GDP}{Monthly real GDP, level (interpolated to monthly)}
#'   \item{GDP Deflator}{Monthly GDP deflator, level}
#'   \item{Commodity Price Idx.}{Commodity price index, level}
#'   \item{Total Reserves}{Total reserves, level}
#'   \item{Non-Borrowed Reserves}{Non-borrowed reserves, level}
#'   \item{Fed. Funds Rate}{Federal funds rate, in percent}
#'   \item{RR}{Romer-Romer narrative monetary policy shock, in percentage
#'     points. \code{NA} before January 1969, where the series does not exist.}
#' }
#'
#' @details
#' The first six series are the same variables as \code{\link{Uhlig2005}} and
#' \code{\link{ADRR2018}}, which makes the three directly comparable. They are
#' stored as levels rather than logs: the replication code takes
#' \code{100 * log()} of all but the funds rate before estimation, and that
#' transformation is left explicit rather than baked into the object.
#'
#' \code{RR} is a narrative instrument, so the 118 exact zeros inside its
#' observed span are genuine months without an identified policy change, not
#' missing data. The 48 leading \code{NA} values form one contiguous block, so
#' the instrument aligns with a single sub-sample of the VAR.
#'
#' Obtained from the replication package of Braun and Bruggemann (2023)
#' (\url{https://github.com/r-a-braun/SVAR-IVSR}), file
#' \code{dataset_ACR.mat}, which carries the Arias, Caldara and Rubio-Ramirez
#' dataset. Note that Braun and Bruggemann identify the model with zero
#' restrictions and sign restrictions on both \eqn{B} and \eqn{A_0 = B^{-1}} in
#' a proxy-augmented Bayesian setup; \code{\link{fSignRestr}} implements the
#' different scheme of Cesa-Bianchi and Sokol (2022), so this object supports
#' their data, not a replication of their results.
#'
#' @references
#' Arias, J. E., Caldara, D., & Rubio-Ramirez, J. F. (2019). The systematic
#' component of monetary policy in SVARs: An agnostic identification procedure.
#' \emph{Journal of Monetary Economics}, 101, 1-13.
#' \doi{10.1016/j.jmoneco.2018.07.011}
#'
#' Braun, R., & Bruggemann, R. (2023). Identification of SVAR Models by
#' Combining Sign Restrictions With External Instruments. \emph{Journal of
#' Business & Economic Statistics}, 41(4), 1077-1089.
#' \doi{10.1080/07350015.2022.2104857}
#'
#' Romer, C. D., & Romer, D. H. (2004). A New Measure of Monetary Shocks:
#' Derivation and Implications. \emph{American Economic Review}, 94(4),
#' 1055-1084. \doi{10.1257/0002828042002651}
#'
#' @seealso \code{\link{fSignRestr}}, \code{\link{Uhlig2005}},
#'   \code{\link{ADRR2018}}
#'
#' @examples
#' data("ACR2019")
#'
#' y <- ACR2019 |>
#'   dplyr::mutate(dplyr::across(-c(Date, `Fed. Funds Rate`, RR),
#'                               \(x) 100 * log(x))) |>
#'   dplyr::select(-Date, -RR) |>
#'   as.matrix()
#'
#' # Shock 1 is pinned by the instrument, so `sign` has k - 1 columns.
#' SIGN <- matrix(0, 6, 5)
#' SIGN[, 1] <- c(-1, -1, 0, 0, 0, -1)   # a conventional demand shock
#'
#' fit <- fSignRestr(y, p = 12, c = 0, sign = SIGN, nsteps = 48,
#'                   ndraws = 100, sr_hor = 6,
#'                   instrument = list(Z = as.matrix(ACR2019$RR)), seed = 42)
"ACR2019"
