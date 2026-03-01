#' Oil Market and Macroeconomic Data
#'
#' A dataset containing oil market variables and macroeconomic indicators
#' for Vector Autoregression (VAR) analysis, along with a high-frequency
#' external instrument for identification.
#'
#' @format A data frame with X rows and 7 columns:
#' 
#' @details
#' This dataset is designed for estimating structural Vector Autoregression
#' models of the global oil market. The high-frequency instrument can be used
#' for external instrument identification of oil supply or demand shocks.
#' 
#' @examples
#' # Load the data
#' data(kaenzig_data)
#' 
#' # View the first few rows
#' head(kaenzig_data)
#' 
#' # Summary statistics
#' summary(kaenzig_data)
"kaenzig_data"
