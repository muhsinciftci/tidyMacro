#' Plot Cholesky Impulse Response Functions with Confidence Bands
#'
#' Creates a grid of IRF plots showing responses of all variables to a single shock,
#' with shaded confidence bands using ggplot2. Designed for Cholesky-identified VARs.
#'
#' @param point Array of point estimates (N x N x horizon)
#' @param upper Array of upper confidence bounds (N x N x horizon)
#' @param lower Array of lower confidence bounds (N x N x horizon)
#' @param shock Integer, which shock to plot (1-indexed, typically 1 to N)
#' @param varnames Character vector of variable names (length N)
#' @param prc Numeric, confidence level percentage (e.g., 68 or 95)
#' @param return_data Logical, if TRUE returns tibble instead of plot
#' @param ribbon_fill Color for confidence band fill (default: "#407EC9")
#' @param ribbon_alpha Transparency for confidence band (default: 0.3)
#' @param line_color Color for point estimate line (default: "#910048")
#' @param zero_line_color Color for horizontal zero line (default: "#707372")
#' @param facet_scales Facet scales option: "free", "free_y", "free_x", or "fixed" (default: "free")
#' @param facet_ncol Number of columns in facet grid (default: NULL, auto-calculated)
#'
#' @return A ggplot object (if return_data = FALSE) or tibble (if return_data = TRUE)
#'
#' @details
#' This function extracts the impulse responses of all N variables to a specific shock
#' from 3D arrays (typically output from fbootstrapChol). The shock parameter selects
#' which column (shock) to visualize across all response variables.
#'
#' The plot uses the currently active ggplot2 theme. For a consistent look with other
#' tidyMacro plots, use \code{ftheme_tidyMacro()} or set it globally with
#' \code{set_tidyMacro_theme()}.
#'
#' @examples
#' \dontrun{
#' # Estimate VAR and compute Cholesky IRF
#' var_result <- fVAR(y, p = 2, c = 1)
#' wold <- fwoldIRF(var_result, horizon = 20)
#' S <- t(chol(var_result$sigma_full))
#' point_irf <- fcholeskyIRF(wold, S)
#' 
#' # Bootstrap confidence bands
#' boot_result <- fbootstrapChol(y, var_result, nboot = 1000, 
#'                               horizon = 20, prc = 68, bootscheme = "residual")
#' 
#' # Plot responses to first shock
#' fplotirf_chol(point_irf, boot_result$upper, boot_result$lower, 
#'               shock = 1, varnames = c("GDP", "Inflation", "Interest Rate"),
#'               prc = 68)
#' 
#' # Use custom theme
#' fplotirf_chol(point_irf, boot_result$upper, boot_result$lower, 
#'               shock = 1, varnames = c("GDP", "Inflation", "Interest Rate"),
#'               prc = 68) +
#'   ftheme_tidyMacro()
#' }
#'
#' @export
fplotirf_chol <- function(point, upper, lower, shock, varnames, prc,
                          return_data = FALSE,
                          ribbon_fill = "#407EC9",
                          ribbon_alpha = 0.3,
                          line_color = "#910048",
                          zero_line_color = "#707372",
                          facet_scales = "free",
                          facet_ncol = NULL) {
    
    # Check if required packages are loaded
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("tibble", quietly = TRUE)) {
        stop("Package 'tibble' is required. Please install it.", call. = FALSE)
    }
    
    # Get dimensions
    dims <- dim(point)
    N <- dims[1]
    horizon <- dims[3]
    
    # Validate shock index
    if (shock < 1 || shock > N) {
        stop(sprintf("shock must be between 1 and %d", N))
    }
    
    # Validate variable names
    if (length(varnames) != N) {
        stop(sprintf("varnames must have length %d (number of variables)", N))
    }
    
    # Calculate auto layout if not specified
    if (is.null(facet_ncol)) {
        facet_ncol <- ceiling(sqrt(N))
    }
    
    # Extract data for the specified shock across all response variables
    irf_list <- vector("list", N)
    
    for (i in 1:N) {
        irf_list[[i]] <- tibble::tibble(
            variable = varnames[i],
            horizon = 0:(horizon - 1),
            point = as.numeric(point[i, shock, ]),
            upper = as.numeric(upper[i, shock, ]),
            lower = as.numeric(lower[i, shock, ])
        )
    }
    
    irf_data <- do.call(rbind, irf_list)
    
    # Convert variable to factor to preserve order
    irf_data$variable <- factor(irf_data$variable, levels = varnames)
    
    # Return tibble if requested
    if (return_data) {
        return(tibble::as_tibble(irf_data))
    }
    
    # Create the plot
    p <- ggplot2::ggplot(irf_data, ggplot2::aes(x = horizon)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), 
                            fill = ribbon_fill, alpha = ribbon_alpha) +
        ggplot2::geom_line(ggplot2::aes(y = point), 
                          color = line_color, linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = 0, 
                           color = zero_line_color, 
                           linetype = "dashed", 
                           linewidth = 0.6) +
        ggplot2::facet_wrap(~ variable, scales = facet_scales, ncol = facet_ncol)
    
    return(p)
}
