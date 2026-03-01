#' Plot Impulse Response Functions with Confidence Bands
#'
#' Creates a grid of IRF plots with shaded confidence bands using ggplot2
#'
#' @param point Matrix of point estimates (N x horizon)
#' @param upper Matrix of upper confidence bounds (N x horizon)
#' @param lower Matrix of lower confidence bounds (N x horizon)
#' @param varnames Character vector of variable names (length N)
#' @param shockname Character string for shock name
#' @param prc Numeric, confidence level percentage (e.g., 68 or 95)
#' @param return_data Logical, if TRUE returns tibble instead of plot
#' @param ribbon_fill Color for confidence band fill (default: "#0072B2")
#' @param ribbon_alpha Transparency for confidence band (default: 0.3)
#' @param line_color Color for point estimate line (default: "black")
#' @param zero_line_color Color for horizontal zero line (default: "red")
#' @param facet_scales Facet scales option: "free", "free_y", "free_x", or "fixed" (default: "free")
#' @param facet_ncol Number of columns in facet grid (default: 2)
#'
#' @return A ggplot object (if return_data = FALSE) or tibble (if return_data = TRUE)
#'
#' @export
fplotirf_iv <- function(point, upper, lower, varnames, shockname, prc, 
                             return_data = FALSE,
                             ribbon_fill = "#407EC9",
                             ribbon_alpha = 0.3,
                             line_color = "#910048",
                             zero_line_color = "#707372",
                             facet_scales = "free",
                             facet_ncol = 2) {
    
    # Get dimensions
    N <- nrow(point)
    horizon <- ncol(point)
    
    # Create tibble for plotting
    irf_list <- vector("list", N)
    
    for (i in 1:N) {
        irf_list[[i]] <- tibble(
            variable = varnames[i],
            horizon = 0:(horizon - 1),
            point = as.numeric(point[i, ]),
            upper = as.numeric(upper[i, ]),
            lower = as.numeric(lower[i, ]),
            shock = shockname
        )
    }
    
    irf_data <- do.call(rbind, irf_list)
    
    # Convert variable to factor to preserve order
    irf_data$variable <- factor(irf_data$variable, levels = varnames)
    
    # Return tibble if requested
    if (return_data) {
        return(as_tibble(irf_data))
    }
    
    # Create the plot
    p <- ggplot(irf_data, aes(x = horizon)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), 
                    fill = ribbon_fill, alpha = ribbon_alpha) +
        geom_line(aes(y = point), color = line_color, linewidth = 0.8) +
        geom_hline(yintercept = 0, color = zero_line_color, linetype = "dashed", linewidth = 0.6) +
        facet_wrap(~ variable, scales = facet_scales, ncol = facet_ncol)
    
    return(p)
}
