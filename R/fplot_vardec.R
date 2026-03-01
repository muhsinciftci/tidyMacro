#' Plot Forecast Error Variance Decomposition
#'
#' Creates area plots showing the forecast error variance decomposition (FEVD)
#' for multiple variables across horizons.
#'
#' @param fevd A numeric array containing forecast error variance shares.
#'   Either a 3-dimensional array (N x N x horizon) for multiple shocks,
#'   or a 2-dimensional matrix (N x horizon) for a single shock.
#' @param varnames Character vector of variable names. Length should equal N
#'   (the number of variables).
#' @param shocknames Character vector of shock names. For 3-dimensional fevd,
#'   length should equal N. For 2-dimensional fevd, should be a single name.
#' @param return_df Logical. If \code{TRUE}, returns a data frame instead of
#'   a plot. Default is \code{FALSE}.
#'
#' @return If \code{return_df = FALSE}, returns a \code{ggplot} object showing
#'   the FEVD for each variable in separate facets. If \code{return_df = TRUE},
#'   returns a data frame in long format with columns: horizon, variable, shock,
#'   and value.
#'
#' @details
#' The function creates stacked area charts showing how forecast error variance
#' is decomposed across different shocks over the forecast horizon. Each variable
#' is displayed in a separate facet panel.
#'
#' For 2-dimensional input (single shock), the function automatically calculates
#' and displays "Other Shocks" as the complement to ensure the decomposition
#' sums to 1.
#'
#' @importFrom ggplot2 ggplot aes geom_area facet_wrap coord_cartesian labs theme_minimal
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr bind_rows
#' @importFrom rlang .data
#'
#' @export
fplot_vardec <- function(fevd, varnames, shocknames, return_df = FALSE) {
    
    # Input validation
    if (!is.numeric(fevd)) {
        stop("fevd must be a numeric array or matrix")
    }
    
    if (!is.character(varnames)) {
        stop("varnames must be a character vector")
    }
    
    if (!is.character(shocknames)) {
        stop("shocknames must be a character vector")
    }
    
    # Determine dimensions
    if (length(dim(fevd)) == 3) {
        dims <- dim(fevd)
        N <- dims[1]
        N_shocks <- dims[2]
        horizon <- dims[3]
        multiple_shocks <- TRUE
        
        if (length(varnames) != N) {
            stop(paste0("Length of varnames (", length(varnames), 
                        ") must equal first dimension of fevd (", N, ")"))
        }
        
        if (length(shocknames) != N_shocks) {
            stop(paste0("Length of shocknames (", length(shocknames), 
                        ") must equal second dimension of fevd (", N_shocks, ")"))
        }
        
    } else if (length(dim(fevd)) == 2) {
        dims <- dim(fevd)
        N <- dims[1]
        horizon <- dims[2]
        multiple_shocks <- FALSE
        
        if (length(varnames) != N) {
            stop(paste0("Length of varnames (", length(varnames), 
                        ") must equal first dimension of fevd (", N, ")"))
        }
        
        if (length(shocknames) != 1) {
            stop("For 2D fevd, shocknames must be a single character string")
        }
        
    } else {
        stop("fevd must be a 2D matrix or 3D array")
    }
    
    # Build dataframe
    df_list <- list()
    
    for (i in seq_len(N)) {
        if (multiple_shocks) {
            # Extract fevd[i, , ] and transpose to horizon x N_shocks
            fevd_var <- t(fevd[i, , ])
        } else {
            # fevd[i, ] plus other shocks
            fevd_var <- cbind(fevd[i, ], 1 - fevd[i, ])
        }
        
        # Create dataframe for this variable
        temp_df <- as.data.frame(fevd_var)
        
        if (multiple_shocks) {
            colnames(temp_df) <- shocknames
        } else {
            colnames(temp_df) <- c(shocknames, "Other Shocks")
        }
        
        temp_df$horizon <- seq(0, horizon - 1)
        temp_df$variable <- varnames[i]
        
        df_list[[i]] <- temp_df
    }
    
    # Combine all dataframes
    df <- dplyr::bind_rows(df_list)
    
    # Reshape to long format
    df_long <- tidyr::pivot_longer(
        df, 
        cols = -c("horizon", "variable"), 
        names_to = "shock", 
        values_to = "value"
    )
    
    # Convert to factor for proper ordering
    df_long$variable <- factor(df_long$variable, levels = varnames)
    
    if (multiple_shocks) {
        df_long$shock <- factor(df_long$shock, levels = shocknames)
    } else {
        df_long$shock <- factor(df_long$shock, levels = c("Other Shocks", shocknames))
    }
    
    if (return_df) {
        return(df_long)
    }
    
    # Create plot
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$horizon, 
                                               y = .data$value, 
                                               fill = .data$shock)) +
        ggplot2::geom_area() +
        ggplot2::facet_wrap(~ .data$variable, nrow = 1) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(y = "% FEVD", x = "Horizon", fill = "Shock")
    
    return(p)
}
