#' Custom ggplot2 Theme for tidyMacro Package
#'
#' A customized ggplot2 theme based on \code{theme_bw()} with Palatino font family,
#' adjusted grid lines, and centered plot elements. Designed for professional
#' publications and presentations of VAR model results.
#'
#' @param base_size Numeric. Base font size in points. Default is 14.
#'   All other text sizes are scaled relative to this value. The plot title
#'   will be 1.5 times the base size.
#'
#' @return A ggplot2 theme object that can be added to ggplot2 plots.
#'
#' The theme inherits from \code{theme_bw()} and modifies specific elements
#' to create a consistent visual style suitable for economic time series
#' and impulse response function plots.
#'
#' @note
#' The Palatino font family must be available on your system for proper rendering.
#' If Palatino is not available, ggplot2 will fall back to the default font family.
#' On most systems, you can check available fonts using \code{extrafont::fonts()}
#' after loading the extrafont package.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(tidyMacro)
#'
#' # Basic usage with default base size
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   labs(title = "Miles per Gallon vs Weight") +
#'   ftheme_tidyMacro()
#'
#' # Customize base size
#' ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
#'   geom_point(size = 3) +
#'   labs(title = "MPG vs Weight by Cylinders",
#'        x = "Weight (1000 lbs)",
#'        y = "Miles per Gallon") +
#'   scale_color_tidyMacro() +
#'   ftheme_tidyMacro(base_size = 16)
#' }
#'
#' @seealso
#' \code{\link[ggplot2]{theme_bw}}, \code{\link[ggplot2]{theme}},
#' \code{\link{scale_color_tidyMacro}}, \code{\link{set_tidyMacro_theme}}
#'
#' @importFrom ggplot2 theme_bw theme element_text element_line element_rect element_blank rel
#'
#' @export
ftheme_tidyMacro <- function(base_size = 14) {
    theme_bw() +
        theme(
            plot.title = element_text(
                hjust = 0.5,
                vjust = 2,
                size = base_size * 1.5,
                family = "Palatino"
            ),
            axis.title = element_text(family = "Palatino"),
            axis.title.x = element_text(
                hjust = 0.5,
                size = base_size,
                family = "Palatino"
            ),
            axis.title.y = element_text(
                hjust = 0.5,
                size = base_size,
                family = "Palatino"
            ),
            axis.line = element_line(colour = 'black', linewidth = 0.35),
            axis.text = element_text(family = "Palatino", size = base_size),
            
            # Facet strip (title) modifications
            strip.text = element_text(
                hjust = 0.5,
                colour = 'black',
                size = base_size
            ),
            strip.background = element_rect(fill = 'gray98', color = NA),
            panel.border = element_rect(colour = NA),
            panel.grid.major = element_line(
                linewidth = rel(0.15),
                colour = 'gray90'
            ),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            legend.position = 'bottom'
        )
}
