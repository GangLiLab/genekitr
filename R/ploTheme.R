#' Themes for all plots
#'
#' Change ggplot text, font, legend and border
#'
#' @param theme_type ggplot theme,
#' @param main_text_size Numeric, main text size
#' @param legend_text_size Numeric, legend text size
#' @param font_type Character, specify the plot text font family, default is "sans".
#' @param border_thick Numeric, border thickness, default is 1.
#' @param remove_grid Logical, remove background grid lines, default is FALSE.
#' @param remove_border Logical, remove border line, default is FALSE.
#' @param remove_text Logical, remove all text, default is FALSE.
#' @param remove_legend Logical, remove legend, default is FALSE.
#' @return ggplot theme
#' @importFrom ggplot2 theme_bw theme margin unit element_text element_rect element_line
#' @importFrom rlang .data
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(x=wt, y=mpg))+ geom_point()+
#'   plot_theme(theme_type = 'bw', font_type = 'Times', border_thick = 2)
#' @export
plot_theme <- function(theme_type = c('bw','classic'),
                       main_text_size = 8,
                       legend_text_size = 6,
                       font_type = "sans",
                       border_thick = 1.5,
                       remove_grid = TRUE,
                       remove_border = FALSE,
                       remove_text = FALSE,
                       remove_legend = FALSE
) {

  theme_type <- match.arg(theme_type)

  ptheme <- if(theme_type == 'bw'){
    theme_bw()
  }else{
    theme_classic()
  }

  font_theme <- theme(
    axis.text.x = element_text(
      color = "black", size = main_text_size,
      vjust = 1, family = font_type
    ),
    axis.text.y = element_text(
      color = "black", size = main_text_size,
      hjust = 1, family = font_type
    ),
    axis.title.x = element_text(
      color = "black", size = main_text_size,
      margin = margin(8, 6, 0, 0),
      family = font_type
    ),
    axis.title.y = element_text(angle = 90, size = main_text_size, family = font_type),
    legend.title = element_text(
      color = "black",
      size = legend_text_size, family = font_type
    ),
    legend.text = element_text(
      color = "black",
      size = legend_text_size, family = font_type
    ),
    panel.border = element_rect(colour = "black", size = border_thick),
    axis.ticks = element_line(colour = "black", size = as.numeric(border_thick) / 3),
    axis.ticks.length = unit(.1, "cm"),
  )

  # remove background grid line
  if (remove_grid) {
    bkg_theme <- theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }else{
    bkg_theme = NULL
  }

  # remove border line
  if (remove_border) {
    bod_theme <- theme(
      panel.border = element_blank(),
    )
  }else{
    bod_theme = NULL
  }

  # remove axis text
  if (remove_text) {
    txt_theme <- theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
    )
  }else{
    txt_theme = NULL
  }

  # remove legend
  if (remove_legend) {
    leg_theme <- theme(
      legend.position = "none"
    )
  }else{
    leg_theme <- NULL
  }

  all_theme <- ptheme + font_theme + bkg_theme + bod_theme + txt_theme + leg_theme

  return(all_theme)
}
