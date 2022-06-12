#' Themes for all plots
#'
#' Change ggplot text, font, legend and border
#'
#' @param main_text_size Numeric, main text size
#' @param legend_text_size Numeric, legend text size
#' @param font_type Character, specify the plot text font family, default is "sans".
#' @param border_thick Numeric, border thickness, default is 1.
#' If set 0, remove both border and ticks.
#' @param remove_grid Logical, remove background grid lines, default is FALSE.
#' @param remove_border Logical, remove border line, default is FALSE.
#' @param remove_main_text Logical, remove all axis text, default is FALSE.
#' @param remove_legend_text Logical, remove all legend text, default is FALSE.
#' @param remove_legend Logical, remove entire legend, default is FALSE.
#' @return ggplot theme
#' @importFrom ggplot2 theme_bw theme margin unit element_text element_rect element_line
#' @importFrom rlang .data
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   plot_theme(font_type = "Times", border_thick = 2)
#' @export
plot_theme <- function(main_text_size = 8,
                       legend_text_size = 6,
                       font_type = "sans",
                       border_thick = 1.5,
                       remove_grid = TRUE,
                       remove_border = FALSE,
                       remove_main_text = FALSE,
                       remove_legend_text = FALSE,
                       remove_legend = FALSE) {
  ptheme <- theme_bw()

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
    axis.ticks.length = unit(.1, "cm")
  )

  # remove background grid line
  if (remove_grid) {
    bkg_theme <- theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  } else {
    bkg_theme <- NULL
  }

  # remove border line
  if (remove_border) {
    bod_theme <- theme(
      panel.border = element_blank(),
    )
  } else {
    bod_theme <- NULL
  }

  # remove axis text
  if (remove_main_text) {
    main_txt_theme <- theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
    )
  } else {
    main_txt_theme <- NULL
  }

  # remove legend text
  if (remove_legend_text) {
    legend_txt_theme <- theme(
      legend.title = element_blank(),
      legend.text = element_blank()
    )
  } else {
    legend_txt_theme <- NULL
  }

  # remove legend
  if (remove_legend) {
    leg_theme <- theme(
      legend.position = "none"
    )
  } else {
    leg_theme <- NULL
  }

  all_theme <- ptheme + font_theme + bkg_theme + bod_theme +
    main_txt_theme + legend_txt_theme + leg_theme

  return(all_theme)
}
