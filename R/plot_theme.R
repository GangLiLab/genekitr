##' ggplot theme
##'
##' @title plot_theme
##' @param main_text_size numeric, main text font size.
##' @param legend_text_size numeric, legend text font size.
##' @param font_type font family.
##' @return ggplot theme
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 ylim
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @examples
##' library(ggplot2)
##' ggplot(mtcars,aes(x=mpg)) + geom_histogram(binwidth=5) + plot_theme()
##' @export
plot_theme <- function(main_text_size=14,
                       legend_text_size=10,
                       font_type = 'Arial',
                       ...) {
  theme_bw()+
    theme(axis.text.x = element_text(color = "black", size = main_text_size,
                                     vjust =1 ,family = font_type),
          axis.text.y = element_text(color = "black", size = main_text_size,
                                     hjust =1 ,family = font_type),
          axis.title.x = element_text(color = "black", size = main_text_size, margin=margin(8, 6, 0, 0),
                                      family = font_type),
          axis.title.y = element_text(angle=90, size = main_text_size, family = font_type),
          legend.title=element_text(color = "black", size = legend_text_size, family = font_type),
          legend.text=element_text(color = "black", size = legend_text_size, family = font_type)
    )
}

