##' ggplot theme
##'
##' @title plot_theme
##' @param font.size font size
##' @param xleft x axis left limit
##' @param xright x axis right limit
##' @param yleft y axis left limit
##' @param yright y axis right limit
##' @return ggplot theme
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 ylim
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @examples
##' library(ggplot2)
##' qplot(1:10) + plot_theme()
##' @export
plot_theme <- function(text_size=14) {
  theme_bw()+
    theme(axis.text.x = element_text(colour = "black",
                                     size = text_size, vjust =1 ),
          axis.text.y = element_text(colour = "black",
                                     size = text_size, hjust =1 ),
          axis.title = element_text(margin=margin(10, 5, 0, 0),
                                    color = "black",
                                    size = text_size),
          axis.title.y = element_text(angle=90)
    )
}

