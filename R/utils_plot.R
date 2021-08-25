#--- ggplot theme settting ---#
plot_theme <- function(main_text_size=14,
                       legend_text_size=10,
                       font_type = 'Arial',
                       border_thick = 1,
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
          legend.text=element_text(color = "black", size = legend_text_size, family = font_type),
          panel.border = element_rect(colour = "black", size= border_thick),
          axis.ticks = element_line(colour = "black", size = border_thick/3),
          axis.ticks.length=unit(.1, "cm"),
    )
}

#--- wrap text if too long ---#
text_wraper <- function(width) {
  function(x) {
    lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse="\n")
  }
}


