#' Venn plot for groups of genes
#'
#' If gene group over 4, plot will be visulized using UpSet plot.
#'
#' @param venn_list A list of gene id.
#' @param color Colors for gene lists, default is NULL.
#' @param alpha_degree Alpha transparency of each circle's area, default is 0.3.
#' @param text_size Text size, default is 1.
#' @param remove_grid Logical, remove circle or grid lines, default is `FALSE`.
#' @param use_venn Logical, use venn to plot, default is `TRUE`, the other
#'   option is upsetplot for large list.
#' @inheritParams plot_theme
#' @return  A ggplot object
#' @importFrom VennDiagram venn.diagram
#' @importFrom dplyr as_tibble filter select group_by summarize %>%
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot geom_bar aes geom_text after_stat theme
#'   element_blank scale_y_continuous
#' @export
#' @examples
#' library(ggplot2)
#' set1 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
#' set2 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
#' set3 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
#' set4 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
#' set5 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
#' sm_gene_list = list(gset1 = set1, gset2 = set2, gset3 = set3)
#' la_gene_list = list(gset1 = set1, gset2 = set2, gset3 = set3,
#'   gset4 = set4, gset5 = set5 )
#' plotVenn(sm_gene_list,text_size = 1.5,alpha_degree = 1,
#'   remove_grid = TRUE,color = ggsci::pal_lancet()(3))
#' plotVenn(la_gene_list,text_size = 15,alpha_degree = 0.2,border_thick = 2,
#' remove_grid = TRUE, use_venn = FALSE)
#'

plotVenn <- function(venn_list,
                     color = NULL,
                     alpha_degree = 0.3,
                     text_size = 1,
                     remove_grid = FALSE,
                     use_venn = TRUE,
                     main_text_size = 10,
                     legend_text_size = 8,
                     font_type = 'Arial',
                     border_thick = 1) {

  #--- args ---#
  stopifnot(is.list(venn_list))
  # if gene list too long, use upset plot
  # use_venn <- ifelse(length(venn_list) <= 4, TRUE, FALSE)

  if (!requireNamespace("futile.logger", quietly = TRUE)) {
    stop("Package futile.logger needed for this function to work. Please install it.",
         call. = FALSE)
  }

  #--- codes ---#
  if (use_venn) {
    # suppress venn.diagram log
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

    # choose color
    if (is.null(color) | length(color) != length(venn_list)) {
      message("Color length should be same with venn_list, auto assign colors...")
      if(length(venn_list) !=2){
        color <- RColorBrewer::brewer.pal(length(venn_list), "Set1")
      }else{
        color <- ggsci::pal_lancet()(2)
      }

    }

    # hide background grid line
    if (remove_grid) border_thick <- 0

    p <- VennDiagram::venn.diagram(
      x = venn_list,
      filename = NULL,
      category.names = names(venn_list),
      output = TRUE,
      cat.cex = text_size,
      main.cex = text_size,
      cex = text_size,
      lwd = border_thick,
      col = color,
      cat.col = color,
      fill = sapply(color, function(x) scales::alpha(x, alpha_degree))
    ) %>%
      cowplot::as_grob() %>%
      ggplotify::as.ggplot()
  } else {
    # use ggupset
    # https://github.com/const-ae/ggupset

    if (text_size < 10) {
      message("Text size is too low, auto set: text_size=10 ...")
      text_size <- 10
    }
    p <- sapply(venn_list, function(x) unique(unlist(venn_list)) %in% x) %>%
      t() %>%
      as.data.frame() %>%
      stats::setNames(., unique(unlist(venn_list))) %>%
      as.matrix() %>%
      dplyr::as_tibble(rownames = "sets") %>%
      tidyr::gather(item, type, -sets) %>%
      dplyr::filter(type) %>%
      dplyr::select(-type) %>%
      dplyr::group_by(item) %>%
      dplyr::summarize(sets = list(sets)) %>%
      ggplot(aes(x = sets)) +
      geom_bar() +
      geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1, size = 3) +
      ggupset::scale_x_upset(name = "") +
      ggplot2::scale_y_continuous(name = "") +
      plot_theme()

    # hide background grid line
    if (remove_grid) {
      p <- p + theme(
        # panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    }
  }

  return(p)
}
