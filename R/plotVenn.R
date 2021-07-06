##' Venn plot for list of genes
##' If venn list length is over 4, use UpSet plot
##'
##' @param venn_list a group list of genes or others to plot venn
##' @param color colors for venn_list, default is NULL.
##' @param alpha_degree alpha transparency of each circle's area, default is 0.3.
##' @param border_thick thickness of each circle, default is 1.
##' @param text_size text size, default is 1.
##' @param remove_grid logical, remove circle or grid lines, default is FALSE.
##' @return ggplot object
##' @importFrom futile.logger flog.threshold ERROR
##' @importFrom VennDiagram venn.diagram
##' @importFrom RColorBrewerbrewer.pal
##' @importFrom cowplot as_grob
##' @importFrom ggplotify as.ggplot
##' @importFrom scales alpha
##' @importFrom stats setNames
##' @importFrom dplyr as_tibble filter select group_by summarize
##' @importFrom tidyr gather
##' @importFrom ggplot2 ggplot geom_bar aes geom_text after_stat theme element_blank
##' @importFrom ggupset scale_x_upset  scale_y_continuous
##' @export
##' @examples
##' \dontrun{
##' set1 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
##' set2 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
##' set3 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
##' set4 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
##' set5 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
##' sm_gene_list = list(gset1 = set1, gset2 = set2, gset3 = set3)
##' la_gene_list = list(gset1 = set1, gset2 = set2, gset3 = set3, gset4 = set4, gset5 = set5 )
##' plotVenn(sm_gene_list,text_size = 1.5,alpha_degree = 1,remove_grid = T,color = ggsci::pal_lancet()(3))
##' plotVenn(la_gene_list,text_size = 15,alpha_degree = 0.2,border_thick = 2,remove_grid = T)
##'
##' }

plotVenn <- function(venn_list,
                     color = NULL,
                     alpha_degree = 0.3,
                     border_thick = 1,
                     text_size = 1,
                     remove_grid = FALSE,
                     ...) {

  #--- args ---#
  stopifnot(is.list(venn_list))
  # if gene list too long, use upset plot
  use_venn <- ifelse(length(venn_list) <= 4, TRUE, FALSE)

  #--- codes ---#
  if (use_venn) {
    # suppress venn.diagram log
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

    # choose color
    if (is.null(color) | length(color) != length(venn_list)) {
      message("Color length should be same with venn_list, auto assign colors...")
      color <- RColorBrewer::brewer.pal(length(venn_list), "Set1")
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
      ggupset::scale_y_continuous(name = "") +
      plot_theme(border_thick = border_thick, main_text_size = text_size)

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
