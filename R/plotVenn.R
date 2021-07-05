##' Venn plot for gene id list
##'
##' @param gene_list a group list of genes to compare.
##' @param color colors for genelist, default is NULL.
##' @param alpha_degree alpha transparency of each circle's area, default is 0.3.
##' @param circle_thick thickness of each circle, default is 1.
##' @param text_size text size, default is 1.
##' @return ggplot object
##' @importFrom futile.logger  flog.threshold ERROR
##' @importFrom VennDiagram  venn.diagram
##' @importFrom RColorBrewer  brewer.pal
##' @importFrom cowplot  as_grob
##' @importFrom ggplotify  as.ggplot
##' @importFrom ggupset  as.ggplot
##' @export
##' @examples
##' \dontrun{
##'set1 <- paste(rep("gene" , 50) , sample(c(1:500) , 50 , replace=F) , sep="")
##'set2 <- paste(rep("gene" , 50) , sample(c(1:500) , 50 , replace=F) , sep="")
##'set3 <- paste(rep("gene" , 50) , sample(c(1:500) , 50 , replace=F) , sep="")
##'gene_list = list(gset1 = set1, gset2 = set2, gset3 = set3)
##'plotVenn(gene_list,text_size = 1,alpha_degree = 0,circle_thick = 2)
##' }

plotVenn <- function(gene_list,
                     color = NULL,
                     alpha_degree = 0.3,
                     circle_thick = 1,
                     text_size =1,
                     ...){

  #--- args ---#
  stopifnot(is.list(gene_list))
  # if gene list too long, use upset plot
  use_venn <- ifelse(length(gene_list) <= 4, TRUE,FALSE)

  #--- codes ---#
  if(use_venn){
    # suppress venn.diagram log
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    # choose color
    if( is.null(color) | length(color) != length(gene_list) ){
      color <- RColorBrewer::brewer.pal(length(gene_list), "Set1")
    }

    p = VennDiagram::venn.diagram(
      x = gene_list,
      filename=NULL,
      category.names = names(gene_list),
      output=TRUE,
      cat.cex = text_size,
      main.cex = text_size,
      cex = text_size,
      lwd = circle_thick,
      col=color,
      cat.col=color,
      fill = sapply(color, function(x) alpha(x,alpha_degree))
    ) %>%
      cowplot::as_grob() %>%
      ggplotify::as.ggplot()

  }else{
  # use ggupset
  # https://github.com/const-ae/ggupset






  }



  return(p)
}

