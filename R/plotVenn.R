#' Venn plot for groups of genes
#'
#' If gene group over 4, plot will be visulized using UpSet plot.
#'
#' @param venn_list A list of gene id.
#' @param use_venn Logical, use venn to plot, default is `TRUE`, the other
#'   option is upsetplot for large list.
#' @param color Colors for gene lists, default is NULL.
#' @param alpha_degree Alpha transparency of each circle's area, default is 0.3.
#' @param venn_percent Logical to show both number and percentage in venn plot.
#' @param ... other arguments transfer to `plot_theme` function
#' @return  A ggplot object
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 ggplot geom_bar aes geom_text after_stat theme
#'   element_blank scale_y_continuous geom_segment geom_point element_text
#' @importFrom rlang .data
#' @export
#' @examples
#' k1 = requireNamespace("ComplexUpset",quietly = TRUE)
#' k2 = requireNamespace("futile.logger",quietly = TRUE)
#' k3 = requireNamespace("ggsci",quietly = TRUE)
#' k4 = requireNamespace("RColorBrewer",quietly = TRUE)
#' if(k1&k2&k3&k4){
#' library(ggplot2)
#' set1 <- paste0(rep("gene", 30), sample(1:1000, 30))
#' set2 <- paste0(rep("gene", 40), sample(1:1000, 40))
#' set3 <- paste0(rep("gene", 50), sample(1:1000, 50))
#' set4 <- paste0(rep("gene", 60), sample(1:1000, 60))
#' set5 <- paste0(rep("gene", 70), sample(1:1000, 70))
#' sm_gene_list <- list(gset1 = set1, gset2 = set2, gset3 = set3)
#' la_gene_list <- list(
#'   gset1 = set1, gset2 = set2, gset3 = set3,
#'   gset4 = set4, gset5 = set5
#' )
#' plotVenn(sm_gene_list,
#'   use_venn = TRUE,
#'   alpha_degree = 0.5,
#'   main_text_size = 3,
#'   border_thick = 0,
#'   venn_percent = TRUE
#' )
#' plotVenn(la_gene_list,
#'   use_venn = FALSE,
#'   main_text_size = 15,
#'   legend_text_size = 8,
#'   legend_position = 'left'
#' )
#' }
#'
plotVenn <- function(venn_list,
                     use_venn = TRUE,
                     color = NULL,
                     alpha_degree = 0.3,
                     venn_percent = FALSE,
                     ...) {

  #--- args ---#
  lst <- list(...) # store outside arguments in list

  if (!requireNamespace("futile.logger", quietly = TRUE)) {
    warning("Package futile.logger needed for this function to work. Install first...",
      call. = FALSE
    )
    # utils::install.packages("futile.logger")
  }

  #--- codes ---#
  ## Venn Diagram
  if (use_venn) {
    # if (!"label_text_size" %in% names(lst)) lst$label_text_size <- 6
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 3
    if (!"border_thick" %in% names(lst)) lst$border_thick <- 1
    # if (!"digits" %in% names(lst)) lst$digits <- 2
    # suppress venn.diagram log
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

    # choose color
    if (is.null(color) | length(color) != length(venn_list)) {
      message("Color length should be same with venn_list, auto assign colors...")
      if (length(venn_list) != 2) {
        color <- RColorBrewer::brewer.pal(length(venn_list), "Set1")
      } else {
        color <- ggsci::pal_lancet()(2)
      }
    }

    if(venn_percent){
      print_mode <- c('raw','percent')
    }else{
      print_mode <- 'raw'
    }
    if(is.null(names(venn_list))) names(venn_list) = paste0('group',1:length(venn_list))

    p <- ggvenn(venn_list,
           show_elements = F,
           show_percentage = venn_percent,
           digits = 2,
           label_sep = "\n",
           fill_color = color,
           # stroke_color = color,
           # set_name_color = color,
           fill_alpha = alpha_degree,
           text_size = lst$main_text_size ,
           set_name_size = lst$main_text_size,
           stroke_size = lst$border_thick
    )

    # p <- VennDiagram::venn.diagram(
    #   x = venn_list,
    #   filename = NULL,
    #   category.names = names(venn_list),
    #   output = TRUE,
    #   cat.cex = lst$main_text_size,
    #   main.cex = lst$main_text_size,
    #   cex = lst$main_text_size,
    #   lwd =  lst$border_thick,
    #   fontfamily = 'Arial',
    #   cat.fontfamily = 'Arial',
    #   col = color,
    #   cat.col = color,
    #   fill = sapply(color, function(x) scales::alpha(x, alpha_degree)),
    #   ext.text = F,
    #   cat.pos = 0,
    #   print.mode= print_mode,
    #   sigdigs = 2,
    #   disable.logging = TRUE,
    #   scaled = FALSE
    # ) %>%
    #   cowplot::as_grob() %>%
    #   ggplotify::as.ggplot()

    } else {
     ## ComplexUpset Diagram
    # if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    #   utils::install.packages("ComplexUpset")
    # }

    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 10
    if (!"legend_text_size" %in% names(lst)) lst$legend_text_size <- 8
    if (!"border_thick" %in% names(lst)) lst$border_thick <- 2
    if (!"legend_position" %in% names(lst)) lst$legend_position <- 'left'

    if(is.null(color))
      color <- c(
        "#B2DF8A", "#FB9A99", "#E31A1C", "#B15928", "#6A3D9A", "#CAB2D6",
        "#A6CEE3", "#1F78B4", "#FDBF6F", "#999999", "#FF7F00",
        "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D", "#58CDD9",
        "#7A142C", "#5D90BA", "#029149", "#431A3D", "#91612D", "#6E568C",
        "#E0367A", "#D8D155", "#64495D", "#7CC767"
      )

    if(is.list(venn_list)){
      dat <- venn_list %>% as.upset() %>%
        do.call(cbind,.) %>% as.data.frame()
    }else{
      dat <- venn_list
    }


    p <- ComplexUpset::upset(dat, colnames(dat), name = '', width_ratio=0.1,
          queries=lapply(colnames(dat), function(x){
            ComplexUpset::upset_query(set = x,
                        fill = color[which(colnames(dat)%in%x)])
          }),
          base_annotations=list(
            'Intersection size'=(
              ComplexUpset::intersection_size(
                # label height great than this will be placed inside bar
                bar_number_threshold=0.95,
                width=0.5,
                text = list(size = lst$legend_text_size/2 )
              )+ scale_y_continuous(expand=expansion(mult=c(0, 0.05)))+
                plot_theme()+
                theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank(),
                      axis.text.y = element_text(size = (lst$legend_text_size + 3)),
                      panel.border =element_rect(colour = "black", size = lst$border_thick))
            )
          ),
          matrix=ComplexUpset::intersection_matrix(
            geom=geom_point(
              shape='circle filled',
              size=3,
              stroke=0.5
            ),
            segment = geom_segment(size = 0.7,color = 'grey46')
          ),
          set_sizes=(
            ComplexUpset::upset_set_size(geom=geom_bar(width=0.5),
                           position = lst$legend_position)+
            theme(
              axis.line.x=element_line(colour='black'),
              axis.ticks.x=element_line(),
              axis.title.x=element_blank(),
              axis.text.x=element_text(size = lst$legend_text_size)
            )
          ),
          stripes=ComplexUpset::upset_stripes(
            geom=ggplot2::geom_segment(size=5),
            colors=c('grey95', 'white')
          ),
          sort_sets='ascending',
          sort_intersections='ascending',
          themes = ComplexUpset::upset_modify_themes(
            list(
              'intersections_matrix'=theme(text=element_text(size=lst$main_text_size))
            )
          )

    )

    # p <- sapply(venn_list, function(x) unique(unlist(venn_list)) %in% x) %>%
    #   t() %>%
    #   as.data.frame() %>%
    #   stats::setNames(., unique(unlist(venn_list))) %>%
    #   as.matrix() %>%
    #   dplyr::as_tibble(rownames = "sets") %>%
    #   tidyr::gather(item, type, -sets) %>%
    #   dplyr::filter(type) %>%
    #   dplyr::select(-type) %>%
    #   dplyr::group_by(item) %>%
    #   dplyr::summarize(sets = list(sets)) %>%
    #   ggplot(aes(x = sets)) +
    #   geom_bar() +
    #   geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1, size = 3) +
    #   ggupset::scale_x_upset(name = "") +
    #   ggplot2::scale_y_continuous(name = "") +
    #   plot_theme(
    #     main_text_size = text_size,
    #     remove_grid = remove_grid,
    #     border_thick = border_thick,
    #     ...
    #   )
  }

  return(p)
}


# obj: dataframe or list
# inter_by: select column to intersect
# name_by: select column as shown names in plot
as.upset <- function(obj,
                     inter_by = "geneID",
                     name_by = "Description"){
  if(is.data.frame(obj)){
    allg_lst <- obj[,inter_by] %>%
      stringr::str_split(.,'\\/')
    names(allg_lst) = obj[,name_by]

    allg <- allg_lst %>% unlist() %>% unique()
    res <- lapply(allg_lst, function(x) allg%in%x)
    return(res)
  }else if(is.list(obj)){
    allg <- obj %>% unlist() %>% unique()
    res <- lapply(obj, function(x) allg%in%x)
    return(res)
  }

}



