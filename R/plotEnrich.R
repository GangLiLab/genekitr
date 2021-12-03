#' Dotplot for GO and KEGG enrichment analysis
#'
#' @param enrich_df `data.frame` of enrichment analysis result.
#' @param plot_type One of "bar" and "dot"
#' @param xlab_type X-axis label type, one of 'GeneRatio','Count','FoldEnrich'.
#' @param legend_type Stats legend type, one of "pvalue", "p.adjust", "qvalue".
#' @param low_color Legend color for low pvalue or qvalue, default is "red".
#' @param high_color Legend color for high pvalue or qvalue, default is "blue".
#' @param show_item Numeric, select top N rows to show, default is 10.
#' @param xlim_left X-axis left limit, default is 0.
#' @param xlim_right X-axis right limit, default is NA.
#' @param wrap_length Numeric, wrap text if longer than this length, default is NULL.
#' @param ... other arguments transfer to `plot_theme` function
#'
#' @importFrom dplyr pull %>% arrange mutate slice_head
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_continuous theme
#'   guide_colorbar scale_y_discrete element_blank xlab labs xlim scale_x_discrete
#'   facet_grid
#' @importFrom stringr str_to_title
#' @importFrom rlang .data
#'
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' data(geneList, package = "genekitr")
#' id <- names(geneList)[1:100]
#' ego <- genGO(id,
#'   org = "human", ont = "bp", pvalueCutoff = 0.05,
#'   qvalueCutoff = 0.05, use_symbol = TRUE
#' )
#' plotEnrich(ego,plot_type = "dot")
#'
#' plotEnrich(ego,plot_type = "bar")
#' }
#'
plotEnrich <- function(enrich_df,
                       plot_type = c('bar','dot'),
                       xlab_type = c("FoldEnrich", "GeneRatio", "Count"),
                       legend_type = c("p.adjust", "pvalue", "qvalue"),
                       low_color = "red",
                       high_color = "blue",
                       show_item = 10,
                       xlim_left = 0,
                       xlim_right = NA,
                       wrap_length = NULL,
                       ...) {
  #--- args ---#
  stopifnot(is.numeric(show_item))
  plot_type <- match.arg(plot_type)
  xlab_type <- match.arg(xlab_type)
  legend_type <- match.arg(legend_type)
  compare_group <- any(grepl('cluster',colnames(enrich_df),ignore.case = T))
  all_go <- any(grepl('ONTOLOGY',colnames(enrich_df),ignore.case = T))

  if(compare_group) plot_type = 'dot'
  if(all_go) plot_type = 'bar'

  types <- c("GeneRatio", "Count", "FoldEnrich")
  legends <- c("p.adjust", "pvalue", "qvalue")

  if (!xlab_type %in% colnames(enrich_df)) {
    stop(
      xlab_type, " not included in this dataframe, try: ",
      paste(intersect(colnames(enrich_df), types), collapse = " | ")
    )
  }
  if (!legend_type %in% colnames(enrich_df)) {
    stop(
      legend_type, " not included in this dataframe, try: ",
      paste(intersect(colnames(enrich_df), legends), collapse = " | ")
    )
  }

  #--- codes ---#
  xlab_title <- ifelse(xlab_type == "FoldEnrich", "Fold Enrichment",
    ifelse(xlab_type == "GeneRatio", "Gene Ratio", "Count")
  )
  legend_title <- ifelse(legend_type == "pvalue", "Pvalue",
    ifelse(legend_type == "p.adjust", "P.adjust", "FDR")
  )

  if(!compare_group & !all_go){
    if (show_item <= nrow(enrich_df)) {
      enrich_df <- enrich_df %>%
        dplyr::arrange(eval(parse(text = xlab_type))) %>%
        dplyr::mutate(Description = factor(.$Description, levels = .$Description, ordered = T)) %>%
        dplyr::slice_head(., n = show_item)
    } else {
      enrich_df <- enrich_df %>%
        dplyr::arrange(eval(parse(text = xlab_type))) %>%
        dplyr::mutate(Description = factor(.$Description, levels = .$Description, ordered = T))
    }
  }else if (!compare_group & all_go){
    if (show_item <= nrow(enrich_df)) {
      enrich_df <- enrich_df %>%
        dplyr::mutate(Description = factor(.$Description, levels = .$Description, ordered = T)) %>%
        dplyr::group_by(ONTOLOGY) %>%
        dplyr::arrange(eval(parse(text = xlab_type)),.by_group = T) %>%
        dplyr::slice_head(.,n=show_item)
    } else {
      enrich_df <- enrich_df %>%
        dplyr::arrange(eval(parse(text = xlab_type))) %>%
        dplyr::mutate(Description = factor(.$Description, levels = .$Description, ordered = T))
    }
  }

  #--- plot ---#
  if(plot_type == 'dot'){
    if(!compare_group){
      p <- ggplot(enrich_df, aes_string(x = xlab_type, y = "Description")) +
        geom_point(aes_string(
          color = legend_type,
          size = "Count"
        )) +
        scale_color_continuous(
          low = low_color, high = high_color, name = legend_title,
          guide = guide_colorbar(reverse = TRUE),
          labels = function(x) format(x, scientific = T)
        ) +
        xlab(xlab_title) +
        labs(color = legend_type)+
        xlim(xlim_left,xlim_right)+
        plot_theme(...)
    }else{
      xtick_lab <- paste0(enrich_df$Cluster,'\n(',enrich_df$Count,')')
      p <- ggplot(enrich_df, aes_string(x = 'Cluster', y = "Description")) +
        geom_point(aes_string(
          color = legend_type,
          size = xlab_type
        )) +
        scale_color_continuous(
          low = low_color, high = high_color, name = legend_title,
          guide = guide_colorbar(reverse = TRUE),
          labels = function(x) format(x, scientific = T)
        ) +
        scale_x_discrete(labels= xtick_lab)+
        xlab('Group') +
        labs(color = legend_type)+
        plot_theme(...)+
        theme(axis.text.x = element_text(
          angle = 45,
          vjust = 0.5, hjust = 0.5
        ))
    }

  }

  if(plot_type == 'bar'){
    p <- ggplot(data=enrich_df, aes_string(x = xlab_type, y = 'Description', fill = legend_type)) +
      geom_bar(stat="identity")+
      ggplot2::scale_fill_continuous(
        low = low_color, high = high_color, name = legend_title,
        guide = guide_colorbar(reverse = TRUE),
        labels = function(x) format(x, scientific = T))+
      xlab(xlab_title)+
      labs(color = legend_type)+
      xlim(xlim_left,xlim_right)+
      plot_theme(...)

    if(all_go) p <- p + facet_grid(ONTOLOGY~., scales = "free")+
        plot_theme(...)

  }


  # wrap long text
  if (!is.null(wrap_length) & is.numeric(wrap_length)) {
    p <- p + scale_y_discrete(labels = text_wraper(wrap_length))
  }

  return(p)
}


#--- wrap text if too long ---#
text_wraper <- function(width) {
  function(x) {
    lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n")
  }
}


