#' Plot for GO and KEGG enrichment analysis
#'
#' @param enrich_df `data.frame` of enrichment analysis result.
#' @param plot_type Choose from "bar", "dot", "heat", "chord"
#' @param xlab_type X-axis label type, one of 'GeneRatio','Count','FoldEnrich'.
#' @param legend_type Stats legend type, one of "pvalue", "p.adjust", "qvalue".
#' @param low_color Legend color for low pvalue or qvalue, default is "red".
#' @param high_color Legend color for high pvalue or qvalue, default is "blue".
#' @param show_item Numeric, select top N rows to show, default is 10.
#' @param xlim_left X-axis left limit, default is 0.
#' @param xlim_right X-axis right limit, default is NA.
#' @param wrap_length Numeric, wrap text if longer than this length, default is NULL.
#' @param show_genes Select genes in heat plot. Default is "all".
#' @param ... other arguments transfer to `plot_theme` function
#'
#' @importFrom dplyr pull %>% arrange mutate slice_head
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_continuous theme
#'   guide_colorbar scale_y_discrete element_blank xlab labs xlim scale_x_discrete
#'   facet_grid
#' @importFrom stringr str_to_title
#' @importFrom rlang .data
#' @importFrom stats setNames
#'
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
#'
#' plotEnrich(ego,plot_type = "heat",
#'   show_genes = c('LRRTM2','LRRN1','SLC14A1','AQP3','GRIN2A'))
#'
#' plotEnrich(ego,plot_type = "chord", show_genes = "all")
#' }
#'
plotEnrich <- function(enrich_df,
                       plot_type = c('bar','dot','heat','chord'),
                       xlab_type = c("FoldEnrich", "GeneRatio", "Count"),
                       legend_type = c("p.adjust", "pvalue", "qvalue"),
                       low_color = "red",
                       high_color = "blue",
                       show_item = 10,
                       xlim_left = 0,
                       xlim_right = NA,
                       wrap_length = NULL,
                       show_genes = 'all',
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
  ## dot plot
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

  ## bar plot
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

  ## heat plot
  if(plot_type == 'heat'){
    if(all(show_genes == 'all')){
      plot_df = enrich_df %>% dplyr::select(Description,geneID) %>%
        tidyr::separate_rows(geneID,sep = '\\/')
    }else{
      plot_df = enrich_df %>% dplyr::select(Description,geneID) %>%
        tidyr::separate_rows(geneID,sep = '\\/') %>%
        dplyr::filter(geneID %in% show_genes)
    }

    p <- ggplot(plot_df, aes_(~geneID, ~Description)) +
      geom_tile(color = 'white')+
      xlab(NULL) + ylab(NULL) +
      plot_theme(theme_type = 'bw',
                 border_thick = 0,...)+
      theme(panel.grid.major = element_blank(),
            axis.text.x=element_text(angle = 50, hjust = 1))
  }

  ## chord plot
  if(plot_type == 'chord'){
    if(all(show_genes == 'all')){
      plot_df = enrich_df %>% dplyr::select(Description,geneID) %>%
        tidyr::separate_rows(geneID,sep = '\\/')
    }else{
      plot_df = enrich_df %>% dplyr::select(Description,geneID) %>%
        tidyr::separate_rows(geneID,sep = '\\/') %>%
        dplyr::filter(geneID %in% show_genes)
    }

    term = unique(plot_df$Description)
    id = unique(plot_df$geneID)
    dat = sapply(id, function(x){
      check = plot_df %>% dplyr::filter(geneID %in% x) %>%
        dplyr::pull(Description) %>%
        unique()
      ifelse(term%in%check,1,0)
    }) %>% t() %>%
      as.data.frame() %>%
      stats::setNames(term)

    my_cols= c("#B2DF8A","#FB9A99","#E31A1C","#B15928","#6A3D9A","#CAB2D6",
               "#A6CEE3","#1F78B4","#FDBF6F","#999999","#FF7F00")
    backup_cols = c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9",
                    "#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C",
                    "#E0367A","#D8D155","#64495D","#7CC767")
    if(length(my_cols) < length(term)){
      my_cols = c(my_cols,backup_cols)
    }
    cols = sample(my_cols,length(term),replace = F)

    p = suppressWarnings(
      GOplot::GOChord(dat,
              space = 0.02,
              gene.order = 'none',
              gene.size = 3,
              border.size = 0.1,
              process.label = 8,
              ribbon.col = cols)
    )
  }


  # wrap long text
  if (!is.null(wrap_length) & is.numeric(wrap_length)) {
    p <- p + scale_y_discrete(labels = text_wraper(wrap_length))
  }

  suppressWarnings(print(p))

}


#--- sub-function: wrap text if too long ---#
text_wraper <- function(width) {
  function(x) {
    lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n")
  }
}


