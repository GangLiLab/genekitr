#' Dotplot for GO and KEGG enrichment analysis
#'
#' @param enrich_df `data.frame` of enrichment analysis result .
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
#'   guide_colorbar scale_y_discrete element_blank xlab labs xlim
#' @importFrom stringr str_to_title
#'
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' data(geneList, package = "genekitr")
#' id <- names(geneList)[1:100]
#' ego <- genGO(id,
#'   org = "human", ont = "mf", pvalueCutoff = 0.05,
#'   qvalueCutoff = 0.1, use_symbol = FALSE
#' )
#' plotEnrichDot(ego,remove_grid = T, main_text_size = 8,
#'   legend_text_size = 6,border_thick = 1.5)
#' }
#'
plotEnrichDot <- function(enrich_df,
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
  xlab_type <- match.arg(xlab_type)
  legend_type <- match.arg(legend_type)

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

  #--- plot ---#
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


