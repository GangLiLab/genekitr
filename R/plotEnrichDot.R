##' Dotplot for enrichment analysis
##'
##' @param enrich_df dataframe of enrichment analysis result .
##' @param xlab_type x-axis label type, one of 'GeneRatio','Count','FoldEnrich'.
##' @param legend_by stats legend type, one of "pvalue", "p.adjust", "qvalue".
##' @param show_item numeric, select top N rows to show.
##' @param xleft numeric, specify the x-axis left limit, default is 0.
##' @param xright numeric, specify the x-axis right limit, default is NA.
##' @param main_text_size numeric, specify the plot text size.
##' @param font_type character, specify the plot text font family, example "Times New Roman", "Arial".
##' @param remove_grid logical, remove background grid lines, default is FALSE.
##' @param wrap_width numeric, wrap text longer than this number.
##' @param border_thick numeric, border thickness in mm.
##' @return ggplot object
##' @importFrom dplyr pull
##' @importFrom ggplot2 ggplot aes geom_point scale_color_continuous theme guide_colorbar scale_y_discrete element_blank
##' @importFrom stringr str_to_title
##' @importFrom clusterProfiler enrichGO
##' @importFrom DOSE setReadable
##' @export
##' @author Yunze Liu
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' id = names(geneList)[1:100]
##' ego <- genGO(id, org = 'human',ont = 'mf',pvalueCutoff = 0.05,qvalueCutoff = 0.1 ,use_symbol = T)
##' plotEnrichDot(ego,xlab_type =  'FoldEnrich', legend_by = 'qvalue',show_item = 10,remove_grid = T)
##' }


plotEnrichDot <- function(enrich_df,
                          xlab_type = c("GeneRatio", "Count", "FoldEnrich"),
                          legend_by = c("p.adjust", "pvalue", "qvalue"),
                          low_color = "red",
                          high_color = "blue",
                          font_type = "Arial",
                          show_item = 10,
                          xleft = 0, xright = NA,
                          main_text_size = 10,
                          legend_text_size = 8,
                          border_thick = 1,
                          remove_grid = FALSE,
                          wrap_width = NULL,
                          ...) {
  #--- args ---#
  stopifnot(is.numeric(show_item))
  xlab_type <- match.arg(xlab_type)
  legend_by <- match.arg(legend_by)

  enrich_df_bk1 <- .check_colname(enrich_df, xlab_type)
  enrich_df <- .check_colname(enrich_df_bk1, legend_by)


  #--- codes ---#
  xlab_title <- ifelse(xlab_type == "FoldEnrich", "Fold Enrichment",
    ifelse(xlab_type == "GeneRatio", "Gene Ratio", "Count")
  )
  legend_title <- ifelse(legend_by == "pvalue", "Pvalue",
    ifelse(legend_by == "p.adjust", "P.adjust", "FDR")
  )

  # Panther GO result
  check_panther <- enrich_df %>%
    dplyr::pull(1) %>%
    stringr::str_detect(".*\\(GO")
  if (any(check_panther)) {
    enrich_df <- enrich_df %>%
      dplyr::rename("Description" = names(.)[1]) %>%
      dplyr::rename("Count" = names(.)[3]) %>%
      dplyr::mutate(Description = stringr::str_to_sentence(.$Description) %>%
        stringr::str_replace_all(., "go:", "GO:"))
  } else {
    # clusterP  GO result
    enrich_df <- enrich_df %>%
      dplyr::mutate(Description = stringr::str_to_sentence(.$Description)) %>%
      dplyr::mutate(GeneRatio = sapply(.$GeneRatio, function(x) eval(parse(text = x))))
  }

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
  p <- ggplot(enrich_df, aes(x = eval(parse(text = xlab_type)), y = Description)) +
    geom_point(aes(
      color = eval(parse(text = legend_by)),
      size = Count
    )) +
    scale_color_continuous(
      low = low_color, high = high_color, name = legend_title,
      guide = guide_colorbar(reverse = TRUE),
      labels = function(x) format(x, scientific = T)
    ) +
    xlab(xlab_title) +
    plot_theme(main_text_size, legend_text_size, font_type, border_thick) +
    xlim(xleft, xright) +
    labs(color = legend_by)

  # hide background grid line
  if (remove_grid) {
    p <- p + theme(
      # panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }
  # wrap long text
  if (!is.null(wrap_width) & is.numeric(wrap_width)) {
    p <- p + scale_y_discrete(labels = text_wraper(wrap_width))
  }

  return(p)
}
