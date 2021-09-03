#' Dotplot for GO and KEGG enrichment analysis
#'
#' @param enrich_df `data.frame` of enrichment analysis result .
#' @param xlab_type X-axis label type, one of 'GeneRatio','Count','FoldEnrich'.
#' @param legend_by Stats legend type, one of "pvalue", "p.adjust", "qvalue".
#' @param remove_grid Logical, remove background grid lines, default is FALSE.
#' @param remove_text Logical, remove all text, default is FALSE.
#' @param remove_legend Logical, remove legend, default is FALSE.
#' @param low_color Legend color for low pvalue or qvalue, default is "red".
#' @param high_color Legend color for high pvalue or qvalue, default is "blue".
#' @param font_type Character, specify the plot text font family, example "Times
#'   New Roman", "Arial".
#' @param show_item Numeric, select top N rows to show, default is 10.
#' @inheritParams plot_theme
#' @param wrap_width Numeric, wrap text if longer than this number, default is NULL.
#'
#' @importFrom dplyr pull %>% arrange mutate slice_head
#' @importFrom ggplot2 ggplot aes geom_point scale_color_continuous theme
#'   guide_colorbar scale_y_discrete element_blank xlab labs
#' @importFrom stringr str_to_title
#'
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' data(geneList, package="DOSE")
#' id = names(geneList)[1:100]
#' ego = genGO(id, org = 'human',ont = 'mf',pvalueCutoff = 0.05,
#'   qvalueCutoff = 0.1 ,use_symbol = FALSE)
#' ego = as.enrichdat(ego)
#' plotEnrichDot(ego)
#' }


plotEnrichDot <- function(enrich_df,
                          xlab_type = c("FoldEnrich", "GeneRatio", "Count"),
                          legend_by = c("p.adjust", "pvalue", "qvalue"),
                          remove_grid = FALSE,
                          remove_text = FALSE,
                          remove_legend = FALSE,
                          low_color = "red",
                          high_color = "blue",
                          show_item = 10,
                          # xleft = 0, xright = NA,
                          main_text_size = 10,
                          legend_text_size = 8,
                          font_type = 'Arial',
                          border_thick = 1,
                          wrap_width = NULL) {
  #--- args ---#
  stopifnot(is.numeric(show_item))
  xlab_type <- match.arg(xlab_type)
  legend_by <- match.arg(legend_by)

  types = c("GeneRatio", "Count", "FoldEnrich")
  legends = c("p.adjust", "pvalue", "qvalue")

  if(! xlab_type %in% colnames(enrich_df)){
    stop(xlab_type,' not included in this dataframe, try: ',
         paste(intersect(colnames(enrich_df), types),collapse = ' | '))
  }
  if(! legend_by %in% colnames(enrich_df)){
    stop(legend_by,' not included in this dataframe, try: ',
         paste(intersect(colnames(enrich_df), legends),collapse = ' | '))
  }

  #--- codes ---#
  xlab_title <- ifelse(xlab_type == "FoldEnrich", "Fold Enrichment",
    ifelse(xlab_type == "GeneRatio", "Gene Ratio", "Count")
  )
  legend_title <- ifelse(legend_by == "pvalue", "Pvalue",
    ifelse(legend_by == "p.adjust", "P.adjust", "FDR")
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
    # xlim(xleft, xright) +
    labs(color = legend_by)

  # hide background grid line
  if (remove_grid) {
    p <- p + theme(
      # panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }

  # hide axis text
  if (remove_text) {
    p <- p + theme(
      # panel.border = element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
    )
  }

  # hide legend
  if (remove_legend) {
    p <- p + theme(
      legend.position = "none"
    )
  }

  # wrap long text
  if (!is.null(wrap_width) & is.numeric(wrap_width)) {
    p <- p + scale_y_discrete(labels = text_wraper(wrap_width))
  }

  return(p)
}




##' Adjust dataframe for enrichment plot
##'
##' make sure colname contains Description, Count, FoldEnrich/GeneRatio, pvalue/qvalue/p.adjust
##'
##' @param enrich_df dataframe of enrichment analysis result .
##'
##' @importFrom stringr str_remove_all str_split str_remove
##' @importFrom dplyr mutate pull
##' @return A `data.frame`.
##' @export

as.enrichdat <- function(enrich_df){

  ## get lower case colnames
  remove <- c("\\(", "\\)", " ",'-','_')
  to_check = stringr::str_remove_all(tolower(colnames(enrich_df)), paste(remove, collapse = "|"))

  ## find description col
  if(!any(grepl('description', to_check))){
    check_description = apply(enrich_df, 2, function(x) all(grepl('[A-Za-z]{3,}',x)) )
    check2 = which(check_description)
    if(any(check2 <  (ncol(enrich_df)/2))){
      colnames(enrich_df)[check_description] [check2 < 5] = 'Description'
    }else{
      stop('Not found description column!')
    }
  }else{
    colnames(enrich_df)[grepl('description', to_check)] = 'Description'
  }

  ## find count col
  # if finds genes col, calc gene num as count; else find another col as count
  if(!any(grepl('count', to_check))){
    check_gene = apply(enrich_df, 2, function(x) all(grepl('[A-Za-z]{3,}|\\/|,',x) & !grepl('tags|list',x)) )
    check2 = which(check_gene)
    if(any(check2 > (ncol(enrich_df)/2))){
      colnames(enrich_df)[check_gene] [check2 > (ncol(enrich_df)/2)] = 'geneID'
      gen_num  = stringr::str_split(enrich_df$geneID,',|\\/') %>% lapply(.,length) %>% unlist()
      enrich_df = enrich_df %>% dplyr::mutate(Count = gen_num)
    }else if(any(grepl('\\([1-9]{,4}\\)', colnames(enrich_df)))){
      gen_num = enrich_df[grepl('\\([1-9]{,4}\\)', colnames(enrich_df))] %>% dplyr::pull(1) %>% as.numeric()
      enrich_df = enrich_df %>% dplyr::mutate(Count = gen_num)
    } else{
      stop('Please rename the gene count column as "Count"!')
      # head(enrich_df[1:2,])
      # message("Cannot auto-select count column...","\n","Please specify the column number which includes gene count...")
      # answer <- scan(what = "character", n =1,quiet =T)
      # message('Choose the No. ',answer,' column as gene count...')
      # enrich_df = enrich_df %>% dplyr::rename(Count = eval(parse(text = answer)))
    }
  }else{
    colnames(enrich_df)[grepl('count', to_check)] = 'Count'
  }

  ## find FoldEnrich col
  # GSEA result has no FoldEnrich, need to exclude
  if(! (any(grepl('\\benrichmentscore\\b', to_check)) & any(grepl('\\bleadingedge\\b', to_check)))){
    if(any(grepl('foldenrich|enrichment', to_check))){
      colnames(enrich_df)[grepl('foldenrich|enrichment', to_check)] = 'FoldEnrich'
    }else{
      stop('Please rename the fold enrichment column as "FoldEnrich"!')
      # head(enrich_df[1:2,])
      # message("Cannot auto-select fold enrichment column...","\n",
      #         "Please specify the column number which includes fold enrichment...")
      # answer <- scan(what = "character", n =1,quiet =T)
      # message('Choose the No. ',answer,' column as fold enrichment...')
      # enrich_df = enrich_df %>% dplyr::rename(FoldEnrich = eval(parse(text = answer)))
    }
  }


  ## find GeneRatio col
  if(any(grepl('generatio', to_check))){
    colnames(enrich_df)[grepl('generatio', to_check)] = 'GeneRatio'
    enrich_df = enrich_df %>% dplyr::mutate(GeneRatio = sapply(.$GeneRatio, function(x) eval(parse(text = x))))
  }else{
    # gsea
    if(any(grepl('setsize', to_check))){
      colnames(enrich_df)[grepl('setsize', to_check)] = 'setSize'
      enrich_df = enrich_df %>% dplyr::mutate(GeneRatio = as.numeric(Count)/as.numeric(setSize))
    }else if(any(grepl('\\([1-9]{,4}\\)', colnames(enrich_df)))){
      # panther result
     setsize = colnames(enrich_df)[grepl('\\([1-9]{,4}\\)', colnames(enrich_df))] %>%
       stringr::str_remove(.,'.*\\(') %>% stringr::str_remove(.,'\\)') %>% as.numeric()
     enrich_df = enrich_df %>% dplyr::mutate(GeneRatio = as.numeric(Count)/setsize)
    }else if( apply(enrich_df, 2, function(x) length( unique(x)) == 1)){
      setsize = enrich_df[1,apply(enrich_df, 2, function(x) length( unique(x)) == 1)] %>%
        stringr::str_remove('0') %>% as.numeric() %>% sort() %>% .[1]
      enrich_df = enrich_df %>% dplyr::mutate(GeneRatio = as.numeric(Count)/setsize)
    }
  }

  ## find pvalue col
  if(any(grepl('pvalue', to_check))){
    colnames(enrich_df)[grepl('pvalue', to_check)] = 'pvalue'
  }else if(any(grepl('\\buncorrectedpvalue\\b', to_check))){
    colnames(enrich_df)[grepl('\\buncorrectedpvalue\\b', to_check)] = 'pvalue'
  }

  ## find p.adjust col
  if(any(grepl('p.adjust', to_check))){
    colnames(enrich_df)[grepl('p.adjust', to_check)] = 'p.adjust'
  }else if(any(grepl('\\bcorrectedpvalue\\b', to_check))){
    colnames(enrich_df)[grepl('\\bcorrectedpvalue\\b', to_check)] = 'p.adjust'
  }

  ## find qvalue col
  if(any(grepl('qvalue', to_check))){
    colnames(enrich_df)[grepl('qvalue', to_check)] = 'qvalue'
  }else if(any(grepl('fdr', to_check) & !grepl('fdrrate', to_check))){
    colnames(enrich_df)[grepl('fdr', to_check)] = 'qvalue'
  }

  return(enrich_df)
}



