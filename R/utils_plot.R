# AnnoGenes utilities for plotting

#---GO & KEGG enrichment  dotplot--#
plotEnrichDot <- function(enrich_df,
                          xlab_type = c('GeneRatio','Count','FoldEnrich'),
                          legend_by = c("pvalue", "p.adjust", "qvalue"),
                          show_item = 10,
                          ...){
  #--- args ---#
  stopifnot(is.numeric(show_item))
  xlab_type = match.arg(xlab_type)
  legend_by = match.arg(legend_by)

  enrich_df_bk1 = .check_colname(enrich_df, xlab_type)
  enrich_df = .check_colname(enrich_df_bk1, legend_by)


  #--- codes ---#
  xlab = ifelse(xlab_type == 'FoldEnrich', "Fold Enrichment",
                ifelse(xlab_type == 'GeneRatio', "Gene Ratio", "Count"))

  # Panther GO result
  check = enrich_df %>% dplyr::pull(1) %>% stringr::str_detect('.*\\(GO')
  if(any(check)){
    enrich_df <- enrich_df %>%
      dplyr::rename('Description' = names(.)[1]) %>%
      dplyr::rename('Count' = names(.)[3]) %>%
      dplyr::mutate(Description = stringr::str_to_sentence(.$Description) %>%
                      stringr::str_replace_all(., 'go:','GO:'))

  }else{
    # clusterP  GO result
    enrich_df <- enrich_df %>%
      dplyr::mutate(Description = stringr::str_to_sentence(.$Description)) %>%
      dplyr::mutate(GeneRatio = sapply(.$GeneRatio, function(x)eval(parse(text = x))))
  }

  if(show_item <= nrow(enrich_df)){
    enrich_df <- enrich_df %>% dplyr::arrange(eval(parse(text = xlab_type))) %>%
      dplyr::mutate(Description = factor(.$Description,levels = .$Description,ordered = T)) %>%
      dplyr::slice_head(.,n=show_item)
  }else{
    enrich_df <- enrich_df %>% dplyr::arrange(eval(parse(text = xlab_type))) %>%
      dplyr::mutate(Description = factor(.$Description,levels = .$Description,ordered = T))
  }

  ggplot(enrich_df,aes(x = eval(parse(text = xlab_type)),y = Description))+
    geom_point(aes(color =  eval(parse(text = legend_by)),
                   size = Count))+
    scale_color_gradient(low = "red", high = "green")+
    xlab(xlab)+
    theme_bw()+
    guides( color = guide_colorbar(reverse = TRUE))+
    labs(color = legend_by)

}
