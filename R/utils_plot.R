# AnnoGenes utilities for plotting

#---GO & KEGG enrichment barplot  ---#
plotEnrichBar <- function(enrich_df,
                          xlab_type = c('GeneRatio','Count','FoldEnrich'),
                          legend_by = c("pvalue", "p.adjust", "qvalue"),
                          show_item = 10,
                          ...){
  #--- args ---#
  xlab_type = match.arg(xlab_type)
  legend_by = match.arg(legend_by)

  x = check_type(enrich_df, xlab_type)

  #--- codes ---#







}
