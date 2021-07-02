##' KEGG enrichment analysis for gene id
##'
##' @param id a vector of entrez gene.
##' @param org  character of organism name which can test by `mapKeggOrg()`.
##' @param readable logical to output as gene symbol, default is TRUE.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
##' @param pvalueCutoff numberic of adjusted pvalue cutoff, default is 0.05.
##' @param qvalueCutoff numberic of adjusted pvalue cutoff, default is 0.1.
##' @param minGSSize numberic of minimal size of each geneSet for analyzing, default is 10.
##' @param maxGSSize numberic of maximal size of each geneSet for analyzing, default is 500.
##' @param universe background genes. If missing, the orgdb all gene list will be used as background.
##'
##'
genKEGG <- function(id,
                    org,
                    # readable = TRUE,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1,
                    minGSSize = 10,
                    maxGSSize = 500,
                    universe,
                    ...){

  #--- args ---#
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)
  options(warn=-1)
  stopifnot(is.character(id))
  if (missing(universe)) universe <- NULL

  org = mapKeggOrg(tolower(org))
  keyType = .gentype(id, org)
  if(! keyType %in% c('ENTREZID') ) {
    message(paste0(keyType), ' gene will be mapped to entrez id')
    trans_id = mapId(id,tolower(keyType),'entrezid',org)
  }else{
    trans_id = id
  }


  #--- codes ---#
  keg <- suppressMessages(
    clusterProfiler::enrichKEGG(gene = trans_id, organism = org, keyType = 'kegg',
                                pvalueCutoff = pvalueCutoff,
                                pAdjustMethod = pAdjustMethod,
                                qvalueCutoff = qvalueCutoff,
                                universe  = universe,
                                minGSSize = minGSSize,
                                maxGSSize = maxGSSize )
  )

  return(keg)

}
