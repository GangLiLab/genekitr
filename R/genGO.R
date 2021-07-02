##' GO enrichment analysis for gene id
##'
##' @param id a gene vector which can be entrez, ensembl or symbol.
##' @param org  organism name from `biocOrg_name()`.
##' @param ont  One of "bp", "mf", and "cc" subontologies, or "all" for all three.
##' @param readable logical to output as gene symbol, default is TRUE.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
##' @param pvalueCutoff adjusted pvalue cutoff, default is 0.05.
##' @param qvalueCutoff adjusted pvalue cutoff, default is 0.1.
##' @param minGSSize minimal size of each geneSet for analyzing, default is 10.
##' @param maxGSSize maximal size of each geneSet for analyzing, default is 500.
##' @param universe background genes. If missing, the orgdb all gene list will be used as background.
##' @return a dataframe of gene info.
##' @importFrom dplyr pull
##' @importFrom stringr str_to_title
##' @importFrom clusterProfiler enrichGO
##' @importFrom DOSE setReadable
##' @export
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' id = names(geneList)[1:100]
##' ego <- genGO(id, org = 'human',ont = 'mf',pvalueCutoff = 0.05,qvalueCutoff = 0.1 ,readable = T)
##' head(ego)
##' }
genGO <- function(id,
                  org,
                  ont,
                  readable = TRUE,
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

  org_bk = org
  org = mapBiocOrg(tolower(org))
  .load_orgdb(org)

  keyType = .gentype(id, org)
  if(! keyType %in% c('SYMBOL','ENSEMBL','ENTREZID')) {
    stop('Gene id type should be one of: SYMBOL, ENSEMBL and ENTREZID')
  }

  if(! .genInorg(id,org) ) stop('Gene id is not matched with organism: ', org_bk, ' !')

  #--- codes ---#
  ego <- suppressMessages(
    clusterProfiler::enrichGO(gene = id, OrgDb = pkg, keyType = keyType, ont = toupper(ont),
                              pvalueCutoff = pvalueCutoff,
                              pAdjustMethod = pAdjustMethod,
                              universe = universe,
                              qvalueCutoff = qvalueCutoff,
                              minGSSize = minGSSize,
                              maxGSSize  = maxGSSize)
  )
  if( readable | keyType != 'SYMBOL'){
    ego <- DOSE::setReadable(ego, OrgDb = pkg)
  }

  return(ego)

}
