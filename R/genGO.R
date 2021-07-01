##' GO enrichment analysis for a genelist
##'
##' @param id a gene vector which can be entrez, ensembl or symbol.
##' @param org  organism name from `biocOrg_data()`.
##' @param ont  One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
##' @param readable logical to output as gene symbol, default is TRUE.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
##' @param pvalueCutoff adjusted pvalue cutoff, default is 0.05.
##' @param qvalueCutoff adjusted pvalue cutoff, default is 0.1.
##' @param minGSSize minimal size of each geneSet for analyzing, default is 10.
##' @param maxGSSize maximal size of each geneSet for analyzing, default is 500.
##' @param universe background genes. If missing, the orgdb all gene list will be used as background.
##' @return a dataframe of gene info.
##' @importFrom dplyr select
##' @importFrom clusterProfiler GSEA
##' @export
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' id = names(geneList)[1:100]
##' ego <- genGO(id, org = 'human',ont = 'CC',pvalueCutoff = 0.05,qvalueCutoff = 0.2)
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
  stopifnot(is.character(id))
  if (missing(universe)) universe <- NULL

  org = mapBiocOrg(org)
  if(!org %in% (biocOrg_data() %>% dplyr::pull(short_name))){
    stop('Check organism name! \n USE FULL NAME: ',
         paste0(biocOrg_data() %>% dplyr::pull(full_name),' | '),
         '\n OR USE SHORT NAME: ',
         paste0(biocOrg_data() %>% dplyr::pull(short_name),' | '))
  }
  org <- stringr::str_to_title(org)
  keyType = gentype(id, org)

  #--- codes ---#
  pkg=paste0("org.", org, ".eg.db")
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package ",pkg," is required!")
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))

  ego <- enrichGO(gene = id, OrgDb = pkg, keyType = keyType, ont = ont,
                  pvalueCutoff, pAdjustMethod,universe, qvalueCutoff,
                  minGSSize,maxGSSize )
  if( readable | keyType != 'SYMBOL'){
    ego <- setReadable(ego, OrgDb = pkg)
  }

  return(ego)

}
