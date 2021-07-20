#' Gene enrichment of GO analysis
#'
#' @param id a gene vector which can be entrez, ensembl or symbol.
#' @param org  organism name from `biocOrg_name`.
#' @param ont  One of "bp", "mf", and "cc" subontologies, or "all" for all three.
#' @param use_symbol logical to output as gene symbol, default is TRUE.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param pvalueCutoff adjusted pvalue cutoff, default is 0.05.
#' @param qvalueCutoff adjusted pvalue cutoff, default is 0.1.
#' @param minGSSize minimal size of each geneSet for analyzing, default is 10.
#' @param maxGSSize maximal size of each geneSet for analyzing, default is 500.
#' @param universe background genes. If missing, the orgdb all gene list will be used as background.
#' @return a dataframe of gene info.
#' @importFrom dplyr  %>% mutate
#' @importFrom clusterProfiler enrichGO
#' @importFrom stringr str_split
#' @export
#' @examples
#' \dontrun{
#' data(geneList, package="DOSE")
#' id = names(geneList)[1:100]
#' ego <- genGO(id, org = 'human',ont = 'mf',pvalueCutoff = 0.05,qvalueCutoff = 0.1 ,use_symbol = T)
#' head(ego)
#' }
genGO <- function(id,
                  org,
                  ont,
                  use_symbol = TRUE,
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
  pkg=paste0("org.", org, ".eg.db")
  keyType = .gentype(id, org)

  #--- codes ---#
  ego <- suppressMessages(
    clusterProfiler::enrichGO(gene = id, OrgDb = pkg, keyType = keyType, ont = toupper(ont),
                              pvalueCutoff = pvalueCutoff,
                              pAdjustMethod = pAdjustMethod,
                              universe = universe,
                              qvalueCutoff = qvalueCutoff,
                              minGSSize = minGSSize,
                              maxGSSize  = maxGSSize))

  if(nrow(as.data.frame(ego)) == 0){
    stop('No GO terms enriched ...')
  }

  if( use_symbol){
    info = genInfo(id,org)
    new_geneID = stringr::str_split(ego$geneID,'\\/') %>%
      lapply(., function(x) {
        info[x,'symbol']
      }) %>% sapply(., paste0, collapse = "/")
    new_ego =  ego %>% as.data.frame() %>%
      dplyr::mutate(geneID = new_geneID) %>% calcFoldEnrich()

  }else{
    new_ego = ego %>% as.data.frame() %>% calcFoldEnrich()
  }

  return(new_ego)
}
