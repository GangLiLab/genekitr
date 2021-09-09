#' Gene GO enrichment analysis
#'
#' @param id A vector of gene id which can be entrez, ensembl or symbol.
#' @param org  Organism name from `biocOrg_name`.
#' @param ont  One of "bp", "mf", and "cc" subontologies, or "all" for all
#'   three.
#' @param use_symbol Logical to set result gene id as gene symbol, default is TRUE.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH",
#'   "BY", "fdr", "none".
#' @param pvalueCutoff Adjusted pvalue cutoff, default is 0.05.
#' @param qvalueCutoff Adjusted pvalue cutoff, default is 0.1.
#' @param minGSSize Minimal size of each gene set for analyzing, default is 10.
#' @param maxGSSize Maximal size of each gene set for analyzing, default is 500.
#' @param universe Background genes. If missing, then all gene list in
#'   orgdb will be used as background.
#' @param ... other argument to `enrichGO` function
#' @importFrom dplyr select filter pull mutate %>%
#' @importFrom stringr str_split
#' @importFrom clusterProfiler enrichGO
#'
#' @return A `data.frame` contains gene ratio and fold enrichment.
#' @export
#'
#' @examples
#' \donttest{
#' data(geneList, package="genekitr")
#' id <- names(geneList)[1:100]
#' ego <- genGO(id, org = 'human',ont = 'cc',pvalueCutoff = 0.01,
#'   qvalueCutoff = 0.1 ,use_symbol = FALSE)
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
  stopifnot(is.character(id))
  if (missing(universe)) universe <- NULL

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
                              maxGSSize  = maxGSSize,
                              ...))

  if(nrow(as.data.frame(ego)) == 0){
    stop('No GO terms enriched ...')
  }

  if( use_symbol){
    info = genInfo(id,org,unique = T) %>% na.omit()
    new_geneID = stringr::str_split(ego$geneID,'\\/') %>%
      lapply(., function(x) {
        info %>% dplyr::filter(input_id %in% x) %>% dplyr::pull(symbol)
      }) %>% sapply(., paste0, collapse = "/")
    new_ego =  ego %>% as.data.frame() %>%
      dplyr::mutate(geneID = new_geneID) %>% calcFoldEnrich()

  }else{
    new_ego = ego %>% as.data.frame() %>% calcFoldEnrich()
  }

  return(new_ego)
}
