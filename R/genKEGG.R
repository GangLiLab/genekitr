#' Gene enrichment of KEGG analysis
#'
#' @param id A vector of entrez gene.
#' @param org  KEGG organism name from `keggOrg_name`.
#' @param use_symbol Logical to set result gene id as gene symbol, default is TRUE.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH",
#'   "BY", "fdr", "none".
#' @param pvalueCutoff Numberic of adjusted pvalue cutoff, default is 0.05.
#' @param qvalueCutoff Numberic of adjusted pvalue cutoff, default is 0.1.
#' @param minGSSize Numberic of minimal size of each geneSet for analyzing,
#'   default is 10.
#' @param maxGSSize Numberic of maximal size of each geneSet for analyzing,
#'   default is 500.
#' @param universe Background genes. If missing, the orgdb all gene list will be
#'   used as background.
#' @param ... Other argument to `enrichKEGG` function
#' @importFrom dplyr  %>% mutate filter
#' @importFrom stringr  str_split
#' @importFrom stringi stri_omit_na
#' @importFrom clusterProfiler enrichKEGG
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' \donttest{
#' data(geneList, package="genekitr")
#' id <- names(geneList)[1:100]
#' keg <- genKEGG(id, org = 'human')
#' }

genKEGG <- function(id,
                    org,
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

  org = mapKeggOrg(tolower(org))
  keyType = .gentype(id, org)

  if(! keyType %in% c('ENTREZID') ) {
    message(paste0(keyType), ' gene will be mapped to entrez id')
    trans_id = suppressMessages(transId(id,'entrezid',org)) %>% stringi::stri_remove_na()
  }else{
    trans_id = id
  }

  info = genInfo(trans_id,org,unique = T) %>% na.omit()

  #--- codes ---#
  keg <- suppressMessages(
    clusterProfiler::enrichKEGG(gene = trans_id, organism = org, keyType = 'kegg',
                                pvalueCutoff = pvalueCutoff,
                                pAdjustMethod = pAdjustMethod,
                                qvalueCutoff = qvalueCutoff,
                                universe  = universe,
                                minGSSize = minGSSize,
                                maxGSSize = maxGSSize,
                                ...))

  if(nrow(as.data.frame(keg)) == 0){
    stop('No KEGG terms enriched ...')
  }

  if( use_symbol ){
    new_geneID = stringr::str_split(keg$geneID,'\\/') %>%
      lapply(., function(x) {
        info %>% dplyr::filter(input_id %in% x) %>% dplyr::pull(symbol)
      }) %>% sapply(., paste0, collapse = "/")
    new_keg =  keg %>% as.data.frame() %>%
      dplyr::mutate(geneID = new_geneID) %>% calcFoldEnrich()

  }else{
    new_keg =  keg %>% as.data.frame()  %>% calcFoldEnrich()
  }

  return(new_keg)

}

