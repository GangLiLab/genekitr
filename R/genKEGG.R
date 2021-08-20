#' Gene enrichment of KEGG analysis
#'
#' @param id a vector of entrez gene.
#' @param org  character of organism name which can test by `mapKeggOrg()`.
#' @param use_symbol logical to output as gene symbol, default is TRUE.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH",
#'   "BY", "fdr", "none".
#' @param pvalueCutoff numberic of adjusted pvalue cutoff, default is 0.05.
#' @param qvalueCutoff numberic of adjusted pvalue cutoff, default is 0.1.
#' @param minGSSize numberic of minimal size of each geneSet for analyzing,
#'   default is 10.
#' @param maxGSSize numberic of maximal size of each geneSet for analyzing,
#'   default is 500.
#' @param universe background genes. If missing, the orgdb all gene list will be
#'   used as background.
#' @return a dataframe of gene info.
#' @importFrom dplyr  %>%
#' @importFrom stringi stri_omit_na
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom DOSE setReadable
#' @export
#' @examples
#' \dontrun{
#' data(geneList, package="DOSE")
#' id = names(geneList)[1:100]
#' keg <- genKEGG(id, org = 'human')
#' head(keg)
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
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)
  options(warn=-1)
  stopifnot(is.character(id))
  if (missing(universe)) universe <- NULL

  org.bk = org
  org = mapKeggOrg(tolower(org))
  keyType = .gentype(id, org)

  if(! keyType %in% c('ENTREZID') ) {
    message(paste0(keyType), ' gene will be mapped to entrez id')
    trans_id = suppressMessages(transId(id,'entrezid',org.bk)) %>% stringi::stri_remove_na()
  }else{
    trans_id = id
  }

  info = genInfo(trans_id,org) %>% dplyr::mutate(entrezid := rownames(.))

  #--- codes ---#
  keg <- suppressMessages(
    clusterProfiler::enrichKEGG(gene = trans_id, organism = org, keyType = 'kegg',
                                pvalueCutoff = pvalueCutoff,
                                pAdjustMethod = pAdjustMethod,
                                qvalueCutoff = qvalueCutoff,
                                universe  = universe,
                                minGSSize = minGSSize,
                                maxGSSize = maxGSSize ))

  if( use_symbol ){
    new_geneID = stringr::str_split(keg$geneID,'\\/') %>%
      lapply(., function(x) {
        info[x,'symbol']
      }) %>% sapply(., paste0, collapse = "/")
    new_keg =  keg %>% as.data.frame() %>%
      dplyr::mutate(geneID = new_geneID) %>% calcFoldEnrich()

  }else{
    new_geneID = stringr::str_split(keg$geneID,'\\/') %>%
      lapply(., function(x) {
        info %>%
          dplyr::filter(entrezid %in% x) %>% dplyr::pull(tolower(keyType))
      }) %>% sapply(., paste0, collapse = "/")
    new_keg =  keg %>% as.data.frame() %>%
      dplyr::mutate(geneID = new_geneID) %>% calcFoldEnrich()
  }

  return(new_keg)

}

