#' GSEA for a gene list with decreasing logFC value
#'
#' @param genelist Order ranked genelist in decreasing order, gene can be
#'   entrez, ensembl or symbol.
#' @param org Organism name from `msig_org`.
#' @param category MSigDB collection abbreviation, one of C1','C2','C3',
#'   'C4','C5','C6','C7','C8','H'.
#' @param subcategory MSigDB sub-collection abbreviation, choose from
#'   `msig_category`.
#' @param minGSSize Minimal size of each geneSet for analyzing, default is 10.
#' @param maxGSSize Maximal size of each geneSet for analyzing, default is 500.
#' @param pvalueCutoff Adjusted pvalue cutoff, default is 0.05.
#' @param ... Other argument to `GSEA` function
#' @importFrom dplyr select filter pull mutate %>%
#' @importFrom stringr str_split
#' @importFrom clusterProfiler GSEA
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' \donttest{
#' data(geneList, package="genekitr")
#' gse = genGSEA(genelist = geneList,org = 'human', category='H')
#' }


genGSEA <- function(genelist,
                    org,
                    category = c('C1','C2','C3','C4','C5','C6','C7','C8','H'),
                    subcategory = NULL,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    ...){

  #--- args ---#
  category = match.arg(category)

  stopifnot(
    is.numeric(minGSSize),
    is.numeric(maxGSSize),
    is.numeric(pvalueCutoff)
  )

  if (is.unsorted(rev(genelist))) stop("genelist should be a decreasing sorted vector...")

  #--- codes ---#
  geneset <- getMsigdb(org, category, subcategory)

  # use entrez id or symbol
  if (any(names(genelist) %in% geneset$gene_symbol)) {
    geneset = geneset %>%
      dplyr::select(gs_name,gene_symbol)
  }else if (any(names(genelist) %in% geneset$entrez_gene)) {
    geneset = geneset %>%
      dplyr::select(gs_name,entrez_gene)
  }else{
    names(genelist) =transId(names(genelist),trans_to = 'entrez',org)
    geneset = geneset %>%
      dplyr::select(gs_name,entrez_gene)
  }

  egmt <- suppressWarnings(clusterProfiler::GSEA(genelist, TERM2GENE=geneset,
                                                 pvalueCutoff, verbose=F,
                                                 ...))

  # if( use_symbol){
  #   info = genInfo(names(genelist),org,unique = T) %>% na.omit()
  #   new_geneID = stringr::str_split(egmt$geneID,'\\/') %>%
  #     lapply(., function(x) {
  #       info %>% dplyr::filter(input_id %in% x) %>% dplyr::pull(symbol)
  #     }) %>% sapply(., paste0, collapse = "/")
  #   new_egmt =  egmt %>% as.data.frame() %>%
  #     dplyr::mutate(core_enrichment = new_geneID)
  #
  # }else{
  #   new_egmt = egmt %>% as.data.frame()
  # }
  #
  new_egmt = egmt %>% as.data.frame()
  return(new_egmt)

}
