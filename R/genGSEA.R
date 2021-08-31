#' GSEA for a gene list with decreasing logFC value
#'
#' @param genelist Order ranked genelist in decreasing order, gene can be
#'   entrez, ensembl or symbol.
#' @param org Organism name from `msig_org`.
#' @param category MSigDB collection abbreviation, one of C1','C2','C3',
#'   'C4','C5','C6','C7','C8','H'.
#' @param subcategory MSigDB sub-collection abbreviation, choose from
#'   `msig_category`.
#' @param use_symbol Logical to set result gene id as gene symbol, default is TRUE.
#' @param minGSSize Minimal size of each geneSet for analyzing, default is 10.
#' @param maxGSSize Maximal size of each geneSet for analyzing, default is 500.
#' @param pvalueCutoff Adjusted pvalue cutoff, default is 0.05.
#' @param ... Other argument to `GSEA` function
#' @importFrom dplyr select
#' @importFrom clusterProfiler GSEA
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' data(geneList, package="DOSE")
#' genGSEA(genelist = geneList,org = 'human', category='C3',
#'   subcategory = 'TFT:GTRD',use_symbol = FALSE)

genGSEA <- function(genelist,
                    org,
                    category = c('C1','C2','C3','C4','C5','C6','C7','C8','H'),
                    subcategory = NULL,
                    use_symbol = TRUE,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    ...){

  #--- args ---#
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)
  category = match.arg(category)

  stopifnot(
    is.numeric(minGSSize),
    is.numeric(maxGSSize),
    is.numeric(pvalueCutoff)
  )

  if (is.unsorted(rev(genelist))) stop("genelist should be a decreasing sorted vector...")

  #--- codes ---#
  org.bk = org
  geneset <- getMsigdb(org.bk, category, subcategory)

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

  if( use_symbol){
    info = genInfo(names(genelist),org)
    new_geneID = stringr::str_split(egmt$core_enrichment,'\\/') %>%
      lapply(., function(x) {
        info[x,'symbol']
      }) %>% sapply(., paste0, collapse = "/")
    new_egmt =  egmt %>% as.data.frame() %>%
      dplyr::mutate(core_enrichment = new_geneID)

  }else{
    new_egmt = egmt %>% as.data.frame()
  }

  return(new_egmt)

}
