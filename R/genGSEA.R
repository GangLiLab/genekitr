#' GSEA for a genelist with decreasing logFC value
#'
#' @param genelist order ranked genelist in decreasing order, gene can be entrez or symbol.
#' @param org organism name from `msigdb_org_data()`.
#' @param category MSigDB collection abbreviation, C1 to C8 and H.
#' @param subcategory MSigDB sub-collection abbreviation, such as REACTOME or BP.
#' @param use_symbol logical to output as gene symbol, default is TRUE.
#' @param minGSSize minimal size of each geneSet for analyzing, default is 10.
#' @param maxGSSize maximal size of each geneSet for analyzing, default is 500.
#' @param pvalueCutoff adjusted pvalue cutoff, default is 0.05.
#' @return a dataframe of gene info.
#' @importFrom dplyr select
#' @importFrom clusterProfiler GSEA
#' @importFrom DOSE setReadable
#' @export
#' @examples
#' \dontrun{
#' data(geneList, package="DOSE")
#' genGSEA(genelist = geneList,org = 'human', category='C3',subcategory = 'TFT:GTRD',use_symbol = F)
#' }
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
  options(warn=-1)
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

  # gene id or symbol
  if (any(names(genelist) %in% geneset$gene_symbol)) {
    geneset = geneset %>%
      dplyr::select(gs_name,gene_symbol)
  }else{
    geneset = geneset %>%
      dplyr::select(gs_name,entrez_gene)
  }

  egmt <- suppressWarnings(clusterProfiler::GSEA(genelist, TERM2GENE=geneset, pvalueCutoff, verbose=F))

  if(use_symbol){
    biocOrg = mapBiocOrg(tolower(org.bk))
    pkg=paste0("org.", biocOrg, ".eg.db")
    keyType = .gentype(names(genelist), biocOrg)
    egmt <- DOSE::setReadable(egmt, OrgDb = pkg, keyType)
  }

  new_egmt = egmt %>% as.data.frame()

  return(new_egmt)

}
