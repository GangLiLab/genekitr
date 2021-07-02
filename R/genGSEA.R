##' GSEA for a genelist with decreasing logFC value
##'
##' @param genelist order ranked genelist in decreasing order, gene can be entrez or symbol.
##' @param geneset  gene set collections from `getMsigdb()`.
##' @param minGSSize minimal size of each geneSet for analyzing, default is 10.
##' @param maxGSSize maximal size of each geneSet for analyzing, default is 500.
##' @param pvalueCutoff adjusted pvalue cutoff, default is 0.05.
##' @return a dataframe of gene info.
##' @importFrom dplyr select
##' @importFrom clusterProfiler GSEA
##' @export
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' msigdb <- getMsigdb(org='human', category='C3',subcategory = 'TFT:GTRD')
##' egmt <- genGSEA(genelist = geneList,geneset = msigdb)
##' egmt2 <- DOSE::setReadable(egmt, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
##' }
genGSEA <- function(genelist,
                    geneset,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    ...){

  #--- args ---#
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)
  options(warn=-1)
  stopifnot(
    is.data.frame(geneset),
    is.numeric(minGSSize),
    is.numeric(maxGSSize),
    is.numeric(pvalueCutoff)
  )

  if (is.unsorted(rev(genelist)))
    stop("genelist should be a decreasing sorted vector...")

  #--- codes ---#
  # gene id or symbol
  if (any(names(genelist) %in% geneset$gene_symbol)) {
    geneset = geneset %>%
      dplyr::select(gs_name,gene_symbol)
  }else{
    geneset = geneset %>%
      dplyr::select(gs_name,entrez_gene)
  }

  egmt <- suppressWarnings(clusterProfiler::GSEA(genelist, TERM2GENE=geneset, pvalueCutoff, verbose=F))

  return(egmt)

}
