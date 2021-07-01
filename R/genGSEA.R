##' GSEA for a genelist with logFC
##'
##' @param genelist order ranked genelist in decreasing order, gene can be entrez or symbol.
##' @param msigdb  gene set collections from `getMsigdb()`.
##' @param minGSSize minimal size of each geneSet for analyzing.
##' @param maxGSSize maximal size of each geneSet for analyzing.
##' @param pvalueCutoff adjusted pvalue cutoff.
##' @return a dataframe of gene info.
##' @importFrom stringr str_to_title
##' @export
##' @examples
##' \dontrun{
##'
##' }
genGSEA <- function(genelist,
                    msigdb,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    ...){

  #--- args ---#
  options(warn=-1)
  stopifnot(
    is.character(org),
    is.numeric(minGSSize),
    is.numeric(maxGSSize),
    is.numeric(pvalueCutoff)
  )

  if (is.unsorted(rev(genelist)))
    stop("genelist should be a decreasing sorted vector...")

  #--- codes ---#
  msigdb <- getMsigdb(org, category, subcategory)
  # gene id or symbol
  if (any(names(genelist) %in% msigdb$gene_symbol)) {
    msigdb = msigdb %>%
      dplyr::select(gs_name,gene_symbol)
  }else{
    msigdb = msigdb %>%
      dplyr::select(gs_name,entrez_gene)
  }

  egmt <- suppressWarnings(GSEA(genelist, TERM2GENE =msigdb, pvalueCutoff, verbose=T))
  gsea_results_df <- egmt@result


}
