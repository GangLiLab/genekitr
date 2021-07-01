##' GSEA for a genelist with logFC
##'
##' @param genelist order ranked genelist in decreasing order, gene can be entrez or symbol.
##' @param species species name from `msigdb_species_data()`.
##' @param category MSigDB collection abbreviation, C1 to C8 and H.
##' @param subcategory MSigDB sub-collection abbreviation, such as REACTOME or BP.
##' @param minGSSize minimal size of each geneSet for analyzing.
##' @param maxGSSize maximal size of each geneSet for analyzing.
##' @return a dataframe of gene info.
##' @importFrom stringr str_to_title
##' @export
##' @examples
##' \dontrun{
##'
##' }
genGSEA <- function(genelist,
                    species,
                    category = c('C1','C2','C3','C4','C5','C6','C7','C8','H'),
                    subcategory = NULL,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    ...){

  #--- args ---#
  options(warn=-1)
  stopifnot(
    is.character(species),
    is.numeric(minGSSize),
    is.numeric(maxGSSize),
    is.numeric(pvalueCutoff)
  )

  if (is.unsorted(rev(genelist)))
    stop("genelist should be a decreasing sorted vector...")

  # species
  msig_species <- msigdb_species_data()
  all_species = c(msig_species[,1],
                  stringr::str_split(msig_species[,2],', ',simplify = T) %>%
                    as.character() %>%
                    stringi::stri_remove_empty_na())
  if (!species %in% all_species) stop("choose a valid species...")

  # category
  if(! category %in% c('C1','C2','C3','C4','C5','C6','C7','C8','H')){
    stop("Category should choose from: C1, C2, C3, C4, C5, C6, C7, C8, H...")
  }else{
    category <- match.arg(category)
  }

  # subcategory
  msig_category <- msigdb_category_data()
  all_sub <- msig_category[,2] %>%
    stringi::stri_remove_empty_na()

  som_sub <- msig_category %>%
    dplyr::filter(gs_cat==category) %>%
    dplyr::pull(gs_subcat)

  if( is.null(subcategory)){
    if(som_sub == ''){
      message(paste0(category,' has no subcategory...'))
      subcategory = ''
    }else{
      stop("choose a valid subcategory for ",category,"...\n",paste0(som_sub,' | '))
    }
  }else if(! subcategory %in% all_sub){
    if(som_sub == ''){
      message(paste0(category,' has no subcategory...'))
      subcategory = ''
    }else{
      stop("choose a valid subcategory for ",category,"...\n",paste0(som_sub,' | '))
    }
  }

  #--- codes ---#
  msigdb <- getMsigdb(species, category, subcategory)
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
