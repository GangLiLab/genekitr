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
#' @importFrom dplyr select filter arrange pull mutate
#' @importFrom stringr str_split
#' @importFrom stringi stri_remove_empty_na
#' @importFrom clusterProfiler GSEA
#' @importFrom rlang .data
#'
#' @return GSEA list
#' @export
#'
#' @examples
#' \donttest{
#' data(geneList, package = "genekitr")
#' gse <- genGSEA(genelist = geneList, org = "human",
#'   category = "H",use_symbol = TRUE)
#' }
#'
genGSEA <- function(genelist,
                    org,
                    category = c("C1", "C2", "C3", "C4",
                                 "C5", "C6", "C7", "C8", "H"),
                    subcategory = NULL,
                    use_symbol = TRUE,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    ...) {

  #--- args ---#
  category <- match.arg(category)

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
    geneset <- geneset %>%
      dplyr::select(gs_name, gene_symbol)
  } else if (any(names(genelist) %in% geneset$entrez_gene)) {
    geneset <- geneset %>%
      dplyr::select(gs_name, entrez_gene)
  } else {
    names(genelist) <- transId(names(genelist), transTo = "entrez", org)
    geneset <- geneset %>%
      dplyr::select(gs_name, entrez_gene)
  }

  egmt <- suppressWarnings(clusterProfiler::GSEA(genelist,
    TERM2GENE = geneset,
    pvalueCutoff = pvalueCutoff,
    verbose = F,
    ...
  ))

  exponent <-  egmt@params[["exponent"]]

  egmt =  egmt %>% as.data.frame() %>% as.enrichdat()
  if( use_symbol){
    # transform id to symbol
    egmt_id = stringr::str_split(egmt$geneID,'\\/') %>% unlist()
    id_all = suppressMessages(transId(egmt_id,'symbol',org = org))

    new_geneID <- stringr::str_split(egmt$geneID, "\\/") %>%
      lapply(., function(x) {
        id_all %>% dplyr::filter(input_id %in% x) %>%
          dplyr::arrange(match(input_id, x)) %>%
          dplyr::pull(symbol)
      }) %>%
      sapply(., paste0, collapse = "/")

    egmt =  egmt %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol,.after = geneID)
  }

  res <- list(genelist = genelist, geneset = geneset, gsea_df = egmt, exponent = exponent, org = org)

  return(res)
}

getMsigdb <- function(org,
                      category = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "H"),
                      subcategory = NULL) {

  #--- args ---#
  if (!requireNamespace("msigdbr", quietly = TRUE)) utils::install.packages("msigdbr")
  org <- tolower(org)
  if (org == "hg" | org == "hsa" | org == "hs" | org == "homo sapiens") org <- "human"
  if (org == "mm" | org == "mmu") org <- "mouse"

  # org
  msigOrg <- msigdb_org_data()
  rm(msig_org, envir = .GlobalEnv)
  all_org <- c(
    msigOrg[, 1],
    stringr::str_split(msigOrg[, 2], ", ", simplify = T) %>%
      as.character() %>%
      stringi::stri_remove_empty_na()
  )
  if (!org %in% tolower(all_org)) stop("Choose a valid organism!\n\n", paste0(all_org, " | "))

  # category
  if (!category %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "H")) {
    stop("Category should choose from: C1, C2, C3, C4, C5, C6, C7, C8, H...")
  } else {
    category <- match.arg(category)
  }

  # subcategory
  msigCategory <- msigdb_category_data()
  rm(msig_category, envir = .GlobalEnv)
  all_sub <- msigCategory[, 2] %>%
    stringi::stri_remove_empty_na()

  som_sub <- msigCategory %>%
    dplyr::filter(gs_cat == category) %>%
    dplyr::pull(gs_subcat)

  if (is.null(subcategory)) {
    if (som_sub == "") {
      message(paste0(category, " has no subcategory, continue..."))
      subcategory <- ""
    } else {
      stop("choose a valid subcategory for ", category, "...\n", paste0(som_sub, " | "))
    }
  } else if (!subcategory %in% som_sub) {
    stop("choose a valid subcategory for ", category, "...\n", paste0(som_sub, " | "))
  }


  #--- codes ---#
  msigdb <- msigdbr::msigdbr(org, category, subcategory) %>%
    dplyr::select(., c("gs_name", "gene_symbol", "entrez_gene")) %>%
    as.data.frame()

  return(msigdb)
}



utils::globalVariables(c("gs_name","gene_symbol","entrez_gene","input_id","symbol",
                         "msig_org","msig_category","gs_cat","gs_subcat"))





