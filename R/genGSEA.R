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
#' @param pvalueCutoff pvalue cutoff, default is 0.05.
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
#' gse <- genGSEA(
#'   genelist = geneList,
#'   org = "human",
#'   category = "H"
#' )
#' }
#'
genGSEA <- function(genelist,
                    org,
                    category = c(
                      "C1", "C2", "C3", "C4",
                      "C5", "C6", "C7", "C8", "H"
                    ),
                    subcategory = NULL,
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
  ens_org <- msig_org[with(msig_org, grepl(org, paste(species_name, species_common_name),ignore.case = T)),1] %>% mapEnsOrg()
  geneset <- getMsigdb(org, category, subcategory)
  id <- names(genelist)
  keyType <- gentype(id = id, org = org)

  # use entrez id or symbol (if not, trans to entrez)
  id <- names(genelist)
  old_id <- id
  if (keyType == "SYMBOL") {
    id_dat1 <- suppressMessages(transId(id, "symbol", ens_org, unique = T))
    id <- id_dat1 %>% dplyr::pull(symbol)
  }

  if (keyType != "ENTREZID") {
    message(paste0(keyType), " gene will be mapped to entrez id")
    id_dat2 <- suppressMessages(transId(id, "entrezid", ens_org, unique = T))
    entrez_id <- id_dat2 %>% dplyr::pull(entrezid)
  } else {
    entrez_id <- id
  }

  names(genelist) <- entrez_id
  geneset <- geneset %>%
    dplyr::select(gs_name, entrez_gene)

  egmt <- suppressWarnings(clusterProfiler::GSEA(genelist,
    TERM2GENE = geneset,
    pvalueCutoff = pvalueCutoff,
    verbose = F,
    ...
  ))

  exponent <- egmt@params[["exponent"]]

  egmt <- egmt %>%
    as.data.frame() %>%
    as.enrichdat() %>%
    dplyr::select(-GeneRatio)


  #--- get geneID_symbol ---#
  # new id
  new_geneID <- get_symbol(egmt$geneID,ens_org)

  # input SYMBOL and no alias
  if( keyType == "SYMBOL" & identical(old_id, id) ){
    egmt <- egmt %>%
      dplyr::mutate(geneID = new_geneID) %>%
      # dplyr::relocate(geneID, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else if(keyType == "SYMBOL" & !identical(old_id, id)){
    old_geneID <- replace_id(id_dat2,egmt$geneID) %>%
      replace_id(id_dat1,.)

    egmt <- egmt %>%
      dplyr::mutate(geneID = old_geneID) %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else if(keyType == "ENTREZID"){
    egmt <- egmt %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else {
    # input other types (including alias)
    old_geneID <- replace_id(id_dat2,egmt$geneID)

    egmt <- egmt %>%
      dplyr::mutate(geneID = old_geneID) %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()
  }



  #--- save as list ---#
  egmt = egmt %>% dplyr::select(-GeneRatio)
  genelist_df = data.frame(ID = names(genelist), logfc = genelist)
  exponent = data.frame(exponent = exponent)
  org = data.frame(org = org)

  res <- list(gsea_df = egmt, genelist = genelist_df, geneset = geneset,  exponent = exponent, org = org)

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



utils::globalVariables(c(
  "gs_name", "gene_symbol", "entrez_gene", "input_id", "symbol",
  "msig_org", "msig_category", "gs_cat", "gs_subcat"
))
