# utilities to AnnoGenes

##' @title Show NCBI database searchable name
##' @param db a character vector. Can be "pubmed" or one or more of `rentrez::entrez_dbs()` result.
##' @return a dataframe including search keyword and information.
##' @importFrom rentrez entrez_db_searchable
##' @export
##' @examples
##' \donttest{
##' showNCBI("pubmed")
##' }
showNCBI <- function(db = "pubmed") {
  # suppress binding notes
  fields <- rentrez::entrez_db_searchable(db)
  res <- as.data.frame(fields)[1:3]

  if (nrow(res) == 0) { # nocov start
    message("Something is wrong in your input, NULL will be returned, please check.")
    return(NULL)
  } # nocov end
  return(res)
}

##' @title Get Msigdb database term and gene information
##' @param org organism name from `msigdb_org_data()`.
##' @param category MSigDB collection abbreviation, C1 to C8 and H.
##' @param subcategory MSigDB sub-collection abbreviation, such as REACTOME or BP.
##' @return a dataframe of 2 columns with term and gene.
##' @importFrom stringr str_split
##' @importFrom dplyr %>% filter pull
##' @importFrom stringi stri_remove_empty_na
##' @export
##' @examples
##' \donttest{
##' msigdb <- getMsigdb(org='human', category='C5',subcategory='GO:CC')
##' }
getMsigdb <- function(org,
                      category = c('C1','C2','C3','C4','C5','C6','C7','C8','H'),
                      subcategory=NULL,
                      ...) {

  #--- args ---#
  # org
  msig_org <- msigdb_org_data()
  all_org = c(msig_org[,1],
                  stringr::str_split(msig_org[,2],', ',simplify = T) %>%
                    as.character() %>%
                    stringi::stri_remove_empty_na())
  if (!org %in% all_org) stop("choose a valid org...")

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
  msigdb <- msigdbr(org, category, subcategory) %>%
    dplyr::select(., c("gs_name","gene_symbol","entrez_gene")) %>%
    as.data.frame()

  return(msigdb)

}

msigdb_org_data <- function() {
  utils::data(list="msig_org", package="AnnoGenes")
  get("msig_org", envir = .GlobalEnv)
}
msigdb_category_data <- function() {
  utils::data(list="msig_category", package="AnnoGenes")
  get("msig_category", envir = .GlobalEnv)
}


##' export result into different sheets
##' @param wb worksheet from `createWorkbook()`.
##' @param sheet_dat dataframe added to sheet.
##' @param sheet_name name of added dataframe.
##' @return a worksheet including many dataframes.
##' @importFrom stringr str_detect
##' @importFrom openxlsx addWorksheet writeData writeFormula createStyle addStyle setColWidths
##' @export
##' @examples
##' \donttest{
##' expo_sheet(wb, sheet_dat =  mtcars, sheet_name = 'mtcars')
##' }
expo_sheet <- function(wb, sheet_dat, sheet_name) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = sheet_dat)

  ## add hyperlink for gene
  check <- apply(sheet_dat, 2, function(x) {
    any(stringr::str_detect(x, "http"))
  })
  check[is.na(check)] <- FALSE
  if (any(check)) {
    sub_dat <- sheet_dat[, check]

    for (i in 1:nrow(sheet_dat)) {
      for (n in which(check)) {
        if (sheet_dat[i, n] != "no_gene_id" & sheet_dat[i, n] != "no_uniprot_id") {
          sheet_datd_link <- paste0("HYPERLINK(\"", sheet_dat[i, n], "\", \"", sheet_dat[i, n], "\")")
          writeFormula(wb, sheet = sheet_name, startRow = i + 1, startCol = n, x = sheet_datd_link)
        }
      }
    }
  }

  ## add hyperlink for pubmed
  check2 <- any(stringr::str_detect(colnames(sheet_dat), "pmid"))
  if (check2) {
    for (i in seq_len(nrow(sheet_dat))) {
      if (sheet_dat[i, 2] != "NA") {
        sheet_datd_link <- paste0(
          "HYPERLINK(\"", paste0("https://pubmed.ncbi.nlm.nih.gov/", sheet_dat[i, 5]),
          "\", \"", sheet_dat[i, 2], "\")"
        )
        writeFormula(wb, sheet = sheet_name, startRow = i + 1, startCol = 2, x = sheet_datd_link)
      }
    }
  }

  ## styling sheet
  headerStyle <- createStyle(textDecoration = "Bold")
  addStyle(wb, sheet = sheet_name, style = headerStyle, rows = 1, cols = seq_len(ncol(sheet_dat)), gridExpand = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(sheet_dat)), widths = "auto")

  invisible(wb)
}



