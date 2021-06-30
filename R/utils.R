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


##' export result into different sheets
##' @param wb worksheet from `createWorkbook()`.
##' @param sheet_dat dataframe added to sheet.
##' @param sheet_name name of added dataframe.
##' @return a worksheet including many dataframes.
##' @importFrom rentrez entrez_db_searchable
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


# .nm <- function(x) deparse(substitute(x))
