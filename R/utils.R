# utilities to AnnoGenes

##' @title Show NCBI database searchable name
##' @param db a character vector. Can be "pubmed" or one or more of `rentrez::entrez_dbs()` result.
##' @return a dataframe including search keyword and information.
##' @importFrom rentrez entrez_db_searchable
##' @importFrom stringr str_detect
##' @importFrom openxlsx addWorksheet writeData writeFormula createStyle addStyle setColWidths
##' @export
##' @examples
##' \donttest{
##' showNCBI("pubmed")
##' }
showNCBI <- function(db = "pubmed") {
  # suppress binding notes
  fields <- rentrez::entrez_db_searchable(db)
  res=as.data.frame(fields)[1:3]

  if (nrow(res) == 0) { # nocov start
    message("Something is wrong in your input, NULL will be returned, please check.")
    return(NULL)
  } # nocov end
  return(res)
}


# export result into different sheets
expo_sheet <- function(wb, sheet_dat, sheet_name){
  addWorksheet(wb, sheet_name)
  writeData(wb,sheet = sheet_name, x = sheet_dat)

  ## if needs to add hyperlink
  check = apply(sheet_dat, 2, function(x){ any(stringr::str_detect(x,'http'))})
  check[is.na(check)]=FALSE
  if(any(check)){
    sub_dat = sheet_dat[,check]

    for(i in 1:nrow(sheet_dat)) {
      for(n in which(check)){
        if(sheet_dat[i,n] != 'NA' & sheet_dat[i,n] != 'no_uniprot_id'){
          sheet_datd_link = paste0("HYPERLINK(\"",sheet_dat[i,n],"\", \"",sheet_dat[i,n],"\")")
          writeFormula(wb, sheet =sheet_name, startRow = i+1, startCol = n, x = sheet_datd_link)
        }
      }
    }
  }

  ## styling sheet
  headerStyle <- createStyle(textDecoration = "Bold")
  addStyle(wb, sheet= sheet_name,style = headerStyle, rows = 1, cols = seq_len(ncol(sheet_dat)), gridExpand = TRUE)
  setColWidths(wb, sheet= sheet_name, cols = seq_len(ncol(sheet_dat)), widths = "auto")

  invisible(wb)
}


.nm <- function(x) deparse(substitute(x))















