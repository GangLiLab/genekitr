#' Export result into different Excel sheets
#'
#' @param wb workbook created from `createWorkbook`.
#' @param sheet_dat `data.frame` to export.
#' @param sheet_name sheet name.
#' @importFrom stringr str_detect
#'
#' @return An Excel file.
#' @export
#' @examples
#' \dontrun{
#' library(openxlsx)
#' wb <- createWorkbook()
#' wb <- expo_sheet(wb, sheet_dat =  mtcars, sheet_name = 'mtcars')
#' saveWorkbook(wb, "./test.xlsx", overwrite = TRUE)
#' }
#'
expo_sheet <- function(wb, sheet_dat, sheet_name) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package openxlsx needed for this function to work. Please install it.",
         call. = FALSE)
  }
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = sheet_dat)

  ## add hyperlink for gene
  check <- apply(sheet_dat, 2, function(x) {
    any(stringr::str_detect(x, "http"))
  })
  # check[is.na(check)] <- FALSE
  if (any(check)) {
    sub_dat <- sheet_dat[, check]

    for (i in 1:nrow(sheet_dat)) {
      for (n in which(check)) {
        if (! is.na(sheet_dat[i, n]) & !is.na(sheet_dat[i, n] )) {
          sheet_datd_link <- paste0("HYPERLINK(\"", sheet_dat[i, n],
                                    "\", \"", sheet_dat[i, n], "\")")
          openxlsx::writeFormula(wb, sheet = sheet_name, startRow = i + 1,
                       startCol = n, x = sheet_datd_link)
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
        openxlsx::writeFormula(wb, sheet = sheet_name, startRow = i + 1,
                     startCol = 2, x = sheet_datd_link)
      }
    }
  }

  ## styling sheet
  headerStyle <- openxlsx::createStyle(textDecoration = "Bold")
  openxlsx::addStyle(wb, sheet = sheet_name, style = headerStyle,
           rows = 1, cols = seq_len(ncol(sheet_dat)), gridExpand = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(sheet_dat)), widths = "auto")

  invisible(wb)
}
