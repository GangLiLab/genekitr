#' Export list of datasets into different Excel sheets
#'
#' @param data_list List of datasets.
#' @param name_list List of data names.
#' @param filename A character string naming an xlsx file.
#' @param dir A character string naming output directory.
#' @param overwrite If TRUE, overwrite any existing file.
#' @importFrom rlang .data
#' @importFrom openxlsx createWorkbook addWorksheet writeData createStyle addStyle setColWidths
#' saveWorkbook
#'
#' @return An Excel file.
#' @export
#' @examples
#' \dontrun{
#' library(openxlsx)
#' expoSheet(
#'   data_list = list(mtcars, ToothGrowth),
#'   name_list = list("mtcars", "tooth"),
#'   filename = "test.xlsx", dir = tempdir()
#' )
#' }
expoSheet <- function(data_list,
                      name_list,
                      filename = NULL,
                      dir = tempdir(),
                      overwrite = TRUE) {

  #--- args ---#
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    utils::install.packages("openxlsx")
  }
  if (!"wb" %in% ls()) wb <- openxlsx::createWorkbook()

  if (length(data_list) != length(name_list)) {
    stop("Datasets number is not equal with names!")
  }

  if (is.null(filename)) filename <- "test.xlsx"

  #--- codes ---#
  invisible(lapply(seq_along(data_list), function(i) {
    openxlsx::addWorksheet(wb, name_list[[i]])
    openxlsx::writeData(wb, sheet = name_list[[i]], x = data_list[[i]])

    ## styling sheet
    headerStyle <- openxlsx::createStyle(textDecoration = "Bold")
    openxlsx::addStyle(wb,
      sheet = name_list[[i]], style = headerStyle,
      rows = 1, cols = seq_len(ncol(data_list[[i]])), gridExpand = TRUE
    )
    openxlsx::setColWidths(wb, sheet = name_list[[i]],
                           cols = seq_len(ncol(data_list[[i]])), widths = "auto")
  }))

  openxlsx::saveWorkbook(wb, paste0(dir, filename), overwrite)
}


