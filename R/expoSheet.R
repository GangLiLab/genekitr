#' Export list of data sets into different 'Excel' sheets
#'
#' @param data_list List of datasets.
#' @param data_name Character of data names.
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
#' library(openxlsx)
#' expoSheet(
#'   data_list = list(mtcars, ToothGrowth),
#'   data_name = c("mtcars", "tooth"),
#'   filename = "test.xlsx", dir = tempfile()
#' )
expoSheet <- function(data_list,
                      data_name,
                      filename = NULL,
                      dir = tempdir(),
                      overwrite = TRUE) {

  #--- args ---#
  if (!"wb" %in% ls()) wb <- openxlsx::createWorkbook()

  if (length(data_list) != length(data_name)) {
    stop("Datasets number is not equal with names!")
  }

  if (is.null(filename)) filename <- "test.xlsx"

  #--- codes ---#
  invisible(lapply(seq_along(data_list), function(i) {
    # i = 1
    openxlsx::addWorksheet(wb, data_name[i])
    openxlsx::writeData(wb, sheet = data_name[i], x = data_list[[i]])

    ## styling sheet
    headerStyle <- openxlsx::createStyle(textDecoration = "Bold")
    if(!is.null(ncol(data_list[[i]]))){
      openxlsx::addStyle(wb,
                         sheet = data_name[i], style = headerStyle,
                         rows = 1, cols = seq_len(ncol(data_list[[i]])), gridExpand = TRUE
      )
    }

    if(!is.null(ncol(data_list[[i]]))){
      openxlsx::setColWidths(wb,
                             sheet = data_name[i],
                             cols = seq_len(ncol(data_list[[i]])), widths = "auto"
      )
    }

  }))

  openxlsx::saveWorkbook(wb, paste0(dir, filename), overwrite)
}
