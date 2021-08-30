#' Export list of datasets into different Excel sheets
#'
#' @param dat_list List of datasets.
#' @param name_list List of data names.
#' @param filename A character string naming an xlsx file.
#' @param dir A character string naming output directory.
#' @param overwrite If TRUE, overwrite any existing file.
#' @importFrom stringr str_detect
#'
#' @return An Excel file.
#' @export
#' @examples
#' \donttest{
#' library(openxlsx)
#' expoSheet(dat_list =  list(mtcars,ToothGrowth), name_list = list('mtcars','tooth'),
#'   filename = 'test.xlsx', dir = tempdir())
#' }
#'

expoSheet <- function(dat_list,
                      name_list,
                      filename = NULL,
                      dir = tempdir(),
                      overwrite = TRUE) {

  #--- args ---#
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package openxlsx needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(!"wb" %in% ls()) wb <- createWorkbook()

  if(length(dat_list) != length(name_list)){
    stop('Datasets number is not equal with names!')
  }

  if(is.null(filename)) filename='test.xlsx'

  #--- codes ---#
  invisible(lapply(seq_along(dat_list), function(i){
    openxlsx::addWorksheet(wb, name_list[[i]])
    openxlsx::writeData(wb, sheet = name_list[[i]], x = dat_list[[i]])

    ## styling sheet
    headerStyle <- openxlsx::createStyle(textDecoration = "Bold")
    openxlsx::addStyle(wb, sheet = name_list[[i]], style = headerStyle,
                       rows = 1, cols = seq_len(ncol(dat_list[[i]])), gridExpand = TRUE)
    openxlsx::setColWidths(wb, sheet =name_list[[i]], cols = seq_len(ncol(dat_list[[i]])), widths = "auto")
  }))

  saveWorkbook(wb, paste0(dir,filename), overwrite)
}
