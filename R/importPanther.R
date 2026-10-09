#' Import 'Panther' web result
#'
#' @param panther_file Panther result file.
#'
#' @importFrom stringr str_remove
#' @importFrom dplyr mutate rename relocate everything
#' @importFrom rlang .data
#' @return  `data.frame`
#' @export

importPanther <- function(panther_file) {
 if (!requireNamespace("rio", quietly = TRUE)) {
   stop("Package \"rio\" needed for this function to work.
         Please install it by install.packages('rio')",call. = FALSE)

  }

  # unmodified panther result (should be txt file)
  if (getExtension(panther_file) == "txt") {
    # remove blank lines and header lines in base R (shell commands are not available on Windows)
    panther_lines <- readLines(panther_file, warn = FALSE)
    panther_lines <- panther_lines[grepl("[^[:space:]]", panther_lines)]
    panther_lines <- panther_lines[!grepl("Analysis Type|Annotation Version|Analyzed List|Reference List|Test Type|Correction", panther_lines)]
    panther_tmp <- file.path(tempdir(), "panther_tmp.txt")
    writeLines(panther_lines, panther_tmp)

    dat <- rio::import(panther_tmp) %>%
      as.enrichdat()
  } else {
    dat <- rio::import(panther_file, col_names = F)
    omit_rows <- grepl("Analysis Type|Annotation Version|Analyzed List|Reference List|Test Type|Correction", dat[, 1], ignore.case = T)
    dat <- dat[!omit_rows, ]
    colnames(dat) <- dat[1, ]
    dat <- dat[-1, ]
    dat <- as.enrichdat(dat)
  }

  bgsize <- colnames(dat)[grepl("\\([0-9]{,4}\\)", colnames(dat))] %>%
    stringr::str_remove(., ".*\\(") %>%
    stringr::str_remove(., "\\)") %>%
    as.numeric() %>%
    max()
  check_reflist <- which(grepl("reflist", tolower(names(dat))))
  dat[, check_reflist] <- as.numeric(dat[, check_reflist])

  dat <- dat %>%
    dplyr::mutate(.[check_reflist] / bgsize) %>%
    dplyr::rename(BgRatio = check_reflist) %>%
    dplyr::mutate(RichFactor = Count / bgsize) %>%
    dplyr::relocate(GeneRatio, .before = BgRatio) %>%
    dplyr::mutate(ID = Description %>%
      stringr::str_remove(".*\\(") %>%
      stringr::str_remove("\\)")) %>%
    dplyr::relocate(ID, .before = dplyr::everything()) %>%
    dplyr::mutate(Description = Description %>%
      stringr::str_remove("\\(.*\\)")) %>%
    dplyr::mutate(
      FoldEnrich = GeneRatio / BgRatio,
      pvalue = as.numeric(pvalue),
      qvalue = as.numeric(qvalue)
    )

  return(dat)
}


getExtension <- function(file) {
  # only keep the last extension (e.g. "result.v2.txt" -> "txt")
  ex <- tolower(tools::file_ext(file))
  return(ex)
}
