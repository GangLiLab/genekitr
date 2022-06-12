#' Import Panther web result
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
    utils::install.packages("rio")
  }

  # unmodified panther result (should be txt file)
  if (getExtension(panther_file) == "txt") {
    system(paste0(
      "cat ", panther_file, ' | grep "\\S" | ',
      "awk '!/Analysis Type|Annotation Version|Analyzed List|Reference List|Test Type|Correction/' > ", tempdir(), "/panther_tmp.txt"
    ))

    dat <- rio::import(paste0(tempdir(), "/panther_tmp.txt")) %>%
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
  ex <- strsplit(basename(file), split = "\\.")[[1]]
  return(ex[-1])
}
