##' @importFrom utils packageDescription

.onLoad <- function(...){
  invisible(suppressPackageStartupMessages(
    sapply(c("stringi", "stringr",
             "ggplot2", "dplyr","devtools"),
           requireNamespace, quietly = TRUE)))

}

.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("Welcome to use ",pkgname, "!\n")

  # citation <- paste0("If you use ", pkgname, " in published research, please acknowledgements:\n",
  #                    "We thank Dr.Jianming Zeng(University of Macau), and all the members of his bioinformatics team, biotrainee, for generously sharing their experience and codes.")

  packageStartupMessage(paste0(msg))
}
