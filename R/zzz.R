##' @importFrom utils packageDescription

.onLoad <- function(...) {
  invisible(suppressPackageStartupMessages(
    sapply(c(
      "stringi", "stringr",
      "ggplot2", "dplyr", "devtools"
    ),
    requireNamespace,
    quietly = TRUE
    )
  ))
}

.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  welcome_message <- paste0("Welcome to use ", pkgname, "! (",
                "Vignette: https://www.genekitr.fun)\n")
  citation <- paste0("Citation for ", pkgname, ":\n",
                     "Liu, Y., Li, G. ",
                     "Empowering biologists to decode omics data: the Genekitr R package and web server. ",
                     "BMC Bioinformatics 24, 214 (2023). https://doi.org/10.1186/s12859-023-05342-9")

  packageStartupMessage(paste0(welcome_message,citation))
}
