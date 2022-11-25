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
  msg <- paste0("Welcome to use ", pkgname, "!\n",
                "Vignette: https://www.genekitr.fun")

  packageStartupMessage(paste0(msg))
}
