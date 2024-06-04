#' Import 'shinyGO' web result
#'
#' @param shinygo_file ShinyGO result file.
#'
#' @importFrom dplyr select mutate relocate
#' @importFrom stats setNames
#' @importFrom rlang .data
#' @return  `data.frame`
#' @export

importShinygo <- function(shinygo_file) {
  if (!requireNamespace("rio", quietly = TRUE)) {
    stop("Package \"rio\" needed for this function to work.
         Please install it by install.packages('rio')",call. = FALSE)
  }

  shinygo <- rio::import(shinygo_file)

  dat <- shinygo %>%
    dplyr::select(c('Ontology Code','Pathway','nGenes','Pathway Genes','Enrichment FDR','Genes','Fold Enrichment')) %>%
    stats::setNames(c('ID','Description','Count','Pathway_count','qvalue','geneID','FoldEnrich')) %>%
    dplyr::mutate(RichFactor = Count/Pathway_count) %>%
    dplyr::relocate(RichFactor,.after = FoldEnrich) %>%
    dplyr::relocate(Count,.after = geneID) %>%
    dplyr::select(-Pathway_count) %>%
    dplyr::mutate(geneID = gsub('\\s+','/',geneID) ) %>%
    dplyr::mutate(FoldEnrich = as.numeric(FoldEnrich)) %>%
    dplyr::mutate(qvalue = as.numeric(qvalue))

  return(dat)
}

utils::globalVariables(c(
  "Pathway_count", "RichFactor"
))
