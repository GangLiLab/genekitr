#' Transform gene id among symbol, entrezid,  ensembl and uniprot
#' @param id Gene ids.
#' @param trans_to Transform to which type, one of "symbol", "entrezid",
#'   "ensembl" and "uniprot."
#' @param org Organism name from `biocOrg_name()`, both full name and short name are fine.
#' @param simple Logical to keep only one matched ID, default is FALSE.
#' @importFrom dplyr %>% filter pull select distinct arrange all_of
#' @importFrom AnnotationDbi toTable
#' @importFrom tibble add_row
#'
#' @return A character of transformed ids.
#' @export
#'
#' @examples
#' transId(
#'   id = c("Cyp2c23", "Fhit", "Gal3st2b", "Trp53", "Tp53"),
#'   trans_to = "ensembl", org = "mouse", simple = TRUE
#' )
#' # input id contains duplicates,fake id and one-to-many match id
#' transId(
#' id = c("MMD2", "HBD", "TP53", "RNR1", "TEC", "BCC7", "FAKEID", "TP53"),
#' trans_to = "entrez", org = "hg", simple = FALSE
#' )
transId <- function(id, trans_to, org, simple = TRUE) {

  #--- args ---#
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)
  org <- mapBiocOrg(tolower(org))
  keytype <- .gentype(id, org)
  from <- tolower(keytype)

  if (grepl(tolower(trans_to), "entrezid")) trans_to <- "entrezid"
  if (grepl(tolower(trans_to), "ensemblid")) trans_to <- "ensembl"
  if (grepl(tolower(trans_to), "symbolid")) trans_to <- "symbol"
  if (grepl(tolower(trans_to), "uniprot")) trans_to <- "uniprot"

  if (!tolower(trans_to) %in% c("symbol", "entrezid", "ensembl", "uniprot")) {
    stop("\nChoose 'trans_to' argument from: \nsymbol | entrezid | ensembl | uniprot")
  }

  #--- codes ---#
  new_id <- genInfo(id, org, simple) %>% dplyr::pull(trans_to)

  if (simple) {
    percent <- paste(round(100 * length(as.character(na.omit(new_id))) / length(id), 2), "%", sep = "")
    message(percent, " genes are mapped from ", from, " to ", trans_to)
    if (any(is.na(new_id))) {
      message(paste0(
        "\nSome ID could not match ", trans_to, ", return NA",
        '...\nMaybe use "na.omit()" for downstream analysis'
      ))
    }
  }

  return(new_id)
}
