#' Transform gene id among symbol, entrezid,  ensembl and uniprot
#' @param id Gene ids.
#' @param trans_to Transform to which type, one of "symbol", "entrezid",
#'   "ensembl" and "uniprot."
#' @param org Short latin name from `ensOrg_name_data`.
#' @param unique If keep only one unique mapped ID when one-to-many gene occurs, default is FALSE.
#' @importFrom dplyr %>% filter pull select distinct arrange all_of
#' @importFrom tibble add_row
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
#' @return A character of transformed ids.
#' @export
#'
#' @examples
#' \dontrun{
#' transId(
#'   id = c("Cyp2c23", "Fhit", "Gal3st2b", "Trp53", "Tp53"),
#'   trans_to = "ensembl", org = "mouse", unique = TRUE
#' )
#' # input id contains fake id and one-to-many match id
#' transId(
#'   id = c("MMD2", "HBD", "RNR1", "TEC", "BCC7", "FAKEID", "TP53"),
#'   trans_to = "entrez", org = "hg", unique = FALSE
#' )
#' }
#'
transId <- function(id, trans_to, org, unique = TRUE) {

  #--- args ---#
  org <- mapEnsOrg(tolower(org))
  keytype <- gentype(id, org)
  from <- tolower(keytype)

  if (grepl(tolower(trans_to), "entrezid")) trans_to <- "entrezid"
  if (grepl(tolower(trans_to), "ensemblid")) trans_to <- "ensembl"
  if (grepl(tolower(trans_to), "symbolid")) trans_to <- "symbol"
  if (grepl(tolower(trans_to), "uniprot")) trans_to <- "uniprot"

  if (!tolower(trans_to) %in% c("symbol", "entrezid", "ensembl", "uniprot")) {
    stop("\nChoose 'trans_to' argument from: \nsymbol | entrezid | ensembl | uniprot")
  }

  #--- codes ---#
  if (unique) {
    new_id <- genInfo(id, org, unique) %>% dplyr::pull(trans_to)
    n_new <- na.omit(new_id) %>%
      as.character() %>%
      length()
  } else {
    new_id <- genInfo(id, org, unique) %>% dplyr::select(input_id, all_of(trans_to))
    n_new <- na.omit(new_id) %>%
      .[1] %>%
      unique() %>%
      nrow()
  }

  percent <- paste(round(100 * n_new / length(id), 2), "%", sep = "")
  message("\n", percent, " genes are mapped from ", from, " to ", trans_to)
  if (n_new != length(id)) {
    message(paste0(
      "Non-matched ID are marked as NA",
      '...\nMaybe use "na.omit()" for downstream analysis'
    ))
  }

  return(new_id)
}
