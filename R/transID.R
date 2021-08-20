#' transform gene id among symbol, entrezid,  ensembl and uniprot
#' @param id character of gene ids.
#' @param trans_to transform to which type, one of symbol, entrezid, ensembl and
#'   uniprot.
#' @param org organism name from `biocOrg_name`.
#' @param keep_unique logical to keep output id is unique, default is `TRUE`.
#' @return a character of transformed ids.
#' @importFrom dplyr %>% filter pull select distinct arrange all_of
#' @importFrom AnnotationDbi toTable
#' @importFrom tibble add_row
#' @export
#' @examples
#' \donttest{
#' transId(id= c("Cyp2c23","Fhit","Gal3st2b","Trp53","Tp53"),
#'   trans_to = 'ensembl', org = 'mouse')
#' }
transId <- function(id, trans_to, org, keep_unique = TRUE){

  #--- args ---#
  options(warn = -1)
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)
  org = mapBiocOrg(tolower(org))
  keytype = .gentype(id, org)
  from = tolower(keytype)

  if (grepl(tolower(trans_to),'entrezid')) trans_to = 'entrezid'
  if (grepl(tolower(trans_to),'ensemblid')) trans_to = 'ensembl'
  if (grepl(tolower(trans_to),'symbolid')) trans_to = 'symbol'
  if (grepl(tolower(trans_to),'uniprot')) trans_to = 'uniprot'

  if(!tolower(trans_to) %in% c('symbol','entrezid','ensembl','uniprot')){
    stop('\nChoose trans_to argument from: \nsymbol | entrezid | ensembl | uniprot')
  }

  # if(tolower(keytype) == trans_to) stop('The argument "trans_to" cannot be the same type with gene id!')

  #--- codes ---#
  new_id = genInfo(id,org) %>%  dplyr::pull(trans_to)

  if(keep_unique){
    # remove other ids behind ";"
    res = stringr::str_remove_all(new_id,';.*')
  }else{
    res = new_id
  }
  percent = paste(round(100*length(as.character(na.omit(new_id)))/length(unique(id)), 2), "%", sep="")
  message(percent,' genes are mapped from ',from, ' to ', trans_to)

  return(res)
}

