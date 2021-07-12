##' transform gene id among symbpl, entrezid and ensembl
##' @param id character of gene ids.
##' @param trans_to transform to which type, one of symbol, entrezid and ensembl.
##' @param org organism name from `biocOrg_name`.
##' @param return_dat logical, default is `TRUE`, return a dataframe of transformed ids.
##' @return a character or dataframe of transformed ids.
##' @importFrom dplyr %>% filter pull select distinct arrange all_of
##' @importFrom AnnotationDbi toTable
##' @export
##' @examples
##' \donttest{
##' msigdb <- getMsigdb(org='human', category='C5',subcategory='GO:CC')
##' }
transId <- function(id, trans_to, org, return_dat = FALSE){

  #--- args ---#
  options(warn = -1)
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)
  org.bk = org
  org = mapBiocOrg(tolower(org))
  keytype = .gentype(id, org)

  if (tolower(trans_to) == "entrez" | tolower(trans_to) == "etrez") trans_to = 'entrezid'
  if (tolower(trans_to) == "ensemblid" | tolower(trans_to) == "ens" |
      tolower(trans_to) == "ensemb" ) trans_to = 'ensembl'
  if (tolower(trans_to) == "sym" ) trans_to = 'symbol'

  if(!tolower(trans_to) %in% c('symbol','entrezid','ensembl')){
    stop('\nChoose trans_to argument from: \nsymbol | entrezid | ensembl ')
  }

  if(tolower(keytype) == trans_to) stop('The argument "trans_to" cannot be the same type with gene id!')

  #--- codes ---#
  # ispecies of human, mouse, rat could use ensembl gtf result
  if(org %in% c('Hs', 'Mm' ,'Rn')){
    res = trans_hsmmrn(id, trans_to, org, return_dat)
  }else{
    # others species could use the orgdb result
    res = trans_others(id, trans_to, org, return_dat)
  }

  return(res)
}


trans_hsmmrn <- function(id, trans_to, org, return_dat = FALSE){
  #--- args ---#
  options(warn = -1)
  if (org == "Hs" ) org = 'human'
  if (org == "Mm" ) org = 'mouse'
  if (org == "Rn" ) org = 'rat'

  #--- codes ---#
  gtf = eval(parse(text = paste0(org,'_gtf')))
  newdat = gtf %>%
    dplyr::filter(eval(parse(text = from)) %in% id  )
    dplyr::arrange(match(.[,from], id)) %>%
    dplyr::relocate(all_of(from),all_of(trans_to) , .before = everything())
  new_id = newdat %>% dplyr::pull(eval(parse(text = tolower(trans_to))))

  percen = paste(round(100*length(as.character(na.omit(new_id)))/length(unique(id)), 2), "%", sep="")
  if(length(new_id) > length(unique(id))){
    message(percen,' genes are mapped from ',from, ' to ', trans_to,'\n',
            'maybe one ', from, ' gene mapps many ', trans_to)
  } else {
    message(percen,' genes are mapped from ',from, ' to ', trans_to)
  }

  if(return_dat){
    res1 = newdat
  }else{
    res1 = new_id
  }

  invisible(res1)
}

trans_others <- function(id, trans_to, org, return_dat = FALSE){
  #--- args ---#
  options(warn = -1)
  .load_orgdb(org)
  keytype = .gentype(id, org)

  symbol_dat = AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egSYMBOL"))))
  colnames(symbol_dat) = c('entrezid','symbol')
  ensem_dat = AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egENSEMBL"))))
  colnames(ensem_dat) = c('entrezid','ensembl')

  from = tolower(keytype)
  if(from == 'ensembl'){
    if(!any(id %in% ensem_dat[,from])){
      stop('Input gene ids do not belong to: ',org.bk)
    }
  }else{
    if(!any(id %in% symbol_dat[,from])){
      stop('Input gene ids do not belong to: ',org.bk)
    }
  }

  #--- codes ---#
  merge_dat = merge(symbol_dat, ensem_dat ,by = 'entrezid', all  = T)
  newdat <- merge_dat %>%
    dplyr::select(c(all_of(from),all_of(trans_to))) %>%
    dplyr::filter(.[,1] %in% id) %>%
    dplyr::distinct() %>%
    dplyr::arrange(match(.[,1], id))

  new_id =  newdat %>% dplyr::pull(2) %>% as.character()

  percen = paste(round(100*length(as.character(na.omit(new_id)))/length(unique(id)), 2), "%", sep="")
  if(length(new_id) > length(unique(id))){
    message(percen,' genes are mapped from ',from, ' to ', trans_to,'\n',
            'maybe one ', from, ' gene mapps many ', trans_to)
  } else {
    message(percen,' genes are mapped from ',from, ' to ', trans_to)
  }

  if(return_dat){
    res2 = newdat
  }else{
    res2 = new_id
  }

  invisible(res2)
}
