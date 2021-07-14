##' Get gene info using org.db
##'
##' @param id gene id.
##' @param org species "mm" or "hs".
##' @return a dataframe of gene info.
##' @importFrom stringr str_to_title
##' @importFrom AnnotationDbi toTable
##' @export
##' @examples
##' \dontrun{
##' x = genInfo(id= c("Cyp2c23","Fhit","Gal3st2b","Gbp4"), org = 'mm')
##' }

genInfo <- function(id,
                    org,
                    ...) {
  #--- args ---#
  options(warn = -1)
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)

  org = mapBiocOrg(tolower(org))
  .load_orgdb(org)
  keytype = .gentype(id, org) %>% tolower()

  if(any(duplicated(id))){
    message('Input gene id has duplicates, keep it unique...')
    id = id[!duplicated(id)]
  }

  #--- code ---#
  all = biocAnno(org)
  gene_info <- all %>%
    dplyr::filter(eval(parse(text = keytype)) %in% id) %>%
    dplyr::arrange(match(eval(parse(text = keytype)), id)) %>%
    dplyr::relocate(all_of(keytype),.before = everything())

  # if some ids are spelled wrong or the name is an alias
  all_alias = all %>% dplyr::pull(gene_alias) %>% stringr::str_split(.,';') %>%
    unlist() %>% unique()

  # first to decide: if input name is alias
  if( any( !id %in% all[,keytype] & id %in% all_alias) ){
    check = which( !id %in% all[,keytype] & id %in% all_alias)
    i = 1
    while (i <= length(check) ) {
      gene_info = gene_info %>%
        tibble::add_row(all[mapply(grepl, paste0('\\b',id[check],'\\b'), all$gene_alias),],
                        .before = check[i])
      gene_info[check,1] = id[check]
      i = i+1
    }

    # then to decide: if input id is wrong
    if (any( !id %in% all[,keytype] & !id %in% all_alias) ){
      check = which(  !id %in% all[,keytype] & !id %in% all_alias )
      tmp = data.frame(rep(NA,ncol(gene_info))) %>% t() %>% as.data.frame() %>%
        stats::setNames(.,colnames(gene_info))
      rownames(tmp) = '1'

      i = 1
      while (i <= length(check) ) {
        gene_info = gene_info %>%
          tibble::add_row(tmp,
                          .before = check[i])
        gene_info[check,1] = id[check]

        i = i+1
      }
    }
  }

  return(gene_info)
}
