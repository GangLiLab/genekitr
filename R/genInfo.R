#' Get gene related information
#'
#' @param id Gene id.
#' @param org Species name from `biocOrg_name`.
#' @importFrom stringr str_detect
#' @importFrom magrittr set_rownames
#' @importFrom dplyr  %>% filter relocate select
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' x = genInfo(id= c("Cyp2c23","Fhit","Gal3st2b","Trp53","Tp53"), org = 'mm')
#' head(x)

genInfo <- function(id,
                    org,
                    ...) {
  #--- args ---#
  options(warn = -1)
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)

  org.bk = org
  org = mapBiocOrg(tolower(org))
  # .load_orgdb(org)
  keytype = .gentype(id, org) %>% tolower()

  if(any(duplicated(id))){
    message('Input gene id has duplicates, keep it unique...')
    id = id[!duplicated(id)]
  }

  #--- code ---#
  all = biocAnno(org) %>% dplyr::relocate(all_of(keytype),.before = everything())
  rm(list = paste0(org,'_anno'), envir = .GlobalEnv)

  tmp1 = data.frame(input_id = id)
  tmp2 <- all %>% dplyr::filter(eval(parse(text = keytype)) %in% id)

  if(keytype != 'symbol'){
    gene_info <- merge(tmp1,tmp2,by.x = 'input_id', by.y = keytype, all.x=T) %>%
      magrittr::set_rownames(.$input_id) %>%
      dplyr::select(-input_id)
    gene_info[ gene_info == "NA" ] <- NA
  }else{
    gene_info <- merge(tmp1,tmp2,by.x = 'input_id', by.y = keytype, all.x=T) %>%
      dplyr::mutate(symbol = dplyr::case_when(input_id%in%all$symbol ~ input_id)) %>%
      dplyr::relocate(symbol, .after = input_id) %>%
      split(., .$input_id) %>%
      lapply(., function(x){
        x = x[order(as.numeric(x$entrezid)),]
        apply(x, 2, function(y){
          if(!any(duplicated(y))) {paste0(y, collapse = "; ") }else{ y[1] }
        })
      }) %>%
      do.call(rbind,.) %>% as.data.frame() %>%
      apply(., 2, function(x) gsub('^NA; ','',x) %>%  gsub('; NA$','',.) %>%
              gsub('^; ','',.) %>% gsub('; $','',.))%>%
      as.data.frame() %>%
      mutate_all(., list(~na_if(.,""))) %>%
      dplyr::select(-input_id)

    gene_info[ gene_info == "NA" ] <- NA
    # check if symbol in alias
    all_alias = data.frame(all_alias = paste(all$ncbi_alias, all$ensembl_alias,sep = '; '))
    check_row = which(is.na(gene_info$symbol))
    for(i in check_row){
      # i = check_row[1]
      alias_row = which(stringr::str_detect(all_alias[,1], paste0('\\b',rownames(gene_info)[i],'\\b')))
      if(length(alias_row) != 0){
        # not match symbol but match alias
        gene_info[i,] = all[alias_row,]
      }
    }
  }

  gene_info = gene_info[id,]
  return(gene_info)
}

