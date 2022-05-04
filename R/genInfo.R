#' Get gene related information
#'
#' @param id Gene id (symbol, ensembl or entrez id) or uniprot id. If this argument is NULL, return all gene info.
#' @param org Latin organism shortname from `ensOrg_name_data`. Default is human.
#' @importFrom stringr str_detect
#' @importFrom dplyr %>% filter relocate select mutate mutate_all na_if
#' @importFrom tidyr unnest
#' @importFrom rlang .data
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' \donttest{
#' # example1: input list with fake id and one-to-many mapping id
#' x = genInfo(id = c(
#'   "MCM10", "CDC20", "S100A9", "MMP1", "BCC7",
#'   "FAKEID", "TP53", "HBD", "NUDT10"))
#'
#' # example2: statistics of human gene biotypes
#' genInfo(org = 'hs') %>% {table(.$gene_biotype)}
#'
#' }
genInfo <- function(id = NULL,
                    org = 'hs') {
  #--- args ---#
  org <- mapEnsOrg(org)

  #--- code ---#
  if(is.null(id)){
    gene_info <- ensAnno(org)
  }else{
    all <- ensAnno(org)
    keytype <- gentype(id = id,data = all, org = org) %>% tolower()

    ## get ensembl/entrez/uniprot/symbol order
    order_dat = getOrder(org, keytype) %>%
      filter(eval(parse(text = keytype))%in%id) %>%
      mutate(!! keytype :=  factor(.[[keytype]], levels = unique(id))) %>%
      arrange(.[[keytype]])

    ## extract info from all data frame
    gene_info  =all[order_dat$rnum,] %>%
      mutate(input_id=order_dat[[keytype]]) %>%
      relocate('input_id',.before = everything()) %>%
      merge(.,as.data.frame(id),
            by.x = 'input_id', by.y = 'id',
            all.y = T) %>%
      filter(!duplicated(cbind(input_id, symbol,chr,start,end))) %>%
      arrange(input_id)

    ## check one-to-many match
    tomany_id <- names(table(gene_info$input_id))[table(gene_info$input_id) > 1]
    tomany_id <- tomany_id[!tomany_id %in% id[duplicated(id)]]
    if (length(tomany_id) > 0 & length(tomany_id) < 3 ) {
      message(paste0(
        'Some ID occurs one-to-many match, like "', paste0(tomany_id, collapse = ", "), '"\n'
      ))
    } else if (length(tomany_id) > 3 ) {
      message(paste0(
        'Some ID occurs one-to-many match, like "', paste0(tomany_id[1:3], collapse = ", "), '"...\n'
      ))
    }
  }

  #   # if only keep one, choose the one with minimum NA
  #   if (!keepAll & length(tomany_id) != 0) {
  #     sub <- gene_info %>% dplyr::filter(input_id %in% tomany_id)
  #     other <- gene_info %>% dplyr::filter(!input_id %in% tomany_id)
  #
  #     uniq_order <- sapply(tomany_id, function(x) {
  #       check <- which(sub$input_id == x)
  #       n_na <- apply(sub[check, ], 1, function(x) sum(is.na(x)))
  #       if (min(n_na) == max(n_na) & keytype != "entrezid") {
  #         check[order(as.numeric(sub$entrezid[check])) == 1]
  #       } else if (min(n_na) == max(n_na) & keytype == "entrezid") {
  #         check[order(as.numeric(sub$input_id[check])) == 1]
  #       } else if (min(n_na) != max(n_na)) {
  #         check[which.min(n_na)]
  #       }
  #     })
  #     gene_info <- rbind(other, sub[uniq_order, ])
  #     gene_info <- gene_info[match(id, gene_info$input_id), ]
  #   } else {
  #     id <- factor(id, ordered = T, levels = unique(id))
  #     gene_info$input_id <- factor(gene_info$input_id, ordered = T, levels = unique(id))
  #     gene_info <- gene_info[order(gene_info$input_id), ]
  #   }

  # reorder column
  if(!is.null(id)){
    if(keytype %in% c('ensembl','entrezid','uniprot')){
      gene_info <- gene_info %>% dplyr::select(!keytype)
    }else{
      gene_info <- gene_info %>% dplyr::relocate(symbol,.after = input_id)
    }
  }


  return(gene_info)
}

utils::globalVariables(c(":=","symbol","uniprot","input_id","symbol","chr","start","end"))
