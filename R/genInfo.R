#' Get gene related information
#'
#' @param id Gene id (symbol, ensembl or entrez id) or uniprot id.
#' @param org Species name from `biocOrg_name`, both full name and short name are fine.
#' @param unique Logical to keep only one unique mapped ID, default is FALSE.
#' @importFrom stringr str_detect
#' @importFrom dplyr %>% filter relocate select mutate mutate_all na_if slice add_row
#' @importFrom tidyr unnest
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' x <- genInfo(id = c(
#'   "MCM10", "CDC20", "S100A9", "MMP1", "BCC7",
#'   "FAKEID", "TP53", "HBD", "NUDT10"),
#'   org = "hg", unique = FALSE)
#' head(x)
genInfo <- function(id,
                    org,
                    unique = FALSE) {
  #--- args ---#
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)

  org <- mapBiocOrg(tolower(org))
  keytype <- .gentype(id, org) %>% tolower()

  #--- code ---#
  all <- biocAnno(org) %>%
    dplyr::relocate(all_of(keytype), .before = everything()) %>%
    dplyr::mutate(!!keytype := strsplit(.[, 1], "; ")) %>%
    tidyr::unnest(cols = all_of(keytype)) %>%
    as.data.frame()
  rm(list = paste0(org, "_anno"), envir = .GlobalEnv)

  tmp1 <- data.frame(input_id = id)
  tmp2 <- all %>% dplyr::filter(eval(parse(text = keytype)) %in% id)

  ## keep each id even has no info
  # only symbol id needs to consider alias
  if (keytype != "symbol") {
    gene_info <- merge(tmp1, tmp2, by.x = "input_id", by.y = keytype, all.x = T)
  } else {
    tmp2 <- tmp2 %>% dplyr::filter_at(dplyr::vars(-symbol, -uniprot), dplyr::all_vars(!is.na(.)))
    gene_info <- merge(tmp1, tmp2, by.x = "input_id", by.y = keytype, all.x = T) %>%
      dplyr::arrange(id) %>%
      dplyr::mutate(symbol = dplyr::case_when(input_id %in% tmp2$symbol ~ input_id)) %>%
      dplyr::relocate(symbol, .after = input_id)

    # check if symbol in alias (only check input ids without matched)
    all_alias <- data.frame(all_alias = paste(all$ncbi_alias,all$ensembl_alias, sep = "; "))
    check_row <- which(is.na(gene_info$symbol))
    for (i in check_row) {
      # i  = check_row[1]
      alias_row <- which(stringr::str_detect(all_alias[, 1], paste0("\\b", gene_info[i, 1], "\\b")))
      if (length(alias_row) != 0) {
        # not match symbol but match alias
        gene_info = dplyr::add_row(gene_info,all[alias_row,],.after = i)
        gene_info$input_id[(i+1):(i+length(alias_row))] = gene_info$input_id[i]
        gene_info = gene_info %>% dplyr::slice(-i)
      }
    }
  }

  # one-to-many match
  check_n <- table(gene_info$input_id)
  tomany_id <- names(check_n)[check_n > 1]
  tomany_id <- tomany_id[!tomany_id %in% id[duplicated(id)]]
  if (length(tomany_id) > 0 & length(tomany_id) < 3 & !unique) {
    message(paste0(
      'Some ID occurs one-to-many match, like "', tomany_id, '"\n',
      'If you want to get one-to-one match, please set "unique=TRUE"'
    ))
  } else if (length(tomany_id) > 3 & !unique) {
    message(paste0(
      'Some ID occurs one-to-many match, like "', paste0(tomany_id[1:3], collapse = ", "), '"...\n',
      'If you want to get one-to-one match, please set "unique=TRUE"'
    ))
  }

  # if keep unique, choose row with minimum NA;
  # if NA number is identical, then choose the smallest entrezid
  if (unique) {
    sub = gene_info %>% dplyr::filter(input_id %in% tomany_id)
    other = gene_info %>% dplyr::filter(!input_id %in% tomany_id)

    uniq_order = sapply(tomany_id, function(x){
      check = which(sub$input_id == x)
      n_na = apply(sub[check,], 1, function(x)sum(is.na(x)))
      if(min(n_na) == max(n_na) & keytype != 'entrezid'){
        check[order(as.numeric(sub$entrezid[check]))==1]
      }else if (min(n_na) == max(n_na) & keytype == 'entrezid'){
        check[order(as.numeric(sub$input_id[check]))==1]
      }else if (min(n_na) != max(n_na) ){
        check[which.min(n_na)]
      }
    })
    gene_info = rbind(other, sub[uniq_order,])
    gene_info = gene_info[match(id,gene_info$input_id),]

  } else {
    id = factor(id,ordered = T,levels = unique(id))
    gene_info$input_id = factor(gene_info$input_id,ordered = T,levels = unique(id))
    gene_info = gene_info[order(gene_info$input_id),]
  }

  return(gene_info)
}
