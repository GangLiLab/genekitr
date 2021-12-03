#' Get gene related information
#'
#' @param id Gene id (symbol, ensembl or entrez id) or uniprot id. If this argument is NULL, return all gene info.
#' @param org Short latin name from `ensOrg_name_data`.
#' @param unique Logical to keep only one matched ID, default is FALSE.
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
#' # input id contains fake id and one-to-many match id
#' x <- genInfo(id = c(
#'   "MCM10", "CDC20", "S100A9", "MMP1", "BCC7",
#'   "FAKEID", "TP53", "HBD", "NUDT10"
#' ), org = "hg", unique = TRUE)
#' head(x)
#'
#' }
genInfo <- function(id = NULL,
                    org,
                    unique = FALSE) {
  #--- args ---#
  org <- mapEnsOrg(org)
  if(is.null(id)){
    gene_info <- ensAnno(org)
    rm(list = paste0(org, "_anno"), envir = .GlobalEnv)
  }else{
    keytype <- gentype(id, org) %>% tolower()

  #--- code ---#
    all <- ensAnno(org) %>%
      dplyr::relocate(all_of(keytype), .before = everything()) %>%
      dplyr::mutate(!!keytype := strsplit(.[, 1], "; ")) %>%
      tidyr::unnest(cols = all_of(keytype)) %>%
      as.data.frame()
    rm(list = paste0(org, "_anno"), envir = .GlobalEnv)

    tmp1 <- data.frame(input_id = id)
    tmp2 <- all %>% dplyr::filter(eval(parse(text = keytype)) %in% id)

    if(any(c("symbol", "uniprot") %in% colnames(all))){
      tmp3 <- tmp2 %>%
        dplyr::select(-c("symbol", "uniprot")) %>%
        apply(., 1, is.na) %>%
        as.data.frame()
    }

    if(keytype != "symbol"){
      gene_info <- merge(tmp1, tmp2, by.x = "input_id", by.y = keytype, all.x = T)
    }else if( any(apply(tmp3, 2, sum) == nrow(tmp3)) ){
      tmp2 <- tmp2[-which(apply(tmp3, 2, sum) == nrow(tmp3)),]
      gene_info <- merge(tmp1, tmp2, by.x = "input_id", by.y = keytype, all.x = T) %>%
        dplyr::left_join(data.frame(input_id=id),.,by="input_id") %>%
        dplyr::mutate(symbol = dplyr::case_when(input_id %in% tmp2$symbol ~ input_id)) %>%
        dplyr::relocate(symbol, .after = input_id)
      all_alias <- data.frame(all_alias = paste(all$ncbi_alias, all$ensembl_alias, sep = "; "))
      check_id <- gene_info$input_i[is.na(gene_info$symbol)]
      for (i in check_id) {
        # i =check_id[3]
        alias_row <- which(stringr::str_detect(all_alias[, 1], paste0("\\b", i, "\\b")))
        loc = match(i,gene_info$input_id)
        if ( length(alias_row) != 0 ) {
          gene_info <- dplyr::add_row(gene_info, all[alias_row, ], .after = loc)
          gene_info$input_id[(loc + 1):(loc + length(alias_row))] <- i
          gene_info <- gene_info %>% dplyr::slice(-loc)
        }
      }
    }else{
      gene_info <- merge(tmp1, tmp2, by.x = "input_id", by.y = keytype, all.x = T) %>%
        dplyr::left_join(data.frame(input_id=id),.,by="input_id") %>%
        dplyr::mutate(symbol = dplyr::case_when(input_id %in% tmp2$symbol ~ input_id)) %>%
        dplyr::relocate(symbol, .after = input_id)
      all_alias <- data.frame(all_alias = paste(all$ncbi_alias, all$ensembl_alias, sep = "; "))
      check_id <- gene_info$input_i[is.na(gene_info$symbol)]
      for (i in check_id) {
        # i =check_id[3]
        alias_row <- which(stringr::str_detect(all_alias[, 1], paste0("\\b", i, "\\b")))
        loc = match(i,gene_info$input_id)
        if ( length(alias_row) != 0 ) {
          gene_info <- dplyr::add_row(gene_info, all[alias_row, ], .after = loc)
          gene_info$input_id[(loc + 1):(loc + length(alias_row))] <- i
          gene_info <- gene_info %>% dplyr::slice(-loc)
        }
      }
    }

    # one-to-many match
    tomany_id <- names(table(gene_info$input_id))[table(gene_info$input_id) > 1]
    tomany_id <- tomany_id[!tomany_id %in% id[duplicated(id)]]
    if (length(tomany_id) > 0 & length(tomany_id) < 3 & !unique) {
      message(paste0(
        'Some ID occurs one-to-many match, like "', paste0(tomany_id, collapse = ", "), '"\n',
        'If you want to get one-to-one match, please set "unique=TRUE"'
      ))
    } else if (length(tomany_id) > 3 & !unique) {
      message(paste0(
        'Some ID occurs one-to-many match, like "', paste0(tomany_id[1:3], collapse = ", "), '"...\n',
        'If you want to get one-to-one match, please set "unique=TRUE"'
      ))
    }

    # if keep unique, choose row with minimum NA
    if (unique & length(tomany_id) != 0) {
      sub <- gene_info %>% dplyr::filter(input_id %in% tomany_id)
      other <- gene_info %>% dplyr::filter(!input_id %in% tomany_id)

      uniq_order <- sapply(tomany_id, function(x) {
        check <- which(sub$input_id == x)
        n_na <- apply(sub[check, ], 1, function(x) sum(is.na(x)))
        if (min(n_na) == max(n_na) & keytype != "entrezid") {
          check[order(as.numeric(sub$entrezid[check])) == 1]
        } else if (min(n_na) == max(n_na) & keytype == "entrezid") {
          check[order(as.numeric(sub$input_id[check])) == 1]
        } else if (min(n_na) != max(n_na)) {
          check[which.min(n_na)]
        }
      })
      gene_info <- rbind(other, sub[uniq_order, ])
      gene_info <- gene_info[match(id, gene_info$input_id), ]
    } else {
      id <- factor(id, ordered = T, levels = unique(id))
      gene_info$input_id <- factor(gene_info$input_id, ordered = T, levels = unique(id))
      gene_info <- gene_info[order(gene_info$input_id), ]
    }
  }

  return(gene_info)
}

utils::globalVariables(c(":=","symbol","uniprot","input_id","symbol"))
