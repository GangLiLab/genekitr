#' Get gene related information
#'
#' @param id Gene id (symbol, ensembl or entrez id) or uniprot id.
#' @param org Species name from `biocOrg_name`, both full name and short name are fine.
#' @param simple Logical to keep only one matched ID, default is FALSE.
#' @importFrom stringr str_detect
#' @importFrom dplyr %>% filter relocate select mutate mutate_all na_if
#' @importFrom tidyr unnest
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' x <- genInfo(id = c(
#'   "MCM10", "CDC20", "S100A9", "MMP1", "BCC7",
#'   "FAKEID", "TP53", "HBD", "TP53", "NUDT10"
#' ), org = "hg", simple = FALSE)
#' head(x)
genInfo <- function(id,
                    org,
                    simple = FALSE) {
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
    gene_info <- merge(tmp1, tmp2, by.x = "input_id", by.y = keytype, all.x = T) %>%
      dplyr::mutate(symbol = dplyr::case_when(input_id %in% all$symbol ~ input_id)) %>%
      dplyr::relocate(symbol, .after = input_id)

    # check if symbol in alias (only check input ids without matched)
    all_alias <- data.frame(all_alias = paste(all$ncbi_alias, all$ensembl_alias, sep = "; "))
    check_row <- which(is.na(gene_info$symbol))
    for (i in check_row) {
      alias_row <- which(stringr::str_detect(all_alias[, 1], paste0("\\b", gene_info[i, 1], "\\b")))
      if (length(alias_row) != 0) {
        # not match symbol but match alias
        gene_info[i, 2:ncol(gene_info)] <- all[alias_row, ]
      }
    }
  }

  # one-to-many match
  tomany_id <- names(table(gene_info$input_id))[table(gene_info$input_id) > 1]
  tomany_id <- tomany_id[!tomany_id %in% id[duplicated(id)]]
  if (length(tomany_id) > 0 & length(tomany_id) < 3 & !simple) {
    message(paste0(
      'Some ID occurs one-to-many match, like "', tomany_id, '"\n',
      'If you want to get one-to-one match, please set "simple=TRUE"'
    ))
  } else if (length(tomany_id) > 3 & !simple) {
    message(paste0(
      'Some ID occurs one-to-many match, like "', paste0(tomany_id[1:3], collapse = ", "), '"...\n',
      'If you want to get one-to-one match, please set "simple=TRUE"'
    ))
  }

  # if keep simple, choose row with minimum NA
  if (simple) {
    gene_info <- gene_info[unlist(lapply(id, function(x) {
      if (x %in% tomany_id) {
        check = which(gene_info$input_id == x)
        apply(gene_info[check,], 1, function(x)sum(is.na(x))) %>%
          which.min() %>% names() %>% as.numeric()
      } else {
        which(gene_info$input_id == x)[1]
      }})), ]
  } else {
    gene_info <- gene_info[unlist(lapply(id, function(x) {
      if (x %in% tomany_id) {
        which(gene_info$input_id == x)
      } else {
        which(gene_info$input_id == x)[1]
      }})), ]
  }

  return(gene_info)
}
