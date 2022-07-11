#' Get gene related information
#'
#' @param id Gene id (symbol, ensembl or entrez id) or uniprot id. If this argument is NULL, return all gene info.
#' @param org Latin organism shortname from `ensOrg_name`. Default is human.
#' @param unique Logical, if one-to-many mapping occurs, only keep one record with fewest NA. Default is FALSE.
#' @param keepNA If some id has no match at all, keep it or not. Default is TRUE.
#' @importFrom dplyr filter mutate arrange relocate select filter_at vars any_vars
#' @importFrom rlang .data
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' \donttest{
#' # example1: input list with fake id and one-to-many mapping id
#' x <- genInfo(id = c(
#'   "MCM10", "CDC20", "S100A9", "MMP1", "BCC7",
#'   "FAKEID", "TP53", "HBD", "NUDT10"
#' ))
#'
#' # example2: statistics of human gene biotypes
#' genInfo(org = "hs") %>%
#'   {
#'     table(.$gene_biotype)
#'   }
#' }
genInfo <- function(id = NULL,
                    org = "hs",
                    unique = FALSE,
                    keepNA = TRUE) {
  #--- args ---#
  org <- mapEnsOrg(org)

  #--- code ---#
  if (is.null(id)) {
    gene_info <- ensAnno(org)
  } else {
    all <- ensAnno(org)

    id <- replace_greek(id)
    keytype <- gentype(id = id, data = all, org = org) %>% tolower()

    ## get ensembl/entrez/uniprot/symbol order
    order_dat <- getOrder(org, all_of(keytype)) %>%
      dplyr::filter(eval(parse(text = keytype)) %in% id) %>%
      dplyr::mutate(!!keytype := factor(.[[keytype]], levels = unique(id))) %>%
      dplyr::arrange(.[[keytype]])

    ## extract info from all data frame
    gene_info <- all[order_dat$rnum, ] %>%
      dplyr::mutate(input_id = order_dat[[keytype]]) %>%
      dplyr::relocate("input_id", .before = everything()) %>%
      merge(., as.data.frame(id),
        by.x = "input_id", by.y = "id",
        all.y = T
      ) %>%
      # dplyr::filter(!duplicated(cbind(input_id, symbol,chr,start,end))) %>%
      dplyr::arrange(input_id)

    ## check one-to-many match
    tomany_id <- names(table(gene_info$input_id))[table(gene_info$input_id) > 1]
    tomany_id <- tomany_id[!tomany_id %in% id[duplicated(id)]]
    if (length(tomany_id) > 0 & length(tomany_id) < 3) {
      message(paste0(
        'Some ID occurs one-to-many match, like "', paste0(tomany_id, collapse = ", "), '"\n'
      ))
    } else if (length(tomany_id) > 3) {
      message(paste0(
        'Some ID occurs one-to-many match, like "', paste0(tomany_id[1:3], collapse = ", "), '"...\n'
      ))
    }
  }

  if (!is.null(id)) {
    # if only keep one, choose the one with minimum NA
    if (unique & length(tomany_id) != 0) {
      sub <- gene_info %>% dplyr::filter(input_id %in% tomany_id)
      other <- gene_info %>% dplyr::filter(!input_id %in% tomany_id)

      # if has entrezid and ensembl column: with exact symbol > minimal entrezid > assembled chr
      if (all(c("entrezid", "chr") %in% colnames(gene_info))) {
        uniq_order <- sapply(tomany_id, function(x) {
          # x = 'GPR1-AS'
          # print(x)
          res <- c()
          check <- which(sub$input_id == x)

          if( keytype == 'symbol' && 'symbol' %in% colnames(gene_info) ){
            sym <- sub[check, "symbol"]
            if(any(sym%in%x)){
              res <- check[which(sym%in%x)]
              if(length(res) >1) res <- c()
            }
          }

          if(is.null(res)){
            n_ent <- as.numeric(sub[check, "entrezid"]) %>% stats::na.omit() %>% as.numeric()
            if (!max(n_ent) == min(n_ent)) {
              min_n <- which(n_ent %in% min(n_ent))
              res <- check[min_n]
              if (length(min_n) > 1) {
                # if entrez is same, then check chr
                if (any(grepl("^[0-9].*$", sub[check, "chr"]))) {
                  real_chr <- which(grepl("^[0-9].*$", sub[check, "chr"]))
                  if (length(real_chr) > 1) real_chr <- real_chr[1]
                } else {
                  real_chr <- check[1]
                }
                res <- check[real_chr]
              }
            } else {
              if (any(grepl("^[0-9].*$", sub[check, "chr"]))) {
                real_chr <- which(grepl("^[0-9].*$", sub[check, "chr"]))
                if (length(real_chr) > 1) real_chr <- real_chr[1]
                res <- check[real_chr]
              } else {
                res <- check[1]
              }
            }
          }

          return(res)
        }) %>% as.character()
      } else {
        # if no entrez or ensembl, then check minimal NA
        uniq_order <- sapply(tomany_id, function(x) {
          res <- c()
          check <- which(sub$input_id == x)


          n_na <- apply(sub[check, ], 1, function(x) sum(is.na(x)))
          if (min(n_na) == max(n_na) & keytype != "entrezid") {
            res <- check[order(as.numeric(sub$entrezid[check])) == 1]
          } else if (min(n_na) == max(n_na) & keytype == "entrezid") {
            res <- check[order(as.numeric(sub$input_id[check])) == 1]
          } else if (min(n_na) != max(n_na)) {
            res <- check[which.min(n_na)]
          }
          return(res)
        }) %>% as.character()
      }

      gene_info <- rbind(other, sub[uniq_order, ])
      gene_info <- gene_info[match(id, gene_info$input_id), ]
    } else {
      id <- factor(id, ordered = T, levels = unique(id))
      gene_info$input_id <- factor(gene_info$input_id, ordered = T, levels = unique(id))
      gene_info <- gene_info[order(gene_info$input_id), ]
    }

    # reorder column
    if (!is.null(id)) {
      if (keytype %in% c("ensembl", "entrezid", "uniprot")) {
        gene_info <- gene_info %>% dplyr::select(!all_of(keytype))
      } else {
        gene_info <- gene_info %>% dplyr::relocate(symbol, .after = input_id)
      }
    }
  }

  if (!keepNA) {
    gene_info <- gene_info %>%
      filter_at(vars(!input_id), any_vars(!is.na(.)))
  }

  # replace back greek letter
  if(!is.null(id)){
    gene_info$input_id = replace_back(gene_info$input_id)
  }

  # convert factor to character
  gene_info[] <- lapply(gene_info, as.character)

  return(gene_info)
}

replace_greek <-function(id){
  id <- stringr::str_replace_all(id, "\\u03b1", 'alpha')
  id <- stringr::str_replace_all(id, "\\u03b2", 'beta')
  id <- stringr::str_replace_all(id, "\\u03b3", 'gamma')
  id <- stringr::str_replace_all(id, "\\u03b4", 'delta')
  id <- stringr::str_replace_all(id, "\\u03b5", 'epsilon')
  id <- stringr::str_replace_all(id, "\\u03bb", 'lambda')
  id <- stringr::str_replace_all(id, "\\u03ba", 'kappa')
  id <- stringr::str_replace_all(id, "\\u03c3", 'sigma')
}

replace_back <-function(id){
  id <- stringr::str_replace_all(id, 'alpha',  "\\u03b1")
  id <- stringr::str_replace_all(id, 'beta',   "\\u03b2")
  id <- stringr::str_replace_all(id, 'gamma',  "\\u03b3" )
  id <- stringr::str_replace_all(id, 'delta',  "\\u03b4" )
  id <- stringr::str_replace_all(id, 'epsilon', "\\u03b5"  )
  id <- stringr::str_replace_all(id, 'lambda', "\\u03bb"  )
  id <- stringr::str_replace_all(id, 'kappa',  "\\u03ba" )
  id <- stringr::str_replace_all(id, 'sigma',  "\\u03c3" )
}


utils::globalVariables(c(":=", "symbol", "uniprot", "input_id", "symbol", "chr", "start", "end"))
