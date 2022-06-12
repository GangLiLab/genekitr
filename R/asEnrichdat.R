#' Modify dataframe for enrichment plot
#'
#' To make sure colname contains Description, Count, FoldEnrich/GeneRatio, pvalue/qvalue/p.adjust
#'
#' @param enrich_df Enrichment analysis `data.frame` result.
#'
#' @importFrom stringr str_remove_all str_split str_remove
#' @importFrom dplyr mutate pull
#' @importFrom rlang .data
#' @return  `data.frame`
#' @export

as.enrichdat <- function(enrich_df) {

  ## get lower case colnames
  remove <- c("\\(", "\\)", " ", "-", "_")
  to_check <- stringr::str_remove_all(tolower(colnames(enrich_df)), paste(remove, collapse = "|"))

  ## find description col
  if (!any(grepl("description", to_check))) {
    check_description <- apply(enrich_df, 2, function(x) all(grepl("[A-Za-z]{3,}", x)))
    check2 <- which(check_description)
    if (any(check2 < (ncol(enrich_df) / 2))) {
      colnames(enrich_df)[check_description][check2 < 5] <- "Description"
    } else {
      stop("Not found description column!")
    }
  } else {
    colnames(enrich_df)[grepl("description", to_check)] <- "Description"
  }

  ## find count col
  # if finds genes col, calc gene num as count; else find another col as count
  if (!any(grepl("count", to_check))) {
    check_gene <- apply(enrich_df, 2, function(x) all(any(grepl("[A-Za-z]{3,}|\\/|,", x)) & !grepl("tags|list", x)))
    check2 <- which(check_gene)
    if (any(check2 > (ncol(enrich_df) / 2))) {
      colnames(enrich_df)[check_gene][check2 > (ncol(enrich_df) / 2)] <- "geneID"
      gen_num <- stringr::str_split(enrich_df$geneID, ",|\\/") %>%
        lapply(., length) %>%
        unlist()
      enrich_df <- enrich_df %>% dplyr::mutate(Count = gen_num)
    } else if (any(grepl("\\([0-9]{,4}\\)", colnames(enrich_df)))) {
      gen_num <- enrich_df[grepl("\\([0-9]{,4}\\)", colnames(enrich_df))] %>%
        dplyr::pull(2) %>%
        as.numeric()
      enrich_df <- enrich_df %>% dplyr::mutate(Count = gen_num)
    } else {
      stop('Please rename the gene count column as "Count"!')
      # head(enrich_df[1:2,])
      # message("Cannot auto-select count column...","\n","Please specify the column number which includes gene count...")
      # answer <- scan(what = "character", n =1,quiet =T)
      # message('Choose the No. ',answer,' column as gene count...')
      # enrich_df = enrich_df %>% dplyr::rename(Count = eval(parse(text = answer)))
    }
  } else {
    colnames(enrich_df)[grepl("count", to_check)] <- "Count"
  }

  ## find FoldEnrich col
  # GSEA result has no FoldEnrich, need to exclude
  if (!(any(grepl("\\benrichmentscore\\b", to_check)) & any(grepl("\\bleadingedge\\b", to_check)))) {
    if (any(grepl("foldenrich|enrichment", to_check))) {
      colnames(enrich_df)[grepl("foldenrich|enrichment", to_check)] <- "FoldEnrich"
    } else {
      stop('Please rename the fold enrichment column as "FoldEnrich"!')
      # head(enrich_df[1:2,])
      # message("Cannot auto-select fold enrichment column...","\n",
      #         "Please specify the column number which includes fold enrichment...")
      # answer <- scan(what = "character", n =1,quiet =T)
      # message('Choose the No. ',answer,' column as fold enrichment...')
      # enrich_df = enrich_df %>% dplyr::rename(FoldEnrich = eval(parse(text = answer)))
    }
  }


  ## find GeneRatio col
  if (any(grepl("generatio", to_check))) {
    colnames(enrich_df)[grepl("generatio", to_check)] <- "GeneRatio"
    enrich_df <- enrich_df %>% dplyr::mutate(GeneRatio = sapply(.$GeneRatio, function(x) eval(parse(text = x))))
  } else {
    # gsea
    if (any(grepl("setsize", to_check))) {
      colnames(enrich_df)[grepl("setsize", to_check)] <- "setSize"
      enrich_df <- enrich_df %>% dplyr::mutate(GeneRatio = as.numeric(Count) / as.numeric(setSize))
    } else if (any(grepl("\\([0-9]{,4}\\)", colnames(enrich_df)))) {
      # panther result
      setsize <- colnames(enrich_df)[grepl("\\([0-9]{,4}\\)", colnames(enrich_df))] %>%
        stringr::str_remove(., ".*\\(") %>%
        stringr::str_remove(., "\\)") %>%
        as.numeric() %>%
        min()
      enrich_df <- enrich_df %>% dplyr::mutate(GeneRatio = as.numeric(Count) / setsize)
    } else if (apply(enrich_df, 2, function(x) length(unique(x)) == 1)) {
      setsize <- enrich_df[1, apply(enrich_df, 2, function(x) length(unique(x)) == 1)] %>%
        stringr::str_remove("0") %>%
        as.numeric() %>%
        min()
      enrich_df <- enrich_df %>% dplyr::mutate(GeneRatio = as.numeric(Count) / setsize)
    }
  }

  ## find pvalue col
  if (any(grepl("pvalue", to_check))) {
    colnames(enrich_df)[grepl("pvalue", to_check)] <- "pvalue"
  } else if (any(grepl("\\buncorrectedpvalue\\b", to_check))) {
    colnames(enrich_df)[grepl("\\buncorrectedpvalue\\b", to_check)] <- "pvalue"
  }

  ## find p.adjust col
  if (any(grepl("p.adjust", to_check))) {
    colnames(enrich_df)[grepl("p.adjust", to_check)] <- "p.adjust"
  } else if (any(grepl("\\bcorrectedpvalue\\b", to_check))) {
    colnames(enrich_df)[grepl("\\bcorrectedpvalue\\b", to_check)] <- "p.adjust"
  }

  ## find qvalue col
  if (any(grepl("qvalue", to_check))) {
    colnames(enrich_df)[grepl("qvalue", to_check)] <- "qvalue"
  } else if (any(grepl("fdr", to_check) & !grepl("fdrrate", to_check))) {
    colnames(enrich_df)[grepl("fdr", to_check)] <- "qvalue"
  }

  return(enrich_df)
}

utils::globalVariables(c("Count", "setSize"))
