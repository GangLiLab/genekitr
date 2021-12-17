#' Gene enrichment of KEGG analysis
#'
#' @param id A vector of entrez gene.
#' @param group_list A list of gene id groups, default is NULL.
#' @param org  KEGG organism name from `keggOrg_name`.
#' @param use_symbol Logical to set result gene id as gene symbol, default is TRUE.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH",
#'   "BY", "fdr", "none".
#' @param pvalueCutoff Numberic of adjusted pvalue cutoff, default is 0.05.
#' @param qvalueCutoff Numberic of adjusted pvalue cutoff, default is 0.1.
#' @param minGSSize Numberic of minimal size of each geneSet for analyzing,
#'   default is 10.
#' @param maxGSSize Numberic of maximal size of each geneSet for analyzing,
#'   default is 500.
#' @param universe Background genes. If missing, the orgdb all gene list will be
#'   used as background.
#' @param ... Other argument to `enrichKEGG` function
#' @importFrom dplyr  %>% mutate filter
#' @importFrom stringr  str_split
#' @importFrom stringi stri_omit_na
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' \dontrun{
#' # only gene ids
#' data(geneList, package = "genekitr")
#' id <- names(geneList)[1:100]
#' keg <- genKEGG(id, org = "human")
#'
#' # gene id with groups
#' id <- c(head(names(geneList),50),tail(names(geneList),50))
#' group <- list(group1  = c(rep('up',50),rep('down',50)),
#'               group2 = c(rep('A',40),rep('B',60)))
#' gkeg <- genKEGG(id, group_list = group,
#'               org = "human", pvalueCutoff = 0.1,
#'               qvalueCutoff = 0.1, use_symbol = FALSE
#' )
#'
#' }
#'
genKEGG <- function(id,
                    group_list = NULL,
                    org,
                    use_symbol = TRUE,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.1,
                    minGSSize = 10,
                    maxGSSize = 500,
                    universe,
                    ...) {

  #--- args ---#
  stopifnot(is.character(id))
  if (missing(universe)) universe <- NULL

  kegg_org <- mapKeggOrg(org)
  full_latin <- getKeggLatin(kegg_org)
  ensorg <- ensOrg_name_data()
  rm(ensOrg_name, envir = .GlobalEnv)

  if(full_latin %in% ensorg$latin_full_name){
    ens_org <- mapEnsOrg(full_latin)
    keyType <- gentype(id, ens_org)
    if (!keyType %in% c("ENTREZID")) {
      message(paste0(keyType), " gene will be mapped to entrez id")
      trans_id <- suppressMessages(transId(id, "entrezid", ens_org)) %>% stringi::stri_remove_na()
    }else{
      trans_id <- id
    }
    info <- genInfo(trans_id, ens_org, unique = T)

  }else{
    trans_id <- id
    use_symbol <- FALSE
  }

  #--- codes ---#
  ## only gene ids
  if(is.null(group_list)){
    keg <- suppressMessages(
      clusterProfiler::enrichKEGG(
        gene = trans_id, organism = kegg_org, keyType = "kegg",
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        qvalueCutoff = qvalueCutoff,
        universe = universe,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        ...
      )
    )

    if (nrow(as.data.frame(keg)) == 0) {
      stop("No KEGG terms enriched ...")
    }

    if (use_symbol) {
      new_geneID <- stringr::str_split(keg$geneID, "\\/") %>%
        lapply(., function(x) {
          info %>%
            dplyr::filter(input_id %in% x) %>%
            dplyr::pull(symbol)
        }) %>%
        sapply(., paste0, collapse = "/")
      new_keg <- keg %>%
        as.data.frame() %>%
        dplyr::mutate(geneID = new_geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    } else {
      new_keg <- keg %>%
        as.data.frame() %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    }
  }else{
    ## gene id plus groups
    df <- as.data.frame(group_list) %>% dplyr::mutate(id = id)
    lkeg <- clusterProfiler::compareCluster(eval(parse(text =paste0('id~',paste(colnames(df)[-ncol(df)],collapse = '+')))),
                                            data=df,
                                            fun='enrichKEGG', organism = kegg_org,
                                            pvalueCutoff = pvalueCutoff,
                                            pAdjustMethod = pAdjustMethod,
                                            qvalueCutoff = qvalueCutoff,
                                            universe = universe,
                                            minGSSize = minGSSize,
                                            maxGSSize = maxGSSize,
                                            ...)

    if (nrow(as.data.frame(lkeg)) == 0) {
      stop("No KEGG terms enriched ...")
    }
    if (use_symbol) {
      new_geneID <- stringr::str_split(lkeg@compareClusterResult$geneID, "\\/") %>%
        lapply(., function(x) {
          info %>%
            dplyr::filter(input_id %in% x) %>%
            dplyr::pull(symbol)
        }) %>%
        sapply(., paste0, collapse = "/")
      new_keg <- lkeg %>%
        as.data.frame() %>%
        dplyr::mutate(geneID = new_geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    } else {
      new_keg <- lkeg %>%
        as.data.frame() %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    }

  }

  return(new_keg)
}


utils::globalVariables(c("input_id","symbol"))

