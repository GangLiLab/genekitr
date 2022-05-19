#' Gene enrichment of KEGG analysis
#'
#' @param id A vector of entrez gene.
#' @param group_list A list of gene id groups, default is NULL.
#' @param org  KEGG organism name from `keggOrg_name`.
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
#' @importFrom dplyr pull filter arrange mutate relocate
#' @importFrom stringr str_split
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom rlang .data
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' \dontrun{
#' # only gene ids
#' data(geneList, package = "genekitr")
#' id <- names(geneList)[abs(geneList) > 1]
#' keg <- genKEGG(id, org = "human")
#'
#' # gene id with groups
#' id <- c(head(names(geneList),100),tail(names(geneList),100))
#' group <- list(group1  = c(rep('up',100),rep('down',100)),
#'               group2 = c(rep('A',130),rep('B',70)))
#' gkeg <- genKEGG(id, group_list = group,
#'               org = "human", pvalueCutoff = 0.05,
#'               qvalueCutoff = 0.05)
#'
#' }
#'
genKEGG <- function(id,
                    group_list = NULL,
                    org,
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
    keyType <- gentype(id = id, org=ens_org)
    if (!keyType %in% c("ENTREZID")) {
      message(paste0(keyType), " gene will be mapped to entrez id")
      trans_id <- suppressMessages(transId(id, "entrezid", ens_org)) %>%
        dplyr::pull(entrezid)
    }else{
      trans_id <- id
    }
  }else{
    trans_id <- id
  }

  #--- codes ---#
  ## NO GROUP INFO
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

    # transform entrez id to symbol
    keg_id = stringr::str_split(keg$geneID, "\\/") %>% unlist()
    id_all = suppressMessages(transId(keg_id,'symbol',org = org,unique = T))

    new_geneID <- stringr::str_split(keg$geneID, "\\/") %>%
      lapply(., function(x) {
        id_all %>% dplyr::filter(input_id %in% x) %>%
          dplyr::arrange(match(input_id, x)) %>%
          dplyr::pull(symbol)
      }) %>%
      sapply(., paste0, collapse = "/")

    new_keg <- keg %>%
      as.data.frame() %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol,.after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  }else{
    ## WITH GROUP INFO
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
    # transform entrez id to symbol
    keg_id = stringr::str_split(lkeg@compareClusterResult$geneID, "\\/") %>% unlist()
    id_all = suppressMessages(transId(keg_id,'symbol',org = org,unique = T))

    new_geneID <- stringr::str_split(lkeg@compareClusterResult$geneID, "\\/") %>%
      lapply(., function(x) {
        id_all %>% dplyr::filter(input_id %in% x) %>%
          dplyr::arrange(match(input_id, x)) %>%
          dplyr::pull(symbol)
      }) %>%
      sapply(., paste0, collapse = "/")

    new_keg <- lkeg %>%
      as.data.frame() %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol,.after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  }

  # add rich factor
  new_keg <- new_keg %>%
    dplyr::mutate(RichFactor = Count/as.numeric(sub("/\\d+","", BgRatio)))

  return(new_keg)
}


utils::globalVariables(c("input_id","symbol","entrezid"))

