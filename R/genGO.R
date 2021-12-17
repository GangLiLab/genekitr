#' Gene GO enrichment analysis
#'
#' @param id A vector of gene id which can be entrez, ensembl or symbol.
#' @param group_list A list of gene id groups, default is NULL.
#' @param org  Organism name from `biocOrg_name`.
#' @param ont  One of "bp", "mf", and "cc" subontologies, or "all" for all
#'   three.
#' @param use_symbol Logical to set result gene id as gene symbol, default is TRUE.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH",
#'   "BY", "fdr", "none".
#' @param pvalueCutoff Adjusted pvalue cutoff, default is 0.05.
#' @param qvalueCutoff Adjusted pvalue cutoff, default is 0.1.
#' @param minGSSize Minimal size of each gene set for analyzing, default is 10.
#' @param maxGSSize Maximal size of each gene set for analyzing, default is 500.
#' @param universe Background genes. If missing, then all gene list in
#'   orgdb will be used as background.
#' @param ... other argument to `enrichGO` function
#' @importFrom dplyr select filter pull mutate %>%
#' @importFrom stringr str_split
#' @importFrom clusterProfiler enrichGO
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
#' @return A `data.frame` contains gene ratio and fold enrichment.
#' @export
#'
#' @examples
#' \dontrun{
#' data(geneList, package = "genekitr")
#'
#' # only gene ids
#' id <- names(geneList)[1:100]
#' ego <- genGO(id,
#'   org = "human", ont = "cc", pvalueCutoff = 0.01,
#'   qvalueCutoff = 0.1, use_symbol = FALSE
#' )
#' head(ego)
#'
#' # gene id with groups
#' id <- c(head(names(geneList),50),tail(names(geneList),50))
#' group <- list(group1  = c(rep('up',50),rep('down',50)),
#'   group2 = c(rep('A',40),rep('B',60)))
#'
#' gego <- genGO(id, group_list = group,
#'   org = "human", ont = "bp", pvalueCutoff = 0.1,
#'   qvalueCutoff = 1, use_symbol = FALSE
#' )
#'
#' }
genGO <- function(id,
                  group_list = NULL,
                  org,
                  ont,
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
  if(!is.null(group_list) & lapply(group_list, function(x) length(x) == length(id)) %>% unlist() %>% sum() == 0){
    stop('Each element in "group_list" should have same length with gene id!')
  }

  bioc_org <- mapBiocOrg(org)
  org <- mapEnsOrg(org)
  pkg <- paste0("org.", bioc_org, ".eg.db")
  keyType <- gentype(id, org)

  #--- codes ---#
  ## only gene ids
  if(is.null(group_list)){
    ego <- suppressMessages(
      clusterProfiler::enrichGO(
        gene = id, OrgDb = pkg, keyType = keyType, ont = toupper(ont),
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        universe = universe,
        qvalueCutoff = qvalueCutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        ...
      )
    )

    if (nrow(as.data.frame(ego)) == 0) {
      stop("No GO terms enriched ...")
    }

    if (use_symbol) {
      info <- genInfo(id, org, unique = T)
      new_geneID <- stringr::str_split(ego$geneID, "\\/") %>%
        lapply(., function(x) {
          info %>%
            dplyr::filter(input_id %in% x) %>%
            dplyr::pull(symbol)
        }) %>%
        sapply(., paste0, collapse = "/")
      new_ego <- ego %>%
        as.data.frame() %>%
        dplyr::mutate(geneID = new_geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    } else {
      new_ego <- ego %>%
        as.data.frame() %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    }
  }else{
    ## gene id plus groups
    df <- as.data.frame(group_list) %>% dplyr::mutate(id = id)
    lego <- clusterProfiler::compareCluster(eval(parse(text =paste0('id~',paste(colnames(df)[-ncol(df)],collapse = '+')))),
                                            data=df,
                                            fun='enrichGO', OrgDb=pkg,
                                            pvalueCutoff = pvalueCutoff,
                                            pAdjustMethod = pAdjustMethod,
                                            qvalueCutoff = qvalueCutoff,
                                            ...)

    if (nrow(as.data.frame(lego)) == 0) {
      stop("No GO terms enriched ...")
    }

    if (use_symbol) {
      info <- genInfo(id, org, unique = T)
      new_geneID <- stringr::str_split(lego@compareClusterResult$geneID, "\\/") %>%
        lapply(., function(x) {
          info %>%
            dplyr::filter(input_id %in% x) %>%
            dplyr::pull(symbol)
        }) %>%
        sapply(., paste0, collapse = "/")
      new_ego <- lego %>%
        as.data.frame() %>%
        dplyr::mutate(geneID = new_geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    }else{
      new_ego <- lego %>%
        as.data.frame() %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    }
  }


  return(new_ego)
}


utils::globalVariables(c("input_id"))
