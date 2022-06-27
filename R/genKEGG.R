#' Gene ORA-KEGG enrichment analysis
#'
#' @param id A vector of gene id which can be entrezid, ensembl, symbol or uniprot.
#' @param group_list A list of gene id groups, default is NULL.
#' @param org  KEGG organism name from `keggOrg_name`.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH",
#'   "BY", "fdr", "none".
#' @param pvalueCutoff Numberic of pvalue cutoff, default is 0.05.
#' @param qvalueCutoff Numberic of adjusted pvalue cutoff, default is 0.05.
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
#' id <- c(head(names(geneList), 100), tail(names(geneList), 100))
#' group <- list(
#'   group1 = c(rep("up", 100), rep("down", 100)),
#'   group2 = c(rep("A", 130), rep("B", 70))
#' )
#' gkeg <- genKEGG(id,
#'   group_list = group,
#'   org = "human", pvalueCutoff = 0.05,
#'   qvalueCutoff = 0.05
#' )
#' }
#'
genKEGG <- function(id,
                    org,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 10,
                    maxGSSize = 500,
                    group_list = NULL,
                    universe,
                    ...) {

  #--- args ---#
  # stopifnot(is.character(id))
  id <- as.character(id)
  if (missing(universe)) universe <- NULL

  kegg_org <- mapKeggOrg(org)
  full_latin <- getKeggLatin(kegg_org)
  # ensorg_data <- ensOrg_name_data()
  ens_org <- mapEnsOrg(full_latin)
  # rm(ensOrg_name, envir = .GlobalEnv)

  keyType <- gentype(id = id, org = ens_org)

  # here we convert all symbol and alias to symbol
  old_id <- id
  if (keyType == "SYMBOL") {
    id_dat1 <- suppressMessages(transId(id, "symbol", ens_org, unique = T))
    id <- id_dat1 %>% dplyr::pull(symbol)
  }

  if (keyType != "ENTREZID") {
    message(paste0(keyType), " gene will be mapped to entrez id")
    id_dat2 <- suppressMessages(transId(id, "entrezid", ens_org, unique = T))
    entrez_id <- id_dat2 %>% dplyr::pull(entrezid)
  } else {
    entrez_id <- id
  }


  if (!is.null(universe)) {
    universe <- suppressMessages(transId(universe, transTo = "entrezid", ens_org, unique = T)[, 2])
  }

  #--- analyse ---#
  ## NO GROUP INFO
  if (is.null(group_list)) {
    keg <- suppressMessages(
      clusterProfiler::enrichKEGG(
        gene = entrez_id, organism = kegg_org, keyType = "kegg",
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        qvalueCutoff = qvalueCutoff,
        universe = universe,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        ...
      )
    )

  } else {
    ## WITH GROUP INFO
    df <- as.data.frame(group_list) %>% dplyr::mutate(id = entrez_id)
    if(ncol(df) >2 ){
      df <- df %>%
        dplyr::mutate(Cluster = apply(df[,1:(ncol(df)-1)],1,paste,collapse="."))
    }else{
      df <- df %>%
        dplyr::mutate(Cluster = .[[1]])
    }

    keg <- df %>%
      dplyr::select(id, Cluster) %>%
      split(.$Cluster) %>%
      lapply(function(x) x %>% dplyr::pull(id)) %>%
      lapply(function(x)
        suppressMessages(clusterProfiler::enrichKEGG(
          gene = x, organism = kegg_org, keyType = "kegg",
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          qvalueCutoff = qvalueCutoff,
          universe = universe,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          ...
        )) %>%
          as.data.frame()) %>%
      do.call(rbind,.) %>%
      dplyr::mutate(Cluster = gsub("\\.[^\\.]*$", "", rownames(.), perl=TRUE)) %>%
      dplyr::relocate(Cluster,.before = dplyr::everything()) %>%
      `rownames<-`(seq_len(nrow(.)))
  }

  if (nrow(as.data.frame(keg)) == 0) {
    stop("No KEGG terms enriched ...")
  }else{
    keg = as.data.frame(keg)
  }

  #--- get geneID_symbol ---#
  # new id
  new_geneID <- get_symbol(keg$geneID,ens_org)

  # input SYMBOL and no alias
  if( keyType == "SYMBOL" & identical(old_id, id) ){
    new_keg <- keg %>%
      dplyr::mutate(geneID = new_geneID) %>%
      # dplyr::relocate(geneID, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else if(keyType == "SYMBOL" & !identical(old_id, id)){
    old_geneID <- replace_id(id_dat2,keg$geneID) %>%
      replace_id(id_dat1,.)

    new_keg <- keg %>%
      dplyr::mutate(geneID = old_geneID) %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else if(keyType == "ENTREZID"){
    new_keg <- keg %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else {
    # input other types (including alias)
    old_geneID <- replace_id(id_dat2,keg$geneID)

    new_keg <- keg %>%
      dplyr::mutate(geneID = old_geneID) %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()
  }

  #--- add rich factor ---#
  new_keg <- new_keg %>%
    dplyr::mutate(RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

  return(new_keg)
}


get_symbol <- function(id,org){
  ori_id <- stringr::str_split(id, "\\/") %>% unlist() %>% unique()
  id_all <- suppressMessages(transId(ori_id, "symbol", org, unique = T))

  new_geneID <- stringr::str_split(id, "\\/") %>%
    lapply(., function(x) {
      id_all %>%
        dplyr::filter(input_id %in% x) %>%
        dplyr::arrange(match(input_id, x)) %>%
        dplyr::pull(symbol)
    }) %>%
    sapply(., paste0, collapse = "/")
  return(new_geneID)
}

replace_id <- function(dat, id){
  stringr::str_split(id, "\\/") %>%
    lapply(., function(x) {
      dat %>% dplyr::filter(.[[2]]%in%x) %>% dplyr::pull(1) %>%
        paste0(.,collapse = '/')
    }) %>%  do.call(rbind,.) %>% as.character()
}


utils::globalVariables(c("input_id", "symbol", "entrezid"))
