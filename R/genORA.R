#' Gene Over-Representation enrichment Analysis (ORA method)
#'
#' @param id A vector of gene id which can be entrezid, ensembl, symbol or uniprot.
#' @param geneset Gene set is a two-column data.frame with term id and gene id.
#' Please use package `geneset` to select available gene set or make new one.
#' @param group_list A list of gene group information, default is NULL.
#' @param padj_method One of "BH", "BY", "bonferroni","fdr","hochberg",
#' "holm", "hommel", "none"
#' @param p_cutoff Numeric of cutoff for both unadjusted and adjusted pvalue, default is 0.05.
#' @param q_cutoff Numeric of cutoff for qvalue, default is 0.15.
#' @param min_gset_size Numeric of minimal size of each geneset for analyzing,
#'   default is 10.
#' @param max_gset_size Numeric of maximal size of each geneset for analyzing,
#'   default is 500.
#' @param universe Character of background genes. If missing, all genes in
#' geneset will be used as background.
#' @importFrom dplyr pull filter arrange mutate relocate
#' @importFrom stringr str_split
#' @importFrom geneset getEnrichrdb getGO getHgDisease getKEGG getMesh getMsigdb getReactome getWiki
#' @importFrom clusterProfiler enricher
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
#' gs <- geneset::getGO(org = "human",ont = "mf")
#' ora <- genORA(id, geneset = gs)
#'
#' # gene id with groups
#' id <- c(head(names(geneList), 100), tail(names(geneList), 100))
#' group <- list(
#'   group1 = c(rep("up", 100), rep("down", 100)),
#'   group2 = c(rep("A", 130), rep("B", 70))
#' )
#' gora <- genORA(id, geneset = gs, group_list = group)
#'
#' }

genORA <- function(id,
                   geneset,
                   group_list = NULL,
                   padj_method = "BH",
                   p_cutoff = 0.05,
                   q_cutoff = 0.15,
                   min_gset_size = 10,
                   max_gset_size = 500,
                   universe){

  #--- args ---#
  id <- as.character(id)
  if (missing(universe)) universe <- NULL
  if(missing(geneset)) stop('Please provide gene set...\nWe recommend to use package `geneset` to select available gene set or make new one.')

  genesetType <- geneset$type
  transToSym <- ifelse(genesetType %in% c("enrichrdb","bp","mf","cc","covid19"), TRUE, FALSE)

  org <- geneset$organism
  ens_org <- mapEnsOrg(org)
  keyType <- gentype(id = id, org = ens_org)
  if(!is.null(universe)) UnikeyType <- gentype(id = universe, org = ens_org)

  #--- initialize ---#
  # input id must be symbol or entrezid for id_dat gene sets

  if(transToSym){
    id_dat <- suppressMessages(transId(id, "symbol", ens_org, unique = T))
    id <- id_dat$symbol
    if (!is.null(universe) ) {
      if(UnikeyType != "SYMBOL")
        universe <- suppressMessages(transId(universe, transTo = "symbol", ens_org, unique = T)) %>% dplyr::pull(symbol)
    }
  }else if(keyType != "ENTREZID"){
    id_dat <- suppressMessages(transId(id, "entrezid", ens_org, unique = T))
    id <- id_dat$entrezid
    if (!is.null(universe) ) {
      if(UnikeyType != "ENTREZID")
        universe <- suppressMessages(transId(universe, transTo = "entrezid", ens_org, unique = T)) %>% dplyr::pull(entrezid)
    }
  }else if(keyType == "ENTREZID"){
    id_dat <- suppressMessages(transId(id, "symbol", ens_org, unique = T)) %>%
      dplyr::relocate(input_id,.after = symbol)
    id <- id_dat$input_id
    if (!is.null(universe) ) {
      if(UnikeyType != "ENTREZID")
        universe <- suppressMessages(transId(universe, transTo = "entrezid", ens_org, unique = T)) %>% dplyr::pull(entrezid)
    }
  }

  #--- analyse ---#
  if (is.null(group_list)) {
    ## NO GROUP INFO
    ora <- suppressMessages(
      enricher(
        gene = id,
        pvalueCutoff = p_cutoff,
        qvalueCutoff = q_cutoff,
        pAdjustMethod = padj_method,
        universe = universe,
        minGSSize = min_gset_size,
        maxGSSize = max_gset_size,
        TERM2GENE = geneset$geneset,
        TERM2NAME = geneset$geneset_name
      ))
  }else {
    ## WITH GROUP INFO
    df <- as.data.frame(group_list)
    if(nrow(df) != length(id)) stop('Please check "group_list"! It should have the same length of elements as input id.')
    df$id = id

    if(ncol(df) >2 ){
      df <- df %>%
        dplyr::mutate(Cluster = apply(df[,1:(ncol(df)-1)],1,paste,collapse="."))
    }else{
      df <- df %>%
        dplyr::mutate(Cluster = .[[1]])
    }

    ora <- df %>%
      dplyr::select(id, Cluster) %>%
      split(.$Cluster) %>%
      lapply(function(x) x %>% dplyr::pull(id)) %>%
      lapply(function(x)
        suppressMessages(
          clusterProfiler::enricher(
            gene = x,
            pvalueCutoff = p_cutoff,
            qvalueCutoff = q_cutoff,
            pAdjustMethod = padj_method,
            universe = universe,
            minGSSize = min_gset_size,
            maxGSSize = max_gset_size,
            TERM2GENE = geneset$geneset,
            TERM2NAME = geneset$geneset_name
          )) %>%
          as.data.frame()) %>%
      do.call(rbind,.) %>%
      dplyr::mutate(Cluster = gsub("\\.[^\\.]*$", "", rownames(.), perl=TRUE)) %>%
      dplyr::relocate(Cluster,.before = dplyr::everything()) %>%
      `rownames<-`(seq_len(nrow(.)))
  }

  #--- post-process ---#
  if (nrow(as.data.frame(ora)) == 0) {
    stop("No terms enriched ...")
  }else{
    ora = as.data.frame(ora)
  }

  ## transToSym means geneset in "enrichrdb","go" and "covid19"
  if(!transToSym){
    # part 1-1
    if(keyType != "SYMBOL"){
      # part 1-1-1
      if(keyType == 'ENTREZID'){
        new_geneID <- get_symbol(ora$geneID,ens_org)
        new_ora <- ora %>%
          dplyr::mutate(geneID_symbol = new_geneID) %>%
          dplyr::relocate(geneID_symbol, .after = geneID)
      }else{
        # part 1-1-2
        old_geneID <- replace_id(id_dat,ora$geneID)
        new_geneID <- get_symbol(ora$geneID,ens_org)
        new_ora <- ora %>%
          dplyr::mutate(geneID_symbol = new_geneID) %>%
          dplyr::mutate(geneID = old_geneID) %>%
          dplyr::relocate(geneID_symbol, .after = geneID)
      }
    # part 1-2
    }else{
      # new_ora <- ora
      old_geneID <- replace_id(id_dat,ora$geneID)
      # new_geneID <- get_symbol(ora$geneID,ens_org)
      new_ora <- ora %>%
        dplyr::mutate(geneID = old_geneID)
    }

  # part 2
  }else{
    # part 2-1
    if(keyType != "SYMBOL" ){
      old_geneID <- replace_id(id_dat,ora$geneID)
      new_ora <- ora %>%
        dplyr::mutate(geneID_symbol = geneID) %>%
        dplyr::mutate(geneID = old_geneID) %>%
        dplyr::relocate(geneID_symbol, .after = geneID)
    }else{
      # part 2-2
      new_ora <- ora
    }

  }

  ## add fold enrcih/rich factor
  new_ora <- new_ora %>%
    calcFoldEnrich() %>%
    as.enrichdat() %>%
    dplyr::mutate(RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

  ## modify id column name for GO
  bioc_org <- ensOrg_name %>%
    dplyr::filter(tolower(latin_short_name) %in% geneset$organism) %>%
    dplyr::pull(bioc_name) %>%
    stringr::str_to_sentence()

  if(genesetType %in% c('bp','cc','mf')){
    colnames(new_ora)[1] = paste0(bioc_org,'_',toupper(genesetType),'_ID')
  }

  return(new_ora)

}



