#' Gene Set Enrichment Analysis (FCS method)
#'
#' @param genelist Order ranked genelist in decreasing order, gene can be
#'   entrez, ensembl or symbol.
#' @param geneset Gene set is a two-column data.frame with term id and gene id.
#' Please use package `geneset` to select available gene set or make new one.
#' @param padj_method One of "BH", "BY", "bonferroni","fdr","hochberg",
#' "holm", "hommel", "none"
#' @param p_cutoff Numeric of cutoff for both unadjusted and adjusted pvalue, default is 0.05.
#' @param q_cutoff Numeric of cutoff for qvalue, default is 0.05.
#' @param min_gset_size Numeric of minimal size of each geneset for analyzing,
#'   default is 10.
#' @param max_gset_size Numeric of maximal size of each geneset for analyzing,
#'   default is 500.
#' @param set_seed GSEA permutations are performed using random reordering,
#' which causes slightly difference results after every time running.
#' If user want to get same result every time for same input,
#' please set `set_seed = TRUE` or `set.seed()` prior to running.
#' @importFrom dplyr pull filter arrange mutate relocate
#' @importFrom stringr str_split
#' @importFrom geneset getEnrichrdb getGO getHgDisease getKEGG getMesh getMsigdb getReactome getWiki
#' @importFrom clusterProfiler GSEA
#' @importFrom rlang .data
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' \dontrun{
#' # only gene ids
#' data(geneList, package = "genekitr")
#' gs <- geneset::getGO(org = "human",ont = "mf")
#' gse <- genGSEA(genelist = geneList, geneset = gs)
#'
#' }

genGSEA <- function(genelist,
                    geneset,
                    padj_method = "BH",
                    p_cutoff = 0.05,
                    q_cutoff = 0.05,
                    min_gset_size = 10,
                    max_gset_size = 500,
                    set_seed = FALSE){

  #--- args ---#
  id <- as.character(names(genelist))
  if(missing(geneset)) stop('Please provide gene set...\nWe recommend to use package `geneset` to select available gene set or make new one.')

  genesetType <- geneset$type
  transToSym <- ifelse(genesetType %in% c("enrichrdb","bp",'mf','cc',"covid19"), TRUE, FALSE)

  org <- geneset$organism
  ens_org <- tryCatch( mapEnsOrg(org), error = function(e) NULL)
  if(is.null(ens_org)) stop('Please make sure the gene types of genelist and geneset are the same.')
  keyType <- gentype(id = id, org = ens_org)

  #--- initialize ---#
  # input id must be symbol or entrezid for c gene sets

  if(transToSym){
    id_dat <- suppressMessages(transId(id, "symbol", ens_org, unique = T))
    genelist <- genelist[names(genelist)%in%id_dat$input_id]
    names(genelist) <- id_dat$symbol
  }else if(keyType != "ENTREZID"){
    id_dat <- suppressMessages(transId(id, "entrezid", ens_org, unique = T))
    genelist <- genelist[names(genelist)%in%id_dat$input_id]
    names(genelist) <- id_dat$entrezid
  }else if(keyType == "ENTREZID"){
    id_dat <- suppressMessages(transId(id, "symbol", ens_org, unique = T)) %>%
      dplyr::relocate(input_id,.after = symbol)
    genelist <- genelist[names(genelist)%in%id_dat$input_id]
    # id <- id_dat$input_id
  }

  #--- analyse ---#
  fcs <- suppressWarnings(
    GSEA(
      geneList = genelist,
      pvalueCutoff = p_cutoff,
      pAdjustMethod = padj_method,
      minGSSize = min_gset_size,
      maxGSSize = max_gset_size,
      TERM2GENE = geneset$geneset,
      TERM2NAME = geneset$geneset_name,
      exponent = 1,
      eps  = 0,
      verbose = FALSE,
      seed = set_seed,
      by = 'fgsea'
    ))


  #--- post-process ---#
  if (nrow(as.data.frame(fcs)) == 0) {
    stop("No terms enriched ...")
  }else{
    exponent <- fcs@params[["exponent"]]
    fcs <- fcs %>%
      as.data.frame() %>%
      as.enrichdat() %>%
      dplyr::select(-GeneRatio) %>%
      dplyr::filter(qvalue < q_cutoff)
  }

  ## transToSym means geneset in "enrichrdb","go" and "covid19"
  if(!transToSym){
    # part 1-1
    if(keyType != "SYMBOL"){
      # part 1-1-1
      if(keyType == 'ENTREZID'){
        new_geneID <- get_symbol(fcs$geneID,ens_org)
        new_fcs <- fcs %>%
          dplyr::mutate(geneID_symbol = new_geneID) %>%
          dplyr::relocate(geneID_symbol, .after = geneID)
      }else{
        # part 1-1-2
        old_geneID <- replace_id(id_dat,fcs$geneID)
        new_geneID <- get_symbol(fcs$geneID,ens_org)
        new_fcs <- fcs %>%
          dplyr::mutate(geneID_symbol = new_geneID) %>%
          dplyr::mutate(geneID = old_geneID) %>%
          dplyr::relocate(geneID_symbol, .after = geneID)
      }
      # part 1-2
    }else{
      # new_fcs <- fcs
      old_geneID <- replace_id(id_dat,fcs$geneID)
      # new_geneID <- get_symbol(fcs$geneID,ens_org)
      new_fcs <- fcs %>%
        dplyr::mutate(geneID = old_geneID)
    }

    # part 2
  }else{
    # part 2-1
    if(keyType != "SYMBOL" ){
      old_geneID <- replace_id(id_dat,fcs$geneID)
      new_fcs <- fcs %>%
        dplyr::mutate(geneID_symbol = geneID) %>%
        dplyr::mutate(geneID = old_geneID) %>%
        dplyr::relocate(geneID_symbol, .after = geneID)
    }else{
      # part 2-2
      new_fcs <- fcs
    }

  }

  ## modify id column name for GO
  bioc_org <- ensOrg_name %>%
    dplyr::filter(tolower(latin_short_name) %in% geneset$organism) %>%
    dplyr::pull(bioc_name) %>%
    stringr::str_to_sentence()

  if(genesetType %in% c('bp','cc','mf')){
    colnames(new_fcs)[1] = paste0(bioc_org,'_',toupper(genesetType),'_ID')
  }

  ## save as list
  genelist_df = data.frame(ID = names(genelist), logfc = genelist)
  exponent = data.frame(exponent = exponent)
  org = data.frame(org = org)
  new_geneset <- geneset$geneset

  res <- list(gsea_df = new_fcs, genelist = genelist_df, geneset = new_geneset,  exponent = exponent, org = org)


  return(res)

}



