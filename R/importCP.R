#' Import clusterProfiler result
#'
#' @param object clusterProfiler object.
#' @param type object type from "go", "gsea" and "other". "other" includes ORA (over-representation analysis) of KEGG, DOSE,...
#'
#' @importFrom dplyr pull mutate rename relocate select
#' @importFrom rlang .data
#' @return  `data.frame`
#' @export

importCP <- function(object,
                     type = c('go','gsea','other')) {
  type <- match.arg(type)

  #--- IMPORT GO ---#
  if(type == 'go'){
    ont <- object@ontology
    org <- mapEnsOrg(object@organism)
    ensorg_data <- ensOrg_name_data()
    bioc_org <- ensorg_data[ensorg_data$latin_full_name == object@organism, 'bioc_name']

    ego <- as.data.frame(object)
    id <- strsplit(ego$geneID,'\\/') %>% unlist() %>% unique()

    old_id <- id
    keyType <- gentype(id = id, org = org)
    if (keyType == "SYMBOL") {
      id_dat <- suppressMessages(transId(id, "symbol", org, unique = T))
      id <- id_dat %>% dplyr::pull(symbol)
    }

    ### get geneID_symbol
    # new id
    new_geneID <- get_symbol(ego$geneID,org)

    # input SYMBOL and no alias
    if( keyType == "SYMBOL" & identical(old_id, id) ){
      new_ego <- ego %>%
        calcFoldEnrich() %>%
        as.enrichdat()

    } else if(keyType == "SYMBOL" & !identical(old_id, id) ){
      # input SYMBOL and have alias
      old_geneID <- replace_id(id_dat,ego$geneID)

      new_ego <- ego %>%
        dplyr::mutate(geneID = old_geneID) %>%
        dplyr::mutate(geneID_symbol = new_geneID) %>%
        dplyr::relocate(geneID_symbol, .after = geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    } else {
      # input other types
      new_ego <- ego %>%
        dplyr::mutate(geneID_symbol = new_geneID) %>%
        dplyr::relocate(geneID_symbol, .after = geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    }

    ### add rich factor
    new_obj <- new_ego %>%
      dplyr::mutate(RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>%
      dplyr::rename(!!paste(bioc_org, ont, "ID", sep = "_") := ID)
  }

  #--- IMPORT GSEA ---#
  if(type == 'gsea'){
    egmt <- object
    ens_org <- mapEnsOrg(object@organism)
    exponent <- egmt@params[["exponent"]]
    geneset <- egmt@geneSets %>%
      reshape2::melt() %>%
      stats::setNames(c('entrez_gene','gs_name')) %>%
      dplyr::relocate(gs_name,.before = dplyr::everything())

    genelist = egmt@geneList
    egmt <- egmt %>%
      as.data.frame() %>%
      as.enrichdat() %>%
      dplyr::select(-GeneRatio)
    id <- strsplit(egmt$geneID,'\\/') %>% unlist() %>% unique()

    old_id <- id
    keyType <- gentype(id = id, org = ens_org)
    if (keyType == "SYMBOL") {
      id_dat1 <- suppressMessages(transId(id, "symbol", ens_org, unique = T))
      id <- id_dat1 %>% dplyr::pull(symbol)
    }

    if (keyType != "ENTREZID") {
      id_dat2 <- suppressMessages(transId(id, "entrezid", ens_org, unique = T))
    }

    ### get geneID_symbol
    # new id
    new_geneID <- get_symbol(egmt$geneID,ens_org)

    # input SYMBOL and no alias
    if( keyType == "SYMBOL" & identical(old_id, id) ){
      egmt <- egmt %>%
        dplyr::mutate(geneID = new_geneID) %>%
        # dplyr::relocate(geneID, .after = geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()

    } else if(keyType == "SYMBOL" & !identical(old_id, id)){
      old_geneID <- replace_id(id_dat2,egmt$geneID) %>%
        replace_id(id_dat1,.)

      egmt <- egmt %>%
        dplyr::mutate(geneID = old_geneID) %>%
        dplyr::mutate(geneID_symbol = new_geneID) %>%
        dplyr::relocate(geneID_symbol, .after = geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()

    } else if(keyType == "ENTREZID"){
      egmt <- egmt %>%
        dplyr::mutate(geneID_symbol = new_geneID) %>%
        dplyr::relocate(geneID_symbol, .after = geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()

    } else {
      # input other types (including alias)
      old_geneID <- replace_id(id_dat2,egmt$geneID)

      egmt <- egmt %>%
        dplyr::mutate(geneID = old_geneID) %>%
        dplyr::mutate(geneID_symbol = new_geneID) %>%
        dplyr::relocate(geneID_symbol, .after = geneID) %>%
        calcFoldEnrich() %>%
        as.enrichdat()
    }

    ### save as list
    egmt = egmt %>% dplyr::select(-GeneRatio)
    genelist_df = data.frame(ID = names(genelist), logfc = genelist)
    exponent = data.frame(exponent = exponent)
    org = data.frame(org = ens_org)

    new_obj <- list(gsea_df = egmt, genelist = genelist_df, geneset = geneset,  exponent = exponent, org = org)
  }

  #---- IMPORT OTHER ---#
  # e.g. KEGG/DO/WikiPathways
  if(type == 'other'){
    if(is.null(object@organism)) object@organism = 'human'
    ens_org <- mapEnsOrg(object@organism)
    keg <- as.data.frame(object)
    id <- strsplit(keg$geneID,'\\/') %>% unlist() %>% unique()
    keyType <- gentype(id = id, org = ens_org)

    old_id <- id
    if (keyType == "SYMBOL") {
      id_dat1 <- suppressMessages(transId(id, "symbol", ens_org, unique = T))
      id <- id_dat1 %>% dplyr::pull(symbol)
    }

    if (keyType != "ENTREZID") {
      id_dat2 <- suppressMessages(transId(id, "entrezid", ens_org, unique = T))
    }

    ### get geneID_symbol ---#
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
    new_obj <- new_keg %>%
      dplyr::mutate(RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

  }


  return(new_obj)
}

utils::globalVariables(c(
  "genelist", "geneset", "gs_name"
))
