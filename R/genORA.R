#' Gene ORA enrichment analysis for other data sets
#'
#' @param id A vector of gene id which can be entrezid, ensembl, symbol or uniprot.
#' @param type ORA analyse type, choose one from "do" (Disease Ontology), "ncg" (Network of Cancer Gene),
#' "dgn" (DisGeNET) and "wiki" (WikiPathways)
#' @param org human as default
#' @param group_list A list of gene id groups, default is NULL.
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
#' @importFrom dplyr pull filter arrange mutate relocate
#' @importFrom stringr str_split
#' @importFrom DOSE enrichDO enrichNCG enrichDGN
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
#' # Disease ontology (DO)
#' edo <- genORA(id, type = "do")
#' # Network of Cancer Gene (NCG)
#' encg <- genORA(names(geneList), type = "ncg")
#' # DisGeNET (DGN)
#' edgn <- genORA(id, type = "dgn", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
#' # WikiPathways (wiki)
#' ewiki <- genORA(id, type = "wiki", org = "human")
#'
#' # gene id with groups
#' id <- c(head(names(geneList), 100), tail(names(geneList), 100))
#' group <- list(
#'   group1 = c(rep("up", 100), rep("down", 100)),
#'   group2 = c(rep("A", 130), rep("B", 70))
#' )
#' gedo <- genORA(id, type = "do",group_list = group)
#'
#' }

genORA <- function(id,
                   type = c('do','ncg','dgn','wiki'),
                   org = 'human',
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   minGSSize = 10,
                   maxGSSize = 500,
                   group_list = NULL,
                   universe){

  #--- args ---#
  type <- match.arg(type)
  id <- as.character(id)
  if (missing(universe)) universe <- NULL

  # different method has own supported species
  org <- tolower(org)
  if(type %in% c('do','ncg','dgn')){
    if (!org %in% c('hs','human','hsa','hsapiens'))  stop('DO enrichment analysis only supports org = "human"')
  }

  if(type == 'wiki'){
    # rWikiPathways::listOrganisms()
    all_org <- c("Arabidopsis thaliana", "Bos taurus", "Caenorhabditis elegans", "Canis familiaris",
                 "Danio rerio","Drosophila melanogaster", "Equus caballus",
                 "Gallus gallus", "Homo sapiens","Mus musculus", "Pan troglodytes", "Populus trichocarpa",
                 "Rattus norvegicus", "Saccharomyces cerevisiae", "Solanum lycopersicum",
                 "Sus scrofa")
    all_org = all_org[all_org%in% genekitr::ensOrg_name$latin_full_name]
    if (org == "hg" | org == "human" | org == "hsa" | org == "hs") org <- "Homo sapiens"
    if (org == "mm" | org == "mouse") org <- "Mus musculus"
    if (org == "rn" | org == "rat") org <- "Rattus norvegicus"
    if (org == "dm" | org == "fly") org <- "Drosophila melanogaster"
    if (org == "dr" | org == "zebrafish") org <- "Danio rerio"

    if (!org %in% all_org)  stop(paste0('WikiPathways enrichment analysis supports organism: \n',paste(sort(all_org),collapse = ' | ')))
  }


  ens_org <- mapEnsOrg(org)
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

  #--- HUMAN: DO analyse ---#
  if(type == 'do'){
    ## NO GROUP INFO
    if (is.null(group_list)) {
      ora <- suppressMessages(
        DOSE::enrichDO(
          gene = entrez_id, ont = "DO",
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          qvalueCutoff = qvalueCutoff,
          universe = universe,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize
        ))
    }else {
      ## WITH GROUP INFO
      df <- as.data.frame(group_list) %>% dplyr::mutate(id = entrez_id)
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
          suppressMessages(DOSE::enrichDO(
            gene = x, ont = "DO",
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod,
            qvalueCutoff = qvalueCutoff,
            universe = universe,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize
          )) %>%
            as.data.frame()) %>%
        do.call(rbind,.) %>%
        dplyr::mutate(Cluster = gsub("\\.[^\\.]*$", "", rownames(.), perl=TRUE)) %>%
        dplyr::relocate(Cluster,.before = dplyr::everything()) %>%
        `rownames<-`(seq_len(nrow(.)))
    }
  }

  #--- HUMAN: NCG (Network of Cancer Gene) analyse ---#
  if(type == 'ncg'){
    ## NO GROUP INFO
    if (is.null(group_list)) {
      ora <- suppressMessages(
        DOSE::enrichNCG(
          gene = entrez_id,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          qvalueCutoff = qvalueCutoff,
          universe = universe,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize
        ))
    }else {
      ## WITH GROUP INFO
      df <- as.data.frame(group_list) %>% dplyr::mutate(id = entrez_id)
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
            DOSE::enrichNCG(
              gene = x,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              universe = universe,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize
            )) %>%
            as.data.frame()) %>%
        do.call(rbind,.) %>%
        dplyr::mutate(Cluster = gsub("\\.[^\\.]*$", "", rownames(.), perl=TRUE)) %>%
        dplyr::relocate(Cluster,.before = dplyr::everything()) %>%
        `rownames<-`(seq_len(nrow(.)))
    }
  }

  #--- HUMAN: DGN (DisGeNET) analyse ---#
  if(type == 'dgn'){
    ## NO GROUP INFO
    if (is.null(group_list)) {
      ora <- suppressMessages(
        DOSE::enrichDGN(
          gene = entrez_id,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          qvalueCutoff = qvalueCutoff,
          universe = universe,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize
        ))
    }else {
      ## WITH GROUP INFO
      df <- as.data.frame(group_list) %>% dplyr::mutate(id = entrez_id)
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
            DOSE::enrichDGN(
              gene = x,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              universe = universe,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize
            )) %>%
            as.data.frame()) %>%
        do.call(rbind,.) %>%
        dplyr::mutate(Cluster = gsub("\\.[^\\.]*$", "", rownames(.), perl=TRUE)) %>%
        dplyr::relocate(Cluster,.before = dplyr::everything()) %>%
        `rownames<-`(seq_len(nrow(.)))
    }
  }

  #--- 16 orgs: WikiPathways analyse ---#
  if(type == 'wiki'){
    wiki_all <-  getWiki(org)
    geneset <- wiki_all %>% dplyr::select(wpid,gene)
    setname <- wiki_all %>%
      dplyr::select(wpid,name) %>%
      dplyr::distinct() %>%
      split(.$wpid) %>%
      lapply(function(x) x %>% dplyr::pull(name))


    ## NO GROUP INFO
    if (is.null(group_list)) {
      ora <- suppressMessages(
        clusterProfiler::enricher(
          gene = entrez_id,
          TERM2GENE = geneset,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          qvalueCutoff = qvalueCutoff,
          universe = universe,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize
        ))
    }else {
      ## WITH GROUP INFO
      df <- as.data.frame(group_list) %>% dplyr::mutate(id = entrez_id)
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
              TERM2GENE = geneset,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              universe = universe,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize
            )) %>%
            as.data.frame()) %>%
        do.call(rbind,.) %>%
        dplyr::mutate(Cluster = gsub("\\.[^\\.]*$", "", rownames(.), perl=TRUE)) %>%
        dplyr::relocate(Cluster,.before = dplyr::everything()) %>%
        `rownames<-`(seq_len(nrow(.)))
    }
    if (nrow(as.data.frame(ora)) == 0) {
      stop("No terms enriched ...")
    }else{
      ora <- as.data.frame(ora)
      ora$Description <- setname[ora$Description] %>% unlist()
    }
  }

  ################################################
  ### Next, all results are processed equally...
  if (nrow(as.data.frame(ora)) == 0) {
    stop("No terms enriched ...")
  }else{
    ora = as.data.frame(ora)
  }

  #--- get geneID_symbol ---#
  # new id
  new_geneID <- get_symbol(ora$geneID,ens_org)

  # input SYMBOL and no alias
  if( keyType == "SYMBOL" & identical(old_id, id) ){
    new_ora <- ora %>%
      dplyr::mutate(geneID = new_geneID) %>%
      # dplyr::relocate(geneID, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else if(keyType == "SYMBOL" & !identical(old_id, id)){
    old_geneID <- replace_id(id_dat2,ora$geneID) %>%
      replace_id(id_dat1,.)

    new_ora <- ora %>%
      dplyr::mutate(geneID = old_geneID) %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else if(keyType == "ENTREZID"){
    new_ora <- ora %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()

  } else {
    # input other types (including alias)
    old_geneID <- replace_id(id_dat2,ora$geneID)

    new_ora <- ora %>%
      dplyr::mutate(geneID = old_geneID) %>%
      dplyr::mutate(geneID_symbol = new_geneID) %>%
      dplyr::relocate(geneID_symbol, .after = geneID) %>%
      calcFoldEnrich() %>%
      as.enrichdat()
  }


  #--- add rich factor ---#
  new_ora <- new_ora %>%
    dplyr::mutate(RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

  return(new_ora)

}

utils::globalVariables(c("name","wpid"))


