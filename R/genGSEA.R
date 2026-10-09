#' Gene Set Enrichment Analysis
#'
#' @param genelist Pre-ranked genelist with decreasing order, gene can be
#'   entrez, ensembl or symbol. If several ids are mapped to the same gene after
#'   id conversion (e.g. alias symbols), only the first (highest ranked) one is kept.
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
#' \donttest{
#' if(requireNamespace("geneset",quietly = TRUE)){
#' # only gene ids
#' data(geneList, package = "genekitr")
#' gs <- geneset::getGO(org = "human",ont = "mf",data_dir = tempdir())
#' gse <- genGSEA(genelist = geneList, geneset = gs)
#'   }
#'}

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
  ens_org <- NULL

  tryCatch(
    {
      ens_org <- tryCatch( mapEnsOrg(org), error = function(e) NULL)
      if(is.null(ens_org)) stop('Please make sure the gene types of genelist and geneset are the same.')
      keyType <- gentype(id = id, org = ens_org)

      #--- initialize ---#
      # input id must be symbol or entrezid for c gene sets

      if(transToSym){
        id_dat <- suppressMessages(transId(id, "symbol", ens_org, unique = T))
        genelist <- genelist[names(genelist)%in%id_dat$input_id]
        names(genelist) <- id_dat$symbol[match(names(genelist), id_dat$input_id)]
      }else if(keyType != "ENTREZID"){
        id_dat <- suppressMessages(transId(id, "entrezid", ens_org, unique = T))
        genelist <- genelist[names(genelist)%in%id_dat$input_id]
        names(genelist) <- id_dat$entrezid[match(names(genelist), id_dat$input_id)]
      }else if(keyType == "ENTREZID"){
        id_dat <- suppressMessages(transId(id, "symbol", ens_org, unique = T)) %>%
          dplyr::relocate(input_id,.after = symbol)
        genelist <- genelist[names(genelist)%in%id_dat$input_id]
        # id <- id_dat$input_id
      }

    },
    error = function(e) {
      ens_org <- org
      transToSym <- FALSE

    }
  )

  #--- remove duplicated genes after id conversion ---#
  # e.g. alias symbols are mapped to the same entrez id, GSEA does not allow duplicated gene names;
  # genelist is pre-ranked, so the first (highest ranked) one is kept
  genelist <- genelist[!is.na(names(genelist))]
  if(any(duplicated(names(genelist)))){
    message(sum(duplicated(names(genelist))),
            " gene ids are mapped to the same gene after id conversion, only the first one of each is kept...")
    genelist <- genelist[!duplicated(names(genelist))]
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
      seed = set_seed
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
    if (nrow(fcs) == 0) stop("No terms enriched after qvalue filtering ...")
  }

  ## gene ids are kept in the same type as the gene set (e.g. entrez for KEGG/MSigDB, symbol for GO);
  ## for entrez-based gene sets, add "geneID_symbol" for easier interpretation and gene plots in plotEnrich()
  if(!transToSym && !is.null(ens_org)){
    new_fcs <- tryCatch(
      fcs %>%
        dplyr::mutate(geneID_symbol = get_symbol(geneID, ens_org)) %>%
        dplyr::relocate(geneID_symbol, .after = geneID),
      error = function(e) fcs
    )
  }else{
    new_fcs <- fcs
  }

  ## modify id column name for GO (same style as genORA, e.g. "Hs_BP_ID", so that simGO() and GO plots work)
  if(genesetType %in% c('bp','cc','mf')){
    bioc_org <- ensOrg_name_data() %>%
      dplyr::filter(tolower(latin_short_name) %in% tolower(org)) %>%
      dplyr::pull(bioc_name) %>%
      stringr::str_to_sentence()
    if(length(bioc_org) == 1 && !is.na(bioc_org)){
      colnames(new_fcs)[1] = paste0(bioc_org,'_',toupper(genesetType),'_ID')
    }else{
      colnames(new_fcs)[1] = paste0(org,'_',toupper(genesetType),'_ID')
    }
  }

  geneAll <- sapply(seq_len(nrow(new_fcs)),function(x){
    gset_nm <- new_fcs[[1]][x]
    # gset_gene <- geneset$geneset %>% filter(.[[1]] == gset_nm) %>% pull(2)
    gset_gene <- geneset$geneset[geneset$geneset[,1] %in% gset_nm,2]
    ov_genes <- intersect(gset_gene, names(genelist)) %>% paste0(.,collapse = '/')
    ov_genes_len <- length(ov_genes)
    return(ov_genes)
  })

  geneAll_len <- sapply(seq_along(geneAll), function(x){
    gl <- geneAll[x] %>% strsplit('\\/') %>% unlist() %>% length()
  })

  new_fcs <- new_fcs %>% dplyr::rename(core_enriched_geneID = geneID, core_enriched_count = Count) %>%
    mutate(all_enriched_geneID = geneAll, all_enriched_count = geneAll_len)

  ## save as list
  genelist_df = data.frame(ID = names(genelist), logfc = genelist)
  exponent = data.frame(exponent = exponent)
  org = data.frame(org = org)
  new_geneset <- geneset$geneset

  res <- list(gsea_df = new_fcs, genelist = genelist_df, geneset = new_geneset,  exponent = exponent, org = org)


  return(res)

}



