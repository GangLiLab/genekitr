#---  get organism Bioconductor shortname ---#
mapBiocOrg <- function(organism) {
  organism <- tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hsa" | organism == "hs") organism <- "hs"
  if (organism == "mm" | organism == "mouse") organism <- "mm"

  # support organisms: http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
  biocorg <- biocOrg_name_data()
  rm(biocOrg_name, envir = .GlobalEnv)

  if (organism %in% biocorg$short_name) {
    org <- organism
  } else if (organism %in% biocorg$full_name) {
    org <- biocorg %>%
      dplyr::filter(full_name == organism) %>%
      dplyr::pull(short_name)
  } else {
    stop(
      "Check organism name! \n USE FULL NAME: ",
      paste0(biocOrg_name_data() %>% dplyr::pull(full_name), " | "),
      "\n OR USE SHORT NAME: ",
      paste0(biocOrg_name_data() %>% dplyr::pull(short_name), " | ")
    )
  }

  org <- stringr::str_to_title(org)
  return(org)
}

#---  get organism KEGG shortname ---#
mapKeggOrg <- function(organism) {
  organism <- tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hs") organism <- "hsa"
  if (organism == "mm" | organism == "mouse") organism <- "mmu"

  # some common name happens to be same with kegg short name
  is.common <- ifelse(organism %in% c("rat", "cow", "dog"), TRUE, FALSE)

  # support organisms: https://www.genome.jp/kegg/catalog/org_list.html
  kgorg <- keggOrg_name_data()
  rm(keggOrg_name, envir = .GlobalEnv)

  if ((organism %in% kgorg$short_name && !is.common)) {
    org <- organism
  } else if (nchar(organism) >= 3) {
    out <- kgorg[grepl(paste0("(\\b", organism, "\\b)"), kgorg$full_name), ]
    if (nrow(out) > 1) {
      stop(
        "\nPlease choose the SHORT NAME from below: \n",
        paste0(utils::capture.output(print.data.frame(out, row.names = FALSE)), collapse = "\n")
      )
    } else if (nrow(out) == 1) {
      org <- out$short_name
    } else {
      stop("\nCheck the organism name again!")
    }
  } else {
    out <- kgorg[grepl(paste0("^", organism, ""), kgorg$short_name), ]
    if (nrow(out) != 0) {
      stop(
        "\nPlease choose the SHORT NAME from below: \n",
        paste0(utils::capture.output(print.data.frame(out, row.names = FALSE)), collapse = "\n")
      )
    } else if (nrow(out) == 1) {
      org <- out$short_name
    } else {
      stop("\nCheck the organism name again!")
    }
  }

  return(org)
}

#--- get KEGG full latin name ---#
getKeggLatin <- function(organism){
  organism <- tolower(organism)
  kgorg <- keggOrg_name_data()
  rm(keggOrg_name, envir = .GlobalEnv)

  latin <- kgorg[kgorg$short_name %in% organism,'full_name' ] %>%
    gsub(' \\(.*\\)$','',.)

  return(latin)
}

#---  get organism ensembl short latin name ---#
mapEnsOrg <- function(organism) {
  # organism <- tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hsa" | organism == "hs") organism <- "hsapiens"
  if (organism == "mm" | organism == "mouse") organism <- "mmusculus"
  if (organism == "rn" | organism == "rat" ) organism <- "rnorvegicus"
  if (organism == "dm" | organism == "fly" ) organism <- "dmelanogaster"
  if (organism == "dr" | organism == "zebrafish" ) organism <- "drerio"
  if (organism == "bt" | organism == "bovine" ) organism <- "btaurus"
  if (organism == "ce" | organism == "worm" ) organism <- "celegans"
  if (organism == "gg" | organism == "chicken" ) organism <- "ggallus"
  if (organism == "mmu" | organism == "rhesus" ) organism <- "mmulatta"
  if (organism == "pt" | organism == "chipm" ) organism <- "ptroglodytes"
  if (organism == "xenopus" ) organism <- "xtropicalis"

  # ensembl organisms: https://asia.ensembl.org/info/about/species.html
  ensorg <- ensOrg_name_data()
  rm(ensOrg_name, envir = .GlobalEnv)

  check_all = apply(ensorg, 2, function(x) organism %in% x)

  if(any(check_all)){
    org = ensorg %>% dplyr::filter(eval(parse(text = colnames(.)[check_all])) %in% organism) %>%
      dplyr::pull(latin_short_name)
  }else{
    stop("\nCheck the latin_short_name in `genekitr::ensOrg_name`")
  }

  return(org)
}

#---  decide gene id type ---#
gentype <- function(id, data = NULL, org) {
  org <- mapEnsOrg(org)
  if(is.null(data)) data <- ensAnno(org)

  if('symbol' %in% colnames(data)){
    data_symbol <- data$symbol %>% stringi::stri_remove_empty_na()
    n_sym = sum(id %in% data_symbol)
  }else{
    n_sym = 0L
  }

  if('ensembl' %in% colnames(data)){
    data_ensembl <- stringr::str_split(data$ensembl,'; ') %>% unlist() %>% stringi::stri_remove_empty_na()
    n_ens = sum(id %in% data_ensembl)
  }else{
    n_ens = 0L
  }

  if('entrezid' %in% colnames(data)){
    data_entrezid <- data$entrezid %>% stringi::stri_remove_empty_na()
    n_ent = sum(id %in% data_entrezid)
  }else{
    n_ent = 0L
  }

  if('uniprot' %in% colnames(data)){
    data_uniprot <- stringr::str_split(data$uniprot,'; ') %>% unlist()%>% stringi::stri_remove_empty_na()
    n_uni = sum(id %in% data_uniprot)
  }else{
    n_uni = 0L
  }

  if('ncbi_alias' %in% colnames(data) ){
    data_alias <- c(data$ncbi_alias, data$ensembl_alias) %>%
      strsplit("; ") %>%
      unlist() %>%
      stringi::stri_remove_empty_na()
    n_ala = sum(id %in% data_alias)
  }else{
    n_ala = 0L
  }

  # rm(list = paste0(org, "_anno"), envir = .GlobalEnv)

  if(sum(n_sym,n_ens,n_ent,n_uni,n_ala) == 0){
    stop("Wrong organism or input id has no match!")
  }else{
    check = which(c(n_sym,n_ens,n_ent,n_uni,n_ala) %in% max(n_sym,n_ens,n_ent,n_uni,n_ala) )
    typ = c("SYMBOL","ENSEMBL","ENTREZID","UNIPROT","SYMBOL")[check] %>% unique()
  }

  return(typ)

}

#---calc fold enrichment ---#
calcFoldEnrich <- function(df) {
  if (any(grepl("[gene|bg]ratio", tolower(colnames(df))))) {
    check_gr <- which(grepl(".*gene.*ratio", tolower(colnames(df))))
    check_bg <- which(grepl(".*bg*ratio", tolower(colnames(df))))
    to_calc <- paste0("(", df[, check_gr], ")/(", df[, check_bg], ")")

    df <- df %>%
      dplyr::mutate(FoldEnrich = sapply(to_calc, function(x) eval(parse(text = x))))
  }
  return(df)
}

#---  get msigdb org ---#
msigdb_org_data <- function() {
  utils::data(list = "msig_org", package = "genekitr")
  get("msig_org", envir = .GlobalEnv)
}

#---  get msigdb org ---#
msigdb_category_data <- function() {
  utils::data(list = "msig_category", package = "genekitr")
  get("msig_category", envir = .GlobalEnv)
}

#--- bioc org name data ---#
biocOrg_name_data <- function() {
  utils::data(list = "biocOrg_name", package = "genekitr")
  get("biocOrg_name", envir = .GlobalEnv)
}

#---  kegg org name data ---#
keggOrg_name_data <- function() {
  utils::data(list = "keggOrg_name", package = "genekitr")
  get("keggOrg_name", envir = .GlobalEnv)
}

#--- ensembl org name data ---#
ensOrg_name_data <- function(){
  utils::data(list = "ensOrg_name", package = "genekitr")
  get("ensOrg_name", envir = .GlobalEnv)
}

#--- ensembl anno data ---#
ensAnno <- function(org, version) {
  if(missing(version)) version = 106
  org <- mapEnsOrg(tolower(org))
  # data_dir = rappdirs::user_data_dir(appname = 'genekitr')
  data_dir = tools::R_user_dir('genekitr',which = 'data')
  data_dir = paste0(data_dir,'/v',version)
  destfile = paste0(data_dir, "/", org, "_anno.fst")

  if(!dir.exists(data_dir)){
    tryCatch(
      {
        dir.create(data_dir,recursive = TRUE)
      },
      error = function(e) {
        message(paste0("Seems like you cannot access dir: ",data_dir,
                       '\nPlease spefify a valid dir to save data...'))
        data_dir = readline(prompt="Enter directory: ")
        dir.create(data_dir,recursive = TRUE)
      }
    )
  }

  if (!file.exists(destfile)) {
    message('We need to download some data, please wait (just once)...')
    url = paste0("http://112.74.191.19/genekitr/v",version,'/', org, "_anno.fst")
    # web_download(url, paste0(data_dir, "/", org, "_anno.fst"),  mode = "wb", quiet = TRUE)

    tryCatch(
      {
        utils::download.file(url, destfile, quiet = TRUE, mode = 'wb')
      },
      error = function(e) {
        message(paste0('Auto download failed...\nPlease download via: ',url,
                       '\nThen save to: ',data_dir))
      }
    )
  }

  dat = suppressMessages(fst::read.fst(destfile))
  invisible(dat)

#   load(paste0(data_dir, "/", org, "_anno.rda"), envir = .GlobalEnv)
#   get(paste0(org, "_anno"), envir = .GlobalEnv)
}

#--- keytype order data ---#
getOrder <- function(org,keytype,version){
  if(missing(version)) version = 106

  org <- mapEnsOrg(tolower(org))
  data_dir = tools::R_user_dir('genekitr',which = 'data')
  data_dir = paste0(data_dir,'/v',version)
  destfile = paste0(data_dir, "/", org,'_',keytype,'_order.fst')

  if (!file.exists(destfile)) {
    url = paste0("http://112.74.191.19/genekitr/v",version,'/', org,'_',keytype,'_order.fst')
    # web_download(url, paste0(data_dir, "/", org,'_',keytype,'_order.fst'),  mode = "wb", quiet = TRUE)
    tryCatch(
      {
        utils::download.file(url, destfile, quiet = TRUE, mode = 'wb')
      },
      error = function(e) {
       message(paste0('Auto download failed...\nPlease download via: ',url,
                      '\nThen save to: ',data_dir))
      }
    )

  }

  dat = suppressMessages(fst::read.fst(destfile))
  invisible(dat)
  # load(paste0(data_dir, "/",org,'_',keytype,'_order.rda'), envir = .GlobalEnv)
  # get(paste0(org,'_',keytype, "_order"), envir = .GlobalEnv)
}

web_download <- function(url, destfile, try_time = 2L, ...) {
  Sys.sleep(0.01)
  tryCatch(
    {
      if (abs(try_time - 3L) > 1) {
        message(abs(try_time - 3L),' attempt ...')
      }

      utils::download.file(url, destfile, quiet = TRUE,...)

    },
    error = function(e) {
      if (try_time == 0) {
        message("Failed after 2 attempts, please check internet connection!")
        invisible(NULL)
      } else {
        web_download(url, destfile, try_time = try_time - 1L, quiet = TRUE,...)
      }
    }
  )
}


#--- add global variables ---#
utils::globalVariables(c(
  ".", "data_dir","biocOrg_name","full_name","short_name","keggOrg_name","item","type","sets",
  "count","theme_classic","input_id","ensOrg_name","latin_short_name","ES","pathway",
  "plotGseaTable","pval","scale_fill_continuous","scale_x_discrete", "ONTOLOGY","facet_grid",
  "BgRatio", "E", "ID", "V", "delete.edges", "enrichGenes", "geneID.y",
  "geneID_symbol", "geom_edge_link", "geom_node_text", "ggraph", "graph.data.frame",
  "guide_legend", "guides", "logfc", "melt", "new_ego", "scale_size_continuous","E<-", "V<-",
  "method","Term", "arrow", "circle", "geom_node_label", "geom_node_point", "go_id", "gotbl", "parent",
  "FoldEnrich", "GeneRatio", "fct_reorder", "geom_col", "scale_fill_discrete",
  "scale_size", "scale_x_continuous", "sec_axis","everything", "gene","coord_flip",
  "expansion", "index", "nes.group", "padj.group", "change","label", "logFC","stat","pvalue"))
