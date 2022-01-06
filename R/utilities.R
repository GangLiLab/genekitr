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
    stop("\nCheck the latin_short_name in `ensOrg_name_data()`")
  }

  return(org)
}

#---  decide gene id type ---#
gentype <- function(id, org) {
  org <- mapEnsOrg(org)
  all <- ensAnno(org)
  if('symbol' %in% colnames(all)){
    all_symbol <- all$symbol %>% stringi::stri_remove_empty_na()
    n_sym = sum(id %in% all_symbol)
  }else{
    n_sym = 0L
  }

  if('ensembl' %in% colnames(all)){
    all_ensembl <- stringr::str_split(all$ensembl,'; ') %>% unlist() %>% stringi::stri_remove_empty_na()
    n_ens = sum(id %in% all_ensembl)
  }else{
    n_ens = 0L
  }

  if('entrezid' %in% colnames(all)){
    all_entrezid <- all$entrezid %>% stringi::stri_remove_empty_na()
    n_ent = sum(id %in% all_entrezid)
  }else{
    n_ent = 0L
  }

  if('uniprot' %in% colnames(all)){
    all_uniprot <- stringr::str_split(all$uniprot,'; ') %>% unlist()%>% stringi::stri_remove_empty_na()
    n_uni = sum(id %in% all_uniprot)
  }else{
    n_uni = 0L
  }

  if('ncbi_alias' %in% colnames(all) ){
    all_alias <- c(all$ncbi_alias, all$ensembl_alias) %>%
      strsplit("; ") %>%
      unlist() %>%
      stringi::stri_remove_empty_na()
    n_ala = sum(id %in% all_alias)
  }else{
    n_ala = 0L
  }

  rm(list = paste0(org, "_anno"), envir = .GlobalEnv)

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
ensAnno <- function(org, version = NULL) {
  if(is.null(version)) version = 104
  org <- mapEnsOrg(tolower(org))
  # data_dir = rappdirs::user_data_dir(appname = 'genekitr')
  data_dir = tools::R_user_dir('genekitr',which = 'data')

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

  if (!file.exists(paste0(data_dir, "/", org, "_anno.rda"))) {
    url = paste0("http://112.74.191.19/genekitr/v",version,'/', org, "_anno.rda")
    web_download(url, paste0(data_dir, "/", org, "_anno.rda"),  mode = "wb", quiet = TRUE)
  }

  load(paste0(data_dir, "/", org, "_anno.rda"), envir = .GlobalEnv)
  get(paste0(org, "_anno"), envir = .GlobalEnv)
}

web_download <- function(url, destfile, try_time = 2L, ...) {
  Sys.sleep(0.01)
  tryCatch(
    {
      if (abs(try_time - 3L) > 1) {
        message(abs(try_time - 3L),' attempt ...')
      }
      utils::download.file(url, destfile, ...)
    },
    error = function(e) {
      if (try_time == 0) {
        message("Failed after 2 attempts, please check internet connection!")
        invisible(NULL)
      } else {
        web_download(url, destfile, try_time = try_time - 1L, ...)
      }
    }
  )
}

#--- add global variables ---#
utils::globalVariables(c(
  ".", "data_dir","biocOrg_name","full_name","short_name","keggOrg_name","item","type","sets",
  "count","theme_classic","input_id","ensOrg_name","latin_short_name","ES","pathway",
  "plotGseaTable","pval","scale_fill_continuous","scale_x_discrete", "ONTOLOGY","facet_grid"))
