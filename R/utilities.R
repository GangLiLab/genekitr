#############################
### Part I: organism name & pbobe platform name
#############################
.initial <- function() {
  pos <- 1
  envir <- as.environment(pos)
  assign(".genekitrEnv", new.env(), envir = envir)
}

#---  get msigdb org ---#
msigdb_org_data <- function() {
  .initial()
  utils::data(list = "msig_org", package = "genekitr",envir = .genekitrEnv)
  get("msig_org", envir = .genekitrEnv)
}

#---  get msigdb org ---#
msigdb_category_data <- function() {
  .initial()
  utils::data(list = "msig_category", package = "genekitr",envir = .genekitrEnv)
  get("msig_category", envir = .genekitrEnv)
}

#--- bioc org name data ---#
biocOrg_name_data <- function() {
  .initial()
  utils::data(list = "biocOrg_name", package = "genekitr",envir = .genekitrEnv)
  get("biocOrg_name", envir = .genekitrEnv)
}

#---  kegg org name data ---#
keggOrg_name_data <- function() {
  .initial()
  utils::data(list = "keggOrg_name", package = "genekitr",envir = .genekitrEnv)
  get("keggOrg_name", envir = .genekitrEnv)
}

#--- ensembl org name data ---#
ensOrg_name_data <- function() {
  .initial()
  utils::data(list = "ensOrg_name", package = "genekitr",envir = .genekitrEnv)
  get("ensOrg_name", envir = .genekitrEnv)
}

#--- probe platform data ---#
hsapiens_probe_data <- function() {
  .initial()
  utils::data(list = "hsapiens_probe_platform", package = "genekitr",envir = .genekitrEnv)
  get("hsapiens_probe_platform", envir = .genekitrEnv)
}

#############################
### Part II:  match organism
#############################
#---  get organism Bioconductor shortname ---#
mapBiocOrg <- function(organism) {
  # organism <- tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hsa" | organism == "hs") organism <- "hs"
  if (organism == "mm" | organism == "mouse") organism <- "mm"
  if (organism == "rat") organism <- "rn"

  # support organisms: http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
  biocorg <- biocOrg_name_data()
  rm(biocOrg_name, envir = .genekitrEnv)

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
  if (organism == "hg" | organism == "human" | organism == "hs" | organism == "hsapiens" ) organism <- "hsa"
  if (organism == "mm" | organism == "mouse" | organism == "mmusculus") organism <- "mmu"
  if (organism == "rn" | organism == "rat") organism <- "rno"

  # some common name happens to be same with kegg short name
  is.common <- ifelse(organism %in% c("rat", "cow", "dog"), TRUE, FALSE)

  # support organisms: https://www.genome.jp/kegg/catalog/org_list.html
  kgorg <- keggOrg_name_data()
  rm(keggOrg_name, envir = .genekitrEnv)

  if ((organism %in% kgorg$short_name && !is.common)) {
    org <- organism
  } else if (nchar(organism) >= 3) {
    out <- kgorg[grepl(paste0("(\\b", organism, "\\b)"), kgorg$full_name,ignore.case = T), ]
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
    out <- kgorg[grepl(paste0("^", organism, ""), kgorg$short_name,ignore.case = T), ]
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
getKeggLatin <- function(organism) {
  organism <- tolower(organism)
  kgorg <- keggOrg_name_data()
  rm(keggOrg_name, envir = .genekitrEnv)

  latin <- kgorg[kgorg$short_name %in% organism, "full_name"] %>%
    gsub(" \\(.*\\)$", "", .)

  return(latin)
}

#---  get organism ensembl short latin name ---#
mapEnsOrg <- function(organism) {
  # organism <- tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hsa" | organism == "hs") organism <- "hsapiens"
  if (organism == "mm" | organism == "mouse") organism <- "mmusculus"
  if (organism == "rn" | organism == "rat") organism <- "rnorvegicus"
  if (organism == "dm" | organism == "fly") organism <- "dmelanogaster"
  if (organism == "dr" | organism == "zebrafish") organism <- "drerio"
  if (organism == "bt" | organism == "bovine") organism <- "btaurus"
  if (organism == "ce" | organism == "worm") organism <- "celegans"
  if (organism == "gg" | organism == "chicken") organism <- "ggallus"
  if (organism == "mmu" | organism == "rhesus") organism <- "mmulatta"
  if (organism == "pt" | organism == "chipm") organism <- "ptroglodytes"
  if (organism == "xenopus") organism <- "xtropicalis"

  # ensembl organisms: https://asia.ensembl.org/info/about/species.html
  ensorg <- ensOrg_name_data()
  rm(ensOrg_name, envir = .genekitrEnv)

  check_all <- apply(ensorg, 2, function(x) organism %in% x)

  if (any(check_all)) {
    org <- ensorg %>%
      dplyr::filter(eval(parse(text = colnames(.)[check_all])) %in% organism) %>%
      dplyr::pull(latin_short_name)
  } else {
    stop("\nCheck the latin_short_name in `genekitr::ensOrg_name`")
  }

  return(org)
}

#############################
### Part III: data query
#############################
web.url <- "https://genekitr-china.oss-accelerate.aliyuncs.com"

#--- ensembl anno data ---#
# options(geneset.download.method = "wininet")
ensAnno <- function(org, download.method = NULL, hgVersion) {
  if(hgVersion == 'v38'){
    org <- mapEnsOrg(tolower(org))
  }else{
    org <- mapEnsOrg(tolower(org))
    org <- paste0(org,'_v19')
  }

  data_dir <- tools::R_user_dir("genekitr", which = "data")
  sub_dir <- "/info/gene/"
  data_dir <- paste0(data_dir, sub_dir)
  make_dir(data_dir)

  url <- paste0(web.url, sub_dir, org, "_anno.fst")
  destfile <- paste0(data_dir, org, "_anno.fst")
  web_f_size <- check_web_size(url)
  local_f_size <- file.size(destfile)
  if(is.na(local_f_size)) local_f_size = 0

  genekitr_download(url, destfile, method = download.method,
                    data_dir, web_f_size, local_f_size)

  dat <- suppressMessages(fst::read.fst(destfile))
  invisible(dat)

}

#--- probe anno data ---#
probAnno <- function(org, download.method = NULL) {
  if (org == "hg" | org == "human" | org == "hsa" | org == "hs") org <- "hsapiens"
  if (org == "mm" | org == "mouse") org <- "mmusculus"
  if (org == "rn" | org == "rat") org <- "rnorvegicus"
  if (!org %in% c("hsapiens", "mmusculus", "rnorvegicus")) stop('We only support "human", "mouse" and "rat".')

  data_dir <- tools::R_user_dir("genekitr", which = "data")
  sub_dir <- "/info/probe/"
  data_dir <- paste0(data_dir, sub_dir)
  make_dir(data_dir)

  url <- paste0(web.url, sub_dir, org, "_probe.fst")
  destfile <- paste0(data_dir, "/", org, "_anno.fst")
  web_f_size <- check_web_size(url)
  local_f_size <- file.size(destfile)
  if(is.na(local_f_size)) local_f_size = 0

  genekitr_download(url, destfile, method = download.method,
                   data_dir, web_f_size, local_f_size)

  dat <- suppressMessages(fst::read.fst(destfile))
  invisible(dat)
}

#--- keytype order data ---#
getOrder <- function(org, keytype, download.method = NULL, hgVersion) {
  if(hgVersion == 'v38'){
    org <- mapEnsOrg(tolower(org))
  }else{
    org <- mapEnsOrg(tolower(org))
    org <- paste0(org,'_v19')
  }


  data_dir <- tools::R_user_dir("genekitr", which = "data")
  sub_dir <- "/info/gene/"
  data_dir <- paste0(data_dir, sub_dir)
  make_dir(data_dir)

  url <- paste0(web.url, sub_dir, org, "_", keytype, "_order.fst")
  destfile <- paste0(data_dir, "/", org, "_", keytype, "_order.fst")
  web_f_size <- check_web_size(url)
  local_f_size <- file.size(destfile)
  if(is.na(local_f_size)) local_f_size = 0

  genekitr_download(url, destfile, method = download.method,
                    data_dir, web_f_size, local_f_size)

  dat <- suppressMessages(fst::read.fst(destfile))
  invisible(dat)
}

#############################
### Part IV: function of geninfo
#############################
make_dir <- function(data_dir){
  if (!dir.exists(data_dir)) {
    tryCatch(
      {
        dir.create(data_dir, recursive = TRUE)
      },
      error = function(e) {
        message(paste0(
          "Seems like you cannot access dir: ", data_dir,
          "\nPlease spefify a valid dir to save data..."
        ))
        data_dir <- readline(prompt = "Enter directory: ")
        dir.create(data_dir, recursive = TRUE)
      }
    )
  }
}

#--- get web server file size ---#
check_web_size <- function(url){
  web_f_size <- RCurl::getURL(url, nobody = 1L, header = 1L) %>%
    strsplit("\r\n") %>%
    unlist() %>%
    stringr::str_extract("Content-Length.*[0-9]") %>%
    stringr::str_remove_all("Content-Length: ") %>%
    stringi::stri_remove_empty_na() %>%
    as.numeric()
  return(web_f_size)
}

#--- download function ---#
# url: web url
# destfile: local file name
# data_dir: destfile directory
# method: "auto" (default), "wininet" (for windows)
# web_f_size: file size on webserver
# local_f_size: file size in local
genekitr_download <- function(url, destfile,data_dir,
                             method, web_f_size, local_f_size){
  if (!file.exists(destfile)) {
    # message("Initializing, please wait (just once)...")
    if (is.null(method)) method <- getOption("genekitr.download.method")

    if (!is.null(method) && method != "auto") {
      tryCatch(utils::download.file(url, destfile, quiet = TRUE, method = method, mode = "wb"),
               error = function(e) {
                 message(paste0(
                   "Auto download failed...\nPlease download via: ", url,
                   "\nThen save to: ", data_dir, "\n"
                 ))
               })
    }else{
      tryCatch(utils::download.file(url, destfile, quiet = TRUE,method = 'auto', mode = "wb"),
               error = function(e) {
                 message(paste0(
                   "Auto download failed...\nPlease download manually via: ", url,
                   "\nThen save to: ", data_dir, "\n"
                 ))
               }
      )
    }
  } else if (web_f_size != local_f_size) {
    message("Detected new version data, updating...")
    if (!is.null(method) && method != "auto") {
      tryCatch(utils::download.file(url, destfile, quiet = TRUE, method = method, mode = "wb"),
               error = function(e) {
                 message(paste0(
                   "Auto download failed...\nPlease download manually via: ", url,
                   "\nThen save to: ", data_dir, "\n"
                 ))
               })
    }else{
      tryCatch(utils::download.file(url, destfile, quiet = TRUE,method = 'auto', mode = "wb"),
               error = function(e) {
                 message(paste0(
                   "Auto download failed...\nPlease download manually via: ", url,
                   "\nThen save to: ", data_dir, "\n"
                 ))
               }
      )
    }
  }
}

#---  decide gene id type ---#
gentype <- function(id, data = NULL, org, hgVersion='v38') {
  org <- mapEnsOrg(org)
  if (is.null(data)) data <- ensAnno(org,hgVersion = hgVersion)

  if ("symbol" %in% colnames(data)) {
    data_symbol_normal <- data$symbol %>% stringi::stri_remove_empty_na()
    data_symbol_lower <- tolower(data_symbol_normal)
    data_symbol_upper <- toupper(data_symbol_normal)
    n_sym <- sum(id %in% data_symbol_normal |  tolower(id) %in% data_symbol_lower | toupper(id) %in% data_symbol_upper)
  } else {
    n_sym <- 0L
  }

  if ("ncbi_alias" %in% colnames(data)) {
    data_ncbi_alias_normal <- data$ncbi_alias %>%
      stringr::str_split("; ") %>%
      unlist() %>%
      stringi::stri_remove_empty_na()
    data_ncbi_alias_lower <- tolower(data_ncbi_alias_normal)
    data_ncbi_alias_upper <- toupper(data_ncbi_alias_normal)
    n_ala <- sum(id %in% data_ncbi_alias_normal | tolower(id) %in% data_ncbi_alias_lower | toupper(id) %in% data_ncbi_alias_upper)
  } else {
    n_ala <- 0L
  }

  if ("ensembl_alias" %in% colnames(data)) {
    data_ensembl_alias_normal <- data$ensembl_alias %>%
      stringr::str_split("; ") %>%
      unlist() %>%
      stringi::stri_remove_empty_na()
    data_ensembl_alias_lower <- tolower(data_ensembl_alias_normal)
    data_ensembl_alias_upper <- toupper(data_ensembl_alias_normal)
    n_e_ala <- sum(id %in% data_ensembl_alias_normal | tolower(id) %in% data_ensembl_alias_lower | toupper(id) %in% data_ensembl_alias_upper)
  } else {
    n_e_ala <- 0L
  }

  if ("ensembl" %in% colnames(data)) {
    data_ensembl <- stringr::str_split(data$ensembl, "; ") %>%
      unlist() %>%
      stringi::stri_remove_empty_na()
    n_ens <- sum(id %in% data_ensembl)
  } else {
    n_ens <- 0L
  }

  if ("entrezid" %in% colnames(data)) {
    data_entrezid <- data$entrezid %>% stringi::stri_remove_empty_na()
    n_ent <- sum(id %in% data_entrezid)
  } else {
    n_ent <- 0L
  }

  if ("uniprot" %in% colnames(data)) {
    data_uniprot <- stringr::str_split(data$uniprot, "; ") %>%
      unlist() %>%
      stringi::stri_remove_empty_na()
    n_uni <- sum(id %in% data_uniprot)
  } else {
    n_uni <- 0L
  }

  if (sum(n_sym, n_ens, n_ent, n_uni, n_ala, n_e_ala) == 0) {
    stop("Wrong organism or input id has no match!")
  } else {
    check <- which(c(n_sym, n_ens, n_ent, n_uni, n_ala, n_e_ala) %in%
                     max(n_sym, n_ens, n_ent, n_uni, n_ala, n_e_ala))
    typ <- c("SYMBOL", "ENSEMBL", "ENTREZID", "UNIPROT", "SYMBOL","SYMBOL")[check] %>% unique()
  }

  return(typ)
}


#############################
### Part V: function of enrichment analysis
#############################

calcFoldEnrich <- function(df) {
  if (any(grepl("[gene|bg]ratio", tolower(colnames(df))))) {
    check_gr <- which(grepl(".*gene.*ratio", tolower(colnames(df))))
    check_bg <- which(grepl(".*bg*ratio", tolower(colnames(df))))
    to_calc <- paste0("(", df[, check_gr], ")/(", df[, check_bg], ")")

    if(any(grepl("foldenrich", colnames(df),ignore.case = T))){
      colnames(df)[grepl("foldenrich", colnames(df),ignore.case = T)] = 'FoldEnrich'
    }else{
      df <- df %>%
        dplyr::mutate(FoldEnrich = sapply(to_calc, function(x) eval(parse(text = x))))
    }


  }
  return(df)
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
  new_id <- stringr::str_split(id, "\\/") %>% unlist()
  check <- any(new_id%in%dat[,2])
  if(check){
    stringr::str_split(id, "\\/") %>%
      lapply(., function(x) {
        dat %>% dplyr::filter(.[[2]]%in%x) %>% dplyr::pull(1) %>%
          paste0(.,collapse = '/')
      }) %>%  do.call(rbind,.) %>% as.character()
  }else{
    stringr::str_split(id, "\\/") %>%
      lapply(., function(x) {
        dat %>% dplyr::filter(.[[1]]%in%x) %>% dplyr::pull(2) %>%
          paste0(.,collapse = '/')
      }) %>%  do.call(rbind,.) %>% as.character()
  }
}


#############################
### Part VI: add global variables
#############################
utils::globalVariables(c(
  ".", "data_dir", "biocOrg_name", "full_name", "short_name", "keggOrg_name", "item", "type", "sets",
  "count", "theme_classic", "input_id", "ensOrg_name", "latin_short_name", "ES", "pathway",
  "plotGseaTable", "pval", "scale_fill_continuous", "scale_x_discrete", "ONTOLOGY", "facet_grid",
  "BgRatio", "E", "ID", "V", "delete.edges", "enrichGenes", "geneID.y",
  "geneID_symbol", "geom_edge_link", "geom_node_text", "ggraph", "graph.data.frame",
  "guide_legend", "guides", "logfc", "melt", "new_ego", "scale_size_continuous", "E<-", "V<-",
  "method", "Term", "arrow", "circle", "geom_node_label", "geom_node_point", "go_id", "gotbl", "parent",
  "FoldEnrich", "GeneRatio", "fct_reorder", "geom_col", "scale_fill_discrete",
  "scale_size", "scale_x_continuous", "sec_axis", "everything", "gene", "coord_flip",
  "expansion", "index", "nes.group", "padj.group", "change", "label", "logFC", "stat", "pvalue",
  "cluster","Cluster", "go", "Bioc_anno", "Platform", "ensembl", "ensembl_id", "probe_id",
  "hsapiens_probe_platform", "new_x","old_id","entrezid","bioc_name",".genekitrEnv","Ontology","venn_percent",
  "input_id2"
))
