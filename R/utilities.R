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

#---  get organism ensembl short latin name ---#
mapEnsOrg <- function(organism) {
  organism <- tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hsa" | organism == "hs") org <- "hsapiens"
  if (organism == "mm" | organism == "mouse") organism <- "mmusculus"
  if (organism == "rn" | organism == "rat" ) organism <- "rnorvegicus"
  if (organism == "dm" | organism == "fly" ) organism <- "dmelanogaster"
  if (organism == "dr" | organism == "zebrafish" ) organism <- "drerio"

  # ensembl organisms: https://asia.ensembl.org/info/about/species.html
  ensorg <- ensOrg_name_data()
  rm(ensOrg_name, envir = .GlobalEnv)

  check_all = apply(ensorg, 2, function(x) organism %in% x)

  if(any(check_all)){
    org = ensorg %>% dplyr::filter(eval(parse(text = colnames(.)[check_all])) %in% org) %>%
      dplyr::pull(latin_short_name)
  }else{
    stop("\nCheck the organism name with `ensOrg_name_data()`")
  }

  return(org)
}

#---  decide gene id type ---#
gentype <- function(id, org) {
  org <- mapBiocOrg(org)
  if (nchar(org) > 2) {
    org <- substr(org, 1, nchar(org) - 1)
  }
  org <- stringr::str_to_title(org)

  all <- biocAnno(org)
  all_symbol <- all$symbol %>% stringi::stri_remove_empty_na()
  all_ensembl <- stringr::str_split(all$ensembl,'; ') %>% unlist() %>% stringi::stri_remove_empty_na()
  all_entrezid <- all$entrezid %>% stringi::stri_remove_empty_na()
  all_uniprot <- stringr::str_split(all$uniprot,'; ') %>% unlist()%>% stringi::stri_remove_empty_na()
  all_alias <- c(all$ncbi_alias, all$ensembl_alias) %>%
    strsplit("; ") %>%
    unlist() %>%
    stringi::stri_remove_empty_na()

  rm(list = paste0(org, "_anno"), envir = .GlobalEnv)
  n_sym = sum(id %in% all_symbol)
  n_ens = sum(id %in% all_ensembl)
  n_ent = sum(id %in% all_entrezid)
  n_uni = sum(id %in% all_uniprot)
  n_ala = sum(id %in% all_alias)

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

#--- bioconductor anno data ---#
biocAnno <- function(org, version = NULL) {
  if(is.null(version)) version = 104
  org <- mapBiocOrg(tolower(org))
  data_dir = rappdirs::user_data_dir(appname = 'genekitr')

  if(!dir.exists(data_dir)){
    tryCatch(
      {
        dir.create(data_dir)
      },
      error = function(e) {
        message(paste0("Seems like you cannot create dir: ",data_dir,
                       '\nPlease spefify a valid dir to save data...'))
        data_dir = readline(prompt="Enter directory: ")
        dir.create(data_dir)
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

#---  auto-install packages ---#
# auto_install <- function(pkg) {
#
#   # check first time
#   ret <- suppressPackageStartupMessages(
#     sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
#   )
#   missing_pkgs <- names(ret[!ret])
#   if (length(missing_pkgs) > 0) {
#     warning("The following packages are not installed: \n",
#             paste0(sprintf("  - %s", missing_pkgs), collapse = "\n"),
#             immediate. = TRUE
#     )
#     message("\nTry installing via Bioconductor...\n")
#
#     mod <- try(suppressWarnings(BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)), silent = T)
#
#     mod <- try( suppressWarnings(utils::install.packages(missing_pkgs, quiet = TRUE, dependencies = TRUE)), silent = T)
#
#     if (isTRUE(class(mod) == "try-error")) {
#       # check again
#       ret <- suppressPackageStartupMessages(
#         sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
#       )
#       missing_pkgs <- names(ret[!ret])
#       if (length(missing_pkgs) > 0) {
#         message("Try installing via CRAN...\n")
#         suppressWarnings(utils::install.packages(missing_pkgs, quiet = TRUE, dependencies = TRUE))
#
#         # 第三次检查
#         ret <- suppressPackageStartupMessages(
#           sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
#         )
#         missing_pkgs <- names(ret[!ret])
#         if (length(missing_pkgs) > 0) {
#           stop(
#             "Maybe you should check the package name ",
#             paste(missing_pkgs, collapse = ", "),
#             " or try devtools::install_github()"
#           )
#         }
#       }
#     }
#   } else {
#     message(sapply(pkg, function(x) paste0("The package ", x, " exist...\n")))
#   }
# }


#--- add global variables ---#
utils::globalVariables(c(
  ".", "biocOrg_name","full_name","short_name","keggOrg_name","item","type","sets",
  "count","theme_classic","input_id"))
