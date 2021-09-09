#--- add global variables ---#
utils::globalVariables(c(".", ":=", "Count", "download.file","Description", "V1", "chr",
                         "count", "day", "doi", "element_line", "element_rect",
                         "element_text", "end", "ensembl_id", "entrez_gene",
                         "entrezid", "full_name", "gene","gene_id", "gene_symbol",
                         "gs_cat", "gs_name", "gs_subcat", "input_id",
                         "install.packages", "item", "journal", "labs", "margin",
                         "month", "msig_category","msig_org", "na.omit", "pmid",
                         "setSize", "sets", "short_name", "start", "strand", "symbol",
                         "theme_bw", "title", "type", "uniprot", "unit", "width", "xlab", "year",
                         'createWorkbook','saveWorkbook','biocOrg_name', 'keggOrg_name'))

#--- NCBI entrez ---#
showNCBI <- function(db = "pubmed") {
  # suppress binding notes
  fields <- rentrez::entrez_db_searchable(db)
  res <- as.data.frame(fields)[1:3]

  if (nrow(res) == 0) { # nocov start
    message("Something is wrong in your input,
            NULL will be returned, please check.")
    return(NULL)
  } # nocov end
  return(res)
}

#---  get msigdb data ---#
getMsigdb <- function(org,
                      category = c('C1','C2','C3','C4','C5','C6','C7','C8','H'),
                      subcategory=NULL) {

  #--- args ---#
  if (!requireNamespace('msigdbr', quietly = TRUE)) auto_install('msigdbr')
  org = tolower(org)
  if (org == "hg" | org == "hsa" |  org == "hs" | org == 'homo sapiens') org = 'human'
  if (org == "mm" | org == "mmu") org = 'mouse'

  # org
  msigOrg <- msigdb_org_data(); rm(msig_org,envir = .GlobalEnv)
  all_org = c(msigOrg[,1],
              stringr::str_split(msigOrg[,2],', ',simplify = T) %>%
                as.character() %>%
                stringi::stri_remove_empty_na())
  if (!org %in% tolower(all_org)) stop("Choose a valid organism!\n\n",paste0(all_org,' | '))

  # category
  if(! category %in% c('C1','C2','C3','C4','C5','C6','C7','C8','H')){
    stop("Category should choose from: C1, C2, C3, C4, C5, C6, C7, C8, H...")
  }else{
    category <- match.arg(category)
  }

  # subcategory
  msigCategory <- msigdb_category_data() ; rm(msig_category,envir = .GlobalEnv)
  all_sub <- msigCategory[,2] %>%
    stringi::stri_remove_empty_na()

  som_sub <- msigCategory %>%
    dplyr::filter(gs_cat==category) %>%
    dplyr::pull(gs_subcat)

  if( is.null(subcategory) ){
    if(som_sub == ''){
      message(paste0(category,' has no subcategory, continue...'))
      subcategory = ''
    }else{
      stop("choose a valid subcategory for ",category,"...\n",paste0(som_sub,' | '))
    }
  }else if(! subcategory %in% som_sub ){
    stop("choose a valid subcategory for ",category,"...\n",paste0(som_sub,' | '))
  }


  #--- codes ---#
  msigdb <- msigdbr::msigdbr(org, category, subcategory) %>%
    dplyr::select(., c("gs_name","gene_symbol","entrez_gene")) %>%
    as.data.frame()

  return(msigdb)

}

#---  map bioc org fullname to shortname ---#
mapBiocOrg <- function(organism) {
  organism = tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hsa" |  organism == "hs") organism = 'hs'
  if (organism == "mm" | organism == "mouse") organism = 'mm'

  # support organisms: http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
  biocorg = biocOrg_name_data(); rm(biocOrg_name,envir = .GlobalEnv)

  if( organism %in% biocorg$short_name ){
    org = organism
  }else if(organism %in% biocorg$full_name){
    org = biocorg %>% dplyr::filter(full_name == organism) %>% dplyr::pull(short_name)
  }else{
    stop('Check organism name! \n USE FULL NAME: ',
         paste0(biocOrg_name_data() %>% dplyr::pull(full_name),' | '),
         '\n OR USE SHORT NAME: ',
         paste0(biocOrg_name_data() %>% dplyr::pull(short_name),' | '))
  }

  org <- stringr::str_to_title(org)
  return(org)
}

#--- map kegg org fullname to shortname ---#
mapKeggOrg <- function(organism){
  organism = tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hs") organism = 'hsa'
  if (organism == "mm" | organism == "mouse") organism = 'mmu'

  # some common name happens to be same with kegg short name
  is.common = ifelse(organism %in% c("rat","cow","dog"), TRUE, FALSE)

  # support organisms: https://www.genome.jp/kegg/catalog/org_list.html
  kgorg = keggOrg_name_data() ;rm(keggOrg_name,envir = .GlobalEnv)

  if( (organism %in% kgorg$short_name && !is.common ) ){
    org = organism
  }else if ( nchar(organism) >= 3 ){
    out = kgorg[grepl(paste0('(\\b',organism,'\\b)'), kgorg$full_name),]
    if(nrow(out) > 1 ){
      stop('\nPlease choose the SHORT NAME from below: \n',
           paste0(utils::capture.output(print.data.frame(out,row.names = FALSE)), collapse = '\n') )
    }else if ( nrow(out) == 1){
      org = out$short_name
    }else{
      stop('\nCheck the organism name again!')
    }
  }else{
    out = kgorg[grepl(paste0('^',organism,''), kgorg$short_name),]
    if( nrow(out) != 0){
      stop('\nPlease choose the SHORT NAME from below: \n',
           paste0(utils::capture.output(print.data.frame(out,row.names = FALSE)), collapse = '\n') )
    } else if ( nrow(out) == 1){
      org = out$short_name
    } else{
      stop('\nCheck the organism name again!')
    }
  }

  return(org)
}

#--- define gene type: entrezid, ensembl or symbol ---#
.gentype <- function(id, org){
  org = mapBiocOrg(org)
  if(nchar(org) > 2){
    org = substr(org,1,nchar(org)-1)
  }
  org <- stringr::str_to_title(org)

  all = biocAnno(org)
  all_symbol = all$symbol %>% stringi::stri_remove_empty_na()
  all_ensembl = all$ensembl %>% stringi::stri_remove_empty_na()
  all_entrezid = all$entrezid %>% stringi::stri_remove_empty_na()
  all_uniprot = all$uniprot %>% strsplit('; ') %>% unlist() %>% stringi::stri_remove_empty_na()

  rm(list = paste0(org, "_anno"), envir = .GlobalEnv)
  if (any(id %in% all_symbol)) {
    c("SYMBOL")
  } else if(any(id %in% all_ensembl)){
    c("ENSEMBL")
  }else if (any(id %in% all_entrezid)){
    c("ENTREZID")
  }else if (any(id %in% all_uniprot)){
    c("UNIPROT")
  }else{
    stop('Wrong organism or input id has no match!')
  }

}

#---  auto-install packages ---#
auto_install <- function(pkg){

  # check first time
  ret <- suppressPackageStartupMessages(
    sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
  )
  missing_pkgs <- names(ret[!ret])
  if (length(missing_pkgs) > 0) {
    warning("The following packages are not installed: \n",
            paste0(sprintf("  - %s", missing_pkgs), collapse = "\n"),
            immediate. = TRUE
    )
    message("\nTry installing via Bioconductor...\n")

    mod = try(suppressMessages(BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)),silent = T)

    if(isTRUE(class(mod)=="try-error")) {
      # check again
      ret <- suppressPackageStartupMessages(
        sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
      )
      missing_pkgs <- names(ret[!ret])
      if (length(missing_pkgs) > 0) {
        message("Try installing via CRAN...\n")
        suppressWarnings(utils::install.packages(missing_pkgs, quiet = TRUE, dependencies = TRUE))

        # 第三次检查
        ret <- suppressPackageStartupMessages(
          sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
        )
        missing_pkgs <- names(ret[!ret])
        if (length(missing_pkgs) > 0) {
          stop("Maybe you should check the package name ",
               paste(missing_pkgs,collapse = ', '),
               " or try devtools::install_github()")
        }
      }
    }
  }else{
    message(sapply(pkg, function(x) paste0('The package ',x, ' exist...\n')))
  }
}

#---calc fold enrichment ---#
calcFoldEnrich <- function(df){
  if( any(grepl('[gene|bg]ratio',tolower(colnames(df)))) ){
    check_gr = which(grepl('.*gene.*ratio',tolower(colnames(df))))
    check_bg = which(grepl('.*bg*ratio',tolower(colnames(df))))
    to_calc =  paste0('(',df[,check_gr],')/(',df[,check_bg],')')

    df <- df %>%
      dplyr::mutate(FoldEnrich = sapply(to_calc, function(x) eval(parse(text = x)) ))
  }
  return(df)

}



