# AnnoGenes utilities for analysing

##' @title Get Msigdb database term and gene information
##' @param org organism name from `msigdb_org_data()`.
##' @param category MSigDB collection abbreviation, C1 to C8 and H.
##' @param subcategory MSigDB sub-collection abbreviation, such as REACTOME or BP.
##' @return a dataframe of 3 columns with term, entrez and symbol name.
##' @importFrom stringr str_split
##' @importFrom dplyr %>% filter pull
##' @importFrom stringi stri_remove_empty_na
##' @export
##' @examples
##' \donttest{
##' msigdb <- getMsigdb(org='human', category='C5',subcategory='GO:CC')
##' }
getMsigdb <- function(org,
                      category = c('C1','C2','C3','C4','C5','C6','C7','C8','H'),
                      subcategory=NULL,
                      ...) {

  #--- args ---#
  # org
  msig_org <- msigdb_org_data()
  all_org = c(msig_org[,1],
              stringr::str_split(msig_org[,2],', ',simplify = T) %>%
                as.character() %>%
                stringi::stri_remove_empty_na())
  if (!org %in% all_org) stop("Choose a valid organism!\n\n",paste0(all_org,' | '))

  # category
  if(! category %in% c('C1','C2','C3','C4','C5','C6','C7','C8','H')){
    stop("Category should choose from: C1, C2, C3, C4, C5, C6, C7, C8, H...")
  }else{
    category <- match.arg(category)
  }

  # subcategory
  msig_category <- msigdb_category_data()
  all_sub <- msig_category[,2] %>%
    stringi::stri_remove_empty_na()

  som_sub <- msig_category %>%
    dplyr::filter(gs_cat==category) %>%
    dplyr::pull(gs_subcat)

  if( is.null(subcategory)){
    if(som_sub == ''){
      message(paste0(category,' has no subcategory...'))
      subcategory = ''
    }else{
      stop("choose a valid subcategory for ",category,"...\n",paste0(som_sub,' | '))
    }
  }else if(! subcategory %in% all_sub){
    if(som_sub == ''){
      message(paste0(category,' has no subcategory...'))
      subcategory = ''
    }else{
      stop("choose a valid subcategory for ",category,"...\n",paste0(som_sub,' | '))
    }
  }

  #--- codes ---#
  msigdb <- msigdbr(org, category, subcategory) %>%
    dplyr::select(., c("gs_name","gene_symbol","entrez_gene")) %>%
    as.data.frame()

  return(msigdb)

}

msigdb_org_data <- function() {
  utils::data(list="msig_org", package="AnnoGenes")
  get("msig_org", envir = .GlobalEnv)
}
msigdb_category_data <- function() {
  utils::data(list="msig_category", package="AnnoGenes")
  get("msig_category", envir = .GlobalEnv)
}

biocOrg_data <- function() {
  utils::data(list="map_biocOrg", package="AnnoGenes")
  get("map_biocOrg", envir = .GlobalEnv)
}


mapBiocOrg <- function(organism) {
  # support organisms: http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
  # just extract org.db
  if (organism == "anopheles") {
    org <- "ag"
  } else if (organism == "bovine") {
    org <- "bt"
  } else if (organism == "worm") {
    org <- "ce"
  } else if (organism == "canine") {
    org <- "cf"
  } else if (organism == "fly") {
    org <- "dm"
  } else if (organism == "zebrafish") {
    org <- "dr"
  } else if (organism == "ecolik12") {
    org <- "ecK12"
  } else if (organism == "ecolisakai") {
    org <- "ecSakai"
  } else if (organism == "chicken") {
    org <- "gg"
  } else if (organism == "human") {
    org <- "hs"
  } else if (organism == "mouse") {
    org <- "mm"
  } else if (organism == "rhesus") {
    org <- "mmu"
  } else if (organism == "chipm") {
    org <- "pt"
  } else if (organism == "rat") {
    org <- "rn"
  } else if (organism == "pig") {
    org <- "ss"
  } else if (organism == "xenopus") {
    org <- "xl"
  } else {
    org <- organism
  }
  return(org)
}
