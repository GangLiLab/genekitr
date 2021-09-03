#---  get msigdb org ---#
msigdb_org_data <- function() {
  utils::data(list="msig_org", package="genekitr")
  get("msig_org", envir = .GlobalEnv)
}

#---  get msigdb org ---#
msigdb_category_data <- function() {
  utils::data(list="msig_category", package="genekitr")
  get("msig_category", envir = .GlobalEnv)
}

#--- get bioc org name ---#
biocOrg_name_data <- function() {
  utils::data(list="biocOrg_name", package="genekitr")
  get("biocOrg_name", envir = .GlobalEnv)
}

#---  get kegg org name ---#
keggOrg_name_data <- function() {
  utils::data(list="keggOrg_name", package="genekitr")
  get("keggOrg_name", envir = .GlobalEnv)
}

#--- get bioconductor anno ---#
biocAnno <- function(org){
  org = mapBiocOrg(tolower(org))
  if(!file.exists(paste0(tempdir(),'/',org,'_anno.rda'))){
    download.file(paste0("http://112.74.191.19/genekitr/",org,'_anno.rda'),
                  paste0(tempdir(),'/',org,'_anno.rda'),mode = "wb",quiet = TRUE)
  }
  load(paste0(tempdir(),'/',org,'_anno.rda'), envir = .GlobalEnv)
  get(paste0(org,'_anno'), envir = .GlobalEnv)
}






