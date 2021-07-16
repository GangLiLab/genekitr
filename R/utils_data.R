#---  get msigdb org ---#
msigdb_org_data <- function() {
  utils::data(list="msig_org", package="AnnoGenes")
  get("msig_org", envir = .GlobalEnv)
}

#--- get bioc org name ---#
biocOrg_name <- function() {
  utils::data(list="biocOrg_name", package="AnnoGenes")
  get("biocOrg_name", envir = .GlobalEnv)
}

#---  get kegg org name ---#
keggOrg_name <- function() {
  utils::data(list="keggOrg_name", package="AnnoGenes")
  get("keggOrg_name", envir = .GlobalEnv)
}

#--- get ensembl anno ---#
# DEPRECATED NOW!
if(F){
  human_gtf <- function() {
    utils::data(list="data/homo_sapiens_V104_gtf.rda", package="AnnoGenes")
    get("human_gtf", envir = .GlobalEnv)
  }
  mouse_gtf <- function() {
    utils::data(list="data/mus_musculus_V104_gtf.rda", package="AnnoGenes")
    get("mouse_gtf", envir = .GlobalEnv)
  }
  rat_gtf <- function() {
    utils::data(list="data/rattus_norvegicus_V104_gtf.rda", package="AnnoGenes")
    get("rat_gtf", envir = .GlobalEnv)
  }
}


#--- get bioconductor anno ---#
biocAnno <- function(org){
  org = mapBiocOrg(tolower(org))
  utils::data(list=paste0(org,'_anno'), package="AnnoGenes")
  get(paste0(org,'_anno'), envir = .GlobalEnv)
}



