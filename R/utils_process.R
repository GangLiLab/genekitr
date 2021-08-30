# add global variables
utils::globalVariables(c(".", ":=", "Count", "download.file","Description", "V1", "chr",
                         "count", "day", "doi", "element_line", "element_rect",
                         "element_text", "end", "ensembl_id", "entrez_gene",
                         "entrezid", "full_name", "gene","gene_id", "gene_symbol",
                         "gs_cat", "gs_name", "gs_subcat", "input_id",
                         "install.packages", "item", "journal", "labs", "margin",
                         "month", "msig_category","msig_org", "na.omit", "pmid",
                         "setSize", "sets", "short_name", "start", "strand", "symbol",
                         "theme_bw", "title", "type", "uniprot", "unit", "width", "xlab", "year",
                         'createWorkbook','saveWorkbook'))

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


#---   define gene type: entrezid, ensembl or symbol ---#
.gentype <- function(id, org){
  org = mapBiocOrg(org)
  if(nchar(org) > 2){
    org = substr(org,1,nchar(org)-1)
  }
  org <- stringr::str_to_title(org)
  suppressPackageStartupMessages(require(paste0("org.", org, ".eg.db"), character.only = TRUE))
  orgSymbol <- AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egSYMBOL"))))
  orgENSEMBL <- AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egENSEMBL"))))
  if (any(id %in% orgSymbol$symbol)) {
    c("SYMBOL")
  } else if(any(id %in% orgENSEMBL$ensembl_id)){
    c("ENSEMBL")
  }else if (any(id %in% orgENSEMBL$gene_id)){
    c("ENTREZID")
  }else{
    stop('Wrong organism!')
  }
}

#---  gene id in this org or not (return a logical) ---#
# .genInorg <- function(id, org){
#   if(nchar(org) > 2){
#     org = substr(org,1,nchar(org)-1)
#   }
#   org <- stringr::str_to_title(org)
#   suppressPackageStartupMessages(require(paste0("org.", org, ".eg.db"), character.only = TRUE))
#   keytype = .gentype(id,org)
#
#   orgSymbol <- AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egSYMBOL"))))
#   if(keytype == "ENTREZID"){
#     ifelse(any(id%in%orgSymbol$gene_id),TRUE,FALSE)
#   }else if(keytype == "SYMBOL"){
#     ifelse(any(id%in%orgSymbol$symbol),TRUE,FALSE)
#   }else{
#     orgENSEMBL <- AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egENSEMBL"))))
#     ifelse(any(id%in%orgENSEMBL$ensembl_id),TRUE,FALSE)
#   }
# }

#---  auto-install packages ---#
auto_install <- function(pkg){
  options(warn=-1)

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

#--- load org.db ---#
.load_orgdb <- function(org){
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)
  org = mapBiocOrg(tolower(org))
  pkg=paste0("org.", org, ".eg.db")
  if (!requireNamespace(pkg, quietly = TRUE)) auto_install(pkg)
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))
}

#--- get ensembl annotation ---#
# current ensembl version: 104
# DEPRECATED NOW!
if(F){
  .get_ensembl_anno <- function(org, ensembl_version = '104', path = 'data'){
    #--- args ---#

    # For now supports common research species
    if (org == "hg" | org == "human" | org == "hs" | org == "hsa") organism = 'homo_sapiens'
    if (org == "mm" | org == "mouse" ) organism = 'mus_musculus'
    if (org == "rn" | org == "rat" ) organism = 'rattus_norvegicus'
    # if (org == "dm" | org == "fly" ) organism = 'drosophila_melanogaster'
    # if (org == "dre"| org == "dr" | org == "zebrafish" ) organism = 'danio_rerio'

    if(! organism %in% c('homo_sapiens','mus_musculus','rattus_norvegicus','drosophila_melanogaster','danio_rerio')){
      stop("For now we support species from:\n homo_sapiens | mus_musculus | rattus_norvegicus ")
    }

    rda_file = paste0(path,'/',organism,'_V',ensembl_version,'_gtf.rda')
    gtf_fle = list.files(path, pattern = stringr::str_to_title(organism), full.names = T)

    command = paste0("wget -c -r -nd -np -R 'index.html*' -A '",
                     stringr::str_to_title(organism),".*.",ensembl_version,".gtf.gz'",
                     " ftp://ftp.ensembl.org/pub/current_gtf/", organism,"/")

    #--- codes ---#
    if(!file.exists(rda_file)){
      system(command)
      dat = rtracklayer::import(gtf_fle) %>%
        as.data.frame() %>%
        dplyr::filter(type == 'gene')%>%
        dplyr::select(-c(15:ncol(.),ends_with('source'),'score','phase','gene_version')) %>%
        dplyr::rename(ensembl = gene_id) %>%
        dplyr::rename(symbol = gene_name) %>%
        dplyr::rename(chr = seqnames) %>%
        dplyr::relocate(symbol, ensembl,.before =  everything())

      entrz = transId(id = dat$ensembl, trans_to = 'entrez',org, return_dat = T)
      new_dat = merge(entrz, dat, by = 'ensembl', all.y = T)%>%
        dplyr::relocate(entrezid,.before = everything()) %>%
        dplyr::arrange(entrezid)

      assign(paste0(org,'_gtf'), new_dat)
      save(list=paste0(org,'_gtf'), file=rda_file)

    }
  }

}

#--- get biomart alias ---#
# alias info processing is complex than other info
.get_biomrt_alias <- function(org){
  org = tolower(org)
  if (org == "hg" | org == "human" | org == "hsa" |  org == "hs") organism = 'hsapiens'
  if (org == "mm" | org == "mouse" ) organism = 'mmusculus'
  if (org == "rn" | org == "rat" ) organism = 'rnorvegicus'

  # Here, we created dataframe of gene symbol with alias
  biomart_alias = biomaRt::getBM( attributes = c("external_gene_name",'external_synonym','uniprot_gn_symbol'),
                         mart = biomaRt::useMart("ensembl",
                                        dataset = paste0(organism,"_gene_ensembl"),
                                        host = "asia.ensembl.org")) %>%
    split(., .$external_gene_name) %>%
    lapply(., function(x) {
      paste0(x[, 2],collapse = "; ")
    })  %>%
    do.call(rbind,.) %>%
    as.data.frame() %>%
    dplyr::mutate(symbol = rownames(.)) %>%
    dplyr::rename(ensembl_alias = V1) %>%
    dplyr::select(symbol, ensembl_alias) %>%
    within(., ensembl_alias <- lapply(stringr::str_split(.[,'ensembl_alias'],'; '),
                                    function(x){
                                      unique(x) %>% paste0(.,collapse = '; ')}) %>% unlist())

  return(biomart_alias)

}

# support add-in via changing attributes (exp. feature)
.get_biomrt_feature <- function(org){
  org = tolower(org)
  if (org == "hg" | org == "human" | org == "hsa" |  org == "hs") organism = 'hsapiens'
  if (org == "mm" | org == "mouse" ) organism = 'mmusculus'
  if (org == "rn" | org == "rat" ) organism = 'rnorvegicus'

  biomart_other = biomaRt::getBM(attributes = c("ensembl_gene_id",'chromosome_name',
                                                'start_position','end_position','strand',
                                                'percentage_gene_gc_content','gene_biotype','transcript_count'),
                         mart = biomaRt::useMart("ensembl",
                                        dataset = paste0(organism,"_gene_ensembl"),
                                        host = "asia.ensembl.org")) %>%
    data.table::setnames(., old =colnames(.),
                         new = c('ensembl','chr','start','end','strand','gc_content',
                                 'gene_biotype','transcript_count')) %>%
    dplyr::mutate(width = (end - start + 1)) %>%
    dplyr::relocate(width, .after = end)
  return(biomart_other)
}

#--- get organism annotation ---#
if(F){
  gsub('eck12|ecSakai|ss|xl','',biocOrg_name$short_name) %>%
    stringi::stri_remove_empty_na() %>%
    for(i in .){
      # i = biocOrg_name$short_name[1]
      if(!file.exists(paste0('data/',stringr::str_to_title(i),'_anno.rda'))){
        x = .get_org_anno(i)
        assign(paste0(stringr::str_to_title(i),'_anno'), x)
        save(list = paste0(stringr::str_to_title(i),'_anno'),
             file = paste0('data/',stringr::str_to_title(i),'_anno.rda'),
             compress = 'xz')
      }
    }
}

.get_org_anno <- function(org){
  org = mapBiocOrg(tolower(org))
  .load_orgdb(org)

  # keytypes(org.Hs.eg.db)
  # first get main id data (entrez, symbol, ensembl, uniprot)
  symbol_dat = AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egSYMBOL"))))

  ensembl_dat = AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egENSEMBL")))) %>%
    split(., .$gene_id) %>%
    lapply(., function(x) {
      paste0(x[, 2], collapse = "; ")
    })  %>% do.call(rbind,.) %>% as.data.frame() %>%
    dplyr::mutate(gene_id = rownames(.)) %>%
    dplyr::rename( ensembl_id = V1) %>%
    dplyr::select(gene_id,ensembl_id)
  uniprot_dat <- AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egUNIPROT"))))%>%
    split(., .$gene_id) %>%
    lapply(., function(x) {
      paste0(x[, 2], collapse = "; ")
    })  %>% do.call(rbind,.) %>% as.data.frame() %>%
    dplyr::mutate(gene_id = rownames(.)) %>%
    dplyr::rename( uniprot = V1) %>%
    dplyr::select(gene_id,uniprot)

  # then get gene name and alias
  name_dat = AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egGENENAME"))))
  biomart_alias = .get_biomrt_alias(org)
  biomart_feature = .get_biomrt_feature(org)
  ncbi_alias = AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egALIAS2EG")))) %>%
    split(., .$gene_id) %>%
    lapply(., function(x) {
      x = x[-nrow(x),]
      paste0(x[, 2], collapse = "; ")
    })  %>% do.call(rbind,.) %>% as.data.frame() %>%
    dplyr::mutate(gene_id = rownames(.)) %>%
    dplyr::rename( ncbi_alias = V1) %>%
    dplyr::select(gene_id,ncbi_alias)

  all = Reduce(function(x, y)
    merge(x, y, all=TRUE),
    list(symbol_dat, ensembl_dat, uniprot_dat,name_dat,ncbi_alias )) %>%
    dplyr::rename(entrezid = gene_id) %>%
    dplyr::rename(ensembl = ensembl_id) %>%
    merge(., biomart_alias, by='symbol',all.x=TRUE,all.y=FALSE) %>%
    # if add more from biomart, just add here...
    merge(., biomart_feature, by='ensembl',all.x=TRUE,all.y=FALSE) %>%
    dplyr::relocate(entrezid,.before = everything()) %>%
    dplyr::arrange(entrezid) %>%
    dplyr::relocate(chr,start,end,width,strand, .after = uniprot)


  return(all)

}









