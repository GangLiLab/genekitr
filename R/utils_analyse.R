# AnnoGenes utilities for analysing

##' @title Get Msigdb database term and gene information
##' @param org organism name from `msigdb_org_data()`.
##' @param category MSigDB collection abbreviation, C1 to C8 and H.
##' @param subcategory MSigDB sub-collection abbreviation, such as REACTOME or BP.
##' @return a dataframe of 3 columns with term, entrez and symbol name.
##' @importFrom msigdbr msigdbr
##' @importFrom stringr str_split
##' @importFrom utils data
##' @importFrom dplyr %>% filter pull select
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
  options(warn=-1)
  if (!requireNamespace('msigdbr', quietly = TRUE)) auto_install('msigdbr')
  org = tolower(org)
  if (org == "hg" | org == "hsa" |  org == "hs" | org == 'homo sapiens') org = 'human'
  if (org == "mm" | org == "mmu") org = 'mouse'

  # org
  msig_org <- msigdb_org_data()
  all_org = c(msig_org[,1],
              stringr::str_split(msig_org[,2],', ',simplify = T) %>%
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
  msig_category <- msigdb_category_data()
  all_sub <- msig_category[,2] %>%
    stringi::stri_remove_empty_na()

  som_sub <- msig_category %>%
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


#---  get msigdb category ---#
msigdb_category_data <- function() {
  utils::data(list="msig_category", package="AnnoGenes")
  get("msig_category", envir = .GlobalEnv)
}


#---  map bioc org fullname to shortname ---#
mapBiocOrg <- function(organism) {
  organism = tolower(organism)
  if (organism == "hg" | organism == "human" | organism == "hsa" |  organism == "hs") organism = 'hs'
  if (organism == "mm" | organism == "mouse") organism = 'mm'

  # support organisms: http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
  biocorg = biocOrg_name()

  if( organism %in% biocorg$short_name ){
    org = organism
  }else if(organism %in% biocorg$full_name){
    org = biocorg %>% dplyr::filter(full_name == organism) %>% dplyr::pull(short_name)
  }else{
    stop('Check organism name! \n USE FULL NAME: ',
         paste0(biocOrg_name() %>% dplyr::pull(full_name),' | '),
         '\n OR USE SHORT NAME: ',
         paste0(biocOrg_name() %>% dplyr::pull(short_name),' | '))
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
  kgorg = keggOrg_name()
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


#--- transform gene id  ---#
transId <- function(id, trans_to, org, return_dat = FALSE){
  options(warn = -1)
  org.bk = org
  org = mapBiocOrg(tolower(org))
  keytype = .gentype(id, org)

  if (tolower(trans_to) == "entrez" | tolower(trans_to) == "entrezid") trans_to = 'entrezid'
  if (tolower(trans_to) == "ensemblid" | tolower(trans_to) == "ens" |
      tolower(trans_to) == "ensemb" ) trans_to = 'ensembl'
  if (tolower(trans_to) == "sym" ) trans_to = 'symbol'

  if(!tolower(trans_to) %in% c('symbol','entrezid','ensembl')){
    stop('\nChoose trans_to argument from: \nsymbol | entrezid | ensembl ')
  }

  if(tolower(keytype) == trans_to) stop('trans_to argument cannot be the same type with gene id!')
  # if org is one of human, mouse, rat
  # we could use much faster way to get id
  if(org %in% c('Hs', 'Mm' ,'Rn')){
    transId2(id, trans_to, org, return_dat)
  }else{
    .load_orgdb(org)
    symbol_dat = AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egSYMBOL"))))
    colnames(symbol_dat) = c('entrezid','symbol')
    ensem_dat = AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egENSEMBL"))))
    colnames(ensem_dat) = c('entrezid','ensembl')

    from = tolower(keytype)
    if(from == 'ensembl'){
      if(!any(id %in% ensem_dat[,from])){
        stop('IDs are not matched with specified organism: ',org.bk)
      }
    }else{
      if(!any(id %in% symbol_dat[,from])){
        stop('IDs are not matched with specified organism: ',org.bk)
      }
    }

  # others species use the common method (use org.db)
    merge_dat = merge(symbol_dat, ensem_dat ,by = 'entrezid', all.x  = T)
    newdat <- merge_dat %>%
      dplyr::select(c(all_of(from),all_of(trans_to))) %>%
      dplyr::filter(.[,1] %in% id) %>%
      dplyr::distinct()

    new_id =  newdat %>% dplyr::pull(2) %>% unique() %>% na.omit() %>% as.character()

    percen = paste(round(100*length(new_id)/length(unique(id)), 2), "%", sep="")
    if(length(new_id) > length(unique(id))){
      message(percen,' genes are mapped from ',from, ' to ', trans_to,'\n',
              'maybe one ', from, ' gene mapps many ', trans_to)
    } else {
      message(percen,' genes are mapped from ',from, ' to ', trans_to)
    }
  }

  if(return_dat){
    res = newdat
  }else{
    res = new_id
  }
  return(res)
}

transId2 <- function(id, trans_to, org, return_dat){
  options(warn = -1)
  if (org == "Hs" ) org = 'human'
  if (org == "Mm" ) org = 'mouse'
  if (org == "Rn" ) org = 'rat'

  gtf = eval(parse(text = paste0(org,'_gtf')))
  newdat = gtf %>%
    dplyr::filter(eval(parse(text = tolower(keytype))) %in% id)
  new_id = newdat %>% dplyr::pull(eval(parse(text = tolower(trans_to))))

  from = tolower(keytype)
  percen = paste(round(100*length(new_id)/length(unique(id)), 2), "%", sep="")
  if(length(new_id) > length(unique(id))){
    message(percen,' genes are mapped from ',from, ' to ', trans_to,'\n',
            'maybe one ', from, ' gene mapps many ', trans_to)
  } else {
    message(percen,' genes are mapped from ',from, ' to ', trans_to)
  }
}



#---calc fold enrichment ---#
calcFoldEnrich <- function(df){
  if( any(grepl('[gene|bg]ratio',tolower(colnames(df)))) ){
    check_gr = which(grepl('.*gene.*ratio',tolower(colnames(df))))
    check_bg = which(grepl('.*bg*ratio',tolower(colnames(df))))
    to_calc =  paste0('(',df[,check_gr],')/(',df[,check_bg],')')

    df <- df %>%
      dplyr::mutate(FoldEnrich = sapply(to_calc, function(x)eval(parse(text = x)) ))
  }
  return(df)

}




