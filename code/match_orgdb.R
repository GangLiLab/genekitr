match_org <- function(org,
                      dir,
                      ...){
  
  # org <- match.arg(org ,several.ok = T)
  
  if(length(org) >1 ){
    for(x in org){
      .search_db(x,dir)
    }
  }else{
    .search_db(x,dir)
  }

}


.search_db <- function(org, dir ){
  
  eg2symbol=toTable(eval(parse(text = paste0('org.',org,'.egSYMBOL'))))
  eg2name=toTable(eval(parse(text = paste0('org.',org,'.egGENENAME'))))
  eg2alias=toTable(eval(parse(text = paste0('org.',org,'.egALIAS2EG'))))
  eg2alis_list=lapply(split(eg2alias,eg2alias$gene_id),function(x){paste0(x[,2],collapse = ";")})
  eg2uniprot=toTable(eval(parse(text = paste0('org.',org,'.egUNIPROT'))))
  
  
  GeneList=mappedLkeys(eval(parse(text = paste0('org.',org,'.egSYMBOL'))))
  if( GeneList[1] %in% eg2symbol$symbol ){
    symbols=GeneList
    geneIds=eg2symbol[match(symbols,eg2symbol$symbol),'gene_id']
  }else{
    geneIds=GeneList
    symbols=eg2symbol[match(geneIds,eg2symbol$gene_id),'symbol']
  }
  
  geneNames=eg2name[match(geneIds,eg2name$gene_id),'gene_name']
  geneAlias=sapply(geneIds,function(x){ifelse(is.null(eg2alis_list[[x]]),"no_alias",eg2alis_list[[x]])})
  uniprotIds=eg2uniprot[match(geneIds,eg2uniprot$gene_id),'uniprot_id']
  
  createLink <- function(base,val) {
    sprintf('<a href="%s" class="btn btn-link" target="_blank" >%s</a>',base,val) 
    ##  target="_blank" 
  }
  gene_info=data.frame(   symbols=symbols,
                          geneIds=createLink(paste0("http://www.ncbi.nlm.nih.gov/gene/",geneIds),geneIds),
                          uniprotIds = ifelse(is.na(uniprotIds),'no_uniprot_id',
                                              createLink(paste0("https://www.uniprot.org/uniprot/",uniprotIds),uniprotIds)),
                          geneNames=geneNames,
                          geneAlias=geneAlias,
                          stringsAsFactors = F
  ) 
  
  file=paste0(dir,org,'_all_gene_bioconductor.html')
  y <- DT::datatable(gene_info,escape = F,rownames=F)
  DT::saveWidget(y,file)
}
