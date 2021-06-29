##' Get gene info using org.db
##'
##' @param id gene id "symbol" or "entrezid".
##' @param org species "mm" or "hs".
##' @return a dataframe or html of gene info.
##' @importFrom DT datatable saveWidget
##' @importFrom rio export
##' @importFrom stringr str_to_title
##' @export
##' @examples
##' \dontrun{
##' genInfo(id,org = 'mm')
##' }

genInfo <- function(id,
                    org = c("mm", "hs"),
                       # html_result = TRUE,
                       # destdir = tempdir(),
                       ...) {
  #--- args ---#
  stopifnot(
    is.character(id),
    org %in% c("mm", "hs")
  )

  org = str_to_title(org) 
  
  #--- code ---#
  # load org data
  require(paste0("org.", org, ".eg.db"), character.only = TRUE)
  # get org data
  eg2symbol <- toTable(eval(parse(text = paste0("org.", org, ".egSYMBOL"))))
  eg2name <- toTable(eval(parse(text = paste0("org.", org, ".egGENENAME"))))
  eg2alias <- toTable(eval(parse(text = paste0("org.", org, ".egALIAS2EG"))))
  eg2alis_list <- lapply(split(eg2alias, eg2alias$gene_id), function(x) {
    paste0(x[, 2], collapse = ";")
  })
  eg2uniprot <- toTable(eval(parse(text = paste0("org.", org, ".egUNIPROT"))))

  # id could be SYMBOL or EntrezID
  if (id[1] %in% eg2symbol$symbol) {
    symbols <- id
    geneIds <- eg2symbol[match(symbols, eg2symbol$symbol), "gene_id"]
  } else {
    geneIds <- id
    symbols <- eg2symbol[match(geneIds, eg2symbol$gene_id), "symbol"]
  }
  geneIds[which(is.na(geneIds))]='NA' 
  geneNames <- eg2name[match(geneIds, eg2name$gene_id), "gene_name"]
  geneAlias <- sapply(geneIds, function(x) {
    ifelse(is.null(eg2alis_list[[x]]), "no_alias", eg2alis_list[[x]])
  })
  uniprotIds <- eg2uniprot[match(geneIds, eg2uniprot$gene_id), "uniprot_id"]

  gene_info <- data.frame(
    symbols = symbols,
    geneIds = paste0("http://www.ncbi.nlm.nih.gov/gene/", geneIds),
    uniprotIds = ifelse(is.na(uniprotIds), "no_uniprot_id",
                        paste0("https://www.uniprot.org/uniprot/", uniprotIds)
    ),
    geneNames = geneNames,
    geneAlias = geneAlias,
    stringsAsFactors = F
  )
  
  invisible(gene_info)
  
}
