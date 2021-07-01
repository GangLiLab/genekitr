##' Get gene info using org.db
##'
##' @param id gene id "symbol" or "entrezid".
##' @param org species "mm" or "hs".
##' @return a dataframe of gene info.
##' @importFrom stringr str_to_title
##' @export
##' @examples
##' \dontrun{
##' x = genInfo(id= c("Cyp2c23","Fhit","Gal3st2b","Gbp4"), org = 'mm')
##' }

genInfo <- function(id,
                    org = c("mm", "hs"),
                    ...) {
  #--- args ---#
  options(rstudio.connectionObserver.errorsSuppressed = TRUE)

  org = match.arg(org)
  stopifnot(
    is.character(id),
    org %in% c("mm", "hs")
  )
  org <- stringr::str_to_title(org)

  #--- code ---#
  # load org data
  suppressPackageStartupMessages(require(paste0("org.", org, ".eg.db"), character.only = TRUE))
  # get org data
  orgSymbol <- toTable(eval(parse(text = paste0("org.", org, ".egSYMBOL"))))
  orgName <- toTable(eval(parse(text = paste0("org.", org, ".egGENENAME"))))
  orgAlias <- toTable(eval(parse(text = paste0("org.", org, ".egALIAS2EG"))))
  orgAlias_list <- lapply(split(orgAlias, orgAlias$gene_id), function(x) {
    paste0(x[, 2], collapse = ";")
  })
  orgUniprot <- toTable(eval(parse(text = paste0("org.", org, ".egUNIPROT"))))

  # id could be SYMBOL or EntrezID
  if (any(id %in% orgSymbol$symbol)) {
    symbols <- id
    geneIds <- orgSymbol[match(symbols, orgSymbol$symbol), "gene_id"]
  } else {
    geneIds <- id
    symbols <- orgSymbol[match(geneIds, orgSymbol$gene_id), "symbol"]
  }

  if(all(is.na(symbols))){
    stop('No result...\nmaybe choose the wrong species...')
  }
  geneNames <- orgName[match(geneIds, orgName$gene_id), "gene_name"]
  geneAlias <- sapply(geneIds, function(x) {
    ifelse(is.null(orgAlias_list[[x]]), NA, orgAlias_list[[x]])
  })


  uniprotIds <- orgUniprot[match(geneIds, orgUniprot$gene_id), "uniprot_id"]

  gene_info <- data.frame(
    symbols = ifelse(is.na(symbols), NA,symbols),
    geneIds = ifelse(!geneIds%in%orgSymbol$gene_id, NA,
                     paste0("http://www.ncbi.nlm.nih.gov/gene/", geneIds)),
    uniprotIds = ifelse(is.na(uniprotIds), NA,
                        paste0("https://www.uniprot.org/uniprot/", uniprotIds)),
    geneNames = ifelse(is.na(geneNames), NA,geneNames),
    geneAlias = ifelse(is.na(geneAlias), NA,geneAlias),
    stringsAsFactors = F
  )


  return(gene_info)
}
