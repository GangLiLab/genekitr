##' Get gene info using org.db
##'
##' @param id gene id "symbol" or "entrezid".
##' @param org species "mm" or "hs".
##' @param html_result result saved as html or xlsx.
##' @param destdir result directory.
##' @return a dataframe or html of gene info.
##' @importFrom DT datatable saveWidget
##' @importFrom rio export
##'  @importFrom stringr str_to_title
##' @export
##' @examples
##' \dontrun{
##' genOrgInfo(id,org = 'mm',html_result = TRUE, dir = '~/Downloads')
##' }

genOrgInfo <- function(id,
                       org = c("mm", "hs"),
                       html_result = TRUE,
                       destdir = tempdir(),
                       ...) {
  stopifnot(
    is.logical(html_result),
    is.character(id),
    org %in% c("mm", "hs")
  )

  org = str_to_title(org) 
  
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


  if (html_result) {
    createLink <- function(base, val) {
      sprintf('<a href="%s" class="btn btn-link" target="_blank" >%s</a>', base, val)
      ##  target="_blank"
    }

    gene_info <- data.frame(
      symbols = symbols,
      geneIds = createLink(paste0("http://www.ncbi.nlm.nih.gov/gene/", geneIds), geneIds),
      uniprotIds = ifelse(is.na(uniprotIds), "no_uniprot_id",
        createLink(paste0("https://www.uniprot.org/uniprot/", uniprotIds), uniprotIds)
      ),
      geneNames = geneNames,
      geneAlias = geneAlias,
      stringsAsFactors = F
    )

    file <- paste0(destdir, "/", org, "_gene_orgdb.html")
    dt <- DT::datatable(gene_info, escape = F, rownames = F)
    DT::saveWidget(dt, file)
  } else {
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

    rio::export(gene_info, file = paste0(destdir, "/", org, "_gene_orgdb.xlsx"))
  }
}
