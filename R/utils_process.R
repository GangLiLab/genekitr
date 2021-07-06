# AnnoGenes utilities for processing

##' @title Show NCBI database searchable name
##' @param db a character vector. Can be "pubmed" or one or more of `rentrez::entrez_dbs()` result.
##' @return a dataframe including search keyword and information.
##' @importFrom rentrez entrez_db_searchable
##' @export
##' @examples
##' \donttest{
##' showNCBI("pubmed")
##' }
showNCBI <- function(db = "pubmed") {
  # suppress binding notes
  fields <- rentrez::entrez_db_searchable(db)
  res <- as.data.frame(fields)[1:3]

  if (nrow(res) == 0) { # nocov start
    message("Something is wrong in your input, NULL will be returned, please check.")
    return(NULL)
  } # nocov end
  return(res)
}

##' export result into different sheets
##' @param wb worksheet from `createWorkbook()`.
##' @param sheet_dat dataframe added to sheet.
##' @param sheet_name name of added dataframe.
##' @return a worksheet including many dataframes.
##' @importFrom stringr str_detect
##' @importFrom openxlsx addWorksheet writeData writeFormula createStyle addStyle setColWidths
##' @export
##' @examples
##' \donttest{
##' expo_sheet(wb, sheet_dat =  mtcars, sheet_name = 'mtcars')
##' }
expo_sheet <- function(wb, sheet_dat, sheet_name) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = sheet_dat)

  ## add hyperlink for gene
  check <- apply(sheet_dat, 2, function(x) {
    any(stringr::str_detect(x, "http"))
  })
  # check[is.na(check)] <- FALSE
  if (any(check)) {
    sub_dat <- sheet_dat[, check]

    for (i in 1:nrow(sheet_dat)) {
      for (n in which(check)) {
        if (! is.na(sheet_dat[i, n]) & !is.na(sheet_dat[i, n] )) {
          sheet_datd_link <- paste0("HYPERLINK(\"", sheet_dat[i, n], "\", \"", sheet_dat[i, n], "\")")
          writeFormula(wb, sheet = sheet_name, startRow = i + 1, startCol = n, x = sheet_datd_link)
        }
      }
    }
  }

  ## add hyperlink for pubmed
  check2 <- any(stringr::str_detect(colnames(sheet_dat), "pmid"))
  if (check2) {
    for (i in seq_len(nrow(sheet_dat))) {
      if (sheet_dat[i, 2] != "NA") {
        sheet_datd_link <- paste0(
          "HYPERLINK(\"", paste0("https://pubmed.ncbi.nlm.nih.gov/", sheet_dat[i, 5]),
          "\", \"", sheet_dat[i, 2], "\")"
        )
        writeFormula(wb, sheet = sheet_name, startRow = i + 1, startCol = 2, x = sheet_datd_link)
      }
    }
  }

  ## styling sheet
  headerStyle <- createStyle(textDecoration = "Bold")
  addStyle(wb, sheet = sheet_name, style = headerStyle, rows = 1, cols = seq_len(ncol(sheet_dat)), gridExpand = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(sheet_dat)), widths = "auto")

  invisible(wb)
}


#---   define gene type: entrezid, ensembl or symbol ---#
.gentype <- function(id, org){
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
  }else{
    c("ENTREZID")
  }
}

#---  gene id in this org or not (return a logical) ---#
.genInorg <- function(id, org){
  if(nchar(org) > 2){
    org = substr(org,1,nchar(org)-1)
  }
  org <- stringr::str_to_title(org)
  suppressPackageStartupMessages(require(paste0("org.", org, ".eg.db"), character.only = TRUE))
  keytype = .gentype(id,org)

  orgSymbol <- AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egSYMBOL"))))
  if(keytype == "ENTREZID"){
    ifelse(any(id%in%orgSymbol$gene_id),TRUE,FALSE)
  }else if(keytype == "SYMBOL"){
    ifelse(any(id%in%orgSymbol$symbol),TRUE,FALSE)
  }else{
    orgENSEMBL <- AnnotationDbi::toTable(eval(parse(text = paste0("org.", org, ".egENSEMBL"))))
    ifelse(any(id%in%orgENSEMBL$ensembl_id),TRUE,FALSE)
  }
}

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
        suppressWarnings(install.packages(missing_pkgs, quiet = TRUE, dependencies = TRUE))

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
  org = mapBiocOrg(tolower(org))
  pkg=paste0("org.", org, ".eg.db")
  if (!requireNamespace(pkg, quietly = TRUE)) auto_install(pkg)
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))
}

#---check enrichment data colname---#
# 'GeneRatio','Count','FoldEnrich'
.check_colname <- function(enrich_df, x){
  remove <- c("\\(", "\\)", " ",'-')
  to_check = stringr::str_remove_all(tolower(colnames(enrich_df)), paste(remove, collapse = "|"))
  xlabs = c('GeneRatio','Count','FoldEnrich')
  legends = c('pvalue','p.adjust','qvalue')
  # to search xlab type

  if(x %in% xlabs){
    avail_opt = sapply(tolower(xlabs), function(n) any(grepl(n , to_check)))

    if( any(grepl( tolower(x),  to_check)) ){
      colnames(enrich_df)[grepl(tolower(x),  to_check)]  = x
    }else if( any(avail_opt) ){
      stop('\n"',x,'" was not found in enrich_df colnames',
           '\nMaybe could try other options: ', paste(xlabs[avail_opt],collapse = ', '))
    }else{
      stop("Not found any colname of: 'GeneRatio','Count','FoldEnrich'")
    }

    #  to search legend type
    # Q-values are the name given to the adjusted p-values found using an optimised FDR approach.
    # p.adjust has other methods, like "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"
  }else if(x %in% legends){
    if(any(grepl('fdr',  to_check))) colnames(enrich_df)[grepl('fdr',  to_check)]  = 'qvalue'
    to_check = stringr::str_remove_all(tolower(colnames(enrich_df)), paste(remove, collapse = "|"))
    avail_opt2 = sapply(legends, function(n) any(grepl(n , to_check)))

    if( any(grepl( tolower(x),  to_check)) ){
      colnames(enrich_df)[grepl(x,  to_check)]  = x
    }else if( any(avail_opt2) ){
      stop('\n"',x,'" was not found in enrich_df colnames',
           '\nMaybe could try other options: ', paste(legends[avail_opt2],collapse = ', '))
    }else{
      stop("Not found any colname of: 'pvalue','p.adjust','qvalue|fdr'")
    }

  }else{
    stop('Check input argument again!')
  }

  return(enrich_df)
}






