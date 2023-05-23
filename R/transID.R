#' Transform id among symbol, entrezid,  ensembl and uniprot.
#' @param id Gene ids or protein ids.
#' @param transTo Transform to what type. User could select one or more from
#' "symbol", "entrez", "ensembl" or "uniprot."
#' @param org Latin organism shortname from `ensOrg_name`. Default is human.
#' @param unique Logical, if one-to-many mapping occurs, only keep one record with fewest NA. Default is FALSE.
#' @param keepNA If some id has no match at all, keep it or not. Default is FALSE.
#' @param hgVersion Select human genome build version from "v38" (default) and "v19".
#' @importFrom dplyr filter pull select distinct arrange all_of filter_at vars any_vars
#' @importFrom rlang .data
#'
#' @return data frame, first column is input id and others are converted id.
#' @export
#'
#' @examples
#' \donttest{
#' # example1:
#' transId(
#'   id = c("Cyp2c23", "Fhit", "Gal3st2b", "Trp53", "Tp53"),
#'   transTo = "ensembl", org = "mouse", keepNA = FALSE
#' )
#'
#' ## example2: input id with one-to-many mapping and fake one
#' transId(
#'   id = c("MMD2", "HBD", "RNR1", "TEC", "BCC7", "FAKEID", "TP53"),
#'   transTo = c("entrez", "ensembl"), keepNA = TRUE
#' )
#'
#' # example3: auto-recognize ensembl version number
#' transId("ENSG00000141510.11", "symbol")
#'
#' # example4: search genes with case-insensitive
#' transId(c('nc886','ezh2','TP53'),transTo = "ensembl",org = 'hs',unique = TRUE)
#' }
#'
transId <- function(id,
                    transTo,
                    org = "hs",
                    unique = FALSE,
                    keepNA = FALSE,
                    hgVersion = c("v38","v19")) {

  #--- args ---#
  org <- mapEnsOrg(organism = tolower(org))
  hgVersion <- match.arg(hgVersion)
  # if id has ensembl version, remove them
  # if(all(all %>% stringr::str_detect(.,'ENS'))) id <- stringr::str_split(id, "\\.", simplify = T)[, 1]

  transTo <- sapply(tolower(transTo), function(x) {
    if (grepl(x, "entrezid")) x <- "entrezid"
    if (grepl(x, "ensemblid")) x <- "ensembl"
    if (grepl(x, "symbolid")) x <- "symbol"
    if (grepl(x, "uniprot")) x <- "uniprot"
    return(x)
  }) %>% as.character()

  if (!all(tolower(transTo) %in% c("symbol", "entrezid", "ensembl", "uniprot"))) {
    stop("\nChoose 'transTo' argument from: \nsymbol | entrezid | ensembl | uniprot")
  }

  #--- codes ---#
  tryCatch(
    {
      all <- ensAnno(org,hgVersion = hgVersion)
      if(all(id %>% stringr::str_detect(.,'ENS'))) id <- stringr::str_split(id, "\\.", simplify = T)[, 1]
      id <- replace_greek(id)
      keytype <- gentype(id = id, data = all, org = org) %>% tolower()

      if(keytype != 'symbol'){
        res <- genInfo(id, org, unique, keepNA,hgVersion) %>%
          dplyr::mutate( !! keytype := input_id) %>%
          dplyr::select(input_id, all_of(transTo)) %>%
          distinct()
      }else{
        res <- genInfo(id, org, unique, keepNA,hgVersion) %>%
          dplyr::select(input_id, all_of(transTo)) %>%
          distinct()
      }

    },
    error = function(e) {
      message('Maybe your "trans_to" argument is wrong, please check again...')
    }
  )


  if (!keepNA) {
    res <- res %>%
      filter_at(vars(!input_id), any_vars(!is.na(.)))
  }

  # replace back greek letter
  res$input_id = replace_back(res$input_id)

  # convert factor to character
  res[] <- lapply(res, as.character)

  ## calculate percent
  n_new <- sapply(transTo, function(x) {
    # x = transTo[1]
    res %>%
      dplyr::select(input_id, all_of(x)) %>%
      stats::na.omit() %>%
      distinct(input_id) %>%
      nrow()
  })

  percent <- paste(round(100 * n_new / length(unique(id)), 2), "%", sep = "")
  sapply(seq_along(transTo), function(x) {
    # x = 1
    message(percent[x], " genes are mapped to ", transTo[x])
  }) %>% invisible()


  return(res)
}
