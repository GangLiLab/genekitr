#' Gene related pubmed paper records
#'
#' @param id Gene "symbol".
#' @param keywords species "mm" or "hs".
#' @param field pubmed field from `showNCBI('pubmed')`, default is "ALL".
#' @importFrom tidyr unite
#' @importFrom dplyr select mutate relocate everything %>%
#'
#' @return A `data.frame`.
#' @export
#' @examples
#' \donttest{
#' x = genPubmed(id = c("Cyp2c23","Fhit","Gal3st2b","Gbp4"),
#'   keywords = 'stem cell', field = 'tiab')
#' }
genPubmed <- function(id,
                      keywords,
                      field = "ALL") {
  #--- args ---#
  if (!requireNamespace("easyPubMed", quietly = TRUE)) {
    stop("Package easyPubMed needed for this function to work. Please install it.",
         call. = FALSE)
  }

  stopifnot(
    is.character(id),
    is.character(keywords)
  )

  field = toupper(field)

  data <- showNCBI("pubmed")
  if (field %in% data$Name) {
    message(paste0(
      "Search example: ",
      paste0(id[1], " [", field, "] AND ", keywords, " [", field, "] ")
    ))
  } else {
    warning(paste0("Wrong field argument: ", field), "\nUse default: ALL")
    field <- "ALL"
  }


  #--- code ---#
  l <- list()
  l <- sapply(id, function(x) {
    query <- paste0(x, " [", field, "] AND ", keywords, " [", field, "] ")
    entrez_id <- easyPubMed::get_pubmed_ids(query)
    if (entrez_id$Count == "0") {
      pub_df <- data.frame(title = "NA", date = "NA", doi = "NA", pmid = "NA", journal = "NA")
      l[[x]] <- pub_df
    } else {
      myxml <- easyPubMed::fetch_pubmed_data(pubmed_id_list = entrez_id)
      pub_df <- suppressMessages(easyPubMed::table_articles_byAuth(
        pubmed_data = myxml,
        included_authors = "first",
        max_chars = 0,
        encoding = "ASCII"
      ) %>%
        tidyr::unite("date", year, month, day) %>%
        dplyr::select(title, date, doi, pmid, journal))

      l[[x]] <- pub_df
    }
    invisible(l)
  })
  names(l) <- id

  res <- lapply(seq_len(length(l)), function(n) {
    l[[n]] %>%
      dplyr::mutate(gene = rep(names(l)[n], nrow(l[[n]]))) %>%
      dplyr::relocate(gene, .before = dplyr::everything())
  }) %>%
    do.call(rbind, .)

  invisible(res)
}
