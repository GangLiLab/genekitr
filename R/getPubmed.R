#' Get pubmed paper records by searching abstract
#'
#' @param term query terms e.g. gene id, GO/KEGG term or id
#' @param keys other searching keys
#' @importFrom europepmc epmc_search
#' @return A list of `tibble` for pubmed records
#' @export
#' @examples
#' \donttest{
#' term <- c("Tp53","Brca1","Tet2")
#' keys <- c('stem cell','mouse')
#' l <- getPubmed(term,keys)
#' # very easy to output
#' expoSheet(l,name_list = term, filename = 'test.xlsx',dir = tempdir())
#' }
#'
getPubmed <- function(term,keys){
  if (!requireNamespace("europepmc", quietly = TRUE)) {
    utils::install.packages("europepmc")
  }

  supp = paste0('AND ABSTRACT:',keys) %>% paste(collapse = ' ')
  res = lapply(term , function(i){
    europepmc::epmc_search(query = paste0('ABSTRACT:',i," ",supp))
  })
  names(res) = term

  return(res)
}



#---BELOW DEGRADED! ---#
# genPubmed <- function(id,
#                       keywords,
#                       field = "ALL") {
#   #--- args ---#
#   if (!requireNamespace("easyPubMed", quietly = TRUE)) {
#     stop("Package easyPubMed needed for this function to work. Please install it.",
#       call. = FALSE
#     )
#   }
#
#   stopifnot(
#     is.character(id),
#     is.character(keywords)
#   )
#
#   field <- toupper(field)
#
#   data <- showNCBI("pubmed")
#   if (field %in% data$Name) {
#     message(paste0(
#       "Search example: ",
#       paste0(id[1], " [", field, "] AND ", keywords, " [", field, "] ")
#     ))
#   } else {
#     warning(paste0("Wrong field argument: ", field), "\nUse default: ALL")
#     field <- "ALL"
#   }
#
#
#   #--- code ---#
#   l <- list()
#   l <- sapply(id, function(x) {
#     query <- paste0(x, " [", field, "] AND ", keywords, " [", field, "] ")
#     entrez_id <- easyPubMed::get_pubmed_ids(query)
#     if (entrez_id$Count == "0") {
#       pub_df <- data.frame(title = "NA", date = "NA", doi = "NA", pmid = "NA", journal = "NA")
#       l[[x]] <- pub_df
#     } else {
#       myxml <- easyPubMed::fetch_pubmed_data(pubmed_id_list = entrez_id)
#       pub_df <- suppressMessages(easyPubMed::table_articles_byAuth(
#         pubmed_data = myxml,
#         included_authors = "first",
#         max_chars = 0,
#         encoding = "ASCII"
#       ) %>%
#         tidyr::unite("date", year, month, day) %>%
#         dplyr::select(title, date, doi, pmid, journal))
#
#       l[[x]] <- pub_df
#     }
#     invisible(l)
#   })
#   names(l) <- id
#
#   res <- lapply(seq_len(length(l)), function(n) {
#     l[[n]] %>%
#       dplyr::mutate(gene = rep(names(l)[n], nrow(l[[n]]))) %>%
#       dplyr::relocate(gene, .before = dplyr::everything())
#   }) %>%
#     do.call(rbind, .)
#
#   return(res)
# }
#
#
# showNCBI <- function(db = "pubmed") {
#   # suppress binding notes
#   fields <- rentrez::entrez_db_searchable(db)
#   res <- as.data.frame(fields)[1:3]
#
#   if (nrow(res) == 0) { # nocov start
#     message("Something is wrong in your input,
#             NULL will be returned, please check.")
#     return(NULL)
#   } # nocov end
#   return(res)
# }
#
# utils::globalVariables(c("year","month","day","title","doi","pmid","journal","gene"))
