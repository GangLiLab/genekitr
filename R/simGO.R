#' Simplitify GO enrichment result
#'
#' @param enrich_df Enrichment analysis `data.frame` result.
#' @param sim_method Method of calculating the similarity between nodes, one of one of "Resnik",
#' "Lin", "Rel", "Jiang" , "Wang" and "JC" (Jaccard similarity coefficient) methods.
#' @param org  Organism name from `biocOrg_name`.
#' @param ont  One of "bp", "mf", and "cc".
#'
#' @importFrom dplyr arrange distinct pull
#'
#' @return A `data.frame` contains simplified GO terms.
#' @export
#'
simGO <- function(enrich_df,
                  sim_method = c("JC","Resnik", "Lin", "Rel", "Jiang" , "Wang"),
                  org=NULL,
                  ont=NULL){
  sim_method <- match.arg(sim_method)

  simterm = get_sim_data(enrich_df,sim_method = sim_method,
                                   org=org, ont=ont) %>% .[["r"]] %>%
    dplyr::arrange(cluster) %>%
    dplyr::distinct(parent,.keep_all = TRUE) %>%
    dplyr::pull(go)

  egosim <- enrich_df %>% dplyr::filter(.[[1]] %in% simterm)


  return(egosim)
}



