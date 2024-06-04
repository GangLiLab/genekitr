#' Transform probe id to symbol, entrezid,  ensembl or uniprot.
#' @param id probe ids.
#' @param transTo Transform to what type. User could select one or more from
#' "symbol", "entrez", "ensembl" or "uniprot."
#' @param org 'human'.
#' @param platform Probe platform. If NULL, program will detect automatically.
#'
#' @importFrom dplyr all_of filter select relocate everything pull rename group_by full_join distinct
#' @importFrom tidyr separate_rows fill
#' @importFrom rlang .data
#'
#' @return data frame, first column is probe id and others are converted id.
#' @export

transProbe <- function(id,
                       transTo,
                       org = "human",
                       platform = NULL) {

  #--- args ---#
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

  if (org == "hg" | org == "human" | org == "hsa" | org == "hs") org <- "hsapiens"
  if (org == "mm" | org == "mouse") org <- "mmusculus"
  if (org == "rn" | org == "rat") org <- "rnorvegicus"

  #--- codes ---#
  probe_dat <- probAnno(org)

  # get platform
  if (is.null(platform)) {
    platform <- which.max(apply(probe_dat, 2, function(x) {
      sum(id %in% x)
    })) %>% names()
    message(paste0(
      "Genekitr assumes probe ids derive from: ", platform,
      "\nIf not, please specify the platform from `genekitr::", org, "_probe_platform`..."
    ))
  } else if (!tolower(platform) %in% colnames(probe_dat)) {
    stop(paste0("Please select platform from `genekitr::", org, "_probe_platform`"))
  } else {
    platform <- tolower(platform)
  }

  # first get ensembl probe data
  ens_res <- probe_dat[, c("ensembl", dplyr::all_of(platform))] %>%
    tidyr::separate_rows(dplyr::all_of(platform), sep = "; ") %>%
    dplyr::filter(eval(parse(text = dplyr::all_of(platform))) %in% id) %>%
    dplyr::select(ensembl, dplyr::all_of(platform)) %>%
    stats::setNames(c("ensembl", "probe_id")) %>%
    dplyr::relocate(probe_id, .before = dplyr::everything())

  ens_res <- merge(data.frame(probe_id = id), ens_res, by = "probe_id", all.x = T)
  ens_res <- ens_res[match(id, ens_res$probe_id), ]

  # for those nomatched ids
  if (org == "hsapiens") prob_plats <- hsapiens_probe_data()

  bioc_pkg <- prob_plats %>%
    dplyr::filter(Platform %in% platform) %>%
    dplyr::pull(Bioc_anno)

  na_id <- ens_res[is.na(ens_res$ensembl), "probe_id"]
  add_res <- get_bioc_probe(na_id, to_type = "ensembl", bioc_pkg) %>%
    dplyr::rename(ensembl = ensembl_id)

  ens_res <- suppressMessages(dplyr::full_join(ens_res, add_res)) %>%
    dplyr::group_by(probe_id) %>%
    tidyr::fill(-probe_id, .direction = "downup") %>%
    dplyr::distinct()

  if (org == "hsapiens") rm(hsapiens_probe_platform, envir = .genekitrEnv)

  # trans to other types
  if (any(transTo %in% c("entrezid", "symbol", "uniprot"))) {
    type2 <- transTo[transTo %in% c("entrezid", "symbol", "uniprot")]
    new_dat <- suppressMessages(transId(stats::na.omit(ens_res$ensembl),
      transTo = type2,
      org = org, unique = T, keepNA = T
    )) %>%
      dplyr::rename(ensembl = input_id)
    res <- merge(ens_res, new_dat, by = "ensembl", all.x = T, all.y = F) %>%
      dplyr::relocate(probe_id, .before = dplyr::everything())
    res <- res[match(id, res$probe_id), ]
    if (!"ensembl" %in% transTo) {
      res <- res %>% dplyr::select(-ensembl)
    }
  } else {
    res <- ens_res
  }

  return(res)
}


get_bioc_probe <- function(id, to_type, bioc_pkg) {
  to_type <- sapply(tolower(to_type), function(x) {
    if (grepl(x, "entrezid")) x <- "ENTREZID"
    if (grepl(x, "ensemblid")) x <- "ENSEMBL"
    if (grepl(x, "symbolid")) x <- "SYMBOL"
    if (grepl(x, "uniprot")) x <- "UNIPROT"
    return(x)
  }) %>% as.character()

  if (!requireNamespace(bioc_pkg, quietly = TRUE)) {
    # BiocManager::install(bioc_pkg)
    message(paste0('Please firstly install ',bioc_pkg,':\n',
                   'BiocManager::install("',bioc_pkg,'")'))
    # pacman::p_load(bioc_pkg, character.only = TRUE)
  }

  bioc_dat <- AnnotationDbi::toTable(
    eval(parse(text = paste(bioc_pkg, "::", gsub(".db", "", bioc_pkg), to_type, sep = "")))
  ) %>%
    dplyr::filter(probe_id %in% id)

  return(bioc_dat)
}
