#' Plot for gene enrichment analysis of ORA method
#'
#' Over-representation analysis (ORA) is a simple method for objectively deciding whether a set of variables of
#' known or suspected biological relevance, such as a gene set or pathway, is more prevalent in a set of variables
#' of interest than we expect by chance.
#'
#' @param enrich_df Enrichment analysis `data.frame` result.
#' @param fold_change Fold change or logFC values with gene IDs as names. Used in "heat" and "chord"
#' plot.
#' @param plot_type Choose from "bar", "wego","bubble","dot", "lollipop","geneheat", "genechord",
#' "network","gomap","goheat","gotangram","wordcloud","upset".
#' @param term_metric Pathway term metric from one of 'GeneRatio','Count','FoldEnrich' and
#' 'RichFactor'.
#' @param stats_metric Statistic metric from one of "pvalue", "p.adjust", "qvalue".
#' @param sim_method Method of calculating the similarity between nodes, one of one of "Resnik",
#' "Lin", "Rel", "Jiang" , "Wang" or "JC" (Jaccard’s similarity index). Only "JC" supports KEGG data.
#' Used in "map","goheat","gotangram","wordcloud".
#' @param up_color Color of higher statistical power (e.g. Pvalue 0.01) or higher logFC, default is "red".
#' @param down_color Color of lower statistical power (e.g. Pvalue 1) or lower logFC, default is
#' "blue".
#' @param show_gene Select genes to show. Default is "all". Used in "heat" and "chord" plot.
#' @param xlim_left X-axis left limit, default is 0.
#' @param xlim_right X-axis right limit, default is NA.
#' @param wrap_length Numeric, wrap text if longer than this length. Default is NULL.
#' @param org  Organism name from `biocOrg_name`.
#' @param ont  One of "BP", "MF", and "CC".
#' @param scale_ratio Numeric, scale of node and line size.
#' @param layout Grapgh layout in "map" plot, e,g, "circle", "dh", "drl", "fr","graphopt", "grid",
#' "lgl", "kk", "mds", "nicely" (default),"randomly", "star".
#' @param n_term Number of terms (used in WEGO plot)
#' @param ... other arguments from `plot_theme` function
#'
#' @importFrom ggplot2 ggplot aes aes_string aes_ facet_grid element_text element_blank labs geom_point
#' geom_segment geom_bar geom_col geom_tile guide_colorbar guides scale_color_continuous
#' scale_size_continuous scale_x_discrete scale_y_discrete scale_x_continuous scale_size
#' scale_fill_discrete scale_fill_continuous scale_y_continuous sec_axis theme xlab ylab xlim
#' @importFrom dplyr arrange mutate group_by top_n ungroup select case_when distinct rename pull
#' @importFrom stringr str_replace str_split
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom ggraph ggraph geom_node_text geom_edge_link circle geom_node_point
#' geom_node_label
#' @importFrom igraph graph.data.frame delete.edges
#' @importFrom igraph V "V<-"
##' @importFrom igraph E "E<-"
#'
#' @return A ggplot object
#' @export
#' @examples
#' \donttest{
#' ## example data
#' ## More examples please refer to https://www.genekitr.fun/plot-ora-1.html
#' library(ggplot2)
#' library(igraph)
#' library(ggraph)
#' data(geneList, package = "genekitr")
#' id <- names(geneList)[abs(geneList) > 1.5]
#' logfc <- geneList[id]
#'
#' gs <- geneset::getGO(org = "human",ont = "bp")
#' ego <- genORA(id, geneset = gs)
#' ego <- ego[1:10, ]
#'
#' ## example plots
#' plotEnrich(ego, plot_type = "dot")
#'
#' plotEnrich(ego, plot_type = "bubble", scale_ratio = 0.4)
#'
#' plotEnrich(ego, plot_type = "bar")
#'
#' plotEnrich(ego,
#'   plot_type = "lollipop",
#'   down_color = "#325CAC", up_color = "#E69056",
#'   wrap_length = 25, scale_ratio = 0.4
#' )
#'
#' plotEnrich(ego, plot_type = "geneheat")
#'
#' }
#'
plotEnrich <- function(enrich_df,
                       fold_change = NULL,
                       plot_type = c(
                         "bar", "wego", "dot", "bubble", "lollipop", "geneheat", "genechord",
                         "network", "gomap", "goheat", "gotangram", "wordcloud", "upset"
                       ),
                       term_metric = c("FoldEnrich", "GeneRatio", "Count", "RichFactor"),
                       stats_metric = c("p.adjust", "pvalue", "qvalue"),
                       sim_method = c("Resnik", "Lin", "Rel", "Jiang", "Wang", "JC"),
                       up_color = "#E31A1C",
                       down_color = "#1F78B4",
                       show_gene = "all",
                       xlim_left = 0,
                       xlim_right = NA,
                       wrap_length = NULL,
                       org = NULL,
                       ont = NULL,
                       scale_ratio,
                       layout,
                       n_term,
                       ...) {
  #--- args ---#
  lst <- list(...) # store outside arguments in list
  plot_type <- match.arg(plot_type)
  term_metric <- match.arg(term_metric)
  stats_metric <- match.arg(stats_metric)
  sim_method <- match.arg(sim_method)
  if (missing(layout)) layout <- "nicely"

  compare_group <- any(grepl("cluster", colnames(enrich_df), ignore.case = T))
  if (compare_group) plot_type <- "dot"

  all_go <-any(grepl("ONTOLOGY", colnames(enrich_df), ignore.case = T))
  if(all_go){
    colnames(enrich_df)[tolower(colnames(enrich_df))%in%'ontology'] = 'ONTOLOGY'
  }

  if(any(grepl("nes",colnames(enrich_df),ignore.case = T))) term_metric <- "Count"

  # if (all_go & !plot_type %in% c("bar", "wego")) {
  #   warning(paste0(
  #     'If you want to plot all ontologies data, please choose "plot_type" from "bar","wego"',
  #     "\nInstead, bar plot will be plotted..."
  #   ))
  # }

  types <- c("GeneRatio", "Count", "FoldEnrich", "RichFactor")
  legends <- c("p.adjust", "pvalue", "qvalue")

  if (!term_metric %in% colnames(enrich_df)) {
    stop(
      term_metric, " not included in this dataframe, try: ",
      paste(intersect(colnames(enrich_df), types), collapse = " | ")
    )
  }
  if (!stats_metric %in% colnames(enrich_df)) {
    stop(
      stats_metric, " not included in this dataframe, try: ",
      paste(intersect(colnames(enrich_df), legends), collapse = " | ")
    )
  }

  #--- codes ---#
  ## uppercase first letter of description
  tryCatch(
    {
      enrich_df$Description <- stringr::str_replace(enrich_df$Description, "^\\w{1}", toupper)
    },
    error = function(e) {
      message(paste0(
        "We need the 'Description' column which means pathway detailed description", "\n",
        "Maybe you should rename the column name..."
      ))
    }
  )

  ## set labels
  term_metric_label <- ifelse(term_metric == "FoldEnrich", "Fold Enrichment",
    ifelse(term_metric == "GeneRatio", "Gene Ratio",
      ifelse(term_metric == "RichFactor", "Rich Factor", "Count")
    )
  )
  stats_metric_label <- ifelse(stats_metric == "pvalue", "Pvalue",
    ifelse(stats_metric == "p.adjust", "P.adjust", "FDR")
  )

  ## set pathway as factor
  if (!compare_group & !all_go) {
    enrich_df <- enrich_df %>%
      dplyr::arrange(eval(parse(text = term_metric))) %>%
      dplyr::mutate(Description = factor(.$Description, levels = .$Description, ordered = T))
  } else if (!compare_group & all_go) {
    enrich_df <- enrich_df %>%
      dplyr::mutate(Description = factor(.$Description, levels = .$Description, ordered = T)) %>%
      dplyr::group_by(ONTOLOGY) %>%
      dplyr::arrange(eval(parse(text = term_metric)), .by_group = T)
  }

  #--- GO/KEGG: dot plot ---#
  if (plot_type == "dot") {
    ## no group
    if (missing(scale_ratio)) scale_ratio <- 0.3
    if (!compare_group) {
      p <- ggplot(enrich_df, aes_string(x = term_metric, y = "Description")) +
        geom_point(aes_string(
          color = stats_metric,
          size = "Count"
        )) +
        scale_size(range = c(min(enrich_df$Count) / 2, max(enrich_df$Count) / 2) * scale_ratio) +
        scale_color_continuous(
          low = up_color, high = down_color, name = stats_metric_label,
          guide = guide_colorbar(reverse = TRUE),
          labels = function(x) format(x, scientific = T)
        ) +
        xlab(term_metric_label) +
        labs(color = stats_metric) +
        xlim(xlim_left, xlim_right) +
        plot_theme(...)
    } else {
      ## compare groups
      if(missing(n_term)) n_term = 5
      enrich_df <- enrich_df %>%
        dplyr::arrange(desc(FoldEnrich)) %>%
        dplyr::group_by(Cluster) %>%
        dplyr::slice_head(n=n_term)

      calc <- enrich_df%>%
        dplyr::group_by(Cluster) %>%
        dplyr::summarise(sum_count = sum(Count)) %>%
        as.data.frame()

      xtick_lab <- paste0(as.character(calc[,1]), "\n(", calc[,2], ")")
      p <- ggplot(enrich_df, aes_string(x = "Cluster", y = "Description")) +
        geom_point(aes_string(
          color = stats_metric,
          size = term_metric
        )) +
        scale_size(range = c(min(enrich_df$Count) / 2, max(enrich_df$Count) / 2) * scale_ratio) +
        scale_color_continuous(
          low = up_color, high = down_color, name = stats_metric_label,
          guide = guide_colorbar(reverse = TRUE),
          labels = function(x) format(x, scientific = T)
        ) +
        scale_x_discrete(labels = xtick_lab) +
        xlab("Group") +
        labs(color = stats_metric) +
        plot_theme(...) +
        theme(axis.text.x = element_text(
          angle = 45,
          vjust = 0.5, hjust = 0.5
        ))
    }
  }

  #--- GO/KEGG: bubble plot ---#
  ## each bubble is a pathway term
  ## x-axis: stats metric e.g. pvalue/qvalue/p.adjust
  ## y-axis: fold enrichment
  if (plot_type == "bubble") {
    if (missing(scale_ratio)) scale_ratio <- 0.3
    # default main and lengend text size
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 8

    mapping <- aes(
      x = -log10(eval(parse(text = stats_metric))),
      y = FoldEnrich
    )
    p <- ggplot(enrich_df, mapping) +
      geom_point(aes(size = Count, fill = Description), alpha = 2, shape = 21) +
      scale_size(range = c(min(enrich_df$Count) / 2, max(enrich_df$Count) / 2) * scale_ratio) +
      scale_fill_discrete(guide = "none") +
      ggrepel::geom_text_repel(aes(label = Description),
        size = lst$main_text_size / 2.8,
        segment.color = "black"
      ) +
      scale_x_continuous(
        limits = c(0, ceiling(max(-log10(enrich_df[stats_metric]))) + 3),
        breaks = seq(0, ceiling(max(-log10(enrich_df[stats_metric]))) + 3, by = 3)
      ) +
      xlab(paste0("-log10(", stats_metric_label, ")")) +
      ylab("Fold Enrichment") +
      plot_theme( ...)
  }

  #--- GO/KEGG: bar plot ---#
  if (plot_type == "bar") {
    p <- ggplot(data = enrich_df, aes_string(x = term_metric, y = "Description", fill = stats_metric)) +
      geom_bar(stat = "identity") +
      scale_fill_continuous(
        low = up_color, high = down_color, name = stats_metric_label,
        guide = guide_colorbar(reverse = TRUE),
        labels = function(x) format(x, scientific = T)
      ) +
      xlab(term_metric_label) +
      labs(color = stats_metric) +
      xlim(xlim_left, xlim_right) +
      plot_theme(...)

    if (all_go) {
      p <- p + ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free") +
        plot_theme(...)
    }
  }

  #--- GO/KEGG: lollipop plot ---#
  if (plot_type == "lollipop") {
    if (missing(scale_ratio)) scale_ratio <- 0.3
    p <- ggplot(
      data = enrich_df,
      aes(
        eval(parse(text = term_metric)),
        forcats::fct_reorder(Description, eval(parse(text = term_metric)))
      )
    ) +
      geom_segment(aes_string(
        xend = 0, yend = "Description",
        colour = stats_metric, size = 2 * scale_ratio
      ), show.legend = F) +
      geom_point(aes_string(color = stats_metric, size = "Count")) +
      theme(legend.key = element_rect(fill = "transparent")) +
      scale_size_continuous(name = "Count", range = c(min(enrich_df$Count) / 2, max(enrich_df$Count) / 2) * scale_ratio) +
      scale_color_continuous(
        low = up_color, high = down_color, name = stats_metric_label,
        guide = guide_colorbar(reverse = TRUE),
        labels = function(x) format(x, scientific = T)
      ) +
      xlab(term_metric_label) +
      ylab(NULL) +
      labs(color = stats_metric) +
      plot_theme(...)
  }

  #--- GO/KEGG: geneheat plot ---#
  ## heatmap to show interaction between go term and gene id
  if (plot_type == "geneheat") {
    id <- enrich_df %>%
      dplyr::pull(geneID) %>%
      stringr::str_split("\\/") %>%
      unlist()
    id_symbol <- enrich_df %>%
      dplyr::pull(geneID_symbol) %>%
      stringr::str_split("\\/") %>%
      unlist()
    id_df <- data.frame(geneID = id, geneID_symbol = id_symbol) %>%
      dplyr::distinct()

    if (all(show_gene == "all")) {
      plot_df <- enrich_df %>%
        dplyr::select(Description, geneID_symbol) %>%
        tidyr::separate_rows(geneID_symbol, sep = "\\/") %>%
        dplyr::rename(geneID = geneID_symbol)
    } else {
      # if show_gene is not symbol, first extract matching symbol
      if (all(show_gene %in% id)) {
        show_gene <- id_df %>%
          dplyr::filter(geneID %in% show_gene) %>%
          dplyr::pull(geneID_symbol)
      }

      plot_df <- enrich_df %>%
        dplyr::select(Description, geneID_symbol) %>%
        tidyr::separate_rows(geneID_symbol, sep = "\\/") %>%
        dplyr::filter(geneID_symbol %in% show_gene) %>%
        dplyr::rename(geneID = geneID_symbol)
    }


    if (is.null(fold_change)) {
      p <- ggplot(plot_df, aes_(~geneID, ~Description)) +
        geom_tile(color = "white") +
        xlab(NULL) +
        ylab(NULL) +
        plot_theme(border_thick = 0, ...) +
        theme(
          panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 50, hjust = 1)
        )
    } else {
      # add logfc
      fold_change <- data.frame(geneID = names(fold_change), logfc = fold_change)

      if (all(id_df$geneID %in% fold_change$geneID)) {
        m1 <- merge(id_df, fold_change, by = "geneID")
        plot_df <- merge(plot_df, m1, by.x = "geneID", by.y = "geneID_symbol") %>%
          dplyr::select(-geneID.y)
      } else {
        plot_df <- merge(plot_df, fold_change, by.x = "geneID")
      }

      p <- ggplot(plot_df, aes_(~geneID, ~Description)) +
        geom_tile(aes_(fill = ~logfc), color = "white") +
        xlab(NULL) +
        ylab(NULL) +
        plot_theme(border_thick = 0, ...) +
        scale_fill_continuous(
          low = down_color, high = up_color,
          name = "logFC"
        ) +
        theme(
          panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 50, hjust = 1)
        )
    }
  }

  #--- GO/KEGG: genechord plot ---#
  ## chord plot to show interaction between go term and gene id
  if (plot_type == "genechord") {
    if (!requireNamespace("GOplot", quietly = TRUE)) {
      utils::install.packages("GOplot")
    }

    id <- enrich_df %>%
      dplyr::pull(geneID) %>%
      stringr::str_split("\\/") %>%
      unlist()
    id_symbol <- enrich_df %>%
      dplyr::pull(geneID_symbol) %>%
      stringr::str_split("\\/") %>%
      unlist()
    id_df <- data.frame(geneID = id, geneID_symbol = id_symbol) %>%
      dplyr::distinct()

    # if show_gene is not symbol, first extract matching symbol
    if(all(show_gene == 'all')) stop('Please specify gene name to "show_gene"...')
    if (length(show_gene %in% id) > length(show_gene %in% id_symbol)) {
      show_gene <- id_df %>%
        dplyr::filter(geneID %in% show_gene) %>%
        dplyr::pull(geneID_symbol)
    }else{
      show_gene <- id_df %>%
        dplyr::filter(geneID_symbol %in% show_gene) %>%
        dplyr::pull(geneID_symbol)
    }


    plot_df <- enrich_df %>%
      dplyr::select(Description, geneID_symbol) %>%
      tidyr::separate_rows(geneID_symbol, sep = "\\/") %>%
      dplyr::filter(geneID_symbol %in% show_gene) %>%
      dplyr::rename(geneID = geneID_symbol)

    # define color
    my_cols <- c(
      "#B2DF8A", "#FB9A99", "#E31A1C", "#B15928", "#6A3D9A", "#CAB2D6",
      "#A6CEE3", "#1F78B4", "#FDBF6F", "#999999", "#FF7F00"
    )
    backup_cols <- c(
      "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D", "#58CDD9",
      "#7A142C", "#5D90BA", "#029149", "#431A3D", "#91612D", "#6E568C",
      "#E0367A", "#D8D155", "#64495D", "#7CC767"
    )

    term <- unique(plot_df$Description)
    if (length(my_cols) < length(term)) {
      my_cols <- c(my_cols, backup_cols)
      cols <- sample(my_cols, length(term), replace = F)
    }
    cols <- my_cols[1:length(term)]

    # prepare chord data
    id <- unique(plot_df$geneID)
    dat <- sapply(id, function(x) {
      check <- plot_df %>%
        dplyr::filter(geneID %in% x) %>%
        dplyr::pull(Description) %>%
        unique()
      ifelse(term %in% check, 1, 0)
    }) %>%
      t() %>%
      as.data.frame() %>%
      stats::setNames(term)
    dat = dat[show_gene,]

    # default main and lengend text size
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 3
    if (!"legend_text_size" %in% names(lst)) lst$legend_text_size <- 8
    if(!"gene_space" %in% names(lst)) lst$gene_space <- 0.3

    if (is.null(fold_change)) {
      p <- suppressWarnings(
        GOplot::GOChord(dat,
          space = 0.02,
          gene.space = lst$gene_space,
          gene.order = "none",
          gene.size = lst$main_text_size,
          process.label = lst$legend_text_size,
          border.size = 0.1,
          ribbon.col = cols
        )
      )
    } else {
      # add logfc
      fold_change <- data.frame(geneID = names(fold_change), logfc = fold_change)

      if (all(rownames(dat) %in% fold_change$geneID)) {
        m1 <- fold_change %>%
          dplyr::filter(geneID %in% rownames(dat)) %>%
          dplyr::arrange(match(geneID_symbol, rownames(dat))) %>%
          dplyr::pull(logfc)
      } else {
        m1 <- merge(id_df, fold_change, by.x = "geneID") %>%
          dplyr::filter(geneID_symbol %in% rownames(dat)) %>%
          dplyr::arrange(match(geneID_symbol, rownames(dat))) %>%
          dplyr::pull(logfc)
      }
      dat <- dat %>% dplyr::mutate(logFC = m1)
      p <- suppressWarnings(
        GOplot::GOChord(dat,
          space = 0.02,
          gene.space = lst$gene_space,
          gene.order = "logFC",
          gene.size = lst$main_text_size,
          process.label = lst$legend_text_size,
          border.size = 0.1,
          ribbon.col = cols,
          lfc.col = c( up_color,"grey50",down_color )
        )
      )
    }
    # p <- p+ theme(legend.title=element_text(size=lst$legend_text_size))
  }

  #--- GO/KEGG: network plot ---#
  ## plot enriched terms in network with edges connecting overlapping gene sets,
  ## mutually overlapping gene sets are tend to cluster together
  if (plot_type == "network") {
    if (missing(scale_ratio)) scale_ratio <- 1
    pkgs <- c("ggnewscale", "ggraph", "igraph")
    invisible(sapply(pkgs, function(x) {
      if (!requireNamespace(x, quietly = TRUE)) {
        utils::install.packages(x)
      }
    }))


    id <- enrich_df[,grepl("^id",colnames(enrich_df),ignore.case = T)]
    enrichGenes <- strsplit(enrich_df$geneID, "\\/") %>% setNames(id)

    # kegg result just use JC method
    if (!any(grepl("GO:", id))) sim_method <- "JC"

    if (sim_method == "JC") {
      m <- get_JC_data(enrich_df)
    } else {
      m <- get_sim_data(enrich_df, org = NULL, ont = NULL, sim_method)[["m"]]
      rownames(m) <- colnames(m) <- enrich_df %>%
        dplyr::filter(.[[1]] %in% colnames(m)) %>%
        dplyr::pull(Description)
    }

    mm <- reshape2::melt(m)
    mm <- mm[mm[, 1] != mm[, 2], ]
    mm <- mm[!is.na(mm[, 3]), ]
    # construct igraph
    g <- graph.data.frame(mm[, -3], directed = FALSE)
    E(g)$width <- sqrt(mm[, 3] * 5) * scale_ratio
    E(g)$weight <- mm[, 3]
    g <- delete.edges(g, E(g)[mm[, 3] < 0.2])
    id_order <- unlist(sapply(V(g)$name, function(x) which(x == enrich_df$Description)))
    id_genes <- sapply(enrichGenes[id_order], length)
    V(g)$size <- id_genes
    V(g)$color <- enrich_df[id_order, stats_metric]

    # default main and lengend text size
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 3
    if (!"legend_text_size" %in% names(lst)) lst$legend_text_size <- 8

    # igraph to ggplot
    p <- ggraph(g, layout) +
      geom_edge_link(
        alpha = .8, aes_(width = ~ I(width)),
        colour = "darkgrey"
      ) +
      ggnewscale::new_scale_fill() +
      geom_point(shape = 21, aes_(
        x = ~x, y = ~y,
        fill = ~color,
        size = ~size
      )) +
      scale_size_continuous(
        name = "Number of genes",
        guide = "legend",
        range = c(min(V(g)$size) / 2, max(V(g)$size) / 2) * scale_ratio
      ) +
      scale_fill_continuous(
        low = up_color, high = down_color,
        name = stats_metric_label
      ) +
      theme(panel.background = element_blank()) +
      geom_node_text(aes_(label = ~name),
        data = NULL,
        size = lst$main_text_size,
        bg.color = "white",
        repel = TRUE, segment.size = 0.2
      ) +
      guides(fill = guide_colorbar(reverse = TRUE)) +
      plot_theme(
        remove_border = T, remove_main_text = T,
        border_thick = 0, ...
      )
  }

  #--- GO/KEGG: wordcloud ---#
  if (plot_type == 'wordcloud'){
    if (!requireNamespace("wordcloud", quietly = TRUE)) {
      utils::install.packages("wordcloud")
    }

    mypal <- RColorBrewer::brewer.pal(9, "Paired")
    dat <- data.frame(term = enrich_df$Description)
    x <- tm::Corpus(tm::VectorSource(dat$term))
    tdm <- tm::TermDocumentMatrix(x, control = list(removePunctuation = TRUE,
                                                    stopwords = TRUE))
    m <- as.matrix(tdm)
    v <- sort(rowSums(m), decreasing = TRUE)
    d <- data.frame(word = names(v), freq = v)
    p <- invisible(wordcloud::wordcloud(d$word, d$freq, min.freq = 0,colors = mypal))
  }

  #--- GO/KEGG: upset plot ---#
  if (plot_type == "upset") {
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 10
    if (!"legend_text_size" %in% names(lst)) lst$legend_text_size <- 8
    if (!"legend_position" %in% names(lst)) lst$legend_position <- 'left'

    plot_df <-  enrich_df %>% as.upset() %>%
      do.call(cbind,.) %>% as.data.frame()

    p <- plotVenn(enrich_df, use_venn = FALSE,...)
  }

  #--- GO: wego plot ---#
  ## WEGO plot show all ontologies
  if (plot_type == "wego") {
    if (!"ontology" %in% tolower(colnames(enrich_df))) {
      stop("WEGO plot needs a column named 'ontology' which infers ontology type of 'bp', 'cc' or 'mf'...")
    }else{
      colnames(enrich_df)[grep("ontology", colnames(enrich_df), ignore.case = T)] <- "Ontology"
    }

    if(missing(n_term)) n_term = 5
    if(is.null(wrap_length)) wrap_length = 30
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 8

    wego <- enrich_df %>%
      dplyr::select(1,'Description', "Count", "GeneRatio", "Ontology") %>%
      dplyr::mutate(GeneRatio = GeneRatio * 100) %>%
      dplyr::group_by(Ontology) %>%
      dplyr::top_n(n_term, GeneRatio) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Ontology, GeneRatio) %>%
      dplyr::mutate(Position = dplyr::n():1) %>%
      dplyr::mutate(Ontology = dplyr::case_when(
        tolower(Ontology) == "bp" ~ "Biological Process",
        tolower(Ontology) == "cc" ~ "Cellular Component",
        tolower(Ontology) == "mf" ~ "Molecular Function"
      ))

    normalizer <- max(wego$Count) / max(wego$GeneRatio)

    p <- ggplot(data = wego, aes(
      x = forcats::fct_reorder(Description, sort(Position)),
      y = GeneRatio, fill = Ontology
    )) +
      ggsci::scale_fill_nejm() +
      geom_col(data = wego, aes(
        x = forcats::fct_reorder(Description, sort(Position)),
        y = Count / normalizer
      )) +
      scale_y_continuous(sec.axis = sec_axis(
        trans = ~ . * normalizer, name = "Number of genes",
        labels = function(b) {
          round(b, 0)
        }
      )) +
      plot_theme(remove_legend = T, ...) +
      scale_x_discrete(labels = text_wraper(wrap_length)) +
      xlab(NULL) +
      ylab("Gene Ratio(%)") +
      facet_grid(. ~ Ontology, scales = "free") +
      theme(
        axis.text.x = element_text(angle = 70, hjust = 1),
        strip.text.x = element_text(size = lst$main_text_size)
      )
  }

  #--- GO: term map plot ---#
  ## show enriched terms structure (its parents and children terms)
  if (plot_type == "gomap") {
    if (missing(scale_ratio)) scale_ratio <- 1
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      utils::install.packages("ggraph")
    }
    if (!requireNamespace("igraph", quietly = TRUE)) {
      utils::install.packages("igraph")
    }

    requireNamespace("ggraph", quietly = T)
    requireNamespace("igraph", quietly = T)

    rownames(enrich_df) = enrich_df[,1]
    id <- enrich_df[, 1]
    enrichGenes <- strsplit(enrich_df$geneID, "\\/") %>% setNames(id)

    # get all id ancestors
    GOANCESTOR <- get_ancestor_data(enrich_df, ont = NULL)
    id_anc <- AnnotationDbi::mget(id, GOANCESTOR)
    anc1 <- id_anc[[1]]
    for (i in 2:length(id_anc)) {
      anc1 <- intersect(anc1, id_anc[[i]])
    }
    uanc <- unique(unlist(id_anc))
    uanc <- uanc[!uanc %in% anc1]

    # get plot edge and nodes data
    gotbl <- get_gosim_data()
    edge <- gotbl[gotbl$go_id %in% unique(c(id, uanc)), ] %>%
      dplyr::select(c(5, 1, 4))
    node <- gotbl %>%
      dplyr::filter(go_id %in% unique(c(edge[, 1], edge[, 2]))) %>%
      dplyr::select(1:3) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        color = enrich_df[go_id, stats_metric],
        size = sapply(enrichGenes[go_id], length)
      )
    rm(gotbl, envir = .genekitrEnv)

    # only show top nodes/id-related parent and child nodes
    show_node_top <- table(edge$parent) %>%
      sort() %>%
      utils::tail(3) %>%
      names()
    show_node_parent <- edge %>%
      dplyr::filter(go_id %in% id) %>%
      dplyr::pull(parent)
    show_node_child <- edge %>%
      dplyr::filter(parent %in% id) %>%
      dplyr::pull(go_id)
    show_nodes <- unique(c(id, show_node_top, show_node_parent, show_node_child))
    sub_node <- node %>%
      dplyr::mutate(Term = ifelse(go_id %in% show_nodes,
        stringr::str_replace(Term, "^\\w{1}", toupper), NA
      )) %>%
      dplyr::mutate(Term = ifelse(is.na(Term), NA,
        sapply(strwrap(Term, width = wrap_length, simplify = FALSE),
          paste,
          collapse = "\n"
        )
      ))


    g <- igraph::graph.data.frame(edge, directed = TRUE, vertices = sub_node)
    E(g)$Relationship <- edge[, 3]

    # default main and lengend text size
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 3
    if (!"legend_text_size" %in% names(lst)) lst$legend_text_size <- 8

    p <- ggraph(g, layout = "sugiyama") +
      geom_edge_link(aes_(linetype = ~Relationship),
        arrow = grid::arrow(length = unit(1, "mm")),
        end_cap = circle(1, "mm"),
        colour = "darkgrey"
      ) +
      geom_node_point(size = 3 * scale_ratio, aes_(color = ~color)) +
      scale_color_continuous(
        low = up_color, high = down_color, name = stats_metric_label,
        guide = guide_colorbar(reverse = TRUE)
      ) +
      geom_node_label(aes_(label = ~Term, color = ~color),
        size = lst$main_text_size,
        repel = TRUE, segment.size = 0.2,
        max.overlaps = 16
      ) +
      scale_fill_continuous(
        low = up_color, high = down_color, name = stats_metric_label,
        guide = guide_colorbar(reverse = TRUE), na.value = "white"
      ) +
      plot_theme(border_thick = 0, remove_main_text = T, ...)
  }

  #--- GO: term heatmap---#
  if(plot_type == "goheat"){
    if (!sim_method %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
      stop('Please choose "sim_method" from: "Resnik", "Lin", "Rel", "Jiang" , "Wang"!')
    }
    l <- get_sim_data(enrich_df, org = NULL, ont = NULL, sim_method)
    simMatrix <- l[["m"]]
    reducedTerms <- l[["r"]]

    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 10
    if (!"legend_text_size" %in% names(lst)) lst$legend_text_size <- 8
    if (!"remove_legend" %in% names(lst)) lst$remove_legend <- FALSE
    my_cols <- c(
      "#B2DF8A", "#FB9A99", "#E31A1C", "#B15928", "#6A3D9A", "#CAB2D6",
      "#A6CEE3", "#1F78B4", "#FDBF6F", "#999999", "#FF7F00",
      "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D", "#58CDD9",
      "#7A142C", "#5D90BA", "#029149", "#431A3D", "#91612D", "#6E568C",
      "#E0367A", "#D8D155", "#64495D", "#7CC767"
    )

    anno <- data.frame(ParentTerm = factor(reducedTerms[match(rownames(simMatrix),
                                                             reducedTerms$go), "parentTerm"]),
                      row.names = rownames(simMatrix))
    anno_col <- list(ParentTerm = my_cols[1:length(unique(anno$ParentTerm))])
    anno$ParentTerm <- factor(stringr::str_replace(anno$ParentTerm, "^\\w{1}", toupper))
    names(anno_col$ParentTerm) <- levels(anno$ParentTerm)


    p <- pheatmap::pheatmap(simMatrix,
                            # color = colorRampPalette(c("#1F78B4", "white", "#E31A1C"))(50),
                            annotation_row = anno,
                            annotation_colors = anno_col,
                            annotation_names_row = F,
                            show_colnames = F,
                            treeheight_col = 0,
                            fontsize = lst$legend_text_size,
                            fontsize_row = lst$main_text_size,
                            legend = !lst$remove_legend
                            # treeheight_row = 0
                            )
  }

  #--- GO: term tangram ---#
  if (plot_type == "gotangram") {
    if (!sim_method %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
      stop('Please choose "sim_method" from: "Resnik", "Lin", "Rel", "Jiang" , "Wang"!')
    }
    my_cols <- c(
      "#B2DF8A", "#FB9A99", "#E31A1C", "#B15928", "#6A3D9A", "#CAB2D6",
      "#A6CEE3", "#1F78B4", "#FDBF6F", "#999999", "#FF7F00",
      "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D", "#58CDD9",
      "#7A142C", "#5D90BA", "#029149", "#431A3D", "#91612D", "#6E568C",
      "#E0367A", "#D8D155", "#64495D", "#7CC767"
    )

    l <- get_sim_data(enrich_df, org = NULL, ont = NULL, sim_method)
    simMatrix <- l[["m"]]
    reducedTerms <- l[["r"]]

    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 8

    p <- treemap::treemap(reducedTerms,
                          index = c("parentTerm", "term"),
                          vSize = "score", type = "index", title = "",
                          palette = RColorBrewer::brewer.pal(length(unique(reducedTerms$parent)),'Set2'),
                          fontsize.labels = lst$main_text_size,
                          fontcolor.labels = c("#FFFFFFDD", "#00000080"),
                          bg.labels = 0,
                          border.col = "#00000080")


  }


  #--- keggpath ---#
  # if(plot_type == 'keggpath'){
  #   # if (!requireNamespace("pathview", quietly = TRUE)) {
  #   #   warning('Depends on pathview package. Installing...')
  #   #   utils::install.packages('pathview')
  #   # }
  #
  #   if(is.null(fold_change))
  #     stop('Please give fold change or logFC values with gene IDs as names!')
  #
  #   org = substr(enrich_df[1,1],1,3) #e.g. hsa
  #   ids = enrich_df$ID
  #   if(!all(grepl("^[a-z]{3}.*",ids))) stop('Please give a kegg result...')
  #
  #   max_fc <- max(abs(fold_change))
  #   bins <- ceiling(max_fc) * 2
  #   p <- lapply(ids, function(i) {
  #     print(paste0('Now plotting ',which(ids%in%i),'/',length(ids),': ',i))
  #     suppressPackageStartupMessages(pathview::pathview(gene.data=fold_change,
  #                                         pathway.id = i,
  #                                         species = org,
  #                                         limit = list(gene=max_fc, cpd=1),
  #                                         bins = list(gene=bins, cpd=10),
  #                                         low = list(gene="blue", cpd="blue"),
  #                                         high = list(gene="red", cpd="yellow"),
  #                                         out.suffix='genekitr',
  #                                         kegg.native=TRUE,
  #                                         new.signature=FALSE))
  #   })
  #   invisible(p)
  # }


  # wrap long text
  if (!is.null(wrap_length) & is.numeric(wrap_length) & plot_type != 'wego') {
    p <- p + scale_y_discrete(labels = text_wraper(wrap_length))
  }

  if(plot_type != "gotangram") return(p)

}


#--- sub-function---#
text_wraper <- function(width) {
  function(x) {
    lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n")
  }
}

ratio_intersect <- function(n, m) {
  n <- unlist(n)
  m <- unlist(m)
  length(intersect(n, m)) / length(unique(c(n, m)))
}

get_JC_data <- function(enrich_df) {
  id <- enrich_df[, 1]
  enrichGenes <- strsplit(enrich_df$geneID, "\\/") %>% setNames(id)

  n <- nrow(enrich_df)
  m <- matrix(NA, nrow = n, ncol = n)
  colnames(m) <- rownames(m) <- enrich_df$Description
  for (i in seq_len(n - 1)) {
    for (x in (i + 1):n) {
      m[i, x] <- ratio_intersect(enrichGenes[id[i]], enrichGenes[id[x]])
    }
  }
  return(m)
}

get_sim_data <- function(enrich_df, org = NULL, ont = NULL, sim_method) {
  org_name <- genekitr::biocOrg_name$short_name

  if (is.null(org)) {
    tryCatch(
      {
        nm <- strsplit(colnames(enrich_df)[1], "_") %>% unlist()
        org <- nm[1]
      },
      error = function(e) {
        message("Please rename the ID column with organism short name from genekitr::biocOrg_name!\nOR you can use 'JC' method")
      }
    )
  } else if (!tolower(org) %in% org_name) {
    stop("Please rename the ID column with organism short name from genekitr::biocOrg_name!\nOR you can use 'JC' method.")
  }

  if (is.null(ont)) {
    tryCatch(
      {
        nm <- strsplit(colnames(enrich_df)[1], "_") %>% unlist()
        ont <- nm[2]
      },
      error = function(e) {
        message('Please rename the ID column with ontology of "BP", "CC" and "MF"!\nOR you can use "JC" method.')
      }
    )
  } else if (!tolower(ont) %in% c("bp", "cc", "mf")) {
    stop('Please rename the ID column with ontology of "BP", "CC" and "MF"!\nOR you can use "JC" method.')
  }

  ont <- toupper(ont)

  #   if(all(tolower(nm) == 'id')){
  #     stop(paste0('Please give organism name to "org" such as "Hs", "Mm"...\n',
  #                 'also give ontology to "ont" from "BP","CC" and "MF"'))
  #   }else{
  #     org =  nm[1]; ont = nm[2]
  #   }

  if(tolower(org) == 'at'){
    orgdb <- 'org.At.tair.db'
  } else if (tolower(org) == 'sc') {
    orgdb <- 'org.Sc.sgd.db'
  }else{
    orgdb <- paste0("org.", stringr::str_to_title(org), ".eg.db")
  }

  id <- enrich_df[, 1]
  # save godata is saving time
  data_dir <- tools::R_user_dir("genekitr", which = "data")
  data_dir <- paste0(data_dir, "/godata")
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
  destfile <- paste0(data_dir, "/", org, "_", ont, "_godata.rda")

  if (!file.exists(destfile)) {
    orgdb <- loadOrgdb(orgdb)
    semdata <- GOSemSim::godata(orgdb, ont = ont)
    save(semdata, file = destfile)
  } else {
    load(destfile)
  }
  # calculate simMatrix(m)
  id <- unique(id)
  found <- id %in% names(semdata@IC)
  getAncestors <- utils::getFromNamespace("getAncestors", "GOSemSim")
  hasAncestor <- !is.na(sapply(id, function(x) {
    tryCatch(getAncestors(ont)[x],
      error = function(e) NA
    )
  }))
  id <- id[found & hasAncestor]
  m <- suppressMessages(matrix(GOSemSim::goSim(id, id, semData = semdata, measure = sim_method),
    ncol = length(id), dimnames = list(id, id)
  ))
  out <- apply(m, 2, function(x) all(is.na(x)))
  m <- m[!out, !out]
  # reduce redundant terms
  if(all(enrich_df$qvalue%in%NA)){
    scores <- setNames(-log10(enrich_df$p.adjust), enrich_df[, 1])
  }else{
    scores <- setNames(-log10(enrich_df$qvalue), enrich_df[, 1])
  }

  r <- rrvgo::reduceSimMatrix(m, scores, orgdb = orgdb)
  return(list(m = m, r = r))
}

get_ancestor_data <- function(enrich_df, ont = NULL) {
  nm <- strsplit(colnames(enrich_df)[1], "_") %>% unlist()
  if (all(tolower(nm) == "id")) {
    stop(paste0('Please give ontology to "ont" from "BP","CC" and "MF"'))
  } else {
    ont <- nm[2]
  }

  # save ancestor_data is saving time
  data_dir <- tools::R_user_dir("genekitr", which = "data")
  data_dir <- paste0(data_dir, "/ancestor_data")
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
  destfile <- paste0(data_dir, "/", ont, "_ancestor.rda")
  getAncestors <- utils::getFromNamespace("getAncestors", "GOSemSim")
  GOANCESTOR <- getAncestors(ont)
  return(GOANCESTOR)
}

loadOrgdb <- function(orgdb) {
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    stop("Bioconductor orgdb for ", orgdb, " not found. You should install first.",
      call. = FALSE
    )
  }
  eval(parse(text = paste0(orgdb, "::", orgdb)))
}

get_gosim_data <- function(){
  .initial()
  utils::data("gotbl", package = "GOSemSim",envir = .genekitrEnv)
  gotbl <- get("gotbl", envir = .genekitrEnv)
}

