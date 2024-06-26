#' Plot for gene enrichment analysis of GSEA method
#'
#' Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of
#' genes shows statistically significant, concordant differences between two biological states (e.g. phenotypes).
#' @param gsea_list GSEA result from `genGSEA` function
#' @param plot_type GSEA plot type, one of 'volcano', 'classic', 'fgsea', 'ridge' or 'bar'.
#' @param stats_metric Statistic metric from one of "pvalue", "p.adjust", "qvalue".
#' @param show_pathway Select plotting pathways by number (will choose top N pathways)
#' or pathway name (choose from ID column).
#' @param show_gene Select genes to show. Default is "all". Used in "classic" plot.
#' @param colour Colour vector. Deafault is NULL. Used in volcano, ridge and bar plot.
#' @param wrap_length Numeric, wrap text if longer than this length. Default is NULL.
#' @param label_by Select which column as the label. If user wants to modify labels in plot,
#' please modify the "Description" column and set the argument as "description". Default is by 'id'.
#' @param ... other arguments transfer to `plot_theme` function
#'
#' @importFrom ggplot2 ggplot aes aes_ geom_point xlab ylab scale_color_manual geom_hline expansion
#' element_blank element_rect margin geom_linerange scale_y_continuous geom_segment geom_bar
#' scale_fill_manual geom_hline element_line geom_line scale_colour_manual guide_colourbar
#' @importFrom dplyr filter arrange slice_head select rename group_by summarise case_when mutate
#' pull all_of
#'
#' @return A ggplot object
#' @export
#' @examples
#' \donttest{
#' k1 = requireNamespace("cowplot",quietly = TRUE)
#' k2 = requireNamespace("fgsea",quietly = TRUE)
#' k3 = requireNamespace("ggplotify",quietly = TRUE)
#' k4 = requireNamespace("ggridges",quietly = TRUE)
#' if(k1&k2&k3&k4){
#' library(ggplot2)
#' ## get GSEA result
#' data(geneList, package = "genekitr")
#' gs <- geneset::getMsigdb(org = "human",category = "H")
#' gse <- genGSEA(genelist = geneList, geneset = gs)
#'
#' ## volcano plot
#' # get top3 of up and down pathways
#' plotGSEA(gse, plot_type = "volcano", show_pathway = 3)
#' # choose pathway by character
#' pathways <- c('HALLMARK_KRAS_SIGNALING_UP','HALLMARK_P53_PATHWAY','HALLMARK_GLYCOLYSIS')
#' plotGSEA(gse, plot_type = "volcano", show_pathway = pathways)
#'
#' ## classic pathway plot
#' genes <- c('ENG','TP53','MET')
#' plotGSEA(gse, plot_type = "classic", show_pathway = pathways, show_gene = genes)
#'
#' ## fgsea table plot
#' plotGSEA(gse, plot_type = "fgsea", show_pathway = 3)
#'
#' ## ridgeplot
#' plotGSEA(gse,
#'   plot_type = "ridge",
#'   show_pathway = 10, stats_metric = "p.adjust"
#' )
#'
#' ## two-side barplot
#' plotGSEA(gse,
#'   plot_type = "bar", main_text_size = 8,
#'   colour = c("navyblue", "orange")
#' )
#' }
#' }
#'
plotGSEA <- function(gsea_list,
                     plot_type = c("volcano", "classic", "fgsea", "ridge", "bar"),
                     stats_metric = c("p.adjust", "pvalue", "qvalue"),
                     show_pathway = NULL,
                     show_gene = NULL,
                     colour = NULL,
                     wrap_length = NULL,
                     label_by = c('id','description'),
                     ...) {

  #--- args ---#
  lst <- list(...) # store outside arguments in list
  stopifnot(is.list(gsea_list))
  plot_type <- match.arg(plot_type)
  stats_metric <- match.arg(stats_metric)
  label_by <- match.arg(label_by)
  gsea_df <- gsea_list$gsea_df
  colnames(gsea_df)[1] = "ID"
  geneset <- gsea_list$geneset
  gl <- gsea_list$genelist

  legends <- c("p.adjust", "pvalue", "qvalue")
  if (!stats_metric %in% colnames(gsea_list$gsea_df)) {
    stop(
      stats_metric, " not included in this dataframe, try: ",
      paste(intersect(colnames(gsea_list$gsea_df), legends), collapse = " | ")
    )
  }

  if(!'geneID_symbol' %in% colnames(gsea_df)){
    gsea_df <- gsea_df %>% dplyr::mutate(geneID_symbol = geneID)
  }


  #--- codes ---#
  ## set labels
  stats_metric_label <- ifelse(stats_metric == "pvalue", "Pvalue",
    ifelse(stats_metric == "p.adjust", "P.adjust", "FDR")
  )

  #--- volcano plot ---#
  # x-axis: NES
  # y-axis: stats value
  if (plot_type == "volcano") {
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 8
    if(is.null(show_pathway)) show_pathway = 3

    gsea_df <- gsea_df[order(gsea_df$NES, decreasing = TRUE), ]
    plot_df <- gsea_df
    if (is.numeric(show_pathway)) {
      plot_df$group <- c(rep("Up", show_pathway), rep("ignore", nrow(plot_df) - 2 * show_pathway), rep("Down", show_pathway))
    } else {
      nes <- plot_df[plot_df$ID %in% show_pathway, "NES"]
      plot_df$group <- "ignore"
      plot_df[plot_df$ID %in% show_pathway, "group"] <- sapply(nes, function(x) ifelse(x > 0, "Up", "Down"))
    }

    if (is.null(colour)) {
      message('"colour" is NULL, now using "red","grey" and "green"...')
      colour <- c( "red","darkgreen", "grey66")
    }

    if(label_by == 'id'){
      p <- ggplot(plot_df, aes(x = NES, y = -log10(eval(parse(text = stats_metric))), color = group)) +
        geom_point(alpha = 0.6, size = 3.5) +
        xlab("NES") +
        ylab(paste0("-log10(", stats_metric_label, ")")) +
        ggrepel::geom_text_repel(
          data = plot_df[plot_df$group != "ignore", ],
          aes(label = ID),
          size = (lst$main_text_size/4),
          color = "black",
          show.legend = F
        )
    }else{
      p <- ggplot(plot_df, aes(x = NES, y = -log10(eval(parse(text = stats_metric))), color = group)) +
        geom_point(alpha = 0.6, size = 3.5) +
        xlab("NES") +
        ylab(paste0("-log10(", stats_metric_label, ")")) +
        ggrepel::geom_text_repel(
          data = plot_df[plot_df$group != "ignore", ],
          aes(label = Description),
          size = (lst$main_text_size/4),
          color = "black",
          show.legend = F
        )
    }


    p = p +
      xlab("Normalized enrichment score")+
      ylim(0,NA)+
      scale_color_manual(breaks = c("Up", "Down"),values = colour,name = '')+
      plot_theme(...)
  }

  #--- classic plot ---#
  # GSEA rank list
  if (plot_type == "classic") {
    gl1 = gl[,2]
    names(gl1) = gl[,1]
    genelist = gl1
    exponent <- as.numeric(gsea_list$exponent)
    org <- as.character(gsea_list$org)

    if(is.null(show_pathway)) show_pathway = 3

    if (is.numeric(show_pathway)) {
      show_pathway <- gsea_df$ID[show_pathway]
    } else if (any(!show_pathway %in% gsea_df$ID)) {
      stop(paste0(show_pathway[!show_pathway %in% gsea_df$ID], " not in GSEA result!"))
    }

    if (is.null(colour)) {
      colour <- c(
        "\\#5DA5DAFF", "\\#FAA43AFF", "\\#60BD68FF", "\\#F15854FF", "\\#B276B2FF",
        "\\#8D4B08FF", "\\#DECF3FFF", "\\#F17CB0FF", "\\#66E3D9FF", "\\#00FF7FFF",
        "\\#E31A1CFF", "\\#FFFF99FF"
      )
      colour <- stringr::str_remove_all(colour, ".*#") %>% paste0("#", .)
    }

    plot_df <- do.call(rbind, lapply(show_pathway, function(x) {
      calcScore(geneset, genelist, x, exponent, fortify = TRUE, org)
    }))

    if(label_by != 'id'){
      gsea_df <- gsea_list$gsea_df
      colnames(gsea_df)[1] = "ID"

      plot_df = merge(plot_df,gsea_df,by.x = 'Description',by.y = 'ID',all.y = F) %>%
        dplyr::select(2:9) %>%
        dplyr::rename(Description = Description.y)

    }

    description_color <- table(plot_df$Description) %>% names() # match pathway and color
    names(description_color) <- colour[seq_along(description_color)]


    p1 <- ggplot(plot_df, aes_(x = ~x)) +
      xlab(NULL) +
      geom_line(aes_(y = ~runningScore, color = ~Description), linewidth = 1) +
      scale_color_manual(values = names(description_color)) +
      geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
      ylab("Enrichment\n Score") +
      plot_theme(remove_grid = T, ...) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
      )

    # size for three panels
    panel_size <- c(1.5, .5, 1.5)

    i <- 0
    for (term in unique(plot_df$Description)) {
      idx <- which(plot_df$ymin != 0 & plot_df$Description == term)
      plot_df[idx, "ymin"] <- i
      plot_df[idx, "ymax"] <- i + 1
      i <- i + 1
    }

    # p2 <- ggplot(plot_df, aes_(x = ~x)) +
    #   geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    #   xlab(NULL) + ylab(NULL) +
    #   scale_color_manual(values = names(description_color))+
    #   plot_theme(remove_legend = T,remove_grid = T,remove_main_text = T,
    #              remove_border = T)+
    #   scale_y_continuous(expand=c(0,0))+
    #   theme(
    #     plot.margin = margin(t=-.1, b=0,unit="cm"),
    #     axis.ticks = element_blank(),
    #     axis.line.y = element_line(size = 1))+
    #   annotate(geom = 'segment', y = Inf, yend = -Inf,  x = Inf, xend = Inf,
    #            size = lst$border_thick)

    p2 <- ggplot(plot_df, aes_(x = ~x)) +
      geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) +
      xlab(NULL) +
      ylab(NULL) +
      scale_color_manual(values = colour) +
      plot_theme(remove_legend = T, remove_grid = T, remove_main_text = T, ...) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(
        plot.margin = margin(t = -.1, b = 0, unit = "cm"),
        axis.line.x = element_blank(),
        axis.ticks = element_blank()
      )


    plot_df$y <- plot_df$genelist

    if (!is.null(show_gene)) {
      select_genes <- data.frame(gene = show_gene)
      select_genes <- merge(select_genes, plot_df, by = "gene")
      select_genes <- select_genes[select_genes$position == 1, ]
    } else {
      select_genes <- data.frame(gene = "no_select_genes")
      select_genes <- merge(select_genes, plot_df, by = "gene")
      select_genes <- select_genes[select_genes$position == 1, ]
    }

    gene_color <- names(description_color)[description_color %in% select_genes$Description]
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 8

    p3 <- ggplot(plot_df) +
      geom_segment(aes_(x = ~x, xend = ~x, y = ~y, yend = 0),
        color = "grey"
      ) +
      geom_bar(
        data = select_genes,
        aes(x = x, y = y, fill = Description, color = Description),
        position = "dodge", stat = "identity", width = 0.5
      ) +
      scale_fill_manual(values = gene_color, guide = "none") +
      scale_color_manual(values = gene_color, guide = "none") +
      geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
      ylab("Ranked list\n metric") +
      xlab("Rank in ordered dataset") +
      plot_theme(remove_grid = T, ...) +
      # show gene names
      ggrepel::geom_text_repel(
        data = select_genes,
        aes(x = x, y = y, label = gene, color = Description),
        show.legend = FALSE,
        max.overlaps = Inf,
        direction = "x",
        ylim = c(2, NA),
        angle = 90,
        size = lst$main_text_size / 3.5, box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) +
      theme(plot.margin = margin(t = -.1, r = .2, b = .2, l = .2, unit = "cm"))

    pl <- list(p1, p2, p3)
    n <- length(pl)
    pl[[n]] <- pl[[n]] +
      theme(
        axis.line.x = element_line(),
        axis.ticks.x = element_line()
      )

    p <- cowplot::plot_grid(plotlist = pl, ncol = 1, align = "v", rel_heights = panel_size)
  }

  #--- fgsea plot ---#
  if (plot_type == "fgsea") {
    if(is.null(show_pathway)) show_pathway = 3

    geneset_list <- geneset %>%
      split(.[[1]]) %>%
      lapply("[[", 2)
    gl1 = gl[,2]
    names(gl1) = gl[,1]
    genelist = gl1

    fres <- suppressWarnings(fgsea::fgsea(
      pathways = geneset_list,
      stats = genelist
    ))

    if (is.numeric(show_pathway)) {
      up_path <- fres %>%
        dplyr::filter(ES > 0) %>%
        dplyr::arrange(., pval) %>%
        dplyr::slice_head(n = show_pathway) %>%
        dplyr::pull(pathway)
      down_path <- fres %>%
        dplyr::filter(ES < 0) %>%
        dplyr::arrange(., pval) %>%
        dplyr::slice_head(n = show_pathway) %>%
        dplyr::pull(pathway)
      all <- c(up_path, down_path)
    } else if (any(!show_pathway %in% fres$pathway)) {
      stop(paste0(show_pathway[!show_pathway %in% fres$pathway], " not in fgsea result!"))
    } else {
      all <- show_pathway
    }

    p <- fgsea::plotGseaTable(geneset_list[all], genelist, fres, render = F) %>%
      ggplotify::as.ggplot()
  }

  #--- ridge plot ---#
  if (plot_type == "ridge") {
    if(is.null(show_pathway)) show_pathway = 3

    if (is.numeric(show_pathway)) {
      if (length(show_pathway) > 1) {
        show_pathway <- gsea_df$ID[show_pathway]
      } else {
        show_pathway <- gsea_df$ID[1:show_pathway]
      }
    }

    check_gset_type <- suppressWarnings(ifelse(all(is.na(as.numeric(dplyr::slice_head(geneset,n=2) %>% dplyr::pull(2)))),'symbol','entrez'))
    if(check_gset_type == 'symbol'){
      new_gsea_df <- gsea_df %>%
        dplyr::filter(ID %in% show_pathway) %>%
        dplyr::select(ID, dplyr::all_of(stats_metric), geneID_symbol) %>%
        tidyr::separate_rows(geneID_symbol, sep = "\\/") %>%
        dplyr::rename(geneID = geneID_symbol)
    }else{
      new_gsea_df <- gsea_df %>%
        dplyr::filter(ID %in% show_pathway) %>%
        dplyr::select(ID, dplyr::all_of(stats_metric), geneID) %>%
        tidyr::separate_rows(geneID, sep = "\\/")
    }


    # if gsea has no entrezid, transID first
    # if (!all(sapply(new_gsea_df$geneID, function(x) grepl("^[0-9].*[0-9]$", x, perl = T)))) {
    #   id_map <- suppressMessages(transId(new_gsea_df$geneID, "entrezid", gsea_list$org))
    #   new_gsea_df <- merge(new_gsea_df, id_map,
    #     by.x = "geneID", by.y = "input_id",
    #     all.x = T, all.y = F
    #   ) %>% dplyr::select(-geneID) %>% dplyr::rename(geneID = entrezid)
    # }

    logfc <- gsea_list$genelist

    plot_df <- merge(new_gsea_df, logfc, by.x = "geneID", by.y = "ID") %>%
      dplyr::select(-geneID)

    term_order <- plot_df %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(logfc2 = sum(logfc)) %>%
      dplyr::arrange(desc(logfc2))
    plot_df$ID <- factor(plot_df$ID, levels = rev(term_order$ID))

    if(label_by != 'id'){
      gsea_df <- gsea_list$gsea_df
      colnames(gsea_df)[1] = "ID"

      plot_df = merge(plot_df,gsea_df,by = 'ID',all.y = F) %>%
        dplyr::select(4,2,3)
      colnames(plot_df)[1] = 'ID'
      colnames(plot_df)[2] = stats_metric

    }

    if(is.null(colour)){
      colour = c("#E31A1C","#1F78B4")
    }
    up_color = colour[1]; down_color = colour[2]

    p <- ggplot(plot_df, aes_string(x = "logfc", y = "ID", fill = stats_metric)) +
      ggridges::geom_density_ridges() +
      scale_fill_continuous(
        low = up_color, high = down_color, name = stats_metric,
        guide = guide_colorbar(reverse = TRUE)
      ) +
      ylab(NULL) +
      xlab("log2(Fold Change)") +
      guides(fill = guide_colourbar(title = stats_metric_label, reverse = T)) +
      plot_theme(remove_grid = T, ...)
  }

  #--- two-side barplot ---#
  ## bar with p.adjust>0.05 showed grey
  if (plot_type == "bar") {

    if(label_by == 'id'){
      gsea_df <- gsea_df %>%
        dplyr::select(ID, NES, "p.adjust") %>%
        dplyr::mutate(
          padj.group = cut(.$p.adjust, breaks = c(-Inf, 0.05, Inf), labels = c(1, 0)),
          nes.group = cut(.$NES, breaks = c(-Inf, 0, Inf), labels = c(0, 1)),
          comb.group = paste0(padj.group, nes.group)
        ) %>%
        dplyr::mutate(group = dplyr::case_when(
          comb.group == "10" ~ "A", # A is nes<0 & p<0.05
          comb.group == "11" ~ "B", # B is nes>0 & p<0.05
          TRUE ~ "C"
        )) %>% # C is grey
        dplyr::arrange(NES) %>%
        dplyr::mutate(index = 1:dplyr::n())

    }else{
      gsea_df <- gsea_df %>%
        dplyr::select(Description, NES, "p.adjust") %>%
        dplyr::rename(ID = Description) %>%
        dplyr::mutate(
          padj.group = cut(.$p.adjust, breaks = c(-Inf, 0.05, Inf), labels = c(1, 0)),
          nes.group = cut(.$NES, breaks = c(-Inf, 0, Inf), labels = c(0, 1)),
          comb.group = paste0(padj.group, nes.group)
        ) %>%
        dplyr::mutate(group = dplyr::case_when(
          comb.group == "10" ~ "A", # A is nes<0 & p<0.05
          comb.group == "11" ~ "B", # B is nes>0 & p<0.05
          TRUE ~ "C"
        )) %>% # C is grey
        dplyr::arrange(NES) %>%
        dplyr::mutate(index = 1:dplyr::n())
    }


    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 8

    if (is.null(colour)) {
      colour <- c("\\#1F78B4", "\\#E31A1C", "\\#A9A9A9")
      colour <- stringr::str_remove_all(colour, ".*#") %>% paste0("#", .)
    }

    p <- ggplot(gsea_df, aes(x = index, y = NES, fill = group)) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_manual(values = c("A" = colour[1], "B" = colour[2], "C" = colour[3])) +
      scale_x_discrete(expand = expansion(add = .5)) +
      scale_y_continuous(breaks = seq(
        floor(min(gsea_df$NES)), ceiling(max(gsea_df$NES)),
        ceiling((ceiling(max(gsea_df$NES)) - floor(min(gsea_df$NES))) / 6)
      )) +
      coord_flip()

    pos_new <- sum(gsea_df$NES>0); neg_nes <- sum(gsea_df$NES<0)
    if(neg_nes == 0 & pos_new != 0 ){
      p <- p +
        geom_text(
          data = subset(gsea_df, NES > 0),
          aes(x = index, y = 0, label = paste0("  ", ID), color = padj.group),
          size = lst$main_text_size / 3.6, hjust = "outward"
        )
    }else if(neg_nes != 0 & pos_new == 0 ){
      p <- p +
        geom_text(
          data = subset(gsea_df, NES < 0),
          aes(x = index, y = 0, label = paste0("  ", ID), color = padj.group),
          size = lst$main_text_size / 3.6, hjust = "inward"
        )
    }else{

      if(pos_new<=neg_nes){
        p <- p +
          geom_text(
            data = subset(gsea_df, NES > 0),
            aes(x = index, y = 0, label = paste0(ID, "  "), color = padj.group),
            size = lst$main_text_size / 3.6,
            hjust = "inward"
          ) +
          geom_text(
            data = subset(gsea_df, NES < 0),
            aes(x = index, y = 0, label = paste0("  ", ID), color = padj.group),
            size = lst$main_text_size / 3.6, hjust = "outward"
          )
      }else if (pos_new > neg_nes){
        p <- p +
          geom_text(
            data = subset(gsea_df, NES > 0),
            aes(x = index, y = 0, label = paste0(ID, "  "), color = padj.group),
            size = lst$main_text_size / 3.6,
            hjust = "inward"
          ) +
          geom_text(
            data = subset(gsea_df, NES < 0),
            aes(x = index, y = 0, label = paste0("  ", ID), color = padj.group),
            size = lst$main_text_size / 3.6, hjust = "outward"
          )
      }
    }

    p <- p +
      scale_colour_manual(values = c("black", colour[3])) +
      labs(x = "", y = "Normalized enrichment score") +
      plot_theme(remove_grid = T, remove_legend = T, ...)
  }

  suppressMessages(print(p))
}


calcScore <- function(geneset, genelist, item, exponent, fortify = TRUE, org) {
  genelist <- genelist[!duplicated(names(genelist))]
  geneset <- geneset %>%
    dplyr::filter(.[[1]] %in% item) %>%
    dplyr::pull(2) %>%
    as.character()
  geneset <- intersect(geneset, names(genelist))
  L <- length(genelist)
  Ls <- length(geneset)
  Nhit <- Nmiss <- numeric(L)
  hits <- names(genelist) %in% geneset
  Nhit[hits] <- abs(genelist[hits])^exponent
  NR <- sum(Nhit)
  Nhit <- cumsum(Nhit / NR)
  Nmiss[!hits] <- 1 / (L - Ls)
  Nmiss <- cumsum(Nmiss)
  score <- Nhit - Nmiss
  max.ES <- max(score)
  min.ES <- min(score)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  df <- data.frame(
    x = seq_along(score), runningScore = score,
    position = as.integer(hits)
  )


  if(org %in% c('hsapiens','mmusculus','rnorvegicus')){
    df$gene <- suppressMessages(transId(names(genelist), "symbol", org, keepNA = T, unique = T) %>%
                                  dplyr::pull(symbol))
  }else{
    df$gene <- names(genelist)
  }


  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$genelist <- genelist

  df$Description <- item
  return(df)
}


utils::globalVariables(c("NES", "qvalue", "group", "Description", "geom_line", "x", "y", "desc", "geneID", "geom_tile", "logfc2","Description.y"))
