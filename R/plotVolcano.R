#' Volcano plot for differential expression analysis
#'
#' @param deg_df DEG dataframe with gene id, logFC and stat(e.g. pvalue/qvalue).
#' @param stat_metric Statistic metric from "pvalue" or "p.adjust".
#' @param stat_cutoff Statistic cutoff, default is 0.05.
#' @param logFC_cutoff Log2 fold change cutoff, default is 1 which is actually 2 fold change.
#' @param up_color Color of up-regulated genes, default is "dark red".
#' @param down_color Color of down-regulated genes, default is "dark blue".
#' @param show_gene Select genes to show, default is no genes to show.
#' @param dot_size Volcano dot size, default is 1.75.
#' @param ... other arguments from `plot_theme` function
#'
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 ggplot geom_point scale_color_manual geom_hline geom_vline labs  xlim guides
#' guide_legend aes unit
#' @return A ggplot object
#' @export
#' @examples
#' \donttest{
#' library(ggplot2)
#' data(deg, package = "genekitr")
#' plotVolcano(deg, "p.adjust", remove_legend = TRUE, dot_size = 3)
#'
#' # show some genes
#' plotVolcano(deg, "p.adjust",
#'   remove_legend = TRUE,
#'   show_gene = c("CD36", "DUSP6", "IER3","CDH7")
#' )
#' }
#'
plotVolcano <- function(deg_df,
                        stat_metric = c("p.adjust", "pvalue"),
                        stat_cutoff = 0.05,
                        logFC_cutoff = 1,
                        up_color = "#E31A1C",
                        down_color = "#1F78B4",
                        show_gene = NULL,
                        dot_size = 1.75,
                        ...) {

  #--- args ---#
  lst <- list(...) # store outside arguments in list
  stat_metric <- match.arg(stat_metric)

  #--- codes ---#
  ## set labels
  stat_metric_label <- ifelse(stat_metric == "pvalue", "P-value", "adjusted P-value")

  ## check column name
  # check adj.P.Val
  check_adjp <- which(grepl("^adj|qval|fdr", tolower(colnames(deg_df))))
  if (length(check_adjp) > 1) check_adjp <- check_adjp[1]

  # check pvalue
  check_p <- which(grepl("^p.*val.*", tolower(colnames(deg_df))))
  if (length(check_adjp) + length(check_p) == 0) {
    stop('Seems input data missing statistical column e.g. "pvalue"...')
  }

  # check logfc
  check_logfc <- which(grepl("^log.*|fold", tolower(colnames(deg_df))))
  if (length(check_logfc) == 0) {
    stop("Seems input data missing fold change column ...")
  }

  # check gene
  check_gene <- which(grepl("gene|entrezid|symbol|ensembl", tolower(colnames(deg_df))))
  if (length(check_gene) > 1 & !is.null(show_gene)) {
    for(i in check_gene){
      if(any(show_gene%in%deg_df[,i])) check_gene <- seq(ncol(deg_df))[i]
    }
  }else{
    check_gene = check_gene[1]
  }

  ## subset
  if (stat_metric == "pvalue") {
    plot_df <- deg_df[, c(check_gene, check_logfc, check_p)]
  } else {
    plot_df <- deg_df[, c(check_gene, check_logfc, check_adjp)]
  }

  plot_df <- plot_df %>% stats::setNames(c("gene", "logFC", "stat")) %>% # set pvalue/adjP as stat
    dplyr::mutate(change = as.factor(ifelse(stat < stat_cutoff & abs(logFC) > logFC_cutoff,
      ifelse(logFC > logFC_cutoff, "UP", "DOWN"), "NOT"
    )))

  #--- plot ---#
  xlim_range <- ceiling(max(max(plot_df$logFC), abs(min(plot_df$logFC))))
  p <- ggplot(data = plot_df, aes(x = logFC, y = -log10(stat), color = change)) +
    geom_point(alpha = 0.4, size = dot_size) +
    scale_color_manual(values = c(down_color, "black", up_color)) +
    geom_hline(yintercept = -log10(stat_cutoff), lty = 4, lwd = 0.6, alpha = 1) +
    geom_vline(xintercept = c(logFC_cutoff, -logFC_cutoff), lty = 4, lwd = 0.6, alpha = 1) +
    labs(x = "Log2 (fold change)", y = paste0("-Log10 (", stat_metric_label, ")")) +
    xlim(-xlim_range, xlim_range) +
    plot_theme(...)

  # show genes
  if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 5
  if (!is.null(show_gene)) {
    show_gene_df <- plot_df %>%
      dplyr::filter(gene %in% show_gene) %>%
      dplyr::mutate(label = gene)
    p <- p + geom_point(
      data = show_gene_df, alpha = 1, size = dot_size, shape = 1,
      stroke = 1,
      color = "black"
    ) +
      ggrepel::geom_label_repel(
        data = show_gene_df, aes(label = label),
        show.legend = FALSE,
        size = lst$main_text_size / 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) +
      guides(color = guide_legend(title = NULL))
  }

  suppressWarnings(print(p))
}
