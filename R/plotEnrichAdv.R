#' Advanced Plot for GO and KEGG enrichment analysis
#' Both up and down regulated pathways could be plotted in one figure as two-side barplot
#' @param up_enrich_df Enrichment analysis `data.frame` for up-regulated genes.
#' @param down_enrich_df Enrichment analysis `data.frame` for down-regulated genes.
#' @param plot_type Choose from "one" and "two". "One" represents both up and down pathways are plotted
#' together; "two" represents up and down are plotted seperately.
#' @param term_metric Pathway term metric from one of 'GeneRatio','Count','FoldEnrich' and
#' 'RichFactor'.
#' @param stats_metric Statistic metric from one of "pvalue", "p.adjust", "qvalue".
#' @param wrap_length Numeric, wrap text if longer than this length. Default is NULL.
#' @param color Plot colors.
#' @param ... other arguments from `plot_theme` function
#' @importFrom ggplot2 ggplot scale_y_discrete scale_x_reverse theme element_blank geom_bar aes
#' scale_fill_manual coord_flip ylim scale_x_discrete scale_y_continuous element_blank ylab guides guide_legend
#' @importFrom dplyr arrange mutate group_by top_n ungroup select case_when distinct rename pull
#' @importFrom rlang .data
#' @importFrom stringr str_wrap str_replace
#' @return A ggplot object
#' @export
plotEnrichAdv <- function(up_enrich_df,
                          down_enrich_df,
                          plot_type = c('one','two'),
                          term_metric = c("FoldEnrich", "GeneRatio", "Count", "RichFactor"),
                          stats_metric = c("p.adjust", "pvalue", "qvalue"),
                          wrap_length = NULL,
                          color,
                          ...){
  #--- args ---#
  stopifnot("The input enrichment analysis is not data frame!"=
              is.data.frame(up_enrich_df)|is.data.frame(down_enrich_df))
  plot_type <- match.arg(plot_type)

  #--- codes ---#
  tryCatch(
    {
      up_enrich_df$Description <- stringr::str_replace(up_enrich_df$Description, "^\\w{1}", toupper)
      down_enrich_df$Description <- stringr::str_replace(down_enrich_df$Description, "^\\w{1}", toupper)
    },
    error = function(e) {
      message(paste0("We need the 'Description' column which means pathway detailed description","\n",
              "Maybe you should rename the column name..."))
    }
  )

  x_label <- ifelse(stats_metric == "pvalue", "-log10(Pvalue)",
                               ifelse(stats_metric == "p.adjust", "-log10(P.adjust)", "-log10(FDR)")
  )


  #--- plot ---#
  if(plot_type == 'two'){
    if(missing(color)) color = c("#a32a31","#f7dcca","#3665a6","#d5e4ef")
    left = suppressMessages(plotEnrich(up_enrich_df,plot_type = "bar",
                      term_metric = term_metric,
                      stats_metric = stats_metric,
                      up_color = color[1], down_color = color[2],...)+
      scale_y_discrete(limits=rev)+
      scale_x_reverse()+
      theme(axis.title.y = element_blank(),
            legend.position = c(0.2, 0.8)))

    right = suppressMessages(plotEnrich(down_enrich_df,plot_type = "bar",
                       term_metric = term_metric,
                       stats_metric = stats_metric,
                       up_color = color[3],down_color = color[4],...)+
      scale_y_discrete(position = "right")+
      theme(axis.title.y = element_blank(),
            legend.position = c(0.8, 0.2)))

    p <- cowplot::plot_grid(left, right, ncol=2)

  }else{
    if(missing(color)) color = c("#3665a6", "#a32a31")

    up_go = dplyr::mutate(up_enrich_df,change = 'up')
    down_go = dplyr::mutate(down_enrich_df,change = 'down')
    df = rbind(up_go, down_go) %>%
      dplyr::mutate(new_x = ifelse(change == "up", -log10(eval(parse(text = stats_metric))), log10(eval(parse(text = stats_metric))))) %>%
      dplyr::arrange(change,new_x) %>%
      dplyr::mutate(Description = factor(Description, levels = unique(Description),
                                         ordered = TRUE))

    tmp = with(df, labeling::extended(range(new_x)[1], range(new_x)[2], m = 5))
    lm = tmp[c(1, length(tmp))]
    lm = c(floor(min(df$new_x)), ceiling(max(df$new_x)))

    p <- suppressMessages(ggplot(df, aes(x = Description, y = new_x)) +
      geom_bar(stat = "identity",
               aes(fill = change), width = 0.7) +
      scale_fill_manual(values = color,
                        name = "change",
                        labels = c("Down-regulated pathways","Up-regulated pathways")) +
      coord_flip() +  ylim(lm) +
      scale_x_discrete(labels = function(x)
        stringr::str_wrap(x, width = 30)) +
      scale_y_continuous(breaks = tmp, labels = abs(tmp)) +
      ylab(x_label)+
      guides ( fill = guide_legend(reverse = TRUE) ) +
      plot_theme(...)+
      theme(legend.title=element_blank()))
  }

  # wrap long text
  if (!is.null(wrap_length) & is.numeric(wrap_length)) {
    p <- p + scale_y_discrete(labels = text_wraper(wrap_length))
  }

  return(p)

}
