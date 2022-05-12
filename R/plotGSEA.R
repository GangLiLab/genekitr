#' GSEA plot
#'
#' @param gsea_list GSEA result from `genGSEA` function
#' @param plot_type GSEA plot type, one of 'volcano', 'classic', 'fgsea' or 'ridge'.
#' @param show_pathway Select plotting pathways by specifying number or
#'   pathway name.
#' @param show_gene Character to specify gene names included in plot when `plot_type` is "classic".
#' @param colors Character to specify colors when `plot_type` is "classic".
#' @param legend_type Legend type from "p.adjust", "pvalue" or "qvalue" when `plot_type` is "ridge".
#' @param ... other arguments transfer to `plot_theme` function
#'
#' @importFrom ggplot2 ggplot aes aes_ geom_point xlab ylab scale_color_manual geom_hline element_blank
#'   element_rect margin geom_linerange scale_y_continuous geom_segment geom_bar scale_fill_manual
#'   geom_hline element_line geom_line
#' @importFrom dplyr filter pull %>%
#'
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ## get GSEA result
#' data(geneList, package = "genekitr")
#' gse <- genGSEA(genelist = geneList, org = "human",
#'                category = "H",use_symbol = TRUE, pvalueCutoff = 1)
#' ## volcano plot
#' # get top3 of up and down pathways
#' plotGSEA(gse, plot_type = 'volcano', show_pathway = 3)
#' # choose pathway by character
#' plotGSEA(gse, plot_type = 'volcano',
#'   show_pathway = c("HALLMARK_ADIPOGENESIS","HALLMARK_XENOBIOTIC_METABOLISM",
#'   "HALLMARK_INTERFERON_ALPHA_RESPONSE"))
#'
#' ## classic pathway plot
#' plotGSEA(gse, plot_type = 'classic', show_pathway = 1:2)
#'
#' ## fgsea for multiple pathway
#' plotGSEA(gse, plot_type = 'fgsea', show_pathway = 3)
#'
#' ## ridgeplot
#' plotGSEA(gse, plot_type = 'ridge',
#'   show_pathway = 10, legend_type = 'p.adjust')
#'
#' }
#'

plotGSEA <- function(gsea_list,
                     plot_type = c('volcano','classic','fgsea','ridge'),
                     show_pathway = 3,
                     show_gene = NULL,
                     colors = NULL,
                     legend_type = c("p.adjust", "pvalue", "qvalue"),
                     ...){

  #--- args ---#
  stopifnot(is.list(gsea_list))
  plot_type <- match.arg(plot_type)
  legend_type <- match.arg(legend_type)

  if(plot_type == 'volcano'){
    gsea_df <- gsea_list$gsea_df
    gsea_df <- gsea_df[order(gsea_df$NES,decreasing = TRUE),]
    if(is.numeric(show_pathway)){
      if(length(show_pathway) > 1){
        gsea_df = gsea_df[show_pathway,]
      }
      gsea_df$group <- c(rep('high',show_pathway),rep('ignore',nrow(gsea_df)-2*show_pathway),rep('low',show_pathway))
    }else{
      nes <- gsea_df[gsea_df$Description %in% show_pathway,'NES']
      gsea_df$group <- 'ignore'
      gsea_df[gsea_df$Description %in% show_pathway,'group'] <- sapply(nes, function(x) ifelse(x>0,'high','low'))
    }

    p <- ggplot(gsea_df, aes(x = NES,y=-log10(qvalue),color=group)) +
      geom_point(alpha=0.6,size=3.5) +
      xlab("NES") + ylab("-log10(pvalue)") +
      scale_color_manual(values = c("red","grey66","darkgreen"))+
      ggrepel::geom_text_repel(data = gsea_df[gsea_df$group != 'ignore',],
                      aes(label = Description),
                      size =3,
                      color = "black",
                      show.legend = F)+
      plot_theme(...)
  }

  if(plot_type == 'classic'){
    gsea_df <- gsea_list$gsea_df
    geneset <- gsea_list$geneset
    genelist <- gsea_list$genelist
    exponent <- gsea_list$exponent
    org <- gsea_list$org

    if(is.numeric(show_pathway)){
      show_pathway = gsea_df$Description[show_pathway]

    }else if(any(!show_pathway %in% gsea_df$Description)){
      stop(paste0(show_pathway[!show_pathway%in%gsea_df$Description],' not in GSEA result!'))
    }

    df <- do.call(rbind, lapply(show_pathway, function(x) calcScore(geneset,genelist,x, exponent,fortify = TRUE, org)))

    if(is.null(colors)){
      colors <- c("\\#5DA5DAFF","\\#FAA43AFF","\\#60BD68FF","\\#F15854FF","\\#B276B2FF",
                  "\\#8D4B08FF","\\#DECF3FFF","\\#F17CB0FF","\\#66E3D9FF","\\#00FF7FFF",
                  "\\#E31A1CFF","\\#FFFF99FF")
      colors <- stringr::str_remove_all(colors,'.*#') %>% paste0('#',.)
    }

    p1 <- ggplot(df, aes_(x = ~x)) + xlab(NULL) +
      geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
      ylab("Enrichment\n Score") +
      plot_theme(remove_grid = T,...)+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line.x=element_blank(),
            legend.position = "top", legend.title = element_blank(),
            legend.background = element_rect(fill = "transparent"),
            plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))

    # size for three panels
    panel_size <- c(1.5, .5, 1.5)

    i <- 0
    for (term in unique(df$Description)) {
      idx <- which(df$ymin != 0 & df$Description == term)
      df[idx, "ymin"] <- i
      df[idx, "ymax"] <- i + 1
      i <- i + 1
    }

    p2 <- ggplot(df, aes_(x = ~x)) +
      geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
      xlab(NULL) + ylab(NULL) +
      scale_color_manual(values = colors)+
      plot_theme(remove_legend = T,remove_grid = T,remove_text = T,...)+
      scale_y_continuous(expand=c(0,0))+
      theme(
        plot.margin = margin(t=-.1, b=0,unit="cm"),
            axis.line.x = element_blank(),
            axis.ticks = element_blank())

    df$y <- df$genelist

    if(!is.null(show_gene)){
      select_genes <- data.frame(gene =  show_gene)
      select_genes <- merge(select_genes, df, by = "gene")
      select_genes <- select_genes[select_genes$position == 1,]
    }else{
      select_genes <- data.frame(gene =  'no_select_genes')
      select_genes <- merge(select_genes, df, by = "gene")
      select_genes <- select_genes[select_genes$position == 1,]
    }

    p3 <- ggplot(select_genes, aes(x, y, fill = Description, color = Description, label = gene)) +
      geom_segment(data=df, aes_(x=~x, xend=~x, y=~y, yend=0),
                   color = "grey") +
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values = colors, guide="none") +
      scale_color_manual(values = colors, guide="none") +
      geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
      ylab("Ranked list\n metric") +
      xlab("Rank in ordered dataset") +
      plot_theme(remove_grid = T,...)+
      # show gene names
      ggrepel::geom_text_repel(data = select_genes,
                      show.legend = FALSE,
                      direction = "x",
                      ylim = c(2, NA),
                      angle = 90,
                      size = 3, box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.3, "lines")) +
      theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

    pl <- list(p1, p2, p3)
    n <- length(pl)
    pl[[n]] <- pl[[n]] +
      theme(axis.line.x = element_line(),
            axis.ticks.x = element_line())

    p = cowplot::plot_grid(plotlist = pl, ncol = 1, align="v", rel_heights = panel_size)

  }

  if(plot_type == 'fgsea'){
    geneset_list <- gsea_list$geneset %>% split(.$gs_name) %>% lapply( "[[", 2)
    genelist <- gsea_list$genelist

    fres <- suppressWarnings(fgsea::fgsea(pathways = geneset_list,
                                          stats = genelist))

    if(is.numeric(show_pathway)){
      up_path = fres %>% dplyr::filter(ES > 0 ) %>% dplyr::arrange(.,pval) %>%
        dplyr::slice_head(n =show_pathway) %>% dplyr::pull(pathway)
      down_path = fres %>% dplyr::filter(ES < 0 ) %>% dplyr::arrange(.,pval) %>%
        dplyr::slice_head(n =show_pathway) %>% dplyr::pull(pathway)
      all = c(up_path,down_path)
    }else if(any(!show_pathway %in% fres$pathway)){
      stop(paste0(show_pathway[!show_pathway%in%fres$pathway],' not in fgsea result!'))
    }else{
      all = show_pathway
    }

    p = fgsea::plotGseaTable(geneset_list[all], genelist, fres,render = F) %>%
      ggplotify::as.ggplot()

  }

  if(plot_type == 'ridge'){
    gsea_df <- gsea_list$gsea_df
    if(is.numeric(show_pathway)){
      if(length(show_pathway) > 1){
        show_pathway = gsea_df$Description[show_pathway,]
      }else{
        show_pathway = gsea_df$Description[1:show_pathway]
      }

    }
    new_gsea_df <- gsea_df %>%
      dplyr::filter(Description%in%show_pathway) %>%
      dplyr::select(Description,all_of(legend_type),geneID) %>%
      tidyr::separate_rows(geneID,sep='\\/')

    # if gsea has no entrezid, transID first
    if(!all(sapply(new_gsea_df$geneID, function(x) grepl("^[0-9].*[0-9]$", x, perl = T)))){
      id_map = suppressMessages(transId(new_gsea_df$geneID,'entrezid',gsea_list$org))
      new_gsea_df = merge(new_gsea_df,id_map,by.x = 'geneID',by.y = 'input_id',
                 all.x = T, all.y=F)
    }else{
      new_gsea_df = new_gsea_df %>% dplyr::rename(entrezid = geneID)
    }

    logfc = data.frame(entrezid = names(gsea_list$genelist),
                       logfc = gsea_list$genelist)
    plot_df = merge(new_gsea_df,logfc, by = 'entrezid') %>%
      dplyr::select(-entrezid)

    term_order = plot_df %>%
      dplyr::group_by(Description) %>%
      dplyr::summarise(logfc2 = sum(logfc)) %>%
      dplyr::arrange(desc(logfc2))
    plot_df$Description = factor(plot_df$Description,levels = rev(term_order$Description))

    p = ggplot(plot_df, aes_string(x="logfc", y="Description", fill=legend_type)) +
      ggridges::geom_density_ridges()  +
      scale_fill_continuous(low="red", high="blue", name = legend_type,
                            guide=guide_colorbar(reverse=TRUE))+
      ylab(NULL) + xlab('logFC')+
      plot_theme(remove_grid = T,...)

  }

  suppressMessages(print(p))

}



calcScore <- function(geneset,genelist,item, exponent, fortify = TRUE, org) {
  genelist <- genelist[!duplicated(names(genelist))]
  geneset <- geneset %>% dplyr::filter(gs_name %in% item) %>% dplyr::pull('entrez_gene') %>% as.character()
  geneset <- intersect(geneset, names(genelist))
  L <- length(genelist)
  Ls <- length(geneset)
  Nhit <- Nmiss <- numeric(L)
  hits <- names(genelist) %in% geneset
  Nhit[hits] <- abs(genelist[hits])^exponent
  NR <- sum(Nhit)
  Nhit <- cumsum(Nhit/NR)
  Nmiss[!hits] <- 1/(L - Ls)
  Nmiss <- cumsum(Nmiss)
  score <- Nhit - Nmiss
  max.ES <- max(score)
  min.ES <- min(score)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  }else {
    ES <- min.ES
  }
  df <- data.frame(x = seq_along(score), runningScore = score,
                   position = as.integer(hits))

  df$gene = suppressMessages(transId(names(genelist),'symbol',org,keepNA = T) %>%
                               dplyr::pull(symbol))

  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$genelist <- genelist

  df$Description <- item
  return(df)
}


utils::globalVariables(c("NES","qvalue","group","Description","geom_line","x","y","desc", "geneID" ,"geom_tile" ,"logfc2"))






