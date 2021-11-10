# 不指定show_pathway名称时，默认使用前后5个；
# 当然可以自己指定可视化的名称，自己挑选出来，然后自动判断NES属于左边还是右边
# gsea_df =

plotGSEA <- function(gsea_list,
                     plot_type = c('volcano','pathway'),
                     show_pathway = 5,
                     show_genes = NULL,
                     ...){

  #--- args ---#
  stopifnot(is.list(gsea_list))
  plot_type <- match.arg(plot_type)

  if(plot_type == 'volcano'){
    gsea_df <- gsea_list$gsea_df
    gsea_df <- gsea_df[order(gsea_df$NES,decreasing = TRUE),]
    if(is.numeric(show_pathway)){
      gsea_df$group <- c(rep('high',show_pathway),rep('ignore',nrow(gsea_df)-2*show_pathway),rep('low',show_pathway))
    }else{
      nes <- gsea_df[gsea_df$Description %in% show_pathway,'NES']
      gsea_df$group <- 'ignore'
      gsea_df[gsea_df$Description %in% show_pathway,'group'] <- sapply(nes, function(x) ifelse(x>0,'high','low'))
    }

    p <- ggplot(gsea_df, aes(x = NES,y=-log10(qvalue),color=group)) +
      geom_point(alpha=0.6,size=3.5) +
      theme_set(theme_set(theme_bw(base_size = 20))) +
      xlab("NES") + ylab("-log10(pvalue)") +
      scale_color_manual(values = c("red","#00000066","darkgreen"))+
      theme_bw()+theme(legend.position="none")+
      geom_text_repel(data = gsea_df[gsea_df$group != 'ignore',],
                      aes(label = Description),
                      size =3,
                      color = "black",
                      show.legend = F)+
      plot_theme(remove_grid = T,remove_legend = T,border_thick = 1.5, ...)
  }

  if(plot_type == 'pathway'){
    gsea_df <- gsea_list$gsea_df
    geneset <- gsea_list$geneset
    genelist <- gsea_list$genelist
    exponent <- gsea_list$exponent
    org <- gsea_list$org

    if(is.numeric(show_pathway)) show_pathway = gsea_df$Description[show_pathway]

    df <- do.call(rbind, lapply(show_pathway, function(x) calcScore(geneset,genelist,x, exponent,fortify = TRUE, org)))

    colors <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

    p1 <- ggplot(df, aes_(x = ~x)) + xlab(NULL) +
      geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
      ylab("Enrichment\n Score") +
      plot_theme(remove_grid = T)+
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
      plot_theme(remove_legend = T,remove_grid = T,remove_text = T)+
      scale_y_continuous(expand=c(0,0))+
      theme(
        plot.margin = margin(t=-.1, b=0,unit="cm"),
            axis.line.x = element_blank(),
            axis.ticks = element_blank())

    df$y <- df$geneList
    select_genes <- data.frame(gene =  c('MCM2','MCM5','MCM6','PSMB9','PSMD3'))
    select_genes <- merge(select_genes, df, by = "gene")
    select_genes <- select_genes[select_genes$position == 1,]

    p3 <- ggplot(select_genes, aes(x, y, fill = Description, color = Description, label = gene)) +
      geom_segment(data=df, aes_(x=~x, xend=~x, y=~y, yend=0),
                   color = "grey") +
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values = colors, guide="none") +
      scale_color_manual(values = colors, guide="none") +
      geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
      ylab("Ranked list\n metric") +
      xlab("Rank in ordered dataset") +
      plot_theme(remove_grid = T)+
      # show gene names
      geom_text_repel(data = select_genes,
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

  return(p)


}



calcScore <- function(geneset,genelist,item, exponent, fortify = TRUE, org) {
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

  df$gene = suppressMessages(transId(names(genelist),'symbol',org))

  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList

  df$Description <- item
  return(df)
}









