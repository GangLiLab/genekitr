# 不指定show_item名称时，默认使用前后5个；
# 当然可以自己指定可视化的名称，自己挑选出来，然后自动判断NES属于左边还是右边

plotGSEA <- function(gsea_df,
                     plot_type = c('volcano','pathway'),
                     show_item = 5,
                     ...){

  #--- args ---#
  stopifnot(is.data.frame(gsea_df))
  plot_type <- match.arg(plot_type)

  if(plot_type == 'volcano'){
    gsea_df <- gsea_df[order(gsea_df$NES,decreasing = TRUE),]
    if(is.numeric(show_item)){
      gsea_df$group <- c(rep('high',show_item),rep('ignore',nrow(gsea_df)-2*show_item),rep('low',show_item))
    }else{
      nes <- gsea_df[gsea_df$Description %in% show_item,'NES']
      gsea_df$group <- 'ignore'
      gsea_df[gsea_df$Description %in% show_item,'group'] <- sapply(nes, function(x) ifelse(x>0,'high','low'))
    }

    ggplot(gsea_df, aes(x = NES,y=-log10(qvalue),color=group)) +
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
      plot_theme(remove_grid = T,remove_legend = T,border_thick = 1.5)
  }




}














