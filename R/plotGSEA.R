# 不指定show_item名称时，默认使用前后5个；
# 当然可以自己指定可视化的名称，自己挑选出来，然后自动判断NES属于左边还是右边

plotGSEA <- function(gsea_df,
                     plot_type = c('volcano','pathway'),
                     show_item = 5){

  #--- args ---#
  stopifnot(is.data.frame(gsea_df))
  plot_type <- match.arg(plot_type)

  if(plot_type == 'volcano'){
    e_left <- gsea_df$Description[tail(order(gsea_df$NES),show_item)]
    e_right <- gsea_df$Description[head(order(gsea_df$NES),show_item)]
    data$enrich <- ifelse(data$Description %in% e_left ,"left",
                          ifelse(data$Description %in% e_right ,"right","no"))
    ggplot(gsea_df, aes(x = NES,y=-log10(qvalue),color=enrich)) +
      geom_point(alpha=0.6,size=3.5) +
      theme_set(theme_set(theme_bw(base_size = 20))) +
      xlab("NES") + ylab("-log10(pvalue)") +
      scale_color_manual(values = c("red","#00000033","darkgreen"))+
      theme_bw()+theme(legend.position="none")+
      geom_text_repel(data = data[data$Description %in% c(e_left,e_right),],
                      aes(label = Description),
                      size =3,
                      color = "black",
                      show.legend = F)+
      plot_theme(remove_grid = T,remove_legend = T,border_thick = 1.5)
  }




}














