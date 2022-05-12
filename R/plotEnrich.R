#' Plot for GO and KEGG enrichment analysis
#'
#' @param enrich_df Enrichment analysis `data.frame` result.
#' @param logfc_df A two-column data frame, first column is enrichment analysis input id
#' and the second is matching logFC value. Default is NULL. Used in "heat" and "chord" plot.
#' @param plot_type Choose from "bar", "dot", "lollipop","geneHeat", "geneChord",
#' "map","goHeat","tangram","wordcloud".
#' @param xlab_type X-axis label type, one of 'GeneRatio','Count','FoldEnrich'.
#' @param legend_type Stats legend type, one of "pvalue", "p.adjust", "qvalue".
#' @param sim_method Method of calculating the similarity between nodes, one of one of "Resnik",
#' "Lin", "Rel", "Jiang" , "Wang" and "JC" (Jaccard similarity coefficient) methods.
#' Used in "map","goHeat","tangram","wordcloud".
#' @param top_color Legend top color as low pvalue or high logFC, default is "red".
#' @param bottom_color Legend bottom color as high pvalue or low logFC, default is "blue".
#' @param show_gene Select genes to show. Default is "all". Used in "heat" and "chord" plot.
#' @param xlim_left X-axis left limit, default is 0.
#' @param xlim_right X-axis right limit, default is NA.
#' @param wrap_length Numeric, wrap text if longer than this length, default is NULL.
#' @param scale_ratio Numeric, scale of nodes and line thickness. Default is 1. Used in "map" plot.
#' @param layout Grapgh layout in "map" plot, e,g, "circle", "dh", "drl", "fr","graphopt", "grid",
#' "lgl", "kk", "mds", "nicely" (default),"randomly", "star".
#' @param org  Organism name from `biocOrg_name`.
#' @param ont  One of "BP", "MF", and "CC".
#' @param ... other arguments transfer to `plot_theme` function
#'
#' @importFrom dplyr pull %>% arrange mutate slice_head
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_continuous scale_fill_continuous
#'   guide_colorbar scale_y_discrete element_blank xlab labs xlim scale_x_discrete theme
#'   facet_grid
#' @importFrom stringr str_to_sentence
#' @importFrom rlang .data
#' @importFrom stats setNames
#'
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' data(geneList, package = "genekitr")
#' id <- names(geneList)[1:100]
#' logfc <- data.frame(gene = id, logfc = geneList[id])
#' ego <- genGO(id,
#'   org = "human", ont = "bp", pvalueCutoff = 0.05,
#'   qvalueCutoff = 0.05)
#' plotEnrich(ego,plot_type = "dot")
#'
#' plotEnrich(ego,plot_type = "bar")
#'
#' plotEnrich(ego,plot_type = "lollipop",
#' bottom_color = "#325CAC", top_color = "#E69056", wrap_length = 25)
#'
#' plotEnrich(ego,plot_type = "geneHeat")
#' plotEnrich(ego,plot_type = "geneHeat",
#'   show_gene = c('AQP7','ASIC2','GRIN2A','SLITRK6'))
#' plotEnrich(ego,logfc_df = logfc, plot_type = "geneHeat",
#'   wrap_length = 25)
#'
#' plotEnrich(ego,logfc_df = logfc, plot_type = "geneChord")
#'
#' plotEnrich(ego,plot_type = "map", layout = "circle",
#'   scale_ratio = 1, legend_text_size = 8)
#'
#' plotEnrich(ego,plot_type = "goHeat",sim_method = 'Rel')
#'
#' plotEnrich(ego,plot_type = "tangram",sim_method = 'Rel')
#'
#' plotEnrich(ego,plot_type = "wordcloud",sim_method = 'Rel')
#' }
#'
plotEnrich <- function(enrich_df,
                       logfc_df = NULL,
                       plot_type = c('bar','dot','lollipop','geneHeat','geneChord','map',
                                     'goHeat','tangram','wordcloud'),
                       xlab_type = c("FoldEnrich", "GeneRatio", "Count", "RichFactor"),
                       legend_type = c("p.adjust", "pvalue", "qvalue"),
                       sim_method =  c("JC","Resnik", "Lin", "Rel", "Jiang" , "Wang"),
                       top_color = "red",
                       bottom_color = "blue",
                       show_gene = "all",
                       xlim_left = 0,
                       xlim_right = NA,
                       wrap_length = NULL,
                       scale_ratio = 1,
                       layout,
                       org = NULL,
                       ont = NULL,
                       ...) {
  #--- args ---#
  plot_type <- match.arg(plot_type)
  xlab_type <- match.arg(xlab_type)
  legend_type <- match.arg(legend_type)
  sim_method <- match.arg(sim_method)
  if(missing(layout)) layout = 'nicely'

  compare_group <- any(grepl('cluster',colnames(enrich_df),ignore.case = T))
  all_go <- any(grepl('ONTOLOGY',colnames(enrich_df),ignore.case = T))

  if(compare_group) plot_type = 'dot'
  if(all_go) plot_type = 'bar'

  types <- c("GeneRatio", "Count", "FoldEnrich", "RichFactor")
  legends <- c("p.adjust", "pvalue", "qvalue")

  if (!xlab_type %in% colnames(enrich_df)) {
    stop(
      xlab_type, " not included in this dataframe, try: ",
      paste(intersect(colnames(enrich_df), types), collapse = " | ")
    )
  }
  if (!legend_type %in% colnames(enrich_df)) {
    stop(
      legend_type, " not included in this dataframe, try: ",
      paste(intersect(colnames(enrich_df), legends), collapse = " | ")
    )
  }

  #--- codes ---#
  ## uppercase first letter of description
  enrich_df$Description <- stringr::str_to_sentence(enrich_df$Description)

  ## set labels
  xlab_title <- ifelse(xlab_type == "FoldEnrich", "Fold Enrichment",
    ifelse(xlab_type == "GeneRatio", "Gene Ratio",
           ifelse(xlab_type == "RichFactor", "Rich Factor", "Count"))
  )
  legend_title <- ifelse(legend_type == "pvalue", "Pvalue",
    ifelse(legend_type == "p.adjust", "P.adjust", "FDR")
  )

  ## set pathway as factor
  if(!compare_group & !all_go){
    enrich_df <- enrich_df %>%
      dplyr::arrange(eval(parse(text = xlab_type))) %>%
      dplyr::mutate(Description = factor(.$Description, levels = .$Description, ordered = T))
  }else if (!compare_group & all_go){
    enrich_df <- enrich_df %>%
      dplyr::mutate(Description = factor(.$Description, levels = .$Description, ordered = T)) %>%
      dplyr::group_by(ONTOLOGY) %>%
      dplyr::arrange(eval(parse(text = xlab_type)),.by_group = T)
  }

  ## set logfc colnames
  if(!is.null(logfc_df)){
    logfc_df <- logfc_df %>% stats::setNames(c('geneID','logfc'))
  }


  #--- dot plot ---#
  if(plot_type == 'dot'){
    if(!compare_group){
      p <- ggplot(enrich_df, aes_string(x = xlab_type, y = "Description")) +
        geom_point(aes_string(
          color = legend_type,
          size = "Count"
        )) +
        scale_color_continuous(
          low = top_color, high = bottom_color, name = legend_title,
          guide = guide_colorbar(reverse = TRUE),
          labels = function(x) format(x, scientific = T)
        ) +
        xlab(xlab_title) +
        labs(color = legend_type)+
        xlim(xlim_left,xlim_right)+
        plot_theme(...)
    }else{
      xtick_lab <- paste0(enrich_df$Cluster,'\n(',enrich_df$Count,')')
      p <- ggplot(enrich_df, aes_string(x = 'Cluster', y = "Description")) +
        geom_point(aes_string(
          color = legend_type,
          size = xlab_type
        )) +
        scale_color_continuous(
          low = top_color, high = bottom_color, name = legend_title,
          guide = guide_colorbar(reverse = TRUE),
          labels = function(x) format(x, scientific = T)
        ) +
        scale_x_discrete(labels= xtick_lab)+
        xlab('Group') +
        labs(color = legend_type)+
        plot_theme(...)+
        theme(axis.text.x = element_text(
          angle = 45,
          vjust = 0.5, hjust = 0.5
        ))
    }

  }

  #--- bar plot ---#
  if(plot_type == 'bar'){
    p <- ggplot(data=enrich_df, aes_string(x = xlab_type, y = 'Description', fill = legend_type)) +
      geom_bar(stat="identity")+
      scale_fill_continuous(
        low = top_color, high = bottom_color, name = legend_title,
        guide = guide_colorbar(reverse = TRUE),
        labels = function(x) format(x, scientific = T))+
      xlab(xlab_title)+
      labs(color = legend_type)+
      xlim(xlim_left,xlim_right)+
      plot_theme(...)

    if(all_go) p <- p + facet_grid(ONTOLOGY~., scales = "free")+
        plot_theme(...)

  }

  #--- lollipop plot ---#
  if(plot_type == 'lollipop'){
    p <- ggplot(data=enrich_df,
                aes(eval(parse(text = xlab_type)), forcats::fct_reorder(Description,eval(parse(text = xlab_type)))))+
      geom_segment(aes_string(xend = 0, yend = "Description",
                              colour=legend_type))+
      geom_point(aes_string(color = legend_type, size = "Count"))+
      scale_color_continuous(
        low = top_color, high = bottom_color, name = legend_title,
        guide = guide_colorbar(reverse = TRUE),
        labels = function(x) format(x, scientific = T))+
      scale_size_continuous(range = c(min(enrich_df$Count)-1, max(enrich_df$Count)/2))+
      xlab(xlab_title) + ylab(NULL)+
      labs(color = legend_type)+
      plot_theme(...)

  }

  #--- geneHeat plot ---#
  if(plot_type == 'geneHeat'){
    id = enrich_df %>% dplyr::pull(geneID) %>%
      stringr::str_split('\\/') %>% unlist()
    id_symbol = enrich_df %>% dplyr::pull(geneID_symbol) %>%
      stringr::str_split('\\/') %>% unlist()
    id_df = data.frame(geneID = id, geneID_symbol = id_symbol) %>%
      dplyr::distinct()

    if(all(show_gene == 'all')){
      plot_df = enrich_df %>% dplyr::select(Description,geneID_symbol) %>%
        tidyr::separate_rows(geneID_symbol,sep = '\\/') %>%
        dplyr::rename(geneID = geneID_symbol)
    }else{
      # if show_gene is not symbol, first extract matching symbol
      if(all(show_gene%in%id)){
        show_gene = id_df %>% dplyr::filter(geneID %in% show_gene) %>%
          dplyr::pull(geneID_symbol)
      }

      plot_df = enrich_df %>% dplyr::select(Description,geneID_symbol) %>%
        tidyr::separate_rows(geneID_symbol,sep = '\\/') %>%
        dplyr::filter(geneID_symbol %in% show_gene) %>%
        dplyr::rename(geneID = geneID_symbol)
    }


    if(is.null(logfc_df)){
      p <- ggplot(plot_df, aes_(~geneID, ~Description)) +
        geom_tile(color = 'white')+
        xlab(NULL) + ylab(NULL) +
        plot_theme(theme_type = 'bw',
                   border_thick = 0,...)+
        theme(panel.grid.major = element_blank(),
              axis.text.x=element_text(angle = 50, hjust = 1))
    }else{
      # add logfc
      if(all(id_df$geneID %in% logfc_df$geneID)){
        m1 = merge(id_df,logfc_df,by = 'geneID' )
        plot_df = merge(plot_df,m1,by.x = 'geneID',by.y = 'geneID_symbol') %>%
          dplyr::select(-geneID.y)
      }else{
        plot_df = merge(plot_df,logfc_df,by.x = 'geneID')
      }

      p <- ggplot(plot_df, aes_(~geneID, ~Description)) +
        geom_tile(aes_(fill = ~logfc), color = "white") +
        xlab(NULL) + ylab(NULL) +
        plot_theme(theme_type = 'bw',
                   border_thick = 0,...)+
        scale_fill_continuous(low=bottom_color, high=top_color,
                              name = "logFC")+
        theme(panel.grid.major = element_blank(),
              axis.text.x=element_text(angle = 50, hjust = 1))
    }

  }

  #--- geneChord plot ---#
  if(plot_type == 'geneChord'){
    id = enrich_df %>% dplyr::pull(geneID) %>%
      stringr::str_split('\\/') %>% unlist()
    id_symbol = enrich_df %>% dplyr::pull(geneID_symbol) %>%
      stringr::str_split('\\/') %>% unlist()
    id_df = data.frame(geneID = id, geneID_symbol = id_symbol) %>%
      dplyr::distinct()

    if(all(show_gene == 'all')){
      plot_df = enrich_df %>% dplyr::select(Description,geneID_symbol) %>%
        tidyr::separate_rows(geneID_symbol,sep = '\\/') %>%
        dplyr::rename(geneID = geneID_symbol)
    }else{
      # if show_gene is not symbol, first extract matching symbol
      if(all(show_gene%in%id)){
        show_gene = id_df %>% dplyr::filter(geneID %in% show_gene) %>%
          dplyr::pull(geneID_symbol)
      }

      plot_df = enrich_df %>% dplyr::select(Description,geneID_symbol) %>%
        tidyr::separate_rows(geneID_symbol,sep = '\\/') %>%
        dplyr::filter(geneID_symbol %in% show_gene) %>%
        dplyr::rename(geneID = geneID_symbol)
    }

    # define color
    my_cols= c("#B2DF8A","#FB9A99","#E31A1C","#B15928","#6A3D9A","#CAB2D6",
               "#A6CEE3","#1F78B4","#FDBF6F","#999999","#FF7F00")
    backup_cols = c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9",
                    "#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C",
                    "#E0367A","#D8D155","#64495D","#7CC767")

    term = unique(plot_df$Description)
    if(length(my_cols) < length(term)){
      my_cols = c(my_cols,backup_cols)
      cols = sample(my_cols,length(term),replace = F)
    }
    cols = my_cols[1:length(term)]

    # prepare chord data
    id = unique(plot_df$geneID)
    dat = sapply(id, function(x){
      check = plot_df %>% dplyr::filter(geneID %in% x) %>%
        dplyr::pull(Description) %>%
        unique()
      ifelse(term%in%check,1,0)
    }) %>% t() %>%
      as.data.frame() %>%
      stats::setNames(term)

    if(is.null(logfc_df)){
      p = suppressWarnings(
        GOplot::GOChord(dat,
                        space = 0.02,
                        gene.order = 'none',
                        gene.size = 3 ,
                        process.label = 8,
                        border.size = 0.1,
                        ribbon.col = cols)
      )

    }else{
      # add logfc
      if(all(rownames(dat) %in% logfc_df$geneID)){
        m1 = logfc_df %>% dplyr::filter(geneID%in%rownames(dat)) %>%
          dplyr::pull(logfc)
      }else{
        m1 = merge(id_df,logfc_df,by.x = 'geneID') %>%
          dplyr::filter(geneID_symbol %in% rownames(dat)) %>%
          dplyr::pull(logfc)
      }
      dat = dat %>% dplyr::mutate(logFC = m1)
      p = suppressWarnings(
        GOplot::GOChord(dat,
                        space = 0.02,
                        gene.order = 'logFC',
                        gene.size = 3 ,
                        process.label = 8,
                        border.size = 0.1,
                        ribbon.col = cols,
                        lfc.col=c(top_color,'grey50',bottom_color))
      )
    }
  }

  #--- map plot ---#
  if(plot_type == 'map'){
    id <- enrich_df[,1]
    enrichGenes <- strsplit(enrich_df$geneID,'\\/') %>% setNames(id)

    if(sim_method == "JC"){
      m = get_JC_data(enrich_df)
    }else{
      m = get_sim_data(enrich_df,org=NULL,ont=NULL,sim_method)[['m']]
      rownames(m) <- colnames(m) <- enrich_df %>%
        dplyr::filter(.[[1]]%in%colnames(m)) %>%
        dplyr::pull(Description)
    }

    mm <- reshape2::melt(m)
    mm <- mm[mm[,1] != mm[,2],]
    mm <- mm[!is.na(mm[,3]),]
    # construct igraph
    g <- graph.data.frame(mm[, -3], directed=FALSE)
    E(g)$width <- sqrt(mm[, 3] * 5) * scale_ratio
    E(g)$weight <- mm[, 3]
    g <- delete.edges(g, E(g)[mm[, 3] < 0.2])
    id_order <- unlist(sapply(V(g)$name, function(x) which(x == enrich_df$Description)))
    id_genes <- sapply(enrichGenes[id_order], length)
    V(g)$size <- id_genes
    V(g)$color <- enrich_df[id_order, legend_type]
    # igraph to ggplot
    p <- ggraph(g, layout)+
      geom_edge_link(alpha=.8, aes_(width=~I(width)),
                     colour='darkgrey')+
      ggnewscale::new_scale_fill() +
      geom_point(shape = 21, aes_(x =~ x, y =~ y,
                                  fill =~ color,
                                  size =~ size)) +
      scale_size_continuous(name = "Number of genes",
                            guide = "legend",
                            range = c(min(V(g)$size)/2,max(V(g)$size)/2) * scale_ratio) +
      scale_fill_continuous(low = top_color, high = bottom_color,
                            name = legend_type) +
      theme(panel.background = element_blank()) +
      geom_node_text(aes_(label=~name), data = NULL,
                     size = 3 * scale_ratio, bg.color = "white",
                     repel=TRUE, segment.size = 0.2)+
      guides(fill = guide_colorbar(reverse = TRUE))+
      plot_theme(remove_border = T,remove_text = T,border_thick = 0,...)

  }

  #--- goHeatmap/tangram/wordcloud ---#
  if(plot_type %in% c('goHeat','tangram','wordcloud')){
    if(!sim_method %in% c("Resnik", "Lin", "Rel", "Jiang" , "Wang")){
      stop('Please choose "sim_method" from: "Resnik", "Lin", "Rel", "Jiang" , "Wang"!')
    }
    l = get_sim_data(enrich_df,org=NULL,ont=NULL,sim_method)
    simMatrix = l[['m']]; reducedTerms = l[['r']]
    if(plot_type == 'goHeat'){
      p <- rrvgo::heatmapPlot(simMatrix,
                       reducedTerms,
                       annotateParent=TRUE,
                       annotationLabel="parentTerm",
                       fontsize=6)
    }
    if(plot_type == 'tangram'){
      p <- suppressMessages(rrvgo::treemapPlot(reducedTerms))
    }
    if(plot_type == 'wordcloud'){
      p <- suppressWarnings(rrvgo::wordcloudPlot(reducedTerms, min.freq=1, colors="black"))
    }
  }


  # wrap long text
  if (!is.null(wrap_length) & is.numeric(wrap_length)) {
    p <- p + scale_y_discrete(labels = text_wraper(wrap_length))
  }

  suppressWarnings(print(p))

}


#--- sub-function: wrap text if too long ---#
text_wraper <- function(width) {
  function(x) {
    lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n")
  }
}

ratio_intersect <- function(n,m){
  n <- unlist(n); m <- unlist(m)
  length(intersect(n, m))/length(unique(c(n,m)))
}

get_JC_data <- function(enrich_df){
  id <- enrich_df[,1]
  enrichGenes <- strsplit(enrich_df$geneID,'\\/') %>% setNames(id)

  n <- nrow(enrich_df)
  m <- matrix(NA, nrow=n, ncol=n)
  colnames(m) <- rownames(m) <- enrich_df$Description
  for (i in seq_len(n-1)) {
    for (x in (i+1):n) {
      m[i,x] <- ratio_intersect(enrichGenes[id[i]], enrichGenes[id[x]])
    }
  }
  return(m)
}

get_sim_data <- function(enrich_df,
                         org=NULL,
                         ont=NULL,
                         sim_method){

  nm <- strsplit(colnames(enrich_df)[1],"_") %>% unlist()
  if(all(tolower(nm) == 'id')){
    stop(paste0('Please give organism name to "org" such as "Hs", "Mm"...\n',
                'also give ontology to "ont" from "BP","CC" and "MF"'))
  }else{
    org =  nm[1]; ont = nm[2]
  }
  orgdb <- paste0("org.",org, ".eg.db")
  id <- enrich_df[, 1]
  # save godata is saving time
  data_dir = tools::R_user_dir('genekitr',which = 'data')
  data_dir = paste0(data_dir,"/godata")
  if(!dir.exists(data_dir)) dir.create(data_dir,recursive = TRUE)
  destfile =  paste0(data_dir, "/", org,"_",ont, "_godata.rda")

  if (!file.exists(destfile)) {
    orgdb <- loadOrgdb(orgdb)
    semdata = GOSemSim::godata(orgdb, ont = ont)
    save(semdata,file = destfile)
  }else{
    load(destfile)
  }
  # calculate simMatrix(m)
  id <- unique(id)
  found <- id %in% names(semdata@IC)
  getAncestors <- utils::getFromNamespace("getAncestors", "GOSemSim")
  hasAncestor <- !is.na(sapply(id, function(x) tryCatch(getAncestors(ont)[x],
                                                       error = function(e) NA)))
  id <- id[found & hasAncestor]
  m <- suppressMessages(matrix(GOSemSim::goSim(id, id, semData = semdata, measure = sim_method),
              ncol = length(id), dimnames = list(id, id)))
  out <- apply(m, 2, function(x) all(is.na(x)))
  m <- m[!out, !out]
  # reduce redundant terms
  scores <- setNames(-log10(enrich_df$qvalue), enrich_df[,1])
  r <- rrvgo::reduceSimMatrix(m,scores, orgdb=orgdb)
  return(list(m = m, r = r))
}

loadOrgdb <- function(orgdb){
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    stop("Bioconductor orgdb for ", orgdb, " not found. You should install first.",
         call. = FALSE)
  }
  eval(parse(text = paste0(orgdb, "::", orgdb)))
}

