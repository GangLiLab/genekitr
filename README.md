# Annotate genes 

> This tool annotates genes with alias, symbol, full name, function also related papers.
>
> 目的就是：与基因id相关的操作（如转换、可视化交集等）、分析（如富集分析），都可以加进来（后期考虑加入网页版）

## Table of Contents

-   [Installation](#installation)
-   [Features](#features)
-   [Plans](#plans)
-   [Let's mining data!](#lets_mining_data)
-   [Let's plot!](#lets_plot)
-   [Tips](#tips)

## Installation

You can also install devel version of **AnnoGenes** from github with:

``` r
# install.packages("remotes")
remotes::install_github("GangLiLab/AnnoGenes")
```

If you want to build vignette in local, please add two options:

``` r
remotes::install_github("GangLiLab/AnnoGenes", build_vignettes = TRUE, dependencies = TRUE)
```



## Features

- genecards虽然全，但是搜索数量有限制，于是整合了基因信息 =>`genInfo`
- 整合了相关的文献信息，可以自定义搜索关键词 => `genPubmed` 
- 每个操作都能得到一个数据框，可以继续探索，也可以作为不同的sheets导出到同一个excel => `expo_sheet`
- 有了基因的id和对应的logFC（需要排序好），就可以做GSEA => `getMsigdb +  genGSEA`
- 有了基因id，就能做GO分析 => `genGO ` 
- 有了基因id，就能做KEGG分析 => `genKEGG`
  - 默认富集分析`GO & KEGG`的结果为数据框，并且增加一列：`FoldEnrichment`
- **作图函数**

  - 气泡图 => `plotEnrichDot ` 

    



## Plans

- [ ] ~~图片也能导入excel（后期再看看这个有没有意义）~~
- [ ] 增加genVenn，先做成数据框结果。然后如果多于五组比较，就做成usetplot图
- [ ] genInfo增加基因位置，并且支持多个不同版本的基因组；可以整合orgdb中的所有信息；另外biomart也可以考虑整合进来
- [x] 图片的y轴label折叠（比如dotplot的y轴有很多的term，且长度不一，如果出现太长的term，最好可以折叠一下）=> `strwrap()`
- [ ] 设置自己的示例数据，like：`data(geneList, package="DOSE")`





## Let's mining data!

#### example gene id

```R
mm_id =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Bloc1s1') 
mm_id=str_split(mm_id,"\n")[[1]]
```

#### Method1: gene alias, full name

```R
test1 = genInfo(mm_id, org = 'mm')
```

rownames of `test1` are entrez ID

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-29-081721.png)

#### Method2: search pubmed 

```R
test2=genPubmed(mm_id, keywords = 'stem cell', field = 'tiab')
# Search example: Ticam2 [TIAB] AND stem cell [TIAB] 

# or use much specific keyword
genPubmed(mm_id, keywords = 'stem cell AND epithelial', field = 'tiab')
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-29-081925.png)

#### Method3: GSEA

```R
# 加载示例数据
data(geneList, package="DOSE")
# 获得msigdb的gene set
msigdb <- getMsigdb(org='human', category='C3',subcategory = 'TFT:GTRD')
# 直接进行gsea
egmt <- genGSEA(genelist = geneList,geneset = msigdb)
# 如果是extrez id，可以用下面的函数将id变成symbol
egmt2 <- DOSE::setReadable(egmt, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-07-02-100533.png)

#### Method4: GO

函数需要用到物种的`org.db`，如果没有相关物种注释包，函数内部的`auto_install()` 会帮助下载👍

```R
data(geneList, package="DOSE")
id = names(geneList)[1:100]
ego = genGO(id, org = 'human',ont = 'mf',pvalueCutoff = 0.05,qvalueCutoff = 0.1 ,use_symbol = T)
head(ego)
tmp=as.data.frame(ego)
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-07-02-035433.png)

**不知道物种名称？别怕！**

```R
> biocOrg_name()
    full_name short_name
1   anopheles         ag
2      bovine         bt
3        worm         ce
4      canine         cf
5         fly         dm
6   zebrafish         dr
7    ecolik12      eck12
8  ecolisakai    ecSakai
9     chicken         gg
10      human         hs
11      mouse         mm
12     rhesus        mmu
13      chipm         pt
14        rat         rn
15        pig         ss
16    xenopus         xl
```



#### Method5: map gene id

- `org` support many from `biocOrg_name()`

- user can choose output dataframe or not, using `return_dat`

```R
# 以下三种均可
mapId(id = ids, from = 'entrez', to = 'symbol',org='hs', return_dat = T)
mapId(id = ids, from = 'entrez', to = 'symbol',org='hg', return_dat = T)
mapId(id = ids, from = 'entrez', to = 'symbol',org='human', return_dat = T)
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-07-02-083549.png)



#### Method6: KEGG

```R
ids = names(geneList)[1:100]
gkeg <- genKEGG(ids, org = 'human') 
# org可以是common name（如human、mouse），也可以是hg、hs等常见的称呼

# 默认支持readable 参数，结果以symbol name展示
keg_raw <- genKEGG(test, org = 'hs', use_symbol = F)
keg_readable <- genKEGG(test, org = 'hs', use_symbol = T)
# 差别就是：
```

![image-20210702174030869](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-07-02-094031.png)



换个物种试试~理论上，**拿任意物种的symbol、entrez、ensembl基因，给函数投食即可**。不需要再提前进行id转换了

```R
# 小鼠基因为例
head(id)
[1] "Adora1"    "Insl3"    "AF067061"      "Alpk1"         "Arhgap20"      "B020004J07Rik" "Bmp6"
keg <- genKEGG(mm_id, org = 'mouse', use_symbol = T, pvalueCutoff = 1, qvalueCutoff = 1, maxGSSize = 3000)
```



## Let's plot!

#### P1: Enrichment dotplot =>  `plotEnrichDot ` 

- support dataframes with GO term, pvalue/qvalue/p.adjust, GeneRatio/Count/FoldEnrichment 
- Not only for result from R packages like `clusterProfiler` , but also for web analysis result like `panther ` from [Gene Ontology Resource](http://geneontology.org/) 
- 支持定义主图和legend的字体及大小、是否去除网格线、自定义渐变色的顶部和底部颜色、设定x轴起点、折叠y轴title、边框和刻度线宽度

```R
# test dataframe was from GeneOntology web result
p1 = plotEnrichDot(test, xlab_type =  'FoldEnrich', legend_by = 'qvalue',
              show_item = 15, main_text_size = 14,legend_text_size = 10,
              low_color = 'red', high_color = 'blue',
              xleft = 0, font_type = 'Arial', remove_grid = T,
              wrap_width = 30,border_thick = 3 )

# ego dataframe was from clusterP result
p2=plotEnrichDot(ego, xlab_type =  'GeneRatio', legend_by = 'p.adjust',
                 show_item = 10, main_text_size = 14,legend_text_size = 10,
                 low_color = 'orange', high_color = 'green',
                 xleft = 0, font_type = 'Times New Roman', remove_grid = F,
                 wrap_width = NULL ,border_thick = 1)

library(patchwork)
p1+p2
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-07-05-054512.png)









## Tips ~

- support pipe ` %>% ` 

```R
library(openxlsx)
wb <- createWorkbook()
wb <- expo_sheet(wb, sheet_dat = test1, sheet_name = 'genInfo') %>% 
  expo_sheet(., sheet_dat = test2, sheet_name = 'genPub') %>% 
saveWorkbook(wb, "~/Downloads/test.xlsx", overwrite = T)
```

<img src='man/figures/example1.png' align="below" />



- `genInfo` vs `bitr`  （后续`genInfo`可以扩展更多）

  <img src='man/figures/example2.png' align="below" />