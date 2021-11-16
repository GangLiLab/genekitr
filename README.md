---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit Rmd file -->



# genekitr

**genekitr** is an R analysis toolkit based on the gene. It mainly contains five features:

- Search: Gene IDs as input then get gene-related information (exp. location, gene name, gene alias, GC content ...) as well as search gene-related PubMed records

- Transform: Transform gene ID type among "symbol", "entrezid", "ensembl" and "uniprot"

- Analysis: Gene enrichment analysis including ORA (GO and KEGG) and GSEA

- Visualize: Visualization for enrichment analysis and gene overlaps

- Export: Gene IDs and analysis results could be exported as various sheets in one Excel file, which could be easily read and shared with others

**Why develop this R package?**


用户痛点：

- ID转换不够全，很多gene alias 不能够被识别（比如huamn 的BCC7实际对应TP53，mouse的Tp53实际可能对应Trp53, Trp53inp2 & Ano9）;
  另外，像是一些常见的gene symbol（如PD1、PDL1），实际上它们的gene id是PDCD1和CD274，很多表达矩阵中可能因此就忽略掉这些基因

- 基因ID转换支持物种太少，有些工具只支持模式物种。我们可以支持198个物种

- 目前很多函数的参数设置不够友好，太繁琐（我们允许用户只输入基因id，至于它是symbol、entrez还是其他，函数可以自行判断）

- 保持和数据库的同步（ensembl、Uniprot），使用最新的ensemble v104版本（后期还可同步更新）
  然后protein id也随着uniprot数据库更新

- 很多函数导出的结果是一个对象（Object），而不是数据框（dataframe），如果用户想自行挑选terms，就比较麻烦；为此，我们的函数全面支持数据框，直接将数据框给到作图函数，便可轻松绘制，而且绘制的图形也尽可能简洁，比较贴近发表级，后期用户只需要简单AI修改即可


## Table of Contents

* [Installation](#installation)
* [Quick guide](#quick-guide)
* [Vignette](#vignette)
* [Citation](#citation)
* [Welcome to contribute](#welcome-to-contribute)

## Installation

Install CRAN stable version:


```r
install.packages("genekitr")
```

Install GitHub dev version:


```r
# install.packages("remotes")
remotes::install_github("GangLiLab/genekitr")
# To build local vignette:
# remotes::install_github("GangLiLab/genekitr", build_vignettes = TRUE, dependencies = TRUE)
```

## Quick guide

To quickly go through the package usage, we will use built-in gene list from GEO airway

([GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)) DEG analysis.

### Search 

##### Search gene related information

`org` argument could accept fullname or shortname of specific organism. For example, we could 

use common name: `human/hs/hg` to describe human.


```r
library(genekitr)
data(geneList)
id = names(geneList)[1:100]
head(id)
#> [1] "2847"   "148145" "1591"   "2903"   "26045"  "10268"
ginfo = genInfo(id, org = 'human')
dplyr::as_tibble(ginfo)
#> # A tibble: 100 × 21
#>    input_id symbol    ensembl  uniprot  chr   start end   width strand gene_name ncbi_alias ensembl_alias gc_content gene_biotype
#>    <ord>    <chr>     <chr>    <chr>    <chr> <chr> <chr> <chr> <chr>  <chr>     <chr>      <chr>         <chr>      <chr>       
#>  1 2847     MCHR1     ENSG000… Q99705;… 22    4067… 4068… 4065  1      melanin … GPR24; MC… GPR24; MCH1R… 55.82      protein_cod…
#>  2 148145   LINC00906 ENSG000… <NA>     19    2888… 2897… 87445 1      long int… <NA>       <NA>          45.05      lncRNA      
#>  3 1591     CYP24A1   ENSG000… Q07973   20    5415… 5417… 20541 -1     cytochro… CP24; CYP… CP24; CYP24;… 43.37      protein_cod…
#>  4 2903     GRIN2A    ENSG000… Q12879;… 16    9753… 1018… 4295… -1     glutamat… EPND; FES… GluN2A; NMDA… 42.06      protein_cod…
#>  5 26045    LRRTM2    ENSG000… O43300;… 5     1388… 1388… 6448  -1     leucine … <NA>       KIAA0416      38.26      protein_cod…
#>  6 10268    RAMP3     ENSG000… O60896;… 7     4515… 4518… 28512 1      receptor… <NA>       <NA>          51.29      protein_cod…
#>  7 114      ADCY8     ENSG000… P40145;… 8     1307… 1310… 2606… -1     adenylat… AC8; ADCY… AC8; ADCY3; … 40.42      protein_cod…
#>  8 79805    VASH2     ENSG000… Q86V25;… 1     2129… 2129… 41518 1      vasohibi… <NA>       FLJ12505      48.23      protein_cod…
#>  9 6563     SLC14A1   ENSG000… Q13336;… 18    4568… 4575… 65496 1      solute c… HUT11; Hs… HsT1341; JK;… 43.44      protein_cod…
#> 10 6355     CCL8      ENSG000… P80075;… 17    3431… 3432… 1968  1      C-C moti… HC14; MCP… HC14; MCP-2;… 44.16      protein_cod…
#> # … with 90 more rows, and 7 more variables: transcript_count <chr>, hgnc_id <chr>, omim <chr>, ccds <chr>, reactome <chr>,
#> #   ucsc <chr>, mirbase_id <chr>
```

##### Search gene PubMed records


```r
pub = genPubmed(
  id = c("Cyp2c23", "Fhit", "Gal3st2b","Insl3", "Gbp4"),
  keywords = "stem cell", field = "tiab")
#> Search example: Cyp2c23 [TIAB] AND stem cell [TIAB]
head(pub)
#>      gene
#> 1 Cyp2c23
#> 2    Fhit
#> 3    Fhit
#> 4    Fhit
#> 5    Fhit
#> 6    Fhit
#>                                                                                                                                                                          title
#> 1                                                                                                                                                                           NA
#> 2                                                                Changes in Methylation Patterns of Tumor Suppressor Genes during Extended Human Embryonic Stem Cell Cultures.
#> 3                                                                           Methylation status of the <i>FHIT</i> gene in the transformed human mesenchymal F6 stem cell line.
#> 4                                                                Totipotent stem cells bearing del(20q) maintain multipotential differentiation in Shwachman Diamond syndrome.
#> 5                                                                         Fhit-deficient hematopoietic stem cells survive hydroquinone exposure carrying precancerous changes.
#> 6 Induction by 7,12-dimethylbenz(a)anthracene of molecular and biochemical alterations in transformed human mammary epithelial stem cells, and protection by N-acetylcysteine.
#>         date                              doi     pmid                           journal
#> 1         NA                               NA       NA                                NA
#> 2 2021_09_25             10.1155/2021/5575185 34552632          Stem cells international
#> 3 2020_10_01             10.3892/ol.2015.3092 26137124                  Oncology letters
#> 4 2009_02_20 10.1111/j.1365-2141.2008.07448.x 19016724    British journal of haematology
#> 5 2008_06_24    10.1158/0008-5472.CAN-07-5687 18483248                   Cancer research
#> 6 2006_10_05                                  16865267 International journal of oncology
```


### Transform

Regardless of input ID type, function will detect automatically.

User only need to specify which type you want, the left things just give `transId`.

##### Transform ID from entrezid to symbol


```r
id[1:5]
#> [1] "2847"   "148145" "1591"   "2903"   "26045"
transId(id[1:5], trans_to = 'symbol',org='hs')
#> 
#> 100% genes are mapped from entrezid to symbol
#> [1] "MCHR1"     "LINC00906" "CYP24A1"   "GRIN2A"    "LRRTM2"
```

##### Transform gene alias

```r
transId(c('BCC7','PDL1','PD1'), trans_to = 'symbol',org='human')
#> 
#> 100% genes are mapped from symbol to symbol
#> [1] "TP53"  "CD274" "PDCD1"

# When one id have multi alias
transId(c('BCC7','PDL1','PD1'), trans_to = 'symbol',org='hg',unique = F)
#> Some ID occurs one-to-many match, like "PD1"
#> If you want to get one-to-one match, please set "unique=TRUE"
#> 
#> 100% genes are mapped from symbol to symbol
#>   input_id symbol
#> 1     BCC7   TP53
#> 2     PDL1  CD274
#> 3      PD1  PDCD1
#> 4      PD1   SNCA
#> 5      PD1 SPATA2
```

##### Transform protein id to gene id

```r
# to symbol
transId(c('Q12879','Q86V25','Q8N386','Q5T7N3'),'sym','hs')
#> 
#> 100% genes are mapped from uniprot to symbol
#> [1] "GRIN2A" "VASH2"  "LRRC25" "KANK4"
# to ensembl
transId(c('Q12879','Q86V25','Q8N386','Q5T7N3'),'ens','hs')
#> 
#> 100% genes are mapped from uniprot to ensembl
#> [1] "ENSG00000183454" "ENSG00000143494" "ENSG00000175489" "ENSG00000132854"
```

##### Transform ensembl id to symbol 

```r
transId(c('ENSG00000146006','ENSG00000134321','ENSG00000136267','ENSG00000105989'),'sym','hg')
#> 
#> 100% genes are mapped from ensembl to symbol
#> [1] "LRRTM2" "RSAD2"  "DGKB"   "WNT2"
```




### Analysis 

All enrichment analysis just give a gene list (especially GSEA need the gene list with a decreasing fold change)

#### over representation analysis (ORA)

##### GO analysis 

User could choose ontology among "bp", "cc", "mf" or "all".

If set `use_symbol = TRUE`, the result will return gene symbol for easily understand.

If you are not sure the organism name, please type `biocOrg_name` and choose full name or short name


```r
ego_bp = genGO(id[1:100], org = 'human', ont = 'bp',pvalueCutoff = 0.05,qvalueCutoff = 0.05, use_symbol = T)

# user could choose all three ontology
ego_all = genGO(id[1:100], org = 'human', ont = 'all',pvalueCutoff = 0.05,qvalueCutoff = 0.05, use_symbol = T)
```

##### KEGG analysis 

If you are not sure the organism name, please type `keggOrg_name` and choose full name or short name


```r
ekeg = genKEGG(id, org = 'hg',use_symbol = T)

# If we want to compare different groups of genes, we only need to give a gene group list
group = list(groupA = c(rep("up",50),rep("down",50)),
             groupB = c(rep("A",30), rep("B",70)))
ekeg_compare =  genKEGG(c(head(id,50),tail(id,50)), group_list = group,
                    org = 'human',pvalueCutoff = 0.15,qvalueCutoff = 0.15, use_symbol = T)
```

#### gene set enrichment analysis (GSEA) analysis 

`category` argument could choose from 'C1','C2','C3', 'C4','C5','C6','C7','C8' and 'H'

If you are not sure `subcategory`, you can only choose `category` and leave `subcategory` as blank. 

The message will tell what options you could choose in the main `category`.


```r
egsea = genGSEA(genelist = geneList,org = 'hs', 
                category = "H",
                use_symbol = T, pvalueCutoff = 1)
#> H has no subcategory, continue...
```


### Visulaize 

##### Gene overlap

If we only have two or three groups of genes, the function will plot Venn diagram;

If we have at least four groups of genes, the default option will be UpSet diagram.


```r
# if only have three groups
set1 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
set2 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
set3 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
sm_gene_list <- list(gset1 = set1, gset2 = set2, gset3 = set3)
plotVenn(sm_gene_list,
  text_size = 1.5, alpha_degree = 1,
  remove_grid = TRUE)
#> Color length should be same with venn_list, auto assign colors...
```

![plot of chunk venn](man/figures/venn-1.png)

```r

# if only have five groups
set4 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
set5 <- paste0(rep("gene", 100), sample(c(1:1000), 100))
la_gene_list <- list(gset1 = set1, gset2 = set2, gset3 = set3,
  gset4 = set4, gset5 = set5)
plotVenn(la_gene_list,
  text_size = 15, alpha_degree = 0.2, border_thick = 2,
  remove_grid = TRUE, use_venn = FALSE)
```

![plot of chunk venn](man/figures/venn-2.png)

##### Enrichment plot

- Support barplot, dotplot

- x-axis support: GeneRatio/Count/FoldEnrich

- stats support: pvalue/p.adjust/qvalue(FDR)

- easily modify plot line & text pattern


```r
plotEnrich(ego_bp,plot_type = 'bar',remove_grid = T, main_text_size = 8,
  legend_text_size = 6,border_thick = 1.5)
```

![plot of chunk enrich_bar1](man/figures/enrich_bar1-1.png)

```r
plotEnrich(ego_all,plot_type = 'bar', xlab_type = 'GeneRatio',legend_type = 'qvalue',
           remove_grid = T, main_text_size = 8, legend_text_size = 6,border_thick = 1.5)
```

![plot of chunk enrich_bar2](man/figures/enrich_bar2-1.png)

##### Enrichment dotplot


```r
plotEnrich(ego_bp,plot_type = 'dot',xlab_type = 'GeneRatio',legend_type = 'qvalue',
          remove_grid = T, main_text_size = 8, legend_text_size = 6)
```

![plot of chunk enrich_dot](man/figures/enrich_dot-1.png)


```r
plotEnrich(ekeg_compare)
```

![plot of chunk enrich_dot2](man/figures/enrich_dot2-1.png)

##### gsea plot

- support: volcano/classic/multi-pathway(fgsea)


```r
plotGSEA(egsea,plot_type = 'volcano', show_pathway = 3)
```

![plot of chunk gsea](man/figures/gsea-1.png)

```r
plotGSEA(egsea, plot_type = 'classic', show_pathway = c("HALLMARK_UV_RESPONSE_DN","HALLMARK_INTERFERON_ALPHA_RESPONSE"),
         show_genes = c("SELL"))
```

![plot of chunk gsea](man/figures/gsea-2.png)

```r

# default shows top3 up & down pathways
plotGSEA(egsea, plot_type = 'fgsea')
```

![plot of chunk gsea](man/figures/gsea-3.png)


### Export 

If you want to export many data sets in one file, you could use `expoSheet`

For example, since we got 100 genes' GO and KEGG result, then we want to export them with gene information:


```r
expoSheet(
  dat_list = list(ginfo, ego_bp, ego_all,ekeg), name_list = list("GeneInfo","GO-BP", "GO-All","KEGG"),
  filename = "gene_enrich.xlsx", dir = tempdir())
```

The result will be:

<img src='man/figures/exp0.png' height="500" alt="exp0"/>


## Vignette

### English 

- Wait to update...

### 中文区 (推文专区)

- 待更新1：R包搞定你的基因需求


## Citation

Wait for paper ...


## Welcome to contribute

If you are interested in this tool, welcome contribute your ideas as follows:

* Git clone this project
* Double click `genekitr.Rproj` to open RStudio
* Modify source code in `R/` folder
* Run `devtools::check()` to make sure no errors, warnings or notes
* Pull request and describe clearly your changes








