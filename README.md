
<!-- README.md is generated from README.Rmd. Please edit that file -->

# genekitr

**genekitr** is an R analysis toolkit based on the gene. It mainly
contains five features:

-   Search: Gene IDs as input then get gene-related information (exp.
    location, gene name, gene alias, GC content …) as well as search
    gene-related PubMed records

-   Transform: Transform gene ID type among “symbol”, “entrezid”,
    “ensembl” and “uniprot”

-   Analysis: Gene enrichment analysis including ORA (GO and KEGG) and
    GSEA

-   Visualize: Visualization for enrichment analysis and gene overlaps

-   Export: Gene IDs and analysis results could be exported as various
    sheets in one Excel file, which could be easily read and shared with
    others

**Why develop this R package?**

用户痛点：

-   ID转换不够全，很多gene alias 不能够被识别（比如huamn
    的BCC7实际对应TP53，mouse的Tp53实际可能对应Trp53, Trp53inp2 & Ano9）

-   参数设置不够友好，太繁琐（我们允许用户只输入基因id，至于它是symbol、entrez还是其他，函数可以自行判断）

-   保持和数据库的同步（ensembl、Uniprot），很多MassSpec数据中uniprot
    ID对应的symbol存在错误或者不完整，我们可以尽量减少

-   …

## Table of Contents

-   [Installation](#installation)
-   [Quick guide](#quick-guide)
-   [Vignette](#vignette)
-   [Citation](#citation)
-   [Welcome to contribute](#welcome-to-contribute)

## Installation

Install CRAN stable version:

``` r
install.packages("genekitr")
```

Install GitHub dev version:

``` r
# install.packages("remotes")
remotes::install_github("GangLiLab/genekitr")
# To build local vignette:
# remotes::install_github("GangLiLab/genekitr", build_vignettes = TRUE, dependencies = TRUE)
```

## Quick guide

## Vignette

## Citation

Wait for paper …

## Welcome to contribute

If you are interested in this tool, welcome contribute your ideas as
follows:

-   Git clone this project
-   Double click `genekitr.Rproj` to open RStudio
-   Modify source code in `R/` folder
-   Run `devtools::check()` to make sure no errors, warnings or notes
-   Pull request
