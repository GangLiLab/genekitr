<img src="https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2022-05-24-043213.png" align="left" width="200"/>

[![CRANstatus](https://www.r-pkg.org/badges/version/genekitr)](https://cran.r-project.org/package=genekitr) [![](https://cranlogs.r-pkg.org/badges/grand-total/genekitr?color=orange)](https://cran.r-project.org/package=genekitr) [![lifecycle](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) 

![Alt](https://repobeats.axiom.co/api/embed/e42ba06d30de893670c70324f19398ef0a7c26fa.svg "Repobeats analytics image")



## Overview

**Genekitr** is a **gene** analysis tool**kit** based on **R**. It mainly contains five features:

- Search: gene-related information (exp. location, gene name, GC content, gene biotype ...) and PubMed records
- Convert: ID conversion among "Symbol & Alias", "NCBI Entrez", "Ensembl" and "Uniprot"
- Analysis: gene enrichment analysis including ORA (GO and KEGG) and GSEA which also supports multi-group comparison
- Plot: 14 GO plots, 7 KEGG plots, 5 GSEA plots and 2 Venn plots with flexible modification on text, color, border, axis and legend. Feel free to make your own plots.
- Export: easily export multiple data sets as various sheets in one excel file



## Tell a story ~ why develop genekitr?

No matter what omics data you research, **genes are the basic research unit just like cells in our body**. Genes issue is very common and a little tedious.

> Next, I want to tell you a story about Mr. Doodle who is a computational biology student.
> **Now let's welcome our host Mr.Doodle to introduce his daily work...**

### Scene 1: repeat work 

PI gave Doodle 30 gene and let him check their location (better with sequences) and detailed names. Doodle searched on NCBI one by one, copy & paste into excel. 1 hour later, Doodle sent the file to PI and PI smiled, "Well done! Now I have another 50!" 
**Doodle wonder how to avoid this repeat searching work?**

### Scene 2: embarrassing name 
PI gave Doodle a DEG (differential expression analysis) matrix and a target gene list file. PI let him find if target gene is up-regulated after treatment. After a while, Doodle found there is no PDL1 gene in the matrix but indeed exists in the gene list. "Do we actually have PDL1 gene?" he asked PI and PI smiled, "Of course! You need to check gene CD274 instead of PDL1, which is just an alias!" **Doodle was confused: how to distinguish from real gene name and alias?**

### Scene 3: outdated database 
Doodle got the up-regulated gene symbols of the last DEG matrix to do KEGG analysis. KEGG only support entrez id so firstly he needs to convert symbol to entrez. He found some symbol do not have matched entrezid, but NCBI has. Doodle remembered he used org.db v3.12 but current is v3.15. After he updated the annotation package, he finally got all matched id. 
**Doodle wonder if there's any method could help him get updated result instead of self-check every time?**

### Scene 4: imcompatible format 
PI did enrichment analysis alone on [GeneOntology](http://geneontology.org/) website and let Doodle do visulization according to that result. "Could you please help to plot the pathway bubble plot and also I want to show x-axis as FoldEnrichment." PI smiled. Doodle wanted to use clusterProfilter R package to plot but he found it only accepts its own object. So he just bite the bullet and self-coding using ggplot2. 
**Doodle wonder why don't have a tool supports common enrichment data frames?**

### Scene 5: annoying plot theme 

Doodle finished the bubble plot at last and sent to PI. After 15 minutes, PI sent him a message with a smile: "seems plot text size is too small and could you give me a white background with border size 4 pt?" Doodle ajusted `ggplot` theme function and modified 10 minutes. After a while, PI sent a message again, "I saw the second version, maybe border is too thick, could you replot?" 
**Doodle wonder if there's a function could help him process plot theme instead of changing current code again and again?**

### Scene 6: limited plot types

Once Doodle got GO enrichment analysis result, PI let him think how to show them nicely. Doodle found every tool has own specific plot e.g. [WEGO](https://wego.genomics.cn/) could  compare BP, CC and MF terms; [GOplot](https://wencke.github.io/)  has chord plot to show the relationship of gene and GO terms; [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) support enriched map and network which could explore the relationship among enrich terms. One big problem is their input data is not compatible so it is unconvinient to plot WEGO plot using clusterProfiler object. **Doodle wonder if there's any method could involve beautiful plots from different tool with one universal data format?**

### Scene 7: chaotic export files

Doodle has finished differential expression analysis and GO/KEGG enrichment analysis, PI let him sent all result files to him. Doodle firstly save all results into three excel files and named "DEG_data.xlsx", "GO_enrich.xlsx" and "KEGG_enrich.xlsx", then he packed three files into one zipped folder and named with date, finally he sent to PI. After a while, PI sent him a message: "Could you put all three result into one excel file?" 
**Doodle wonder if there's way to save all data into one file without too much manual operation?**

> **If you have ever had one or more similar problems like Mr. Doodle, you may need `genekitr` !**

## 🛠 Installation

#### Install CRAN stable version:

```R
install.packages("genekitr")
```

#### Install GitHub development version:

```R
remotes::install_github("GangLiLab/genekitr")
```

#### Install Gitee (for CHN mainland users):

```R
remotes::install_git("https://gitee.com/genekitr/pacakge_genekitr")
```



## 📚 Vignettes





## ✍️ Authors

[Yunze Liu](https://www.jieandze1314.com/)

[![](https://img.shields.io/badge/follow%20me%20on-WeChat-orange.svg)](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2022-05-24-015641.png)



## 🔖 Citation

Wait to update...



## 💓 Welcome to contribute

If you are interested in `genekitr`, welcome contribute your ideas as follows:

* Git clone this project
* Double click `genekitr.Rproj` to open RStudio
* Modify source code in `R/` folder
* Run `devtools::check()` to make sure no errors, warnings or notes
* Pull request and describe clearly your changes

