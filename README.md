<img src="https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2022-05-24-043213.png" align="left" width="200"/>

[![CRANstatus](https://www.r-pkg.org/badges/version/genekitr)](https://cran.r-project.org/package=genekitr) [![](https://img.shields.io/badge/devel%20version-1.1.0-green.svg)](https://github.com/GangLiLab/genekitr) [![](https://cranlogs.r-pkg.org/badges/grand-total/genekitr?color=orange)](https://cran.r-project.org/package=genekitr) [![lifecycle](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) 

![Alt](https://repobeats.axiom.co/api/embed/e42ba06d30de893670c70324f19398ef0a7c26fa.svg "Repobeats analytics image")



## Overview

> **Genekitr** is a **gene** analysis tool**kit** based on **R**. 

![overview](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2022-07-26-142024.png)

#### Five core features:

- **Search**: gene-related information (exp. location, gene name, GC content, gene biotype ...) and PubMed records
- **Convert**: ID conversion among `Symbol & Alias`, `NCBI Entrez`, `Ensembl` ,`Uniprot` and `  Microarray probe`

- **Analysis**: user could select interested gene set from hundreds of gene sets for both model and non-model species, including [GO](http://geneontology.org/) (BP, CC and MF), [KEGG](https://www.kegg.jp/kegg/) (pathway, module, enzyme, network, drug and disease), [WikiPathways](https://wikipathways.org/), [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/), [EnrichrDB](https://maayanlab.cloud/Enrichr/), [Reactome](https://reactome.org/), [MeSH](https://www.ncbi.nlm.nih.gov/mesh/), [DisGeNET](https://www.disgenet.org/), [Disease Ontology](https://disease-ontology.org/) (DO), [Network of Cancer Gene](http://ncg.kcl.ac.uk/) (NCG) (version 6 and v7) and [COVID-19](https://maayanlab.cloud/covid19/). Gene enrichment analysis (GSA) contains both over representation analysis (ORA) and gene set enrichment analysis (GSEA) methods. ORA could support multi-group comparison.

- **Plot**: 13 ORA plots, 5 GSEA plots, 2 Venn plots  and 1 Volcano plot with flexible modification on text, color, border, axis and legend. All plot function input is dataframe format and supports GeneOntology web result. Feel free to make your own plots.
- **Export**: easily export multiple data sets as various sheets in one excel file

#### Supported organisms:

> For more details, please refer to [this site](https://genekitr.online/docs/species.html).

- Search & ID conversion:  **195** vertebrate species, **120** plant species and **2** bacteria species
- Enrichment analysis: GO supports **143** species, KEGG supports **8213** species, MeSH supports **71** species, MsigDb supports **20** species, WikiPahtwaysupports **16** species, Reactome supports **11** species, EnrichrDB supports 5 species and **human-specific** gene sets (DO, NCG, DisGeNET and COVID-19)



## 🛠 Installation

> New features are available for `version > 1.0.0`
> ```R
> # check current version
> packageVersion('genekitr') 
> ```

#### Install stable version from CRAN:

```R
install.packages("genekitr")
```

#### Install development version from GitHub:

```R
remotes::install_github("GangLiLab/genekitr")
```

#### Install development version from Gitee (for CHN mainland users):

```R
remotes::install_git("https://gitee.com/genekitr/pacakge_genekitr")
```



## 📚 Vignette

https://www.genekitr.fun/



## 🧙🏻‍♂️ Tell a story ~ why develop genekitr?

> **Genes are the basic omics research unit, just like cells in our body**. 

However, the issue of the gene is a little tedious.

Here, I want to tell you a story about Mr. Doodle, a computational biology student. **Now let us welcome our host Mr.Doodle to introduce his daily work with PI...**

### Scene 1: repeat work 

PI gave Doodle 30 genes and let him check their locations (better with sequences) and exact names. Doodle searched on NCBI one by one and copied & paste it into excel. Doodle sent the file to PI one hour later, and PI smiled, "Well done! Now I have another 50!" 

##### Doodle wonders how to avoid this repeat searching work?

### Scene 2: embarrassing name 

PI gave Doodle a DEG (differential expression analysis) matrix and a target gene list file. PI let him find if the target gene is up-regulated after treatment. After a while, Doodle found no PDL1 gene in the matrix but indeed exists in the gene list. "Do we have PDL1 gene?" he asked PI, and PI smiled, "Of course! You need to check gene CD274 instead of PDL1, which is an alias!"

<b style='color:#486CBE'> **Doodle was confused: how to distinguish between a real gene name and an alias?**</b>

### Scene 3: outdated database 

Doodle got the up-regulated gene symbols of the last DEG matrix to analyze KEGG. KEGG only supports Entrez id, so he needs to convert the symbol to Entrez. He found some symbols do not match Entrez id, but NCBI has. Doodle remembered he used org.db v3.12, but the current is v3.15. After he updated the annotation package, he finally got all matched IDs. 

<b style='color:#486CBE'>**Doodle wonders if there is any method to help him get updated results instead of self-check every time?**</b>

### Scene 4: imcompatible format 

PI did enrichment analysis alone on [the GeneOntology](http://geneontology.org/) website and let Doodle do visualization according to that result. "Could you please help plot the pathway bubble plot? Meanwhile, I want to show the x-axis as FoldEnrichment?" PI smiled. Doodle wanted to use the clusterProfiler R package for the plot, but he found it only accepts its object. So he bites the bullet and self-coding using ggplot2. 

<b style='color:#486CBE'>**Doodle wonders why does not have a tool that supports standard enrichment data frames?**</b>

### Scene 5: annoying plot theme 

Doodle finished the bubble plot at last and sent it to PI. After 15 minutes, PI sent him a message with a smile: "seems plot text size is too small, and could you give me a white background with border size 4 pt?" Doodle adjusted the ggplot theme function and modified 10 minutes. After a while, PI sent a message again, "I saw the second version; maybe the border is too thick. Could you replot?"

<b style='color:#486CBE'>**Doodle wonder if there is a function that could help him process the plot theme instead of changing the current code again and again?**</b>

### Scene 6: limited plot types

Once Doodle got the GO enrichment analysis result, PI let him think about how to show them nicely. Doodle found that every tool has its specific plot. For example, [WEGO](https://academic.oup.com/nar/article/46/W1/W71/4999241) could compare BP, CC, and MF terms; [GOplot](https://wencke.github.io/) has a chord plot to show the relationship of gene and GO terms; [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) support enriched map and network, which could explore the relationship among enrich terms. One big problem is that their input data is not compatible, so it is inconvenient to plot WEGO plots using clusterProfiler objects. 

<b style='color:#486CBE'>**Doodle wonder if there is any method that could involve beautiful plots from different tools with one universal data format?**</b>

### Scene 7: chaotic export files

Doodle has finished differential expression analysis and GO/KEGG enrichment analysis; PI let him send all result files to him. Doodle firstly saved all results into three excel files and named "DEG_data.xlsx," "GO_enrich.xlsx," and "KEGG_enrich.xlsx" then, he packed three files into one zipped folder and named them the date, finally, he sent to PI. After a while, PI sent him a message: "Could you put all three results into one excel file?" 

<b style='color:#486CBE'>**Doodle wonders if there is a way to save all data into one file without much manual operation?**</b>

> **If you have ever had one or more similar problems like Mr. Doodle, try `genekitr` !**



## ✍️ Author

[Yunze Liu](https://www.jieandze1314.com/)

[![](https://img.shields.io/badge/follow%20me%20on-WeChat-orange.svg)](https://genekitr.online/img/bioinfoplanet.png)



## 🔖 Citation

Wait to update...



## 💓 Welcome to contribute

If you are interested in `genekitr`, welcome contribute your ideas as follows:

* Git clone this project
* Double click `genekitr.Rproj` to open RStudio
* Modify source code in `R/` folder
* Run `devtools::check()` to make sure no errors, warnings or notes
* Pull request and describe clearly your changes

