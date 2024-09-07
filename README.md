<img src="https://github.com/GangLiLab/genekitr/assets/34302878/b7337c2c-076a-46e8-9777-b2dfb742d613" align="left" width="200"/>

[![CRANstatus](https://www.r-pkg.org/badges/version/genekitr)](https://cran.r-project.org/package=genekitr) [![](https://img.shields.io/badge/devel%20version-1.2.7-green.svg)](https://github.com/GangLiLab/genekitr) [![](https://cranlogs.r-pkg.org/badges/grand-total/genekitr?color=orange)](https://cran.r-project.org/package=genekitr) [![lifecycle](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) 

![Alt](https://repobeats.axiom.co/api/embed/e42ba06d30de893670c70324f19398ef0a7c26fa.svg "Repobeats analytics image")



## Overview

> **Genekitr** is a **gene** analysis tool**kit** based on **R**. 

<img width="1328" alt="image" src="https://github.com/GangLiLab/genekitr/assets/34302878/c60a5fb0-d2da-4a56-8c58-8204b6769adb">

#### Five core features:

- **Search**: gene-related information (exp. gene functional summary, gene name, location,  GC content, gene biotype ...) and PubMed records
- **Convert**: ID conversion among `Symbol & Alias`, `NCBI Entrez`, `Ensembl` ,`Uniprot` and `  human microarray probe`

- **Analysis**: users could select interested gene set from hundreds of gene sets for both model and non-model species, including [GO](https://geneontology.org/) (BP, CC and MF), [KEGG](https://www.kegg.jp/kegg/) (pathway, module, enzyme, network, drug and disease), [WikiPathways](https://www.wikipathways.org/), [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/), [EnrichrDB](https://maayanlab.cloud/Enrichr/), [Reactome](https://reactome.org/), [MeSH](https://www.ncbi.nlm.nih.gov/mesh/), [DisGeNET](https://disgenet.com/), [Disease Ontology](https://www.disease-ontology.org/) (DO), Network of Cancer Gene (NCG) (version 6 and v7) and [COVID-19](https://maayanlab.cloud/covid19/). Gene enrichment analysis (GSA) contains both over representation analysis (ORA) and gene set enrichment analysis (GSEA) methods. ORA is capable of supporting multi-group comparisons.

- **Plot**: easily generate 13 ORA plots, 5 GSEA plots, 2 Venn plots, and 1 Volcano plot with customizable features such as text, color, border, axis, and legend. The function is capable of accepting a dataframe as input and supports GeneOntology website results based on [PantherDB](http://www.pantherdb.org/)..
- **Export**: quickly export numerous datasets as different sheets within a single Excel file.

#### Supported organisms:

> For more details, please refer to [this site](https://genekitr.online/docs/species.html).

- Search & ID conversion:  **195** vertebrate species, **120** plant species and **2** bacteria species
- Enrichment analysis: GO supports **143** species, KEGG supports **8213** species, MeSH supports **71** species, MsigDb supports **20** species, WikiPahtwaysupports **16** species, Reactome supports **11** species, EnrichrDB supports 5 species and **human-specific** gene sets (DO, NCG, DisGeNET and COVID-19)



## 🛠 Installation

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

**ENGLISH:** https://www.genekitr.fun/



## 🧙🏻‍♂️ Tell a story ~ why develop genekitr?

> **Genes, the essence of life's art,** 
> **Omics research's fundamental part,** 
>
> **Like cells in our physical frame,** 
> **Their study reveals life's vibrant flame.**

Let me tell you a story about Mr. Doodle, a computational biology student working with his PI.

### Scene 1: repeat work 

One day, PI gave him 30 genes to check for their locations and exact names, preferably with sequences.

Mr. Doodle searched for each gene on NCBI, copying and pasting the information into an Excel sheet. He sent the file to PI an hour later, and received praise for his work. But just when he thought he was done, PI gave him another 50 genes to check!

Despite feeling a little overwhelmed, Mr. Doodle repeated the same process with determination, determined to complete the task to the best of his abilities.

##### <b style='color:#486CBE'>Doodle wondered how to avoid having to repeatedly search for the same information.?</b>

### Scene 2: embarrassing name 

Once upon a time, PI gave Mr. Doodle a DEG matrix and a target gene list file. The task was to determine if the target gene was up-regulated after treatment.

Mr. Doodle searched the matrix but couldn't find the PDL1 gene, even though it was in the gene list. He asked PI about it, and PI explained that the gene was listed as CD274, which is an alias for PDL1.

This left Mr. Doodle feeling a little confused. He wondered how to distinguish between real gene names and aliases.

<b style='color:#486CBE'> **Doodle wondered  how to differentiate between a real gene name and an alias.**</b>

### Scene 3: outdated database 

Mr. Doodle was analyzing KEGG pathways for the up-regulated genes in the last DEG matrix. However, KEGG only supported Entrez IDs, and the genes were listed by their symbols.

Mr. Doodle needed to convert the gene symbols to Entrez IDs, but he found that some symbols did not match the corresponding Entrez IDs. However, he discovered that NCBI had the correct IDs.

Mr. Doodle realized that he was using an outdated org.Hs.eg.db v3.15 annotation package. After updating to the current version, v3.17, he was finally able to obtain all the matched IDs, and continue his analysis of the KEGG pathways.

<b style='color:#486CBE'>**Doodle wondered if there was a method to help him obtain updated results automatically, instead of having to check them manually every time.**</b>

### Scene 4: imcompatible format 

PI did some fancy enrichment analysis all by himself on a website called [GeneOntology](https://geneontology.org/).  He then asked Mr.Doodle to help him make a pretty picture of the results. . "Can you make a bubble plot for me and show the FoldEnrichment on the x-axis?" he asked with a smile. Doodle tried to use a fancy R package called clusterProfiler, but it wouldn't work with the data. So, he bravely coded it himself using ggplot2.

<b style='color:#486CBE'>**Doodle wondered why there isn't a tool that supports easy data frames.**</b>

### Scene 5: annoying plot theme 

Doodle finally finished making the bubble plot and sent it to PI. After 15 minutes, PI sent him a message with a smile: "the text is too small, and can you make the background white with a border size of 4 points?" Doodle tweaked the ggplot theme and made the changes in 10 minutes. But, after a little while, PI sent another message saying, "The border is too thick in the second version. Can you please redo it?"

<b style='color:#486CBE'>**Doodle wondered if there was a function that could help him process the plot theme instead of having to modify the current code repeatedly.**</b>

### Scene 6: limited plot types

PI gave Doodle the GO enrichment analysis result and asked him to think of a creative way to display it. Doodle found that each tool had its specific plot. For example, [WEGO](https://wego.genomics.cn/) could compare BP, CC, and MF terms; [GOplot](https://wencke.github.io/) had a chord plot to show the relationship between genes and GO terms; and [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)supported enriched map and network, which could explore the relationship among enriched terms. However, there was a big problem - the input data for each tool was not compatible, making it inconvenient to plot WEGO plots using clusterProfiler objects.

<b style='color:#486CBE'>**Doodle wondered if there was a method that could produce beautiful plots from different tools using a universal data format.**</b>

### Scene 7: chaotic export files

Doodle finished conducting differential expression analysis and GO/KEGG enrichment analysis. PI asked him to send over all the result files. Doodle saved the results into three separate excel files, naming them "DEG_data.xlsx," "GO_enrich.xlsx," and "KEGG_enrich.xlsx." He then compressed the three files into a zipped folder, naming it after the date, and sent it to PI. After a while, PI asked him if he could put all three results into a single excel file.

<b style='color:#486CBE'>**Doodle wondered if there is a way to save all data into a single file without having to perform many manual operations.?**</b>

> **If you have encountered similar problems like Mr. Doodle, give <u>genekitr</u> a try!**



## ✍️ Author

[Yunze Liu](https://www.jieandze1314.com/)

[![](https://img.shields.io/badge/follow%20me%20on-WeChat-orange.svg)](https://genekitr.online/img/bioinfoplanet.png)



## 🔖 Citation

> **For now, the paper is published. Please cite:**

Liu, Y., Li, G. Empowering biologists to decode omics data: the Genekitr R package and web server. BMC Bioinformatics 24, 214 (2023). https://doi.org/10.1186/s12859-023-05342-9





## 💓 Welcome to contribute

If you are interested in `genekitr`, welcome contribute your ideas as follows:

* Git clone this project
* Double click `genekitr.Rproj` to open RStudio
* Modify source code in `R/` folder
* Run `devtools::check()` to make sure no errors, warnings or notes
* Pull request and describe clearly your changes

