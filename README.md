<img src="https://github.com/GangLiLab/genekitr/assets/34302878/b7337c2c-076a-46e8-9777-b2dfb742d613" align="right" width="160" alt="Genekitr logo"/>

# genekitr

**A gene analysis toolkit for R: search, convert, analyze, plot, and export.**

[![CRAN version](https://www.r-pkg.org/badges/version/genekitr)](https://cran.r-project.org/package=genekitr)
[![Development version](https://img.shields.io/badge/devel%20version-1.3.0-green.svg)](DESCRIPTION)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/genekitr?color=orange)](https://cran.r-project.org/package=genekitr)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

[Documentation](https://www.genekitr.fun/) · [Quick start](#quick-start) · [Changelog](ChangeLog) · [Report an issue](https://github.com/GangLiLab/genekitr/issues) · [Citation](#citation)

## Overview

**Genekitr** brings gene annotation, ID conversion, enrichment analysis, visualization, and Excel export into one workflow. Enrichment results use data frames, with import helpers for clusterProfiler, PANTHER, and ShinyGO.

| Task | Functions | What you can do |
| --- | --- | --- |
| Search | `genInfo()`, `getPubmed()` | Retrieve gene annotations and search literature abstracts through Europe PMC. |
| Convert | `transId()`, `transProbe()` | Convert gene symbols and aliases, Entrez, Ensembl, UniProt, and supported human microarray probe IDs. |
| Analyze | `genORA()`, `genGSEA()`, `simGO()` | Run over-representation analysis (ORA), compare gene groups, run ranked gene set enrichment analysis (GSEA), and simplify GO terms. |
| Plot | `plotEnrich()`, `plotEnrichAdv()`, `plotGSEA()`, `plotVenn()`, `plotVolcano()`, `plot_theme()` | Draw enrichment plots, compare up/down enrichment, visualize gene overlaps and differential expression, and customize themes. |
| Import | `importCP()`, `importPanther()`, `importShinygo()`, `as.enrichdat()` | Adapt enrichment results from other tools for plotting. |
| Export | `expoSheet()` | Save multiple tables as separate worksheets in one Excel workbook. |

<img width="1328" alt="Overview of the Genekitr workflow" src="https://github.com/GangLiLab/genekitr/assets/34302878/c60a5fb0-d2da-4a56-8c58-8204b6769adb">

## Installation

Requires **R >= 4.0.0**, **dplyr >= 1.1.0**, and **ggplot2 >= 3.5.0**. Install the Bioconductor dependency `clusterProfiler` before installing genekitr; CRAN dependencies are installed automatically.

### Stable version from CRAN

```r
install.packages("BiocManager")
BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
install.packages("genekitr")
```

### Development version from GitHub

After installing the Bioconductor dependency above:

```r
install.packages("remotes")
remotes::install_github("GangLiLab/genekitr")
```

Alternatively, use the Gitee mirror (for users in mainland China):

```r
remotes::install_git("https://gitee.com/genekitr/pacakge_genekitr")
```

Check your installed version with `packageVersion("genekitr")`. See the [changelog](ChangeLog) for development changes.

Some plot types and import helpers need optional packages. For example, `importPanther()` and `importShinygo()` use `rio`; GO semantic similarity uses `GOSemSim`, `rrvgo`, and an organism annotation package such as `org.Hs.eg.db`. Install these as needed:

```r
install.packages("rio")
BiocManager::install(c("GOSemSim", "rrvgo", "org.Hs.eg.db"))
```

## Quick start

### Run ORA and plot the results

This example uses the bundled human `geneList` dataset, a named numeric vector of fold changes with Entrez IDs as names. The cutoff below selects an example gene subset; use your own differential expression criteria for real data.

```r
library(genekitr)

data("geneList", package = "genekitr")
id <- names(geneList)[abs(geneList) > 1.5]

# Retrieve GO molecular function gene sets from the geneset package.
gs <- geneset::getGO(org = "human", ont = "mf")

# Use the genes measured in this example as the background.
ora <- genORA(id = id, geneset = gs, universe = names(geneList))
head(ora)

# Select the terms to display before plotting.
plotEnrich(
  head(ora, 10),
  plot_type = "bubble",
  term_metric = "FoldEnrich",
  stats_metric = "p.adjust",
  main_text_size = 10
)
```

`genORA()` returns a data frame. Set `universe` to the genes eligible for selection in your experiment; if omitted, the default background is the genes in the gene set collection. For group comparisons, supply `group_list`, with one label per input gene in each list element; grouped results are displayed as a dot plot by `plotEnrich()`.

### Run GSEA with the full ranked list

Use a **named numeric vector sorted in decreasing order**, with gene IDs as names and ranking statistics as values. Keep the full ranked list rather than only significant genes.

```r
# Continue with geneList and gs from the ORA example.
ranked_genes <- sort(geneList, decreasing = TRUE)
gsea <- genGSEA(genelist = ranked_genes, geneset = gs, set_seed = TRUE)

head(gsea$gsea_df)
plotGSEA(gsea, plot_type = "bar", label_by = "description")
```

`genGSEA()` returns a **list** containing `gsea_df`, `genelist`, `geneset`, `exponent`, and `org`. Pass the complete list to `plotGSEA()` and use `gsea$gsea_df` for tables and Excel export. The results distinguish leading-edge genes (`core_enriched_geneID`, `core_enriched_count`) from all ranked genes overlapping each set (`all_enriched_geneID`, `all_enriched_count`). If ID conversion maps several entries to the same gene, only the first entry in the ranked list is kept.

### Export analysis tables to Excel

```r
expoSheet(
  data_list = list(ora, gsea$gsea_df),
  data_name = c("ORA", "GSEA"),
  filename = "genekitr_results.xlsx",
  dir = tempdir(),
  overwrite = FALSE
)
```

The workbook is saved to `file.path(tempdir(), "genekitr_results.xlsx")`. Change `dir` to your output folder. Supply one sheet name per table; `expoSheet()` creates the directory if needed.

### Look up genes and convert IDs

```r
ids <- c("TP53", "CD274", "BRCA1")

info <- genInfo(id = ids, org = "human", unique = TRUE)
converted <- transId(
  id = ids,
  transTo = c("entrez", "ensembl"),
  org = "human",
  unique = TRUE,
  keepNA = TRUE
)
```

`genInfo()` and `transId()` recognize gene symbols and aliases without regard to case and accept versioned Ensembl IDs. `unique = FALSE` preserves one-to-many mappings; `unique = TRUE` keeps one record per input ID, preferring records with fewer missing fields. `keepNA = TRUE` retains unmatched IDs. For human annotations, select `hgVersion = "v38"` (default) or `"v19"`.

Annotation and gene set retrieval need an internet connection on first use. Gene annotations are cached under `tools::R_user_dir("genekitr", which = "data")`; `geneset` manages its own cache. Annotation files are reused when the server cannot be reached and local files are available.

## Gene sets and supported organisms

Gene set collections are provided by the companion **[geneset](https://github.com/GangLiLab/geneset)** package:

| Collection | Retrieval function |
| --- | --- |
| Gene Ontology (BP, CC, MF) | `geneset::getGO()` |
| KEGG | `geneset::getKEGG()` |
| WikiPathways | `geneset::getWiki()` |
| MSigDB | `geneset::getMsigdb()` |
| Enrichr libraries | `geneset::getEnrichrdb()` |
| Reactome | `geneset::getReactome()` |
| MeSH | `geneset::getMesh()` |
| Human disease collections | `geneset::getHgDisease()` |

Organism coverage differs between annotation sources and gene set collections. Check the bundled annotation organism table and the relevant `geneset` function's help for accepted organisms, rather than assuming every collection covers the same species:

```r
data("ensOrg_name", package = "genekitr")
head(ensOrg_name)

# Supported human microarray platforms for transProbe().
data("hsapiens_probe_platform", package = "genekitr")
head(hsapiens_probe_platform)
```

## Visualization and external results

- `plotEnrich()` offers 13 plot types: `bar`, `wego`, `dot`, `bubble`, `lollipop`, `geneheat`, `genechord`, `network`, `gomap`, `goheat`, `gotangram`, `wordcloud`, and `upset`.
- `plotGSEA()` supports `volcano`, `classic`, `fgsea`, `ridge`, and `bar`.
- `plotVenn()` draws Venn diagrams; set `use_venn = FALSE` for an UpSet plot when comparing more gene sets.
- `plotVolcano()` visualizes differential expression results. Use `plot_theme()` or plotting function arguments to adjust text, colors, borders, and legends.

For external enrichment results, use `importCP()`, `importPanther()`, or `importShinygo()` before plotting. Select a statistic column present in your table: PANTHER and ShinyGO imports use `qvalue`, so set `stats_metric = "qvalue"`. For other data frames, `as.enrichdat()` can normalize column names; the selected plot still needs its required term, count, statistic, and gene columns. See the [documentation](https://www.genekitr.fun/) and `?plotEnrich` for plot-specific requirements.

For `genGSEA()` results, `plotEnrich(gsea$gsea_df, plot_type = "dot")` also works; it uses the leading-edge gene counts for the term metric.

## Documentation and support

- [User guide and worked examples](https://www.genekitr.fun/)
- R help: `help(package = "genekitr")`, `?genORA`, `?genGSEA`, and `?plotEnrich`
- [Bug reports and feature requests](https://github.com/GangLiLab/genekitr/issues): include a minimal reproducible example and `sessionInfo()`.

## Why genekitr?

<details>
<summary>The story of Mr. Doodle: seven everyday gene analysis problems</summary>

> **Genes, the essence of life's art,**
> **Omics research's fundamental part,**
>
> **Like cells in our physical frame,**
> **Their study reveals life's vibrant flame.**

Let me tell you a story about Mr. Doodle, a computational biology student working with his PI.

### Scene 1: repetitive work

One day, PI gave him 30 genes to check for their locations and exact names, preferably with sequences.

Mr. Doodle searched for each gene on NCBI, copying and pasting the information into an Excel sheet. He sent the file to PI an hour later, and received praise for his work. But just when he thought he was done, PI gave him another 50 genes to check!

Despite feeling a little overwhelmed, Mr. Doodle repeated the same process with determination, determined to complete the task to the best of his abilities.

##### <b style='color:#486CBE'>Doodle wondered how to avoid repeatedly searching for the same information.</b>

### Scene 2: confusing gene names

Once upon a time, PI gave Mr. Doodle a DEG matrix and a target gene list file. The task was to determine if the target gene was up-regulated after treatment.

Mr. Doodle searched the matrix but couldn't find the PDL1 gene, even though it was in the gene list. He asked PI about it, and PI explained that the gene was listed as CD274, the official gene symbol; PD-L1 is a commonly used alias.

This left Mr. Doodle feeling a little confused. He wondered how to distinguish between real gene names and aliases.

<b style='color:#486CBE'> **Doodle wondered  how to differentiate between a real gene name and an alias.**</b>

### Scene 3: outdated database

Mr. Doodle was analyzing KEGG pathways for the up-regulated genes in the last DEG matrix. His analysis workflow required Entrez IDs, but the genes were listed by their symbols.

Mr. Doodle needed to convert the gene symbols to Entrez IDs, but he found that some symbols did not match the corresponding Entrez IDs. However, he discovered that NCBI had the correct IDs.

Mr. Doodle realized that he was using an older release of the org.Hs.eg.db annotation package. After updating his annotations, he could resolve the missing mappings and continue his KEGG pathway analysis.

<b style='color:#486CBE'>**Doodle wondered if there was a method to help him obtain updated results automatically, instead of having to check them manually every time.**</b>

### Scene 4: incompatible formats

PI did some fancy enrichment analysis all by himself on a website called [GeneOntology](https://geneontology.org/).  He then asked Mr.Doodle to help him make a pretty picture of the results. "Can you make a bubble plot for me and show the FoldEnrichment on the x-axis?" he asked with a smile. Doodle tried to use a fancy R package called clusterProfiler, but it wouldn't work with the data. So, he bravely coded it himself using ggplot2.

<b style='color:#486CBE'>**Doodle wondered why there isn't a tool that supports easy data frames.**</b>

### Scene 5: annoying plot theme

Doodle finally finished making the bubble plot and sent it to PI. After 15 minutes, PI sent him a message with a smile: "the text is too small, and can you make the background white with a border size of 4 points?" Doodle tweaked the ggplot theme and made the changes in 10 minutes. But, after a little while, PI sent another message saying, "The border is too thick in the second version. Can you please redo it?"

<b style='color:#486CBE'>**Doodle wondered if there was a function that could help him process the plot theme instead of having to modify the current code repeatedly.**</b>

### Scene 6: limited plot types

PI gave Doodle the GO enrichment analysis result and asked him to think of a creative way to display it. Doodle found that each tool had its specific plot. For example, [WEGO](https://wego.genomics.cn/) could compare BP, CC, and MF terms; [GOplot](https://wencke.github.io/) had a chord plot to show the relationship between genes and GO terms; and [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) supported enriched map and network, which could explore the relationship among enriched terms. However, there was a big problem - the input data for each tool was not compatible, making it inconvenient to plot WEGO plots using clusterProfiler objects.

<b style='color:#486CBE'>**Doodle wondered if there was a method that could produce beautiful plots from different tools using a universal data format.**</b>

### Scene 7: chaotic export files

Doodle finished conducting differential expression analysis and GO/KEGG enrichment analysis. PI asked him to send over all the result files. Doodle saved the results into three separate excel files, naming them "DEG_data.xlsx," "GO_enrich.xlsx," and "KEGG_enrich.xlsx." He then compressed the three files into a zipped folder, naming it after the date, and sent it to PI. After a while, PI asked him if he could put all three results into a single excel file.

<b style='color:#486CBE'>**Doodle wondered if there is a way to save all data into a single file without having to perform many manual operations.**</b>

> **If you have encountered similar problems like Mr. Doodle, give <u>genekitr</u> a try!**

</details>

## Author

[Yunze Liu](https://www.jieandze1314.com/)

[![Follow on WeChat](https://img.shields.io/badge/follow%20me%20on-WeChat-orange.svg)](https://genekitr.online/img/bioinfoplanet.png)

## Citation

If you use genekitr in your research, please cite:

Liu, Y., Li, G. Empowering biologists to decode omics data: the Genekitr R package and web server. *BMC Bioinformatics* **24**, 214 (2023). https://doi.org/10.1186/s12859-023-05342-9

## Contributing

Contributions and suggestions are welcome. For code changes:

1. Fork the repository and clone your fork.
2. Open `genekitr.Rproj` in RStudio and edit the source files in `R/`.
3. If you change public function documentation, regenerate the help files with `devtools::document()`.
4. Run `devtools::check()` and address any new errors, warnings, or notes.
5. Open a pull request describing the change and how you validated it.

Install `devtools` with `install.packages("devtools")` if needed. Genekitr is licensed under [GPL-3](LICENSE.md).

<details>
<summary>Repository activity</summary>

![Repository activity](https://repobeats.axiom.co/api/embed/e42ba06d30de893670c70324f19398ef0a7c26fa.svg "Repobeats analytics")

</details>
