## Resubmission: genekitr 1.3.0

This release fixes a hard error with clusterProfiler >= 4.20 (`GSEA()` no longer
accepts `by =`), removes deprecated ggplot2/igraph/tidyselect usages, stops
writing to the global environment, and drops unused dependencies.
See ChangeLog for details.

## Test environments
* local macOS (x86_64), R 4.6.1, clusterProfiler 4.20.2, ggplot2 4.0.3, dplyr 1.2.1

## R CMD check --as-cran results

0 errors | 0 warnings | 0 notes

* Suggested packages not installed locally (fgsea, GOplot, pheatmap, tm,
  treemap, rrvgo, wordcloud) are only used conditionally.
* `\donttest{}` examples need internet access to download annotation data
  (cached under `tools::R_user_dir("genekitr", "data")`, hence R >= 4.0.0).
