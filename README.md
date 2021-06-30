# Annotate genes 

> This tool annotates genes with alias, symbol, full name, function also related papers.
> 目的就是：将与基因ID相关的内容整合在一起，并且能直接导出为一个excel的多个sheets（后期考虑加入网页版）

### Features

- genecards虽然全，但是搜索数量有限制
- 我们整合了基因信息、相关的文献信息、基因相关的GO、KEGG（如果能提供对应的logFC，还能做GSEA；如果能提供表达矩阵，还能做GSVA）



### Let's do it!

#### example gene id

```R
id =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Bloc1s1') 
id=str_split(id,"\n")[[1]]
```



#### Method1: gene alias, full name

```R
test1 = genInfo(id, org = 'mm')
```

rownames of `test1` are entrez ID

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-29-081721.png)

#### Method2: search pubmed 

```R
test2=genPubmed(id, keywords = 'stem cell', field = 'tiab')
# Search example: Ticam2 [TIAB] AND stem cell [TIAB] 

# or use much specific keyword
genPubmed(id, keywords = 'stem cell AND epithelial', field = 'tiab')
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-29-081925.png)





### Exporting result is very easy!

- support pipe ` %>% ` 

```R
library(openxlsx)
wb <- createWorkbook()
wb <- expo_sheet(wb, sheet_dat = test1, sheet_name = 'genInfo') %>% 
  expo_sheet(., sheet_dat = test2, sheet_name = 'genPub') %>% 
saveWorkbook(wb, "~/Downloads/test.xlsx", overwrite = T)
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-29-082350.png)



### Plans

- plot can also save into excel

