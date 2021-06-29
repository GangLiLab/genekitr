# Annotate genes 

> This tool annotates genes with alias, symbol, full name, function also related papers.

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
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-29-081925.png)





### Exporting result is very easy!

- support pipe ` %>% ` 

```R
library(openxlsx)
wb <- createWorkbook()
expo_sheet(wb, sheet_name = 'genInfo',sheet_dat = genInfo(id, org = 'mm')) %>% 
  expo_sheet(., sheet_name = 'genPub',sheet_dat = pub)
saveWorkbook(wb, "~/Downloads/test.xlsx", overwrite = T)
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-29-081258.png)



### Plans

- plot can also save into excel

