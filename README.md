# Annotate genes 

> This tool annotates genes with alias, symbol, full name and function.

### Method #1

use `org.db`  to annotate genes, including human (hs) and mouse (mm)

```R
library(stringr)
library(dplyr)
# you could paste any ids as you wish
id =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Bloc1s1') 
id=str_split(id,"\n")[[1]]

> id
[1] "Ticam2"     "Arhgap33os" "Insl3"      "Myo15"      "Gal3st2b"   "Bloc1s1" 

genOrgInfo(id, org = 'mm', html_result = T, destdir = '~/Downloads/')
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-15-120520.png)

```R
genOrgInfo(id, org = 'mm', html_result = F, destdir = '~/Downloads/')
```

![](https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2021-06-15-120658.png)

