library(biomaRt)
ensembl <- useMart("ensembl",host = "asia.ensembl.org")
listDatasets(ensembl) %>% dplyr::pull(1) %>% stringr::str_remove_all('_gene_ensembl')

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "asia.ensembl.org")
feature_page <- listAttributes(ensembl) %>%
filter(page == 'feature_page')
sequences <- listAttributes(ensembl) %>%
  filter(page == 'sequences')

View(attributes)

mm_id =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Bloc1s1
Tp53
Gcn5l1')
mm_id=stringr::str_split(mm_id,"\n")[[1]]
attr <- c("external_gene_name",'external_synonym','uniprot_gn_symbol')

biomart_alias = getBM( values = unique(all$symbol),
                       attributes = c("external_gene_name",'external_synonym','uniprot_gn_symbol'),
                       filters = "external_gene_name",
                       mart = useMart("ensembl",
                                      dataset = "mmusculus_gene_ensembl",
                                      host = "asia.ensembl.org")) %>%
  split(., .$external_gene_name) %>%
  lapply(., function(x) {
    paste0(x[, 2],collapse = "; ")
  })  %>%
  do.call(rbind,.) %>%
  as.data.frame() %>%
  dplyr::mutate(symbol = rownames(.)) %>%
  # dplyr::mutate(uniq_symbol = paste(V1,V2,V3)) %>%
  dplyr::rename(uniq_symbol = V1) %>%
  dplyr::select(symbol, uniq_symbol) %>%
  within(., uniq_symbol <- sapply(stringr::str_split(uniq_symbol,'; '), unique)) %>%
  dplyr::arrange(match(symbol, mm_id))


# get all symbols
library(biomaRt)
organism = 'hsapiens'
bmt = getBM( attributes = c("ensembl_gene_id",'chromosome_name','start_position','end_position','strand',
                      'percentage_gene_gc_content','gene_biotype','transcript_count'),
       mart = useMart("ensembl",
                      dataset = paste0(organism,"_gene_ensembl"),
                      host = "asia.ensembl.org")) %>%
  data.table::setnames(., old =colnames(.),
           new = c('ensembl','chr','start','end','strand','gc_content','gene_biotype','transcript_count')) %>%
  dplyr::mutate(width = (end - start + 1))

colnames(bmt)

