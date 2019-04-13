mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
listMarts(mart = mart)
datasets <- listDatasets(mart = mart) %>% as_tibble()


annot <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene"), 
                        filters = "entrezgene", values = egid, mart = mart)


attributes <- listAttributes(mart = mart)
