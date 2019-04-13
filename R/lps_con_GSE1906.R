# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Jan 11 09:22:36 EST 2018

# R version 3.4.3 (2017-11-30)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] BiocInstaller_1.28.0 limma_3.34.9         GEOquery_2.46.15     Biobase_2.38.0      
# [5] BiocGenerics_0.24.0 

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library(DESeq2)


## select il's and ebi3 for heatmap below

# load series and platform data from GEO

gset <- getGEO("GSE48427", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL81", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "014253014253"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G5-G0, G1-G0, G2-G1, G3-G2, G4-G3, G5-G4, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE1906", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL81", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "014253014253"
sml <- c(gsms)
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("0_min_lps","15_min_lps","60_min_lps","240_min_lps","30_min_lps","120_min_lps")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf","#f2cb98","#dff4e4","#f4dff4","#dcdaa5", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE1906", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")


##################################################################################
# create annotated heatmap
##################################################################################


  library(RNAseq123)
  library(lattice)
  library(Biobase)
  #install.packages("pheatmap")
  library(pheatmap)
  ##  install.packages("xlsx")
  #library(xlsx)
  library(gdata)
  
  
  require("rprojroot") || utils::install.packages("rprojroot")
  library(rprojroot)
  root <- find_root_file(criterion = is_rstudio_project)
  
  
  richSetSummarizedExp <- 
    SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(richSet)
  
  richSetSummarizedExp
  
  col_data <- colData(richSetSummarizedExp)
  count_data <- assay(richSetSummarizedExp)
  col_data
  head(count_data)
  str(count_data)
  
  # writes summarizedExperiment to disk as Rda object
  readr::write_rds(richSetSummarizedExp, 
                   path = "./data/summarized_experiment_jp_06022017.Rda")
  
  return(richSetSummarizedExp)
  
}

summarized_experiment <- generate_summarized_experiment()
summarized_experiment@assays[[1]]

###################################################################
## FUNCTION TO GET AN ANNOTATED HEATMAP, YOU CAN CHOOSE THE ATTRIBUTE
## FROM ONE OF THE FOLLOWING VALUES:
# 'go_id',
# 'ensembl_gene_id',                                                             'ensembl_transcript_id',
# 'hgnc_symbol',
# 'hgnc_id', 
# 'external_gene_name',
# 'ensembl_peptide_id',
# 'description'
###################################################################

## setting argumnents for the function outside to check
#number_top_genes = 30000
#fold_change = 1.2
#summarized_experiment = summarized_experiment
#go_term = "GO:0005576"

create_annotated_heatmap <- function(summarized_experiment, 
                                     number_top_genes,
                                     go_term,
                                     fold_change){
  # packages:
  library(biomaRt)
  library(tidyverse)
  library(stringr)
  
  require("rprojroot") || utils::install.packages("rprojroot")
  library(rprojroot)
  root <- find_root_file(criterion = is_rstudio_project)
  
  dds <- DESeqDataSet(summarized_experiment, design = ~ treatment)
  dds_reduced <- dds[rowSums(counts(dds)) > fold_change, ]
  rld <- rlog(dds_reduced, blind=FALSE)
  
  topVarGenes <- head(order(rowVars(assay(rld)), 
                            decreasing=TRUE), 
                      number_top_genes)
  
  counts_matrix <- assay(rld)[ topVarGenes, ]
  
  mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #listDatasets(mart)
  #attributes <- as_tibble(listAttributes(mart))
  
  # help on biomaRt
  #?useMart
  #?listMarts()
  #listMarts()
  #ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  
  annotations <- 
    getBM(attributes=c(
      'go_id',
      'entrezgene',
      'ensembl_gene_id',
      'ensembl_transcript_id',
      'hgnc_id', 
      'external_gene_name',
      'ensembl_peptide_id',
      'description'),
      filters = 'ensembl_gene_id', 
      values = row.names(counts_matrix), 
      mart = mart)
  
  # convert counts matrix to dataframe
  df <- as_tibble(counts_matrix)
  df$gene_id <- row.names(counts_matrix)
  
  attribute_filtered <- annotations %>%
    filter(go_id == go_term)
  
  library(dplyr)
  filtered_df <- df %>%
    filter(gene_id %in% attribute_filtered$ensembl_gene_id)
  
  class(filtered_df$TRA1) 
  
  ## convert to matrix  
  
  ?as.matrix
  
  str(filtered_df)
  
  go_matrix <- as.matrix(filtered_df
                         [,-which(names(filtered_df) == "gene_id")]) 
  
  row.names(go_matrix) <- c(filtered_df$gene_id)
  head(go_matrix)
  
  
  getMatrixWithSelectedIds <- function(df, type, keys){
    require("AnnotationDbi")
    require("org.Mm.eg.db")
    
    geneSymbols <- mapIds(org.Mm.eg.db, 
                          keys=rownames(df), 
                          column=type, 
                          keytype=keys, 
                          multiVals="first")
    
    # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
    inds <- which(!is.na(geneSymbols))
    found_genes <- geneSymbols[inds]
    
    # subset your data frame based on the found_genes
    df2 <- df[names(found_genes), ]
    rownames(df2) <- found_genes
    return(df2)
  }
  
  go_matrix <- getMatrixWithSelectedIds(
    go_matrix, type="SYMBOL", keys="ENSEMBL")
  
  # entrez_id <- getMatrixWithSelectedIds(
  #   go_matrix, type="ENTREZID", keys="ENSEMBL")
  
  
  go_matrix <- go_matrix - rowMeans(go_matrix)
  
  
  
  filename <- str_replace_all(go_term, ":", "_")
  
  png(filename = file.path(root, "images", 
                           paste(filename, "png", sep = ".")), 
      res = 600, width = 15, height = 20, units = "cm")
  pheatmap::pheatmap(go_matrix, 
                     show_rownames = TRUE, 
                     main = go_term)
  dev.off()
  
  heatmap <- pheatmap::pheatmap(go_matrix, 
                                show_rownames = TRUE, 
                                main = go_term)
  
  return(list(go_matrix, annotations)) #, entrez_id))
}

## a heatmap for GO:0005576: extracellular secretion domain
heatmap <- create_annotated_heatmap(
  number_top_genes = 1000,
  fold_change = 1.5,
  summarized_experiment = summarized_experiment,
  go_term = "GO:0005576"
)

## a heatmap for GO:0002534
# https://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0002534 
# Name	cytokine production involved in inflammatory response
# Ontology	Biological Process
# Definition	The synthesis or release of a cytokine following a inflammatory stimulus as part of an inflammatory response, resulting in an increase in its intracellular or extracellular levels.
# GONUTS	GO:0002534 Wiki Page

heatmap <- create_annotated_heatmap(
  number_top_genes = 30000,
  fold_change = 1.5,
  summarized_experiment = summarized_experiment,
  go_term = "GO:0002534"
)


#Accession
#GO:0005125
#Name
#cytokine activity
#Ontology
#molecular_function
#Synonyms
#autocrine activity, paracrine activity
#Alternate IDs
#None
#Definition
#Functions to control the survival, growth, differentiation and effector function of tissues and cells. Source: ISBN:0198599471
#Comment
#Also consider annotating to 'receptor agonist activity ; 
#GO:0048018'.
#History
#See term history for GO:0005125 at QuickGO
#Subset
#goslim_chembl
heatmap <- create_annotated_heatmap(
  number_top_genes = 30000,
  fold_change = 1.5,
  summarized_experiment = summarized_experiment,
  go_term = "GO:0005125"
)

#GO:0004896
#Name	cytokine receptor activity
#Ontology	Molecular Function
#Definition	Combining with a cytokine and transmitting the signal from one side of the membrane to the other to initiate a change in cell activity.
#Secondary IDs	GO:0004907
#GONUTS	GO:0004896 Wiki Page
heatmap <- create_annotated_heatmap(
  number_top_genes = 30000,
  fold_change = 1.5,
  summarized_experiment = summarized_experiment,
  go_term = "GO:0004896"
)


## cell adhesion: GO:0007155
heatmap <- create_annotated_heatmap(
  number_top_genes = 5000,
  fold_change = 1.5,
  summarized_experiment = summarized_experiment,
  go_term = "GO:0007155"
)


# positive meoloid dendritic cell activation GO:0030886 
heatmap <- create_annotated_heatmap(
  number_top_genes = 30000,
  fold_change = 1.5,
  summarized_experiment = summarized_experiment,
  go_term = "GO:0030886"
)


