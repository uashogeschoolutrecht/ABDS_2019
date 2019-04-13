### function to extract .mtx files (Matrix M)

if (!require("rprojroot")) install.packages("rprojroot")
library("rprojroot")
root <- find_root_file(criterion = is_rstudio_project)
root

## packages
library(tidyverse)
library(GEOquery)
library(edgeR)

## example data
## from GEO accession GSE106510
## for details: https://www.ncbi.nlm.nih.gov/pubmed/30392957
# supp_data <- GEOquery::getGEOSuppFiles("GSE106510")
# path_supp <- file.path(root, 
#                       "GSE106510")

## The raw mtx files are in the archive GSE106510_RAW.tar
## untar this archive
# untar(tarfile = file.path(path_supp, "GSE106510_RAW.tar"), 
  #    exdir = path_supp)

## calculated file:
path_supp <- file.path(root, 
                       "GSE106510")

file_list <- list.files(path = path_supp,
                        full.names = TRUE)


mtx_list <- list.files(path = path_supp,
                       full.names = TRUE, 
                       pattern = "*.mtx")
## gunzip all mtx files (multicore use)
library(parallel)
library(doParallel)
detectCores()
registerDoParallel()
makeCluster(7)

# foreach(i=file_list) %dopar% GEOquery::gunzip(i)

## new list (unzipped)
barcodes_list <- list.files(path = path_supp,
                       full.names = TRUE, 
                       pattern = "*_barcodes.tsv")


genes_list <- list.files(path = path_supp,
                           full.names = TRUE, 
                           pattern = "*_genes.tsv")


matrix_list <- list.files(path = path_supp,
                         full.names = TRUE, 
                         pattern = "*_matrix.mtx")

## table with files that 'belong' together
matrix_df <- matrix_list %>% as_tibble() %>%
  arrange()
genes_df <- genes_list %>% as_tibble() %>%
  arrange()
barcodes_df <- barcodes_list %>% as_tibble() %>%
  arrange()

files_data_raw_df <- dplyr::bind_cols(matrix_df,
                                      genes_df,
                                      barcodes_df) %>%
  as_tibble() %>% 
#  print()

  
  dplyr::mutate(basename = basename(value),
                basename1 = basename(value1),
                basename2 = basename(value2)) %>%
  print()


df = files_data_raw_df

    mtx <- df$value %>% as.list
    genes <- df$value1 %>% as.list()
    barcodes <- df$value2 %>% as.list()

    x <- list()
      
    for(i in seq_along(mtx)){
    x[[i]] <- edgeR::read10X(mtx = mtx[[i]], 
               genes = genes[[i]], 
               barcodes = barcodes[[i]])

  }

x[[1]] %>% class
x[[2]]

## Continue workflow with edgeR DGEList >>>
library(DESeq2)
library(SingleCellExperiment)
#library(DESeq2)
#library(edgeR)
library(DEFormats)
library(edgeR)
sce <- SingleCellExperiment(
  assays = x[[1]]$counts,
  rowData = x[[1]]$genes,
  colData = names(x[[1]]$counts %>% as.data.frame())

)
sce
normcounts(sce) <- log2(counts(sce) + 1)
