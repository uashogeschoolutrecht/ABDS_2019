## download from url
library(tidyverse)
library(readr)

url <- c("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76073/suppl/GSE76073_Raw_counts.txt.gz")
download.file(url = url, destfile = "./data-raw/GSE76073_Raw_counts.txt.gz")

?untar
untar(tarfile = "./data-raw/GSE76073_Raw_counts.txt.gz", exdir = "./data/GSE76073_Raw_counts.txt")
data_raw <- read_lines(file = "./data-raw/GSE76073_Raw_counts.txt.gz", ) %>% print()

