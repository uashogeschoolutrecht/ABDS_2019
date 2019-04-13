## get the individual sample files for the "airway2" package

#BiocManager::install("GEOquery", dependencies = TRUE)

library(GEOquery)

## defines the root of the project for later use
require("rprojroot") || utils::install.packages("rprojroot")
library(rprojroot)
root <- find_root_file(criterion = is_rstudio_project)

files_path <- file.path(root, "temp", "airway2_GSE_files")
dir.create(files_path)

getGEOSuppFiles("GS0E52778", baseDir = files_path)
getGEO("GSE52778", destdir = files_path)

gunzip(filename = file.path(files_path, "GSE52778_series_matrix.txt.gz"))
gunzip(filename = file.path(
  files_path, "GSE52778", "GSE52778_All_Sample_FPKM_Matrix.txt.gz")
  )

gunzip(filename = file.path(files_path, ""))



