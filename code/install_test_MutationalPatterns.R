################################################################
## install MutationalPatterns (for Bioconductor Version 3.8) ###
################################################################
library(tidyverse)
## confirm current version
BiocManager::version()
## download dependencies
dir.create("tmp")
## Bioc 3.8 version of rngtools
download.file("https://cran.r-project.org/src/contrib/Archive/rngtools/rngtools_1.2.4.tar.gz",
              destfile = file.path("tmp", "rngtools_1.2.4.tar.gz"))
## install rngtools from local source
install.packages(repos = NULL, pkgs = file.path("tmp", "rngtools_1.2.4.tar.gz"))
## Bioc 3.8 version of NMF
download.file("https://cran.r-project.org/src/contrib/Archive/NMF/NMF_0.20.6.tar.gz",
              destfile = file.path("tmp", "NMF_0.20.6.tar.gz"))
## install NMF from local source
install.packages(repos = NULL, pkgs = file.path("tmp", "NMF_0.20.6.tar.gz"))
## install Bioc 3.8 version of MutationalPatterns package
BiocManager::install("MutationalPatterns", version = "3.8")
library(MutationalPatterns)

#browseVignettes("MutationalPatterns")

