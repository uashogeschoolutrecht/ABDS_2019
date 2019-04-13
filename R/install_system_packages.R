## system- wide install of packages
## start R from terminal by `sudo R`

library(tidyverse)

install_packages_on_system <- function(package_name, lib){

library <- .libPaths()

  pacman::p_load("BiocManager")

  require("pacman") || utils::install.packages("pacman")
  pacman::p_load(package_name)

  
  }

package_df <- readr::read_lines("./data/package_list.txt") %>% 
  tibble::as_tibble() 

ind <- duplicated(package_df)
package_df <- package_df[ind, ]

package_df <- package_df %>%
  mutate(cran = pacman::p_iscran(value))

cran_pkgs <- package_df %>%
  filter(cran == TRUE)

purrr::map(as.list(cran_pkgs$value), install.packages)

non_cran <- package_df %>%
  filter((cran != TRUE))

map(as.list(non_cran$value),BiocManager::install)
##pkgs <- package_df$value %>% as.character() %>% as.list()

library(DESeq2)

browseVignettes("DESeq2")




browseVignettes("RforProteomics")
