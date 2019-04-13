all_rmd <- list.files(path = ".",
                      recursive = TRUE,
                      full.names = TRUE,
                      pattern = ".Rmd")


## dummy:
path = all_rmd[44]

get_rmd_dependencies <- function(path){

  reqlibs = sub(".*library\\(\"(.*?)\"\\).*","\\1",
              grep("library\\(", 
                   readLines(path),
                   value=TRUE))
  
  reqlibs <- str_replace_all(string = reqlibs,
                             pattern = "library\\(", 
                             replacement = "")
  reqlibs <- str_replace_all(string = reqlibs,
                             pattern = "\\)", 
                             replacement = "")
  
  reqlibs
  
 
}

get_rmd_dependencies(path = path)
packages <- map(as.list(all_rmd), get_rmd_dependencies) %>% 
  unlist() %>% 
  unique() %>%
  trimws()

packages


