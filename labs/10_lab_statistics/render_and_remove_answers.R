## renders all Rmd in the current folder
library(tidyverse)
library(filesstrings)

own_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
rmd_files <- list.files(path = own_dir, pattern = "\\.Rmd", 
                        full.names = TRUE) %>%
  as.list()

rmd_files

for(i in seq_along(rmd_files)){
purrr::map(rmd_files[[i]], rmarkdown::render)
}

## remove the Rmd files (exercises only) that contain the 
## answers to the exercises and puts them in the "/answers folder
## put the /answers folder in gitignore 

## TODO: write a function that puts them back in a lab, on the basis 
## of a lab name

#library(tidyverse)
#library(filesstrings)



rmd_files <- list.files(path = own_dir, pattern = "\\.Rmd", 
                        full.names = TRUE)

rmd_files_df <- rmd_files %>%
  enframe(name = NULL)

rmd_files_df <- rmd_files_df %>%
  mutate(file_name = basename(value))

rmd_files_df

ind <- str_detect(string = rmd_files_df$file_name, 
                  pattern = "._exercise_.")

exercises <- rmd_files_df[ind, "value"] %>%
  mutate(file_name = basename(value))

destination <- here::here("answers")

rmd_copied_to <- file.path(destination, exercises$file_name) %>%
  enframe(name = NULL)

## save rmd new locations
write_csv(rmd_copied_to, path = file.path(own_dir, "rmd_copied_to.csv"))


map(exercises, file.move, destinations = destination)

