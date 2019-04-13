## retrieves the Rmd files contining the answers to the exercises 
## from the "/answers" folder

retrieve_rmd_from_answers <- function(){
  
  exercises_back <- read_csv(file = file.path(own_dir, "rmd_copied_to.csv"))
  
  map(exercises_back$value, file.move, destinations = own_dir)
  
}

retrieve_rmd_from_answers()
