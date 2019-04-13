## function to convert XStrings to dataframe 
#' Convert an XString(Set) to a dataframe  
#' @param dss An object of formal calss XString or XStringSet, can be DANAString, DNAStringSet, AAString, AAStringSEt

xss2df <- function(dss){
  
  tbl <- tibble(width=width(dss), 
                seq=as.character(dss), 
                names=names(dss))
  return(tbl)
}