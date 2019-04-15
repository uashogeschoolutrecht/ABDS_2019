options(stringsAsFactors=FALSE)


headers <- read.table(file="SGD_features.headers")
headers <- headers[, 1]
## headers <- headers$V1

SGD <- read.table(file="SGD_features.txt", sep="\t", 
                  header=FALSE, quote="", col.names=headers)
## row.names=1: sgdid
