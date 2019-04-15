# This directory contains the data used throughout the introduction to R course
# and also the scripts that were used to generate them.
# The data contains the expression profiles of GSTFs and chromatin factors
# selected from all deletion mutants profiled for the deleteome data.

# The data files (but not the scripts) are made visible to the students by
# sharing them in the for_students/data folder. Note that any change made
# in masters_intro_R/data/abc.ext will be directly visible in 
# for_students/data/abc.ext.
# However, addition or removal of a file masters_intro_R/data/def.ext is
# not directly visible in the for_students/data folder, you have to 
# explicitly add them (using shift-Z) or delete them (in both folders)


# Retrieve the data from OffBase
source /db/offbase/prod/env.sh
/db/offbase/prod/base/bin/export_dba.pl -hs 1897 -t reporter -ra 'systematicName' -ca 'sampleName' -df 'M,p_value' -o ./deleteome_dba.txt -v

# convert the data to PDL
source /data/www/integromics/env.sh
/data/www/integromics/Methods/perl/offbase-tab_to_pdl.pl -i ./deleteome_dba.txt -nar 9 -rsc -si -ex_row /hpc/dbg_gen/patrick/integromics/data/exp/wt_variable_genes.txt -rc 1.7 -pc 0.05 -c_min 4 -r_min 2

##############################
# within R
##############################

library(integromicsCorePDL)
species <- "cerevisiae"

# Read in the list of functional categories
func.cat <- read.delim("functional_categories.txt", row.names=1, stringsAsFactors=FALSE)
fav.genes <- rownames(func.cat)[func.cat["functional.category"]=="chromatin factor" | func.cat["functional.category"]=="gene-specific transcription factor"]

# read in the full deleteome dataset
deleteome <- new("iCorePDL", file="./deleteome_dba.txt", species="cerevisiae", verbose=T)
stdMean <- rep(TRUE,ncol(dataMatrix(deleteome)))
stdMean[colAnnot(deleteome)[,"dataType"]=="p_value"] <- FALSE
deleteome <- averageDups(deleteome, stdMean=stdMean)

# subselect only for GSTFs and chromatin factors
gstf.cf.pdl <- keepIdentifiers(deleteome, idList=fav.genes, hdrName="systematic_name", dim="col")

# write out tab-delimited files for the M and p-values
# size of each file: 6163 rows x 338 columns
for ( type in c("M", "p_value") ) {
  data <- dataMatrix(gstf.cf.pdl)[ ,colAnnot(gstf.cf.pdl)[,"dataType"]==type]
  rownames(data) <- rowAnnot(gstf.cf.pdl)[ ,"systematic_name"]
  colnames(data) <- colAnnot(gstf.cf.pdl)[colAnnot(gstf.cf.pdl)[ ,"dataType"]==type , "systematicName"]
  colnames(data) <- gsub( " vs. wt.*", "", colnames(data))
  colnames(data) <- gsub( "-matA", "", colnames(data))
  colnames(data) <- gsub( "-", "_", colnames(data))

  write.table(data, file=paste("./exp_profiles_gstf_cf.", type, ".txt", sep=""), sep="\t", quote=FALSE, col.names=NA)
}


simple-IO is the first 100 lines and first 10 columns of exp_profiles_gstf_cf.M.txt


