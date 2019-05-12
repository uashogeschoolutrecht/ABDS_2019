library(pacman)
p_load(GenomeGraphs)
library(GenomeGraphs)
library(biomaRt)



#biomart <- "ensembl"
#id = "ENSG00000140564"
#type = "ensembl_gene_id"
#plot_title = "ENSG00000140564"
#mart <- useMart(biomart = biomart, dataset = "hsapiens_gene_ensembl")
#mart
#marts <- listMarts(mart = mart)

#datasets <- listDatasets(mart = mart) %>% as_tibble()

#dataset <- "hsapiens_gene_ensembl"



#structure = c("gene", "transcript", "snps") 



make_gene_graph <- function(biomart, 
                            dataset, 
                            type, 
                            id){
 

  mart <- useMart(biomart = biomart, 
                  dataset = dataset)

  
  
#listDatasets(mart)
  gene <- makeGene(id = id, 
                 type = type, 
                 biomart = mart)
  
  transcript <- makeTranscript(id = id,
                               type = type, 
                               biomart = mart)
    
 #   } else {
  
#  if(structure == c("transcript", "snps") | 
 #   structure == c("snps", "transcripts")) { 
  
    
#    } 

 #   }  

  my.gene<-getBM(c('chromosome_name',
                 'start_position',
                 'end_position',
                 'external_gene_name'),
               filters = type,
               value = id, mart = mart)

  gdPlot(list(gene, transcript))  
  gdPlot(list(makeTitle(as.character(my.gene$external_gene_name)), 
            makeIdeogram(chromosome = my.gene$chromosome_name), 
            gene, 
            transcript, 
            makeGenomeAxis()))


# get gene start and end pos
# get snp tranck on the plot
#snpmart <- useMart("snp",dataset="hsapiens_snp")

  
  
    snpmart <- useMart("ENSEMBL_MART_SNP",
                   dataset = "hsapiens_snp") #updating database name

  
#snps<-getBM(c('refsnp_id','allele','chrom_start','clinical_significance'),filters = c('chr_name','chrom_start','chrom_end'), values = list(my.gene$chromosome_name,my.gene$start_position,my.gene$end_position), mart = snpmart)

  snps<-getBM(c('refsnp_id','allele','chrom_start'),
            filters = c('chr_name','start','end'), 
            values = list(my.gene$chromosome_name,
                          my.gene$start_position,
                          my.gene$end_position), 
            mart = snpmart)

  snpannot<-makeAnnotationTrack(start=snps$chrom_start,
                              end=snps$chrom_start)


  gdPlot(list(makeTitle
                      (as.character(my.gene$external_gene_name)), 
            makeIdeogram(chromosome = my.gene$chromosome_name), 
            gene, 
            transcript, 
            snpannot,
            makeGenomeAxis()))


  
  
  
  # create custom annotation
#  customann<-makeAnnotationTrack(start=c(90876557,
#                                         90876757,
#                                         90878557),
#                                 end=c(90876597,
#                                       90876897,
#                                       90878697),
#                                 feature=c('bind','bind','del'),
#                                 dp = DisplayPars(bind = 'blue',
#                                                  del='red'))


#  gdPlot(list(makeTitle("Human furin gene"), 
#              makeIdeogram(chromosome = 15), 
#              gene, 
#              transcript, 
#              snpannot,
#              customann,
#              makeGenomeAxis()))

#  return(plot)
}


# make_gene_graph(biomart = biomart,
#                 dataset = dataset,
#                 type = type,
#                 id = id)
# #traceback()
# 
# listDatasets(mart = biomart)
# 
# biomart <- "ensembl"
# id = "ENSG00000140564"
# type = "ensembl_gene_id"
# plot_title = "ENSG00000140564"
# 
# structure = c("gene", "transcript", "snps") 
# 
# 
# biomart <- "ensembl"
# id = "ENSG00000120075"
# type = "ensembl_gene_id"
# plot_title = "ENSG00000140564"
# 
