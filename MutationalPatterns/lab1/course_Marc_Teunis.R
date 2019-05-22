#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(reshape2)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
library("readxl")
library(dplyr)


###output_dir
dirpath <- here::here("MutationalPatterns", "output")
dir.create(dirpath)


###Functions
cos_sim = function(x, y){
  res = x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
  # coerce matrix to numeric
  res = as.numeric(res)
  return(res)
}
colorpalette = c("Signature.1" =  '#8dd3c7',
                 "Signature.2" =  '#ffffb3',
                 "Signature.3" =  '#bebada',
                 "Signature.4" =  '#fb8072',
                 "Signature.5" =  '#80b1d3',
                 "Signature.6" =  '#fdb462',
                 "Signature.7" =  '#b3de69',
                 "Signature.8" =  '#fccde5',
                 "Signature.9" =  '#d9d9d9',
                 "Signature.10" = '#ff1417' ,
                 "Signature.11" = '#ff6611' ,
                 "Signature.12" = '#c4ff00' ,
                 "Signature.13" = '#ff8844' ,
                 "Signature.14" = '#ffee55' ,
                 "Signature.15" = '#ffff99' ,
                 "Signature.16" = '#78FA37' ,
                 "Signature.17" = '#aacc22' ,
                 "Signature.18" = '#bbdd77' ,
                 "Signature.19" = '#c8cf82' ,
                 "Signature.20" = '#92a77e' ,
                 "Signature.21" = '#5599ee' ,
                 "Signature.22" = '#0088cc' ,
                 "Signature.23" = '#226688' ,
                 "Signature.24" = '#175279' ,
                 "Signature.25" = '#557777' ,
                 "Signature.26" = '#ddbb33' ,
                 "Signature.27" = '#d3a76d' ,
                 "Signature.28" = '#a9834b' ,
                 "Signature.29" = '#aa6688',
                 "Signature.30" = '#767676',
                 "Signature.A" = '#458B00' ,
                 "Signature.B" = '#D2691E' ,
                 "Signature.C" = '#6495ED' ,
                 "Signature.D" = '#A2CD5A' ,
                 "Signature.E" = '#CD3333' ,
                 "Signature.F" = '#7AC5CD' ,
                 "Signature.G" = '#009ACD' ,
                 "Signature.H" = '#CD2626' ,
                 "Signature.I" = '#FFB90F' ,
                 "Signature.J" = '#76EEC6' ,
                 "Signature.K" = '#EEB422' ,
                 "Signature.L" = '#97FFFF' ,
                 "Signature.M" = '#E9967A' ,
                 "Signature.N" = '#5F9EA0')




#_________________#
#_________________#

#import SNVs
vcf_organoids <- list.files(here::here("MutationalPatterns", 
                                       "filtered_VAF"), 
                            pattern = ".vcf", 
                            full.names = TRUE)
vcf_files_names <- substr(basename(vcf_organoids), 1, 
                          nchar(basename(vcf_organoids)) - 4) 
vcf_files_names <- sub("_pass_filtered_indel_VAF30_70",
                       "",
                       vcf_files_names)
vcfs_grange <- read_vcfs_as_granges(vcf_organoids, 
                                    vcf_files_names, 
                                    genome = ref_genome)
treatment <- c(rep("Control", 6),
               rep("5-FU", 2))

#total number of SNVs
length(vcfs_grange$"STE072-control-p17_5-FU-2-625-8")
length(vcfs_grange$"STE072-control-p17_5-FU-3-625-7")

#plot type occurences
type_occurrences <- mut_type_occurrences(vcfs_grange, 
                                         ref_genome)
type_occurrences
plot_spectrum(type_occurrences)
plot_spectrum(type_occurrences, 
              CT = TRUE, 
              legend = FALSE)
plot_spectrum(type_occurrences, 
              by = treatment, 
              CT = TRUE, 
              legend = TRUE)

#create 96-mutation matrix
auto <- extractSeqlevelsByGroup(species="Homo_sapiens",
                                style="UCSC",
                                group="auto")

## filter ranges for autosomal 
vcfs <- lapply(vcfs_grange, function(x) keepSeqlevels(
  x, auto, pruning.mode="coarse")
  )

vcfs_mm <- mut_matrix(vcf_list = vcfs, 
                      ref_genome = ref_genome)

colnames(vcfs_mm)
colSums(vcfs_mm)
plot_96_profile(vcfs_mm)

vcfs_mm_treatment <- as.data.frame(vcfs_mm)

vcfs_mm_treatment <- cbind(as.data.frame(
  rowSums(vcfs_mm_treatment[grep("5-FU", 
                                 names(vcfs_mm_treatment),
                                 value = T)])),
                    as.data.frame(rowSums(
                      vcfs_mm_treatment[grep("^STE00", 
                                             names(vcfs_mm_treatment),
                                             value = T)])
                      )
  )

plot_96_profile(vcfs_mm_treatment)

colnames(vcfs_mm_treatment) <- c("5-FU", "Control")
plot_compare_profiles(vcfs_mm_treatment[,1], 
                      vcfs_mm_treatment[,2], 
                      profile_names = c("5-FU", "Control"),
                      condensed = TRUE)

cos_sim(vcfs_mm_treatment[,1], vcfs_mm_treatment[,2])

# generate 5FU specific signature
vcfs_mm_treatment <- vcfs_mm_treatment %>% 
  mutate(relative_5FU = vcfs_mm_treatment[,1]/ 
           sum(vcfs_mm_treatment[,1]),
    relative_control = vcfs_mm_treatment[,2]/ 
      sum(vcfs_mm_treatment[,2]),
    diff = relative_5FU - relative_control,
    diff_5FU = ifelse(diff >0, diff, 0),
    diff_control = ifelse(diff <0, abs(diff), 0))

plot_96_profile(vcfs_mm_treatment[,6:7],ymax = 0.4)
cos_sim(vcfs_mm_treatment[,6], vcfs_mm_treatment[,7])

### compare mutation profiles with COSMIC signatures
#load COSMIC signatures
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]


#method 1: least square regression fitting
fit_res <- fit_to_signatures(vcfs_mm_treatment[,1:2], cancer_signatures)
fit_res$contribution
#select <- which(rowSums(fit_res$contribution) > 200)
select <- which((rowSums(fit_res$contribution)/sum(fit_res$contribution)*100) > 10)
plot_contribution(fit_res$contribution[select,],
                                      cancer_signatures[,select],
                                       coord_flip = FALSE,
                                       mode = "relative")
#method 2: calculate cosine similarity
cos_sim_samples = cos_sim_matrix(as.matrix(vcfs_mm_treatment[,6:7]),cancer_signatures)
plot_cosine_heatmap(cos_sim_samples,
                    col_order = cosmic_order,
                    cluster_rows = TRUE)


#####----#### de novo sigantures
#de novo extraction with breast cancer samples from Norman/Joep et al
#link to VCFs
#You can also read the mut matrix for the breast cancer dataset to gain time
vcf_list <- list.files("~/Documents/courses/Omics_course/VCFs/Exercise_2/", pattern = ".vcf", full.names = TRUE)
#load VCFs in GRangesList (I have downsampled the MSI VCFs to gain time)
#VCFs still contain indels, will not be uploaded in the grange object
vcfs <- read_vcfs_as_granges(vcf_list, vcf_files_names, genome = ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

#optinally
mut_mat = read.table("~/Documents/courses/Omics_course/VCFs/Exercise_2/mut_mat_breast_organoids.txt", sep = "\t")

#extract signatures using the NMF approach
library("NMF")
mut_mat = mut_mat + 0.00001
#estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)
nmf_res <- extract_signatures(mut_mat, rank = 5, nrun = 10)

#Give your babies a name
signatures = c("sign.1","sign.2","sign.3","sign.4","sign.5")
colnames(nmf_res$signatures) = signatures
rownames(nmf_res$contribution) = signatures

plot_96_profile(nmf_res$signatures, condensed = TRUE)
plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "relative")
plot_contribution_heatmap(nmf_res$contribution, cluster_samples=TRUE)

#cluster the COSMICs
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
plot(hclust_cosmic)

#Calculate pairwise cosine similarity between mutational profiles and COSMIC signatures
cos_sim_samples_signatures = cos_sim_matrix(nmf_res$signatures, cancer_signatures)
#and plot
plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order, cluster_rows = TRUE,plot_values = TRUE)


