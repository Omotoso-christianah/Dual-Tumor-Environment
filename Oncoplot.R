
setwd("~/Desktop/Desktop/Xena/Oncoplot_6")
#Running Oncoplot for all clusters#
BiocManager::install("maftools")
library(maftools)
install.packages("tidyverse")
library(tidyverse)



#Running Oncoplot for cluster 1
laml<- read.csv("./Maf_data_analysis/Final_maf_gene_Cluster1.tsv", sep = ",", header = TRUE)

head(laml)

clin <- read.csv("./Clin_Laml_analysis/Clin_Subtype_1.tsv", header = TRUE, sep = ",",row.names = NULL)

head(clin)

laml2 = read.maf(maf = laml,
                 clinicalData = clin,
                 verbose = FALSE,vc_nonSyn = c("frameshift_variant" ,"stop_gained"    ,    "missense_variant"  , "synonymous_variant"))
laml2



oncoplot(maf = laml2, draw_titv = TRUE)
vc_cols = RColorBrewer::brewer.pal(n = 5, name = 'Paired')
names(vc_cols) = c("frameshift_variant" ,"stop_gained"    ,    "missense_variant"  , "synonymous_variant")

print(vc_cols)

oncoplot(maf = laml2, colors = vc_cols, top = 11)


#Running Oncoplot for cluster 2
laml<- read.csv("./Maf_data_analysis/Final_maf_gene_Cluster2.tsv", sep = ",", header = TRUE)

head(laml)

clin <- read.csv("./Clin_Laml_analysis/Clin_Subtype_2.tsv", header = TRUE, sep = ",",row.names = NULL)

head(clin)

laml2 = read.maf(maf = laml,
                 clinicalData = clin,
                 verbose = FALSE,vc_nonSyn = c("frameshift_variant" ,"stop_gained"    ,    "missense_variant"  , "synonymous_variant"))
laml2



oncoplot(maf = laml2, draw_titv = TRUE)
vc_cols = RColorBrewer::brewer.pal(n = 5, name = 'Paired')
names(vc_cols) = c("frameshift_variant" ,"stop_gained"    ,    "missense_variant"  , "synonymous_variant")

print(vc_cols)

oncoplot(maf = laml2, colors = vc_cols, top = 11)














