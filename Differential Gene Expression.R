setwd("~/Desktop/Desktop/Xena/DEG_Michelle")
library(limma)
y <- read.csv2("L_Ordered.tsv", header = TRUE, row.names = 1, sep = ",", dec = ",")
y[1:3,1:3]
# Experimental groups
# Samples are in order of clusters
groups <-  c(rep("C1",16), rep("C2",23))
cluster <- factor(groups, levels = c("C1","C2"))

# Design object - without intercept
design <- model.matrix(~0+cluster)
colnames(design) <- c("C1","C2")

# Contrast matrix
constrasts <- makeContrasts(C1-C2,levels = design)

# Fit with no intercept 
fit <- eBayes(contrasts.fit(lmFit(y, design = design), constrasts))

#topTable
tabl <- topTable(fit, number = Inf)
head(tabl)
write.csv2(tabl, "tabl.csv")

decide<- decideTests(fit,method="separate",adjust.method="BH",p.value=0.05,lfc=1.5)
head(decide)

write.csv2(decide, "decide.csv")

# output filter of sign fold
library(tidyverse)

colnames(decide) <- c("c12")

df <- decide |> 
  as.data.frame() |>
  dplyr::filter(c12 != 0 ) 
dim(df)
write.csv2(df, "decide_genes_upordownonly.csv")

L <- read.csv2("L_Ordered.tsv", header = TRUE, row.names = 1, sep = ",", dec = ",")

exp<- L |>
  filter(row.names(L) %in% row.names(df)) 


write.csv(exp, "L1855genes.csv")




# Extract DE genes 
# Set pval
all_de_genes <- topTable(fit, adjust="fdr", p.value=0.05, number=Inf)

# Output CSV file
write.table(all_de_genes, file="all_de_genes_pval0.05.csv",row.names=T, sep=";")

