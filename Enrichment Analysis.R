library(dplyr)
library(readr)
library(tidyr)
library(biomaRt)
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)


deg <- read_csv("output/c2c1_degs.csv") |>
  dplyr::rename(ENSEMBL = Ensembl_ID)

ens2ntrz <- bitr(deg$ENSEMBL, fromType = "ENSEMBL", 
     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

deg <- deg |>
  left_join(ens2ntrz, by = "ENSEMBL")

deg <- deg |> 
  drop_na(ENTREZID, adj.P.Val)

length(unique(deg$ENTREZID)) 
nrow(deg)

deg <- deg |> 
  arrange(adj.P.Val) |> 
  distinct(ENTREZID, .keep_all = TRUE)

length(unique(deg$ENTREZID)) == nrow(deg)
write.csv(deg, "deg_entrezid.csv")
sig_genes <- deg |> 
  dplyr::filter(adj.P.Val < 0.05) |>
  dplyr::filter(logFC >= 1.5 | logFC <= -1.5) |> 
  dplyr::pull(ENTREZID)

go_ora <- enrichGO(gene = as.character(sig_genes),
                   OrgDb = org.Hs.eg.db,
                   universe = as.character(deg$ENTREZID),
                   ont = "MF",
                   readable = TRUE) 

barplot(go_ora, showCategory = 10)
dotplot(go_ora, showCategory = 10)

write.table(go_ora, file="output/enrich/go_MF.tsv", sep="\t", col.names=NA, quote=F)

mkegg <- enrichMKEGG(gene = sig_genes,
           universe = as.character(deg$ENTREZID)) 

write.table(mkegg, file="output/enrich/kegg/mkegg.tsv", sep="\t", col.names=NA, quote=F)

barplot(mkegg, showCategory = 10)
dotplot(mkegg, showCategory = 10)


  t() |>
  as.data.frame(stringsAsFactors=FALSE) 

colnames(df) <- df[1,]
df <- df[-1,] |>
  tibble::rownames_to_column(var = "Samples")

df <- merge(df, cl, by="Samples", all = FALSE)

df$Group <- factor(df$Group)
df$ENSG00000080839 <- as.numeric(df$ENSG00000080839)
df$ENSG00000110852 <- as.numeric(df$ENSG00000110852)
ggplot(df, aes(x = Group, y = ENSG00000080839, color = Group)) +
  geom_boxplot() +
  ylab("Exp") +
  xlab("Cluster") +
  labs(title = "Expression by Cluster") +
  theme_minimal()


