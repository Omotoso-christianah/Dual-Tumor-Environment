library(dendextend)
library(flashClust)
library(dplyr)
library(data.table)
library(magrittr)
library(RColorBrewer)
library(wordspace)

Tumor <- read.csv("L_new.csv", sep = ",", header = TRUE, row.names = 1, dec = ".")
Tumor.t <- t(Tumor)
Tumor[1:3,1:3]

sampleDist <- Tumor.t %>% dist.matrix(method = "cosine") %>% as.dist

heatDendro <-
  sampleDist %>%
  flashClust(method = "ward") %>%
  as.dendrogram

plot(heatDendro)

heatDendro %>% dendextend::color_branches(.,
                                          h = 400,
                                          col = c("#00AFBB",  "#FC4E07"),
                                          groupLabels = FALSE) %>%
  plot
plot
groups <- heatDendro %>% cutree(k = 2)


# Write samples and cluster number to an excel sheet
groups %>%
  data.frame("Samples" = names(.),
             "Group" = .
  ) %>%
  fwrite("~/Desktop/Desktop/Xena/HC_1/HC_Group_Assignments2.csv",
         sep = ",",
         row.names = FALSE,
         col.names = TRUE
  )
