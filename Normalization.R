setwd(path)

library(tidyverse)
library(magrittr)

early <- read_csv2("Early.csv") |> column_to_rownames(var = "gene")
late <- read_csv2("Late.csv") |> column_to_rownames(var = "gene")

mi <- 1/ncol(early) * rowSums(early)
L <- log(late/mi)

write.csv2(L, file = "L.csv")


