# Run once on the cluster to record the exact package versions:
#   module load R/4.5.0/gnu-12.2.0
#   Rscript environment/capture_sessionInfo.R
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
pkgs <- c("Seurat","SeuratObject","sctransform","DESeq2","fgsea","msigdbr",
          "CellChat","SingleR","SCINA","dplyr","Matrix","ggplot2","ggrepel",
          "httr2","jsonlite")
suppressMessages(invisible(lapply(pkgs, function(p)
  if (requireNamespace(p, quietly=TRUE)) library(p, character.only=TRUE))))
writeLines(capture.output(sessionInfo()), "environment/sessionInfo.txt")
cat("Wrote environment/sessionInfo.txt\n")
