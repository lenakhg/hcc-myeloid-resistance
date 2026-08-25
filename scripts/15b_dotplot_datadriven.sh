#!/bin/bash
#SBATCH --job-name=dotplot_dd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48GB
#SBATCH --time=00:25:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/dotplot_dd_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/dotplot_dd_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)
suppressMessages({library(Seurat); library(ggplot2); library(dplyr)})

DATA <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
FIG  <- "/ibex/user/alghamlk/HCC_P18_excluded/figures"

mye <- readRDS(file.path(DATA,"seurat_P18excl_myeloid_res04_named.rds"))
DefaultAssay(mye) <- "RNA"
mye <- NormalizeData(mye, verbose=FALSE)

# label column
labcol <- NULL
for (cn in c("myeloid_label","cell_type","subset","ident")) {
  if (cn %in% colnames(mye@meta.data)) { labcol <- cn; break }
}
Idents(mye) <- mye@meta.data[[labcol]]

# drop contamination
keep_id <- setdiff(levels(Idents(mye)), grep("ontam", levels(Idents(mye)), value=TRUE))
mye2 <- subset(mye, idents=keep_id); Idents(mye2) <- droplevels(Idents(mye2))

# cluster number -> name map (from run_12)
name_map <- c("1"="C1Q/APOE-TAM","2"="cDC2","3"="FCN1-monocyte","4"="MRC1/CD9-TAM",
              "5"="S100A8-monocyte","6"="Contamination (T/hepatocyte)","7"="cDC2",
              "8"="Inflammatory-TAM","9"="cDC1","10"="Perivascular-TAM","11"="mregDC")

# biological display order
bio_order <- c("C1Q/APOE-TAM","Perivascular-TAM","MRC1/CD9-TAM","Inflammatory-TAM",
               "FCN1-monocyte","S100A8-monocyte","cDC1","cDC2","mregDC")

# read markers file, map cluster->name, keep significant, rank by avg_log2FC
mk <- read.csv(file.path(DATA,"P18excl_myeloid_subcluster_markers_res0.4_final.csv"))
mk$state <- name_map[as.character(mk$cluster)]
mk <- mk[!is.na(mk$state) & mk$state %in% bio_order, ]
mk <- mk[mk$p_val_adj < 0.05, ]

build_genes <- function(topN) {
  genes <- c()
  for (st in bio_order) {
    sub <- mk[mk$state==st, ]
    sub <- sub[order(-sub$avg_log2FC), ]
    # take topN not already used (avoid duplicates across states)
    picks <- setdiff(sub$gene, genes)[1:topN]
    picks <- picks[!is.na(picks)]
    genes <- c(genes, picks)
  }
  unique(genes[genes %in% rownames(mye2)])
}

make_plot <- function(topN, fname) {
  genes <- build_genes(topN)
  cat("top", topN, ":", length(genes), "genes\n")
  # order idents biologically
  mye2@active.ident <- factor(Idents(mye2), levels=rev(bio_order))
  p <- DotPlot(mye2, features=genes, dot.scale=7) +
    scale_color_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0) +
    theme_bw(base_size=13) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=10, face="italic"),
          axis.text.y = element_text(size=12),
          axis.title = element_blank(),
          legend.position = "right",
          panel.grid.major = element_line(color="grey92")) +
    labs(title=sprintf("Myeloid subset marker expression (top %d per state, data-driven)", topN))
  wid <- 6 + length(genes)*0.32
  ggsave(file.path(FIG,fname), p, width=wid, height=6.5, dpi=300, bg="white", limitsize=FALSE)
  cat("Saved", fname, "\n")
}

make_plot(3, "dotplot_myeloid_datadriven_top3.png")
make_plot(4, "dotplot_myeloid_datadriven_top4.png")
make_plot(5, "dotplot_myeloid_datadriven_top5.png")

cat("Done:", as.character(Sys.time()), "\n")
EOF
echo "Job finished: $(date)"
