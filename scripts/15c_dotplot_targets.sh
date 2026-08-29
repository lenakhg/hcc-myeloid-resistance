#!/bin/bash
#SBATCH --job-name=dotplot_targets
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48GB
#SBATCH --time=00:25:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/dotplot_targets_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/dotplot_targets_%J.err

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
cat("Label column:", labcol, "\n")
Idents(mye) <- mye@meta.data[[labcol]]

# drop contamination
keep_id <- setdiff(levels(Idents(mye)), grep("ontam", levels(Idents(mye)), value=TRUE))
mye2 <- subset(mye, idents=keep_id); Idents(mye2) <- droplevels(Idents(mye2))

# biological order
bio_order <- c("C1Q/APOE-TAM","Perivascular-TAM","MRC1/CD9-TAM","Inflammatory-TAM",
               "FCN1-monocyte","S100A8-monocyte","cDC1","cDC2","mregDC")
bio_order <- bio_order[bio_order %in% levels(Idents(mye2))]

# nominated targets + myeloid PD genes (grouped: targets first, then extra PD genes)
target_genes <- c("CD14","FCGR1A","FCGR3A","SYK","MERTK","FLT3","ADORA3")
extra_pd     <- c("SIGLEC7","CD300LF","LRRC25","CPM","FPR1","ALOX15B")
genes <- c(target_genes, extra_pd)
genes <- genes[genes %in% rownames(mye2)]
cat("Genes present:", length(genes), "of", length(target_genes)+length(extra_pd), "\n")
missing <- setdiff(c(target_genes,extra_pd), rownames(mye2))
if(length(missing)) cat("MISSING from data:", paste(missing, collapse=", "), "\n")

# ===== FIGURE 1: targets across the 9 myeloid states =====
mye2@active.ident <- factor(Idents(mye2), levels=rev(bio_order))
p1 <- DotPlot(mye2, features=genes, dot.scale=8) +
  scale_color_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0) +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=12, face="italic"),
        axis.text.y = element_text(size=12),
        axis.title = element_blank(),
        panel.grid.major = element_line(color="grey92"),
        plot.title = element_text(face="bold")) +
  labs(title="Nominated targets & myeloid PD-elevated genes across myeloid states")
# annotate target vs extra with a vertical line after the 7 targets
p1 <- p1 + geom_vline(xintercept=length(target_genes)+0.5, linetype="dashed", color="grey50")
ggsave(file.path(FIG,"dotplot_targets_by_state.png"), p1,
       width=11, height=6, dpi=300, bg="white")
cat("Saved dotplot_targets_by_state.png\n")

# ===== FIGURE 2: split by response (PR vs PD), faceted =====
# find the response column
respcol <- NULL
for (cn in c("recist","response","RECIST","group","condition")) {
  if (cn %in% colnames(mye2@meta.data)) { respcol <- cn; break }
}
cat("Response column:", respcol, "\n")
if (!is.null(respcol)) {
  print(table(mye2@meta.data[[respcol]]))
  # keep only PR and PD
  mye2$resp <- mye2@meta.data[[respcol]]
  mye3 <- subset(mye2, subset = resp %in% c("PR","PD"))
  mye3@active.ident <- factor(Idents(mye3), levels=rev(bio_order))
  p2 <- DotPlot(mye3, features=genes, dot.scale=7, split.by="resp",
                cols=c("#2166AC","#B2182B")) +
    theme_bw(base_size=12) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=11, face="italic"),
          axis.text.y = element_text(size=10),
          axis.title = element_blank(),
          plot.title = element_text(face="bold")) +
    labs(title="Target gene expression by myeloid state, split by response (PR vs PD)")
  ggsave(file.path(FIG,"dotplot_targets_by_state_PRvsPD.png"), p2,
         width=11, height=9, dpi=300, bg="white")
  cat("Saved dotplot_targets_by_state_PRvsPD.png\n")

  # Alternative faceted version (cleaner): two panels side by side
  library(tidyr)
  p2b <- DotPlot(mye3, features=genes, dot.scale=7, group.by=labcol, split.by="resp") 
} else {
  cat("No response column found - skipping Figure 2. Check metadata columns:\n")
  print(colnames(mye2@meta.data))
}

cat("Done:", as.character(Sys.time()), "\n")
EOF
echo "Job finished: $(date)"
