#!/bin/bash
#SBATCH --job-name=P18excl_mye03
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/mye03_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/mye03_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'

.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)
library(Seurat)
library(dplyr)
library(ggplot2)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
FIGDIR  <- "/ibex/user/alghamlk/HCC_P18_excluded/figures"

cat("############################################################\n")
cat("Myeloid res 0.3: markers + contamination + lineage check\n")
cat("(clustering already computed; just re-key idents to 0.3)\n")
cat("############################################################\n")
mye <- readRDS(file.path(DATADIR, "seurat_P18excl_myeloid_subclustered.rds"))
Idents(mye) <- mye@meta.data[["SCT_snn_res.0.3"]]
cat("res 0.3 subclusters:", length(unique(Idents(mye))), "\n")
print(table(Idents(mye)))
cat("\n")

# ── Markers at res 0.3 ──────────────────────────────────────────
mye <- PrepSCTFindMarkers(mye, verbose = FALSE)
mk <- FindAllMarkers(mye, assay = "SCT", only.pos = TRUE,
                     min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
write.csv(mk, file.path(DATADIR, "P18excl_myeloid_subcluster_markers_res0.3.csv"),
          row.names = FALSE)
cat("Top 15 markers per subcluster (res 0.3):\n")
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 15) %>%
  summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop")
print(as.data.frame(top), right = FALSE)
cat("\n")

# ── Targeted lineage/state check via mean expression per subcluster
# (RNA log-norm). Answers: which is contamination, is there SPP1-TAM,
# where are the DC states -- all WITHOUT typed signature scoring,
# just raw mean expression of canonical genes for inspection. ────
DefaultAssay(mye) <- "RNA"
mye <- NormalizeData(mye, verbose = FALSE)
probe <- c(
  # myeloid lineage
  "LYZ","CD68","ITGAM",
  # TAM states
  "C1QC","C1QA","APOE","SPP1","TREM2","GPNMB","MRC1","CD163",
  # monocyte
  "FCN1","S100A8","S100A9","VCAN",
  # DC
  "FCER1A","CD1C","CLEC9A","XCR1","LAMP3","CCR7",
  # contamination sentinels (T / NK / hepatocyte / B)
  "CD3E","CD3D","IL7R","NKG7","ALB","MS4A1"
)
probe <- intersect(probe, rownames(mye))
avg <- AverageExpression(mye, assays = "RNA", features = probe, layer = "data")$RNA
cat("############################################################\n")
cat("Mean expression per subcluster (inspect, don't auto-name):\n")
cat("Look for: high CD3E/IL7R/ALB = contamination to DROP;\n")
cat("SPP1/TREM2 = SPP1-TAM?; LAMP3/CCR7 = mregDC; XCR1 = cDC1.\n")
cat("############################################################\n")
print(round(avg, 2))
cat("\n")

# ── flag likely contamination automatically (high T/hepatocyte,
# low myeloid) so you can SEE it, not hide it ───────────────────
cat("Contamination screen (mean expr): CD3E / IL7R / ALB vs LYZ/CD68\n")
flag <- data.frame(
  subcluster = colnames(avg),
  CD3E = round(avg["CD3E", ], 2),
  IL7R = if ("IL7R" %in% rownames(avg)) round(avg["IL7R", ], 2) else NA,
  ALB  = if ("ALB"  %in% rownames(avg)) round(avg["ALB", ], 2) else NA,
  LYZ  = round(avg["LYZ", ], 2),
  CD68 = round(avg["CD68", ], 2))
print(flag, row.names = FALSE)
cat("\n-> subclusters with high CD3E/IL7R/ALB and low LYZ/CD68 are\n")
cat("   non-myeloid contamination; drop before CellChat.\n\n")

# UMAP at 0.3
p <- DimPlot(mye, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle(sprintf("Myeloid subclusters res 0.3 (n=%d)", ncol(mye)))
ggsave(file.path(FIGDIR, "umap_myeloid_subclusters_res0.3.png"), p,
       width = 9, height = 7, dpi = 300)
cat("Saved umap_myeloid_subclusters_res0.3.png\n")

saveRDS(mye, file.path(DATADIR, "seurat_P18excl_myeloid_res03.rds"))
cat("Saved seurat_P18excl_myeloid_res03.rds (idents = res 0.3)\n\n")
cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
