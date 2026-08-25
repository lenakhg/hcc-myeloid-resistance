#!/bin/bash
#SBATCH --job-name=P18excl_mye04
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/mye04_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/mye04_%J.err
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
RES_COL <- "SCT_snn_res.0.4"

cat("############################################################\n")
cat("Myeloid res 0.4 finalize: data-driven markers only.\n")
cat("No typed gene lists. FindAllMarkers describes each subcluster;\n")
cat("interpretation (vascular / TAM / DC / contamination) is done\n")
cat("by the analyst from the marker output, not by the script.\n")
cat("############################################################\n")
mye <- readRDS(file.path(DATADIR, "seurat_P18excl_myeloid_subclustered.rds"))
Idents(mye) <- mye@meta.data[[RES_COL]]
cat("res 0.4 subclusters:", length(unique(Idents(mye))), "\n")
print(table(Idents(mye)))
cat("\n")

DefaultAssay(mye) <- "SCT"
mye <- PrepSCTFindMarkers(mye, verbose = FALSE)
mk <- FindAllMarkers(mye, assay = "SCT", only.pos = TRUE,
                     min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
write.csv(mk, file.path(DATADIR, "P18excl_myeloid_subcluster_markers_res0.4_final.csv"),
          row.names = FALSE)
cat("Total marker rows:", nrow(mk), "\n\n")

cat("Top 25 markers per subcluster (by log2FC) -- read these to\n")
cat("interpret identity yourself:\n")
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 25) %>%
  summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop")
print(as.data.frame(top), right = FALSE)
cat("\n")

# UMAP colored by subcluster number (no naming)
p <- DimPlot(mye, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle(sprintf("Myeloid subclusters res 0.4 (n=%d)", ncol(mye)))
ggsave(file.path(FIGDIR, "umap_myeloid_subclusters_res0.4_final.png"), p,
       width = 9, height = 7, dpi = 300)
cat("Saved umap_myeloid_subclusters_res0.4_final.png\n")

saveRDS(mye, file.path(DATADIR, "seurat_P18excl_myeloid_res04.rds"))
cat("Saved seurat_P18excl_myeloid_res04.rds (idents = res 0.4)\n\n")
cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
