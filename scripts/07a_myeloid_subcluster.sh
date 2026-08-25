#!/bin/bash
#SBATCH --job-name=P18excl_myeloid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128GB
#SBATCH --time=03:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/myeloid_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/myeloid_%J.err
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
cat("MYELOID SUB-CLUSTERING. Subset Leiden cluster 2 as-is.\n")
cat("REUSE existing SCT normalization; re-run only the\n")
cat("dimensionality reduction + clustering on the subset\n")
cat("(FindVariableFeatures -> PCA -> UMAP -> Neighbors -> Clusters).\n")
cat("Resolution sweep 0.2-0.8, then FindAllMarkers. No naming.\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Whole-tumor cells:", ncol(seurat), "\n")

# ── Subset cluster 2 as-is ──────────────────────────────────────
mye <- subset(seurat, idents = "2")
cat("Myeloid cells (cluster 2):", ncol(mye), "\n")
rm(seurat); gc()

# ── REUSE SCT. Just re-select variable features on the subset and
# re-run PCA/UMAP/clustering so the reduction captures myeloid-
# INTERNAL variation. No re-normalization. ──────────────────────
DefaultAssay(mye) <- "SCT"
cat("\nSCT assay reused (not re-fit). Re-selecting variable features",
    "on myeloid-only cells...\n")
mye <- FindVariableFeatures(mye, assay = "SCT", nfeatures = 2000, verbose = FALSE)
cat("Variable features (myeloid-specific):", length(VariableFeatures(mye)), "\n")

mye <- RunPCA(mye, npcs = 30, verbose = FALSE)
mye <- RunUMAP(mye, dims = 1:20, seed.use = 42, verbose = FALSE)
mye <- FindNeighbors(mye, dims = 1:20, verbose = FALSE)
cat("PCA + UMAP + Neighbors done (20 PCs).\n\n")

# ── Resolution sweep 0.2-0.8 ────────────────────────────────────
cat("############################################################\n")
cat("Resolution sweep (Leiden, algorithm=4)\n")
cat("############################################################\n")
res_vals <- seq(0.2, 0.8, by = 0.1)
for (r in res_vals) {
  mye <- FindClusters(mye, resolution = r, algorithm = 4, verbose = FALSE)
}
cat("Subclusters per resolution:\n")
sweep_tbl <- data.frame(resolution = res_vals,
  n_subclusters = sapply(res_vals, function(r)
    length(unique(mye@meta.data[[paste0("SCT_snn_res.", r)]]))))
print(sweep_tbl, row.names = FALSE)
write.csv(sweep_tbl, file.path(DATADIR, "P18excl_myeloid_resolution_sweep.csv"),
          row.names = FALSE)
cat("\nCluster sizes at each resolution:\n")
for (r in res_vals) {
  col <- paste0("SCT_snn_res.", r)
  cat(sprintf("--- res %.1f ---\n", r))
  print(table(mye@meta.data[[col]]))
}
cat("\n")

# ── Pick a working resolution just to RUN markers on. This is a
# provisional handle, NOT a committed choice -- all sweep columns
# are saved in the object so you can re-pick after seeing markers.
WORK_RES <- 0.4
res_col <- paste0("SCT_snn_res.", WORK_RES)
Idents(mye) <- mye@meta.data[[res_col]]
cat("Running markers at resolution", WORK_RES, "->",
    length(unique(Idents(mye))), "subclusters (provisional).\n")
print(table(Idents(mye)))
cat("\n")

# ── Subcluster markers (SCT assay, PrepSCTFindMarkers first) ────
cat("############################################################\n")
cat("FindAllMarkers on myeloid subclusters (res", WORK_RES, ")\n")
cat("############################################################\n")
mye <- PrepSCTFindMarkers(mye, verbose = FALSE)
mk <- FindAllMarkers(mye, assay = "SCT", only.pos = TRUE,
                     min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
write.csv(mk, file.path(DATADIR, "P18excl_myeloid_subcluster_markers_res0.4.csv"),
          row.names = FALSE)
cat("Total marker rows:", nrow(mk), "\n\n")
cat("Top 12 markers per subcluster (by log2FC):\n")
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12) %>%
  summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop")
print(as.data.frame(top), right = FALSE)
cat("\n")

# ── Save object (all sweep resolutions retained) + UMAP ─────────
saveRDS(mye, file.path(DATADIR, "seurat_P18excl_myeloid_subclustered.rds"))
cat("Saved seurat_P18excl_myeloid_subclustered.rds\n")
cat("(contains SCT_snn_res.0.2 ... 0.8 -- re-pick resolution anytime)\n")

p <- DimPlot(mye, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle(sprintf("Myeloid subclusters (res %.1f, n=%d)", WORK_RES, ncol(mye)))
ggsave(file.path(FIGDIR, "umap_myeloid_subclusters_res0.4.png"), p,
       width = 9, height = 7, dpi = 300)
cat("Saved umap_myeloid_subclusters_res0.4.png\n\n")

cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
