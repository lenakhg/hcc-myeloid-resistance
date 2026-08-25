#!/bin/bash
#SBATCH --job-name=P18excl_celltypist
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=02:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/celltypist_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/celltypist_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'

.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)

library(Seurat)
library(Matrix)
library(dplyr)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("PART 1 (R): Export raw counts for CellTypist\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Cells loaded:", ncol(seurat), "\n")
cat("Clusters present:", length(unique(Idents(seurat))), "\n\n")

counts <- GetAssayData(seurat, assay = "RNA", layer = "counts")
cat("Counts matrix: ", paste(dim(counts), collapse=" x "),
    " nonzero=", length(counts@x), "\n\n")

writeMM(counts, file.path(DATADIR, "celltypist_counts.mtx"))
write.table(rownames(counts), file.path(DATADIR, "celltypist_genes.tsv"),
           row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(counts), file.path(DATADIR, "celltypist_barcodes.tsv"),
           row.names = FALSE, col.names = FALSE, quote = FALSE)
cluster_map <- data.frame(barcode = colnames(seurat), cluster = as.character(Idents(seurat)))
write.csv(cluster_map, file.path(DATADIR, "celltypist_cluster_map.csv"), row.names = FALSE)
cat("Exported matrix, genes, barcodes, cluster map\n\n")

cat("############################################################\n")
cat("PART 2 (Python via reticulate + py_require)\n")
cat("############################################################\n")
library(reticulate)
py_require(c("scanpy", "celltypist", "numpy", "scipy", "pandas"))

sc <- import("scanpy")
celltypist <- import("celltypist")
pd <- import("pandas")

cat("Reading exported matrix into scanpy AnnData...\n")
adata <- sc$read_mtx(file.path(DATADIR, "celltypist_counts.mtx"))
adata <- adata$T

genes <- readLines(file.path(DATADIR, "celltypist_genes.tsv"))
barcodes <- readLines(file.path(DATADIR, "celltypist_barcodes.tsv"))
adata$var_names <- genes
adata$obs_names <- barcodes
cat("AnnData shape (cells x genes):", paste(adata$shape, collapse=" x "), "\n\n")

cat("Normalizing (CP10K + log1p)...\n")
sc$pp$normalize_total(adata, target_sum = 1e4)
sc$pp$log1p(adata)
cat("Normalization complete.\n\n")

cat("Downloading/confirming model...\n")
celltypist$models$download_models(model = "Immune_All_Low.pkl")
cat("Model ready: Immune_All_Low.pkl\n\n")

cat("Running celltypist.annotate()...\n")
predictions <- celltypist$annotate(adata, model = "Immune_All_Low.pkl", majority_voting = TRUE)

cat("Diagnostic: checking how reticulate handed back this object\n")
cat("(discovered last run: reticulate auto-converts pandas\n")
cat("DataFrames to native R data.frames on access, so no further\n")
cat("Python-side .to_csv() call is needed or valid here)...\n")
result_r <- predictions$predicted_labels
cat("Class:", paste(class(result_r), collapse=", "), "\n")
cat("Dimensions:", paste(dim(result_r), collapse=" x "), "\n")
cat("Column names:", paste(colnames(result_r), collapse=", "), "\n")
cat("Row names (first 5):", paste(head(rownames(result_r), 5), collapse=", "), "\n\n")

result_r$barcode <- rownames(result_r)
write.csv(result_r, file.path(DATADIR, "celltypist_predictions.csv"), row.names = FALSE)
cat("Saved celltypist_predictions.csv\n\n")

cat("############################################################\n")
cat("PART 3 (R): Merge predictions back, per-cluster summary,\n")
cat("myeloid cross-validation\n")
cat("############################################################\n")
preds <- read.csv(file.path(DATADIR, "celltypist_predictions.csv"))
cat("Predictions loaded:", nrow(preds), "rows\n")
cat("Columns:", paste(colnames(preds), collapse=", "), "\n\n")

merged <- merge(cluster_map, preds, by = "barcode")
cat("Merged (barcode-matched) rows:", nrow(merged),
    "of", nrow(cluster_map), "expected -- ",
    ifelse(nrow(merged) == nrow(cluster_map), "FULL MATCH", "MISMATCH -- investigate"), "\n\n")

majority_vote_col <- grep("majority_voting", colnames(merged), value = TRUE)
if (length(majority_vote_col) == 0) majority_vote_col <- grep("predicted_labels", colnames(merged), value = TRUE)[1]
cat("Using column for per-cluster summary:", majority_vote_col, "\n\n")

cluster_summary <- merged %>%
  group_by(cluster) %>%
  count(.data[[majority_vote_col]]) %>%
  slice_max(n, n = 1) %>%
  ungroup()
cat("Per-cluster majority CellTypist label:\n")
print(cluster_summary, n = 20)
cat("\n")

write.csv(cluster_summary, file.path(DATADIR, "P18excl_cluster_annotations_CellTypist.csv"), row.names = FALSE)
cat("Saved P18excl_cluster_annotations_CellTypist.csv\n\n")

cluster_labels <- setNames(as.character(cluster_summary[[majority_vote_col]]), cluster_summary$cluster)
seurat@meta.data$auto_annotation_celltypist <- cluster_labels[as.character(Idents(seurat))]

myeloid_hits <- grep("mono|macro|dendritic|myeloid", cluster_labels, value = TRUE, ignore.case = TRUE)
cat("Cluster label(s) matched to myeloid identity:", paste(unique(myeloid_hits), collapse=", "), "\n")
if (length(myeloid_hits) > 0) {
  our_cells <- sum(seurat@meta.data$auto_annotation_celltypist %in% myeloid_hits, na.rm = TRUE)
  overlap <- sum(seurat@meta.data$auto_annotation_celltypist %in% myeloid_hits &
                seurat@meta.data$major_cluster == "Myeloid", na.rm = TRUE)
  author_total <- sum(seurat@meta.data$major_cluster == "Myeloid")
  cat(sprintf("Precision: %.1f%%\n", 100 * overlap / our_cells))
  cat(sprintf("Recall: %.1f%%\n\n", 100 * overlap / author_total))
}

cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
