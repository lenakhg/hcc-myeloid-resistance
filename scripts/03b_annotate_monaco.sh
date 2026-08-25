#!/bin/bash
#SBATCH --job-name=P18excl_monaco
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/monaco_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/monaco_%J.err
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
library(SingleR)
library(celldex)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("MonacoImmuneData -- purpose-built for fine immune/T-cell\n")
cat("subsets. Biologically appropriate here since this dataset\n")
cat("is CD45+ sorted (immune cells only), not general tissue.\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Clusters present:", length(unique(Idents(seurat))), "\n\n")

ref <- tryCatch({
  MonacoImmuneData()
}, error = function(e) {
  stop("FAILED to download MonacoImmuneData: ", conditionMessage(e))
})
cat("Reference loaded. Fine label categories available:\n")
print(table(ref$label.fine))
cat("\n")

query_data <- GetAssayData(seurat, assay = "SCT", layer = "data")
singler_results <- SingleR(test = query_data, ref = ref,
                           labels = ref$label.fine,
                           clusters = Idents(seurat))

cat("SingleR (Monaco, fine) results per cluster:\n")
print(data.frame(cluster = rownames(singler_results),
                 label = singler_results$labels,
                 pruned_label = singler_results$pruned.labels))

annotation_df <- data.frame(
  cluster = rownames(singler_results),
  auto_annotation = ifelse(is.na(singler_results$pruned.labels),
                           paste0(singler_results$labels, " (low confidence)"),
                           singler_results$labels)
)
write.csv(annotation_df, file.path(DATADIR, "P18excl_cluster_annotations_Monaco.csv"), row.names = FALSE)
cat("\nSaved P18excl_cluster_annotations_Monaco.csv\n\n")

# Cross-validation against author's Myeloid label (Monaco calls
# monocytes/myeloid something like "Classical monocytes" etc.)
cluster_labels <- setNames(annotation_df$auto_annotation, annotation_df$cluster)
seurat@meta.data$auto_annotation_monaco <- cluster_labels[as.character(Idents(seurat))]
myeloid_hits <- grep("monocyte|macrophage|dendritic", cluster_labels, value = TRUE, ignore.case = TRUE)
cat("Cluster label(s) matched to myeloid identity:", paste(unique(myeloid_hits), collapse=", "), "\n")
if (length(myeloid_hits) > 0) {
  our_cells <- sum(seurat@meta.data$auto_annotation_monaco %in% myeloid_hits)
  overlap <- sum(seurat@meta.data$auto_annotation_monaco %in% myeloid_hits &
                seurat@meta.data$major_cluster == "Myeloid")
  author_total <- sum(seurat@meta.data$major_cluster == "Myeloid")
  cat(sprintf("Precision: %.1f%%\n", 100 * overlap / our_cells))
  cat(sprintf("Recall: %.1f%%\n\n", 100 * overlap / author_total))
}

cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
