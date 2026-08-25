#!/bin/bash
#SBATCH --job-name=P18excl_silhouette
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128GB
#SBATCH --time=05:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/silhouette_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/silhouette_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'

.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)

library(Seurat)
library(cluster)
library(dplyr)

t_start <- Sys.time()
DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
RESULTS_FILE <- file.path(DATADIR, "silhouette_sweep_P18excl.csv")

cat("############################################################\n")
cat("Loading PCA checkpoint (Job 1) -- pure clustering sweep,\n")
cat("does not require cell-type annotation\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_pca.rds"))
cat("Cells loaded:", ncol(seurat), "\n")
cat("PCA dims available:", ncol(Embeddings(seurat, "pca")), "\n\n")

# ── Fixed subsample, same across all combinations for comparability
# (matches the original P18-included sweep's methodology) ─────
n_sub <- 8000
sub_cells <- sample(colnames(seurat), n_sub)
pca_emb <- Embeddings(seurat, "pca")[sub_cells, 1:30]
cat("Fixed subsample of", n_sub, "cells drawn for silhouette computation\n")
cat("Full clustering (FindNeighbors/FindClusters) still uses ALL",
    ncol(seurat), "cells each time -- only silhouette itself is\n")
cat("computed on this fixed subsample, for tractability.\n\n")

cat("Precomputing distance matrix on the fixed subsample...\n")
dist_mat <- dist(pca_emb)
cat("Distance matrix ready.\n\n")

k_values <- c(10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 100)
res_values <- seq(0.1, 1.0, by = 0.1)
combos <- expand.grid(k = k_values, resolution = res_values)
cat("Total combinations to test:", nrow(combos), "\n\n")

# ── Incremental write: append each result immediately, so a
# walltime timeout does not lose completed work (this is exactly
# what happened in the original sweep -- 106/110 completed before
# hitting the 3hr limit; this version is resumable/partial-safe) ──
if (!file.exists(RESULTS_FILE)) {
  write.csv(data.frame(k=integer(), resolution=numeric(),
                       n_clusters=integer(), mean_silhouette=numeric()),
           RESULTS_FILE, row.names = FALSE)
}
already_done <- read.csv(RESULTS_FILE)
cat("Already-completed combinations found on disk:", nrow(already_done), "\n\n")

for (i in seq_len(nrow(combos))) {
  k_val <- combos$k[i]
  res_val <- combos$resolution[i]

  if (nrow(already_done) > 0 &&
      any(already_done$k == k_val & already_done$resolution == res_val)) {
    next  # skip combinations already completed in a prior partial run
  }

  t_iter <- Sys.time()
  seurat <- FindNeighbors(seurat, dims = 1:30, k.param = k_val, verbose = FALSE)
  seurat <- FindClusters(seurat, resolution = res_val, algorithm = 4, verbose = FALSE)

  cluster_assign <- as.integer(Idents(seurat)[sub_cells])
  n_clust <- length(unique(cluster_assign))

  if (n_clust < 2) {
    mean_sil <- NA
  } else {
    sil <- silhouette(cluster_assign, dist_mat)
    mean_sil <- mean(sil[, "sil_width"])
  }

  result_row <- data.frame(k = k_val, resolution = res_val,
                           n_clusters = n_clust, mean_silhouette = mean_sil)
  write.table(result_row, RESULTS_FILE, sep = ",", append = TRUE,
             row.names = FALSE, col.names = FALSE)

  cat(sprintf("[%d/%d] k=%d, res=%.1f -> %d clusters, silhouette=%.4f (%.1f sec)\n",
             i, nrow(combos), k_val, res_val, n_clust, mean_sil,
             as.numeric(difftime(Sys.time(), t_iter, units="secs"))))
}

cat("\n############################################################\n")
cat("SWEEP COMPLETE (or walltime reached -- check row count below\n")
cat("against total combinations; rerun this same script to resume\n")
cat("if incomplete, already-done rows will be skipped)\n")
cat("############################################################\n")
final_results <- read.csv(RESULTS_FILE)
cat("Total combinations completed:", nrow(final_results), "of", nrow(combos), "\n")
cat("Silhouette score range:", round(min(final_results$mean_silhouette, na.rm=TRUE), 4),
    "to", round(max(final_results$mean_silhouette, na.rm=TRUE), 4), "\n")
cat("\nTop 5 highest silhouette scores:\n")
print(head(final_results[order(-final_results$mean_silhouette), ], 5))

cat("\nTotal elapsed:", round(difftime(Sys.time(), t_start, units="mins"), 1), "min\n")
cat("Done:", as.character(Sys.time()), "\n")

EOF

echo "Job finished: $(date)"
