#!/bin/bash
#SBATCH --job-name=P18excl_azimuth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=02:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/azimuth_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/azimuth_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)

if (!requireNamespace("Azimuth", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("satijalab/azimuth", quiet = TRUE)
}

library(Seurat)
library(dplyr)
library(Azimuth)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("Azimuth (Seurat-native reference mapping), PBMC reference\n")
cat("(multi-level: celltype.l1 broad, l2 medium, l3 finest --\n")
cat("includes CD4 Naive/TCM/TEM, CD8 Naive/TCM/TEM, Treg, MAIT,\n")
cat("gdT, classical/non-classical monocytes, and more)\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Cells loaded:", ncol(seurat), "\n")
cat("Clusters present:", length(unique(Idents(seurat))), "\n\n")

cat("Running RunAzimuth (per-cell reference mapping, requires\n")
cat("network access to download the PBMC reference -- flagging\n")
cat("as a potential failure point, same category of risk as the\n")
cat("earlier celldex download)...\n")

seurat <- tryCatch({
  RunAzimuth(seurat, reference = "pbmcref")
}, error = function(e) {
  cat("FAILED:", conditionMessage(e), "\n")
  stop("Azimuth reference mapping failed -- likely network/proxy ",
       "restriction, or a package version mismatch. The clustering ",
       "checkpoint (seurat_P18excl_clustered.rds) is untouched and ",
       "does not need to be regenerated to retry this.")
})

cat("\nAzimuth prediction columns added:\n")
azimuth_cols <- grep("predicted\\.celltype", colnames(seurat@meta.data), value = TRUE)
print(azimuth_cols)
cat("\n")

# Summarize per-cluster (majority vote at each resolution level),
# for direct comparison with the SingleR-based results
for (col in azimuth_cols) {
  cat("=== Per-cluster majority label:", col, "===\n")
  summary_df <- seurat@meta.data %>%
    group_by(seurat_clusters) %>%
    count(.data[[col]]) %>%
    slice_max(n, n = 1) %>%
    ungroup()
  print(summary_df)
  cat("\n")
}

write.csv(seurat@meta.data %>% select(seurat_clusters, all_of(azimuth_cols)),
    file.path(DATADIR, "P18excl_azimuth_percell_labels.csv"), row.names = TRUE)
cat("Saved P18excl_azimuth_percell_labels.csv\n\n")

# Cross-validation against author's Myeloid label, using the l2
# (medium-resolution) prediction
if ("predicted.celltype.l2" %in% colnames(seurat@meta.data)) {
  myeloid_hits <- grep("Mono|DC|Macrophage", seurat@meta.data$predicted.celltype.l2,
                       value = TRUE, ignore.case = TRUE)
  our_cells <- sum(seurat@meta.data$predicted.celltype.l2 %in% unique(myeloid_hits))
  overlap <- sum(seurat@meta.data$predicted.celltype.l2 %in% unique(myeloid_hits) &
                seurat@meta.data$major_cluster == "Myeloid")
  author_total <- sum(seurat@meta.data$major_cluster == "Myeloid")
  cat(sprintf("Azimuth (l2) myeloid precision: %.1f%%\n", 100 * overlap / our_cells))
  cat(sprintf("Azimuth (l2) myeloid recall: %.1f%%\n\n", 100 * overlap / author_total))
}

saveRDS(seurat, file.path(DATADIR, "seurat_P18excl_azimuth.rds"))
cat("Saved seurat_P18excl_azimuth.rds\n")
cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
