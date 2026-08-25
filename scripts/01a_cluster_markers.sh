#!/bin/bash
#SBATCH --job-name=P18excl_cluster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128GB
#SBATCH --time=03:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/cluster_markers_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/cluster_markers_%J.err
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

t_start <- Sys.time()
DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("Loading checkpoint from Job 1\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_pca.rds"))
cat("Cells loaded:", ncol(seurat), "\n")
cat("P18 confirmed absent:", !("P18" %in% seurat@meta.data$patient), "\n")
cat("PCA dims available:", ncol(Embeddings(seurat, "pca")), "\n\n")

cat("############################################################\n")
cat("STEP 3: UMAP + FindNeighbors + FindClusters (dims=1:30,\n")
cat("k.param=30, resolution=0.3, algorithm=4/Leiden)\n")
cat("############################################################\n")
t3 <- Sys.time()
seurat <- RunUMAP(seurat, dims = 1:30, verbose = TRUE)
seurat <- FindNeighbors(seurat, dims = 1:30, k.param = 30, verbose = TRUE)
seurat <- FindClusters(seurat, resolution = 0.3, algorithm = 4, verbose = TRUE)
cat("Clustering complete. Elapsed:", round(difftime(Sys.time(), t3, units="mins"), 1), "min\n")
cat("Number of clusters found:", length(unique(Idents(seurat))), "\n")
print(table(Idents(seurat)))
cat("\n")

cat("############################################################\n")
cat("STEP 4: PrepSCTFindMarkers + FindAllMarkers\n")
cat("############################################################\n")
t4 <- Sys.time()
seurat <- PrepSCTFindMarkers(seurat, verbose = TRUE)
markers <- FindAllMarkers(seurat, assay = "SCT", only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, verbose = TRUE)
cat("FindAllMarkers complete. Elapsed:", round(difftime(Sys.time(), t4, units="mins"), 1), "min\n")
cat("Total marker rows:", nrow(markers), "\n")
write.csv(markers, file.path(DATADIR, "P18excl_cluster_markers_full.csv"), row.names = FALSE)
cat("Saved P18excl_cluster_markers_full.csv\n\n")

saveRDS(seurat, file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Saved seurat_P18excl_clustered.rds -- checkpoint for Job 2b\n\n")

cat("############################################################\n")
cat("JOB 2a COMPLETE -- clustering and markers only.\n")
cat("Annotation (SingleR) is a SEPARATE job (2b), so a failure\n")
cat("there will not require re-running this step.\n")
cat("############################################################\n")
cat("Total elapsed:", round(difftime(Sys.time(), t_start, units="mins"), 1), "min\n")
cat("Done:", as.character(Sys.time()), "\n")

EOF

echo "Job finished: $(date)"
