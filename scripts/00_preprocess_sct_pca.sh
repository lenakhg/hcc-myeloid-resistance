#!/bin/bash
#SBATCH --job-name=P18excl_sct_pca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128GB
#SBATCH --time=03:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/sct_pca_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/sct_pca_%J.err
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
cat("STEP 0: Load raw object (local, clean copy), exclude P18\n")
cat("BEFORE any processing\n")
cat("############################################################\n")
seurat_raw <- readRDS(file.path(DATADIR, "seurat_tumor_raw.rds"))
cat("Total cells, all patients, before exclusion:", ncol(seurat_raw), "\n")
print(table(seurat_raw@meta.data$patient))

seurat <- subset(seurat_raw, subset = patient != "P18")
cat("\nTotal cells AFTER excluding P18:", ncol(seurat), "\n")
print(table(seurat@meta.data$patient))
cat("P18 confirmed absent:", !("P18" %in% seurat@meta.data$patient), "\n\n")
rm(seurat_raw); gc()

cat("############################################################\n")
cat("DIAGNOSTIC: mitochondrial-PERCENTAGE column, explicitly\n")
cat("validated as bounded 0-100 before use (not the first grep\n")
cat("match blindly -- this exact mistake happened once already)\n")
cat("############################################################\n")
all_cols <- colnames(seurat@meta.data)
cat("All metadata columns present:\n")
print(all_cols)
cat("\n")

if ("percent.mt" %in% all_cols) {
  mx <- max(seurat@meta.data$percent.mt, na.rm = TRUE)
  cat("'percent.mt' exists. Max value:", mx, "\n")
  if (mx > 100) stop("'percent.mt' max=", mx, " > 100 -- not a valid percentage. Stopping.")
  cat("Confirmed valid. Using as-is.\n")
} else if ("pct_counts_mt" %in% all_cols) {
  mx <- max(seurat@meta.data$pct_counts_mt, na.rm = TRUE)
  cat("Found 'pct_counts_mt'. Max value:", mx, "\n")
  if (mx > 100) stop("'pct_counts_mt' max=", mx, " > 100 -- invalid. Stopping rather than proceeding.")
  seurat@meta.data$percent.mt <- seurat@meta.data$pct_counts_mt
  cat("Confirmed valid (bounded 0-100). Assigned pct_counts_mt -> percent.mt.\n")
  cat("Summary:\n"); print(summary(seurat@meta.data$percent.mt))
} else {
  cat("No valid percentage column found. Computing fresh via\n")
  cat("PercentageFeatureSet (pattern='^MT-').\n")
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  print(summary(seurat@meta.data$percent.mt))
}
cat("\n")

cat("############################################################\n")
cat("STEP 1: SCTransform (v2, 5000 features)\n")
cat("############################################################\n")
t1 <- Sys.time()
seurat <- SCTransform(seurat, method = "glmGamPoi", vars.to.regress = "percent.mt",
                      variable.features.n = 5000, vst.flavor = "v2", verbose = TRUE)
cat("SCTransform complete. Elapsed:", round(difftime(Sys.time(), t1, units="mins"), 1), "min\n")
cat("Variable features:", length(VariableFeatures(seurat)), "\n\n")

cat("############################################################\n")
cat("STEP 2: PCA (30 components)\n")
cat("############################################################\n")
t2 <- Sys.time()
seurat <- RunPCA(seurat, npcs = 30, verbose = TRUE)
cat("PCA complete. Elapsed:", round(difftime(Sys.time(), t2, units="mins"), 1), "min\n")
cat("Confirmed npcs:", ncol(Embeddings(seurat, "pca")), "\n\n")

saveRDS(seurat, file.path(DATADIR, "seurat_P18excl_pca.rds"))
cat("Saved seurat_P18excl_pca.rds\n")

cat("\n############################################################\n")
cat("JOB 1 COMPLETE\n")
cat("############################################################\n")
cat("Total elapsed:", round(difftime(Sys.time(), t_start, units="mins"), 1), "min\n")
cat("Done:", as.character(Sys.time()), "\n")

EOF

echo "Job finished: $(date)"
