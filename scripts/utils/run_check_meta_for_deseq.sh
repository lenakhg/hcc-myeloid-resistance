#!/bin/bash
#SBATCH --job-name=P18_metacheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=00:15:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/metacheck_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/metacheck_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0
Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
library(Seurat)
DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
s <- readRDS(file.path(DATADIR,"seurat_P18excl_combined.rds"))
cat("Cells:", ncol(s), "\n\n")
cat("Metadata columns:\n"); print(colnames(s@meta.data)); cat("\n")
cat("=== patient x recist ===\n"); print(table(s$patient, s$recist))
cat("\n=== patient x timepoint ===\n"); print(table(s$patient, s$timepoint_cell))
cat("\n=== assays present ===\n"); print(Assays(s))
cat("\n=== RNA assay layers (need 'counts' for DESeq) ===\n")
print(Layers(s[["RNA"]]))
cat("\n=== recist values ===\n"); print(table(s$recist))
cat("=== timepoint values ===\n"); print(table(s$timepoint_cell))
EOF
echo "Job finished: $(date)"
