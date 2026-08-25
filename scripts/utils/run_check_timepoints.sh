#!/bin/bash
#SBATCH --job-name=P18excl_tpcheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=00:20:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/tpcheck_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/tpcheck_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
library(Seurat); library(dplyr)
DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
s <- readRDS(file.path(DATADIR, "seurat_P18excl_combined.rds"))
md <- s@meta.data

cat("=== Patients x timepoint (does everyone have pre AND post?) ===\n")
print(table(patient = md$patient, timepoint = md$timepoint_cell))

cat("\n=== Patient x response (PR/PD) ===\n")
print(table(patient = md$patient, recist = md$recist))

cat("\n=== PRE-TREATMENT ONLY: cells per group per response arm ===\n")
pre <- md[md$timepoint_cell == "Pre-treatment", ]
print(table(cell_type = pre$cell_type_combined, recist = pre$recist))

cat("\n=== For reference: POOLED (pre+post) cells per group per arm ===\n")
print(table(cell_type = md$cell_type_combined, recist = md$recist))
EOF
echo "Job finished: $(date)"
