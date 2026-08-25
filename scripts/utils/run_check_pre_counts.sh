#!/bin/bash
#SBATCH --job-name=P18excl_precount
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=00:15:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/precount_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/precount_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
library(Seurat); library(dplyr)
DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
s <- readRDS(file.path(DATADIR, "seurat_P18excl_combined.rds"))
md <- s@meta.data

cat("=== WHOLE-COHORT PRE-TREATMENT: cells per group (all patients) ===\n")
pre <- md[md$timepoint_cell == "Pre-treatment", ]
tb <- sort(table(droplevels(factor(pre$cell_type_combined))), decreasing = TRUE)
print(tb)
cat("\nTotal pre-treatment cells:", nrow(pre), "\n")
cat("Groups BELOW 10 cells (would drop from CellChat):\n")
print(names(tb)[tb < 10])
cat("\nGroups BELOW 25 cells (thin, borderline):\n")
print(names(tb)[tb < 25])

cat("\n=== For reference: WHOLE-COHORT POOLED (pre+post) per group ===\n")
print(sort(table(droplevels(factor(md$cell_type_combined))), decreasing = TRUE))
EOF
echo "Job finished: $(date)"
