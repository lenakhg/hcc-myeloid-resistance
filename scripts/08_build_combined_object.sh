#!/bin/bash
#SBATCH --job-name=P18excl_combine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=00:45:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/combine_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/combine_%J.err
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

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("Build COMBINED cell_type on the full 70k object:\n")
cat("myeloid subsets (from subclustering) spliced onto the whole\n")
cat("tumor; non-myeloid = CellTypist per-cluster labels.\n")
cat("Nothing dropped; contamination kept as its own group.\n")
cat("############################################################\n")

# ── full object with committed CellTypist cell_type ─────────────
full <- readRDS(file.path(DATADIR, "seurat_P18excl_final.rds"))
cat("Full object cells:", ncol(full), "\n")
cat("CellTypist cell_type levels:\n"); print(table(full$cell_type))

# ── myeloid object with per-cell myeloid_label (barcodes match) ─
mye <- readRDS(file.path(DATADIR, "seurat_P18excl_myeloid_res04_named.rds"))
cat("\nMyeloid object cells:", ncol(mye), "\n")
cat("Myeloid labels:\n"); print(table(mye$myeloid_label))

mye_lab <- setNames(as.character(mye$myeloid_label), colnames(mye))

# ── Splice: start from CellTypist cell_type; for any cell present
# in the myeloid object, OVERWRITE with its myeloid subset label.
# (These are exactly the old cluster-2 cells; everything else
# keeps its CellTypist label.) ──────────────────────────────────
combined <- as.character(full$cell_type)
names(combined) <- colnames(full)
hit <- intersect(names(mye_lab), names(combined))
cat("\nCells receiving myeloid subset label:", length(hit),
    "(expected ~", ncol(mye), ")\n")
combined[hit] <- mye_lab[hit]

full$cell_type_combined <- factor(combined)
cat("\nCOMBINED cell_type_combined distribution:\n")
print(sort(table(full$cell_type_combined), decreasing = TRUE))
cat("\nTotal groups:", nlevels(full$cell_type_combined), "\n\n")

# sanity: the old myeloid CellTypist label (DC2) should now be gone
cat("Any cells still labelled 'DC2' (should be 0 -- all became",
    "myeloid subsets):", sum(full$cell_type_combined == "DC2"), "\n\n")

saveRDS(full, file.path(DATADIR, "seurat_P18excl_combined.rds"))
cat("Saved seurat_P18excl_combined.rds (cell_type_combined) <- CellChat input\n")
write.csv(data.frame(barcode = colnames(full),
                     cell_type_combined = as.character(full$cell_type_combined)),
          file.path(DATADIR, "P18excl_combined_labels.csv"), row.names = FALSE)
cat("Saved P18excl_combined_labels.csv\n")
cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
