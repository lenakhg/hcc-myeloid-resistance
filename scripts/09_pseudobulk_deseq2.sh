#!/bin/bash
#SBATCH --job-name=P18_deseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/deseq_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/deseq_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'

.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

library(Seurat)
library(DESeq2)
library(dplyr)
library(Matrix)

cat("############################################################\n")
cat("DESeq2 on P18-EXCLUDED object.\n")
cat("Design: PRE-TREATMENT only, PR vs PD (PD reference),\n")
cat("whole-tumor PATIENT-level pseudobulk (sum raw RNA counts).\n")
cat("Filter: gene >=10 counts in >=2 of 5 patients (stricter than\n")
cat("the original total>=10; standard small-n pseudobulk filter).\n")
cat("CR already absent from this cohort.\n")
cat("############################################################\n")
cat("Loading seurat_P18excl_combined.rds:", as.character(Sys.time()), "\n")
seurat <- readRDS("/ibex/user/alghamlk/HCC_P18_excluded/data/seurat_P18excl_combined.rds")
cat("Cells:", ncol(seurat), "\n\n")

# ── Restrict to pre-treatment ONLY (predictive biomarker framing) ──
seurat_pre <- subset(seurat, subset = timepoint_cell == "Pre-treatment")
cat("Pre-treatment cells:", ncol(seurat_pre), "\n")

# Only PR and PD — CR already absent from this cohort
seurat_pre <- subset(seurat_pre, subset = recist %in% c("PR", "PD"))
cat("Cells after PR/PD filter:", ncol(seurat_pre), "\n")
cat("Patients:", length(unique(seurat_pre@meta.data$patient)), "\n")
print(table(seurat_pre@meta.data$patient, seurat_pre@meta.data$recist))
cat("\n")

# ── Pseudobulk: sum RAW counts per patient ──────────────────
raw_counts <- GetAssayData(seurat_pre, assay = "RNA", layer = "counts")

patients <- unique(seurat_pre@meta.data$patient)
pseudobulk_mat <- sapply(patients, function(p) {
  cells_p <- rownames(seurat_pre@meta.data)[seurat_pre@meta.data$patient == p]
  Matrix::rowSums(raw_counts[, cells_p, drop = FALSE])
})
colnames(pseudobulk_mat) <- patients

cat("Pseudobulk matrix dimensions:", dim(pseudobulk_mat), "\n")
cat("(genes x patients)\n\n")

# ── Build sample metadata for DESeq2 ────────────────────────
patient_recist_map <- seurat_pre@meta.data %>%
  distinct(patient, recist) %>%
  arrange(match(patient, colnames(pseudobulk_mat)))

coldata <- data.frame(
  row.names = patient_recist_map$patient,
  recist    = factor(patient_recist_map$recist, levels = c("PD", "PR"))
)
cat("Sample design:\n")
print(coldata)

# ── DESeq2 ────────────────────────────────────────────────
cat("\nRunning DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = round(pseudobulk_mat),
                              colData   = coldata,
                              design    = ~ recist)

# Pre-filter: gene must have >=10 counts in >=2 of 5 patients
# (stricter than original total>=10; standard small-n pseudobulk)
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("recist", "PR", "PD"))
res <- res[order(res$pvalue), ]

cat("\nTop 30 differentially expressed genes (PR vs PD):\n")
print(head(as.data.frame(res), 30))

write.csv(as.data.frame(res),
    "/ibex/user/alghamlk/HCC_P18_excluded/data/pseudobulk_deseq2_PRvsPD_P18excl.csv",
    row.names = TRUE)
cat("\nSaved pseudobulk_deseq2_PRvsPD_P18excl.csv\n")

nsig05 <- sum(!is.na(res$padj) & res$padj < 0.05)
nsig10 <- sum(!is.na(res$padj) & res$padj < 0.1)
cat("\nSignificant genes (padj < 0.05):", nsig05, "\n")
cat("Significant genes (padj < 0.1):", nsig10, "\n")

# ── ranked list for GSEA (by Wald stat) ─────────────────────
res_df <- as.data.frame(res); res_df$gene <- rownames(res_df)
rnk <- res_df[!is.na(res_df$stat), c("gene","stat")]
rnk <- rnk[order(-rnk$stat), ]
write.csv(rnk, "/ibex/user/alghamlk/HCC_P18_excluded/data/deseq2_PRvsPD_P18excl_ranked.csv",
          row.names = FALSE)
cat("Saved deseq2_PRvsPD_P18excl_ranked.csv (for GSEA)\n")

# ── FDR-significant gene list for druggability (padj<0.05) ───
sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
sig_genes <- sig_genes[order(sig_genes$padj), ]
write.csv(sig_genes,
    "/ibex/user/alghamlk/HCC_P18_excluded/data/deseq2_PRvsPD_P18excl_SIG.csv",
    row.names = FALSE)
cat("Saved deseq2_PRvsPD_P18excl_SIG.csv (", nrow(sig_genes), "genes for druggability)\n")

cat("\nDone:", as.character(Sys.time()), "\n")

EOF

echo "Job finished: $(date)"
