#!/bin/bash
#SBATCH --job-name=pr_prepost_P18excl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/pr_prepost_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/pr_prepost_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({library(Seurat); library(DESeq2); library(dplyr); library(Matrix)})

cat("Loading seurat_P18excl_combined.rds:", as.character(Sys.time()), "\n")
seurat <- readRDS("/ibex/user/alghamlk/HCC_P18_excluded/data/seurat_P18excl_combined.rds")

# PR patients with BOTH timepoints
pr_paired_patients <- c("P11", "P27", "P5")
seurat_pr <- subset(seurat,
    subset = patient %in% pr_paired_patients &
             timepoint_cell %in% c("Pre-treatment", "Post-treatment"))
cat("PR paired-only cells:", ncol(seurat_pr), "\n")
print(table(seurat_pr@meta.data$patient, seurat_pr@meta.data$timepoint_cell))
cat("\n")

# Pseudobulk per patient per timepoint
raw_counts <- GetAssayData(seurat_pr, assay = "RNA", layer = "counts")
sample_groups <- seurat_pr@meta.data %>%
  mutate(sample_id = paste(patient, timepoint_cell, sep = "_")) %>%
  distinct(patient, timepoint_cell, sample_id)

pseudobulk_mat <- sapply(sample_groups$sample_id, function(sid) {
  tp  <- ifelse(grepl("Pre-treatment", sid), "Pre-treatment", "Post-treatment")
  pat <- sub(paste0("_", tp), "", sid)
  cells_s <- rownames(seurat_pr@meta.data)[
    seurat_pr@meta.data$patient == pat &
    seurat_pr@meta.data$timepoint_cell == tp ]
  Matrix::rowSums(raw_counts[, cells_s, drop = FALSE])
})
colnames(pseudobulk_mat) <- sample_groups$sample_id
cat("Pseudobulk dims:", dim(pseudobulk_mat), "\n\n")

coldata <- sample_groups %>%
  mutate(timepoint = factor(timepoint_cell, levels=c("Pre-treatment","Post-treatment")),
         patient = factor(patient)) %>%
  select(sample_id, patient, timepoint) %>% as.data.frame()
rownames(coldata) <- coldata$sample_id; coldata$sample_id <- NULL
coldata <- coldata[colnames(pseudobulk_mat), ]
cat("Paired design:\n"); print(coldata)

dds <- DESeqDataSetFromMatrix(round(pseudobulk_mat), coldata, design = ~ patient + timepoint)
keep <- rowSums(counts(dds) >= 5) >= 2
dds <- dds[keep, ]
cat("\nGenes after filtering:", nrow(dds), "\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("timepoint","Post-treatment","Pre-treatment"))
res <- res[order(res$pvalue), ]

write.csv(as.data.frame(res),
  "/ibex/user/alghamlk/HCC_P18_excluded/data/pseudobulk_deseq2_PR_PrePost_P18excl.csv", row.names=TRUE)
cat("\nSig genes (padj<0.05):", sum(!is.na(res$padj) & res$padj<0.05), "\n\n")

# ===== Cross-check against P18-EXCLUDED baseline signature =====
cat("=== Overlap with baseline P18-excluded PR-vs-PD signature ===\n")
baseline <- read.csv("/ibex/user/alghamlk/HCC_P18_excluded/data/pseudobulk_deseq2_PRvsPD_P18excl.csv", row.names=1)

pd_up <- rownames(baseline)[!is.na(baseline$padj) & baseline$padj<0.05 & baseline$log2FoldChange<0]
pr_up <- rownames(baseline)[!is.na(baseline$padj) & baseline$padj<0.05 & baseline$log2FoldChange>0]
cat("Baseline: genes UP in PD:", length(pd_up), " | UP in PR:", length(pr_up), "\n")

post_down <- rownames(res)[!is.na(res$padj) & res$padj<0.05 & res$log2FoldChange<0]
post_up   <- rownames(res)[!is.na(res$padj) & res$padj<0.05 & res$log2FoldChange>0]
cat("On-treatment: genes DOWN post-tx:", length(post_down), " | UP post-tx:", length(post_up), "\n\n")

cat("KEY: PD-resistance genes that DECREASE with successful treatment:\n")
ov1 <- intersect(pd_up, post_down)
cat("  Overlap count:", length(ov1), "\n"); if(length(ov1)>0) print(ov1)

cat("\nOther overlaps (for full transparency):\n")
cat("  PD-up genes that INCREASE post-tx:", length(intersect(pd_up, post_up)), "\n")
cat("  PR-up genes that overlap post-tx-up:", length(intersect(pr_up, post_up)), "\n")

# Is CD14 in the on-treatment sig set, and which direction?
cat("\nCD14 check:\n")
if ("CD14" %in% rownames(res)) {
  r <- res["CD14",]
  cat(sprintf("  CD14 post-vs-pre: log2FC=%+.2f padj=%.2e\n", r$log2FoldChange, ifelse(is.na(r$padj),NA,r$padj)))
} else cat("  CD14 not in filtered on-treatment set\n")

cat("\nDone:", as.character(Sys.time()), "\n")
EOF
echo "Job finished: $(date)"
