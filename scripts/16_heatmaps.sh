#!/bin/bash
#SBATCH --job-name=P18_heatmaps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=00:40:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/heatmaps_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/heatmaps_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)
suppressMessages({library(Seurat); library(dplyr); library(pheatmap); library(Matrix)})

DATA_P <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
DATA_V <- "/ibex/user/alghamlk/HCC_project_v2/tumor/data"
FIG    <- "/ibex/user/alghamlk/HCC_P18_excluded/figures"

# ============================================================
# HEATMAP A: 71 signature genes x 5 patients (pre-treatment)
# ============================================================
cat("=== Heatmap A: signature genes x patients ===\n")
s <- readRDS(file.path(DATA_P,"seurat_P18excl_combined.rds"))
# pre-treatment, PR/PD only
keep <- s$timepoint_cell=="Pre-treatment" & s$recist %in% c("PR","PD")
sp <- s[, keep]
sp$recist <- droplevels(factor(sp$recist))

# significant genes from DESeq
d <- read.csv(file.path(DATA_V,"pseudobulk_deseq2_PRvsPD.csv"), row.names=1)
sig <- rownames(d)[!is.na(d$padj) & d$padj < 0.05]
cat("Significant genes:", length(sig), "\n")

# pseudobulk: mean log-normalized expression per patient
DefaultAssay(sp) <- "RNA"
sp <- NormalizeData(sp, verbose=FALSE)
expr <- GetAssayData(sp, assay="RNA", layer="data")
sig <- sig[sig %in% rownames(expr)]
cat("Signature genes present:", length(sig), "\n")

patients <- unique(sp$patient)
pb <- sapply(patients, function(pt){
  cells <- colnames(sp)[sp$patient==pt]
  Matrix::rowMeans(expr[sig, cells, drop=FALSE])
})
colnames(pb) <- patients

# z-score per gene (row)
pbz <- t(scale(t(pb)))
pbz[is.na(pbz)] <- 0

# column annotation: response
resp <- sapply(patients, function(pt) as.character(sp$recist[sp$patient==pt][1]))
ann_col <- data.frame(Response=resp, row.names=patients)
ann_colors <- list(Response=c(PR="#1F6FA5", PD="#B03A2E"))

# order columns by response (PD first, then PR)
ord <- order(factor(resp, levels=c("PD","PR")))
pbz <- pbz[, ord]; ann_col <- ann_col[ord,,drop=FALSE]

png(file.path(FIG,"heatmap_A_signature_patients.png"), width=1600, height=2600, res=200)
pheatmap(pbz,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  annotation_col = ann_col, annotation_colors = ann_colors,
  cluster_cols = FALSE, cluster_rows = TRUE,
  show_rownames = TRUE, fontsize_row = 6, fontsize_col = 12,
  main = "Signature genes (padj<0.05) x patients, pre-treatment (z-scored)",
  border_color = NA)
dev.off()
cat("Saved heatmap_A_signature_patients.png\n\n")

# ============================================================
# HEATMAP B: 9 myeloid subsets x top marker genes
# ============================================================
cat("=== Heatmap B: myeloid subsets x markers ===\n")
mye <- readRDS(file.path(DATA_P,"seurat_P18excl_myeloid_res04_named.rds"))
cat("Myeloid label column check:\n"); print(table(mye$myeloid_label))

mk <- read.csv(file.path(DATA_P,"P18excl_myeloid_subcluster_markers_res0.4_final.csv"))
cat("Marker file columns:", paste(colnames(mk), collapse=", "), "\n")
# top 5 markers per subcluster by avg_log2FC (exclude contamination)
top <- mk %>% filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  group_by(cluster) %>% slice_max(avg_log2FC, n=5) %>% ungroup()
mgenes <- unique(top$gene)
cat("Marker genes:", length(mgenes), "\n")

# exclude contamination cells for the figure
mye2 <- mye[, mye$myeloid_label != "Contamination (T/hepatocyte)"]
mye2$myeloid_label <- droplevels(factor(mye2$myeloid_label))
DefaultAssay(mye2) <- "RNA"
mye2 <- NormalizeData(mye2, verbose=FALSE)
mexpr <- GetAssayData(mye2, assay="RNA", layer="data")
mgenes <- mgenes[mgenes %in% rownames(mexpr)]

# mean expression per subset
subsets <- levels(mye2$myeloid_label)
mmat <- sapply(subsets, function(ss){
  cells <- colnames(mye2)[mye2$myeloid_label==ss]
  Matrix::rowMeans(mexpr[mgenes, cells, drop=FALSE])
})
colnames(mmat) <- subsets
mmatz <- t(scale(t(mmat))); mmatz[is.na(mmatz)] <- 0

png(file.path(FIG,"heatmap_B_myeloid_markers.png"), width=2200, height=2600, res=200)
pheatmap(mmatz,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  cluster_cols = TRUE, cluster_rows = TRUE,
  show_rownames = TRUE, fontsize_row = 7, fontsize_col = 11,
  main = "Myeloid subsets x top markers (z-scored mean expression)",
  border_color = NA)
dev.off()
cat("Saved heatmap_B_myeloid_markers.png\n")

cat("\nDone:", as.character(Sys.time()), "\n")
EOF
echo "Job finished: $(date)"
