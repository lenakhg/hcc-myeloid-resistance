#!/bin/bash
#SBATCH --job-name=P18excl_scina
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/scina_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/scina_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'

.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)

if (!requireNamespace("SCINA", quietly = TRUE))
  install.packages("SCINA", repos = "https://cloud.r-project.org")

library(Seurat)
library(dplyr)
library(SCINA)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("SCINA -- whole-object LINEAGE-level cross-check (Option A).\n")
cat("Purpose: independent semi-supervised confirmation of the\n")
cat("major-lineage backbone the other methods agreed on. NOT\n")
cat("expected to yield new biology; documentation/reviewer value.\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Cells loaded:", ncol(seurat), "\n")
cat("Clusters present:", length(unique(Idents(seurat))), "\n\n")

# ── RNA assay, log-normalized (SCINA expects normalized data,
# NOT SCT-corrected -- same convention as the marker/CellChat path) ─
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 1e4, verbose = FALSE)
all_genes <- rownames(seurat)

# ── LINEAGE-level signatures (major immune types only). Kept
# deliberately canonical and non-overlapping so lineages separate
# cleanly. This is a lineage cross-check, NOT subset resolution. ─
signatures <- list(
  Macrophage      = c("C1QC","C1QA","C1QB","APOE","CD68","CD163","SELENOP","MRC1"),
  Monocyte        = c("FCN1","S100A8","S100A9","VCAN","LYZ"),
  DC              = c("FCER1A","CD1C","CLEC10A","CLEC9A"),
  pDC             = c("LILRA4","GZMB","IL3RA","CLEC4C"),
  CD8_T           = c("CD3D","CD3E","CD8A","GZMK","NKG7"),
  CD4_T           = c("CD3D","CD3E","CD4","IL7R","CCR7"),
  Treg            = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  NK              = c("KLRD1","NKG7","GNLY","KLRF1","NCAM1"),
  B_cell          = c("MS4A1","CD79A","CD79B","CD19"),
  Plasma          = c("MZB1","XBP1","SDC1","IGHG1")
)

cat("############################################################\n")
cat("Validating signature genes present BEFORE running SCINA\n")
cat("(SCINA errors on genes absent from the matrix; prune first)\n")
cat("############################################################\n")
signatures_present <- list()
for (nm in names(signatures)) {
  present <- intersect(signatures[[nm]], all_genes)
  missing <- setdiff(signatures[[nm]], all_genes)
  cat(sprintf("%-12s %d/%d present", nm, length(present), length(signatures[[nm]])))
  if (length(missing) > 0) cat("  MISSING:", paste(missing, collapse=","))
  cat("\n")
  if (length(present) >= 2) signatures_present[[nm]] <- present
  else cat("   -> dropped (<2 genes present)\n")
}
cat("\nSignatures used:", length(signatures_present), "of", length(signatures), "\n\n")

# ── Restrict expression matrix to signature genes for speed/RAM ──
sig_genes <- unique(unlist(signatures_present))
expr <- as.matrix(GetAssayData(seurat, assay = "RNA", layer = "data")[sig_genes, ])
cat("Expression matrix for SCINA:", paste(dim(expr), collapse=" x "), "\n\n")

cat("############################################################\n")
cat("Running SCINA (allow_unknown=TRUE so genuinely ambiguous\n")
cat("cells are flagged 'unknown' rather than force-assigned)\n")
cat("############################################################\n")
results <- SCINA(exp = expr, signatures = signatures_present,
                 max_iter = 100, convergence_n = 10,
                 convergence_rate = 0.999, sensitivity_cutoff = 0.9,
                 rm_overlap = FALSE, allow_unknown = TRUE)

seurat@meta.data$scina_label <- results$cell_labels
cat("SCINA per-cell label distribution:\n")
print(sort(table(results$cell_labels), decreasing = TRUE))
cat("\n")

unknown_frac <- mean(results$cell_labels == "unknown")
cat(sprintf("Fraction 'unknown': %.1f%%\n\n", 100*unknown_frac))

cat("############################################################\n")
cat("Per-CLUSTER majority SCINA label (this is the cross-check)\n")
cat("############################################################\n")
md <- seurat@meta.data
md$cluster <- as.character(Idents(seurat))
cluster_majority <- md %>%
  group_by(cluster) %>%
  count(scina_label) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(cluster_n = as.integer(table(md$cluster)[cluster]),
         majority_frac = round(n / cluster_n, 3)) %>%
  arrange(as.integer(cluster))
cat("Per-cluster dominant SCINA lineage (with majority fraction):\n")
print(as.data.frame(cluster_majority), row.names = FALSE)
cat("\n")
write.csv(cluster_majority,
          file.path(DATADIR, "P18excl_cluster_annotations_SCINA.csv"),
          row.names = FALSE)
cat("Saved P18excl_cluster_annotations_SCINA.csv\n\n")

cat("############################################################\n")
cat("CROSS-TAB: SCINA per-cell label vs author major_cluster\n")
cat("(independent concordance check on the lineage backbone)\n")
cat("############################################################\n")
print(table(SCINA = md$scina_label, Author = md$major_cluster))
cat("\n")

cat("############################################################\n")
cat("Cluster-2 confirmation (the myeloid cluster of interest)\n")
cat("############################################################\n")
c2 <- md$scina_label[md$cluster == "2"]
cat("SCINA label breakdown within cluster 2:\n")
print(sort(table(c2), decreasing = TRUE))
cat(sprintf("\n-> Cluster 2 dominant SCINA lineage: %s\n",
            names(which.max(table(c2)))))

cat("\nDone:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
