#!/bin/bash
#SBATCH --job-name=P18excl_sctype
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/sctype_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/sctype_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'

.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)

if (!requireNamespace("HGNChelper", quietly = TRUE))
  install.packages("HGNChelper", repos = "https://cloud.r-project.org")
if (!requireNamespace("openxlsx", quietly = TRUE))
  install.packages("openxlsx", repos = "https://cloud.r-project.org")

library(Seurat)
library(dplyr)
library(HGNChelper)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("scType -- CORRECTED INPUT: scale full RNA assay over the\n")
cat("marker-DB genes, not SCT scale.data (which held only the\n")
cat("5000 variable features and produced spurious labels last run)\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Cells loaded:", ncol(seurat), "\n")
cat("Clusters present:", length(unique(Idents(seurat))), "\n\n")

cat("Sourcing scType functions from the official repository...\n")
tryCatch({
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  cat("scType functions loaded successfully.\n\n")
}, error = function(e) {
  stop("FAILED to source scType from GitHub: ", conditionMessage(e),
       " -- likely network/proxy restriction.")
})

cat("Loading scType's Immune-system marker database...\n")
db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
gs_list <- tryCatch({
  gene_sets_prepare(db_url, "Immune system")
}, error = function(e) stop("FAILED to load ScTypeDB: ", conditionMessage(e)))
cat("Marker DB loaded. Cell types available:", length(gs_list$gs_positive), "\n\n")

# ── KEY FIX ────────────────────────────────────────────────────
# scType expects a SCALED matrix over the genes it will score.
# SCT scale.data only contains the 5000 variable features, so many
# marker-DB genes were absent last run -> partial scores -> junk
# labels. Here we normalize + scale the FULL RNA assay, then check
# marker-gene coverage explicitly before scoring.
# ───────────────────────────────────────────────────────────────
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 1e4, verbose = FALSE)

# genes referenced by the marker DB (positive + negative sets)
db_genes <- unique(unlist(c(gs_list$gs_positive, gs_list$gs_negative)))
db_genes <- db_genes[!is.na(db_genes)]
rna_genes <- rownames(seurat)
covered <- intersect(db_genes, rna_genes)
cat("Marker-DB genes total:", length(db_genes),
    "| present in RNA assay:", length(covered),
    sprintf("(%.1f%%)\n", 100*length(covered)/length(db_genes)))
cat("This coverage should be HIGH (>90%). Last run's SCT scale.data\n")
cat("would have covered far fewer -- that was the bug.\n\n")

# scale over union of variable features and all marker-DB genes,
# so every scored gene has a scaled value available
scale_features <- unique(c(VariableFeatures(seurat, assay = "RNA"), covered))
scale_features <- intersect(scale_features, rna_genes)
cat("Scaling RNA assay over", length(scale_features), "features",
    "(variable features + all covered marker genes)...\n")
seurat <- ScaleData(seurat, features = scale_features, verbose = FALSE)

expr_matrix <- as.matrix(GetAssayData(seurat, assay = "RNA", layer = "scale.data"))
cat("Scaled matrix for scoring:", paste(dim(expr_matrix), collapse=" x "), "\n\n")

cat("Scoring cells against marker DB (sctype_score, scaled=TRUE)...\n")
es_max <- sctype_score(scRNAseqData = expr_matrix, scaled = TRUE,
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# ── Per-cluster assignment, with 2nd-place gap for confidence ──
cat("Assigning per-cluster identity (top cumulative score),\n")
cat("reporting runner-up + margin so low-confidence calls show...\n\n")
cluster_results <- lapply(unique(Idents(seurat)), function(cl) {
  cells <- names(Idents(seurat))[Idents(seurat) == cl]
  scores <- sort(rowSums(es_max[, cells, drop = FALSE]), decreasing = TRUE)
  ncell <- length(cells)
  # scType's own low-confidence rule: top score < ncell/4 -> "Unknown"
  conf <- ifelse(scores[1] < ncell/4, "LOW/Unknown", "ok")
  data.frame(cluster = cl,
             type = names(scores)[1], score = round(scores[1],1),
             runner_up = names(scores)[2], runner_score = round(scores[2],1),
             margin = round(scores[1]-scores[2],1),
             n_cells = ncell, confidence = conf)
}) %>% bind_rows() %>% arrange(as.numeric(as.character(cluster)))

cat("Per-cluster scType assignment (corrected):\n")
print(cluster_results, row.names = FALSE)
cat("\n")

write.csv(cluster_results,
          file.path(DATADIR, "P18excl_cluster_annotations_scType_corrected.csv"),
          row.names = FALSE)
cat("Saved P18excl_cluster_annotations_scType_corrected.csv\n\n")

# ── Cluster-2 myeloid focus: show its full myeloid score ranking ─
cat("############################################################\n")
cat("CLUSTER-2 detail: full ranked scores (myeloid adjudication)\n")
cat("############################################################\n")
c2_cells <- names(Idents(seurat))[Idents(seurat) == "2"]
c2_scores <- sort(rowSums(es_max[, c2_cells, drop = FALSE]), decreasing = TRUE)
print(round(head(c2_scores, 8), 1))
cat("\n")

# ── Cross-validation, myeloid grep FIXED to exclude 'dendritic'
# false-positives and to be explicit about what it counts ──────
cluster_labels <- setNames(cluster_results$type, cluster_results$cluster)
seurat@meta.data$auto_annotation_sctype <- cluster_labels[as.character(Idents(seurat))]
myeloid_hits <- grep("monocyte|macrophage", cluster_labels,
                     value = TRUE, ignore.case = TRUE)
cat("Cluster label(s) matched to monocyte/macrophage:",
    paste(unique(myeloid_hits), collapse=", "), "\n")
if (length(myeloid_hits) > 0) {
  our_cells <- sum(seurat@meta.data$auto_annotation_sctype %in% myeloid_hits, na.rm=TRUE)
  overlap <- sum(seurat@meta.data$auto_annotation_sctype %in% myeloid_hits &
                 seurat@meta.data$major_cluster == "Myeloid", na.rm=TRUE)
  author_total <- sum(seurat@meta.data$major_cluster == "Myeloid")
  cat(sprintf("Precision: %.1f%%  Recall: %.1f%%\n",
              100*overlap/our_cells, 100*overlap/author_total))
}

cat("\nDone:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
