#!/bin/bash
#SBATCH --job-name=P18excl_annot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=01:30:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/annot_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/annot_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'

.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)

if (!requireNamespace("SingleR", quietly = TRUE))
  BiocManager::install("SingleR", ask = FALSE, update = FALSE)
if (!requireNamespace("celldex", quietly = TRUE))
  BiocManager::install("celldex", ask = FALSE, update = FALSE)
if (!requireNamespace("scrapper", quietly = TRUE))
  BiocManager::install("scrapper", ask = FALSE, update = FALSE)

library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
library(ggplot2)
library(patchwork)

t_start <- Sys.time()
DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
FIGDIR  <- "/ibex/user/alghamlk/HCC_P18_excluded/figures"
ORIGDIR <- "/ibex/user/alghamlk/HCC_project_v2/tumor/data"  # read-only reference

cat("############################################################\n")
cat("Loading checkpoint from Job 2a\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Cells loaded:", ncol(seurat), "\n")
cat("Clusters present:", length(unique(Idents(seurat))), "\n\n")

cat("############################################################\n")
cat("STEP 5: SingleR annotation (celldex HumanPrimaryCellAtlasData)\n")
cat("NOTE: broad reference -- labels myeloid cells as 'Macrophage'/\n")
cat("'Monocyte', not tissue-specific 'Kupffer cell' terminology.\n")
cat("############################################################\n")
ref <- tryCatch({
  HumanPrimaryCellAtlasData()
}, error = function(e) {
  cat("FAILED to download celldex reference:", conditionMessage(e), "\n")
  stop("celldex download failed -- likely network/proxy restriction. ",
       "Clustering checkpoint (seurat_P18excl_clustered.rds) is safe ",
       "and does not need to be regenerated to retry this step.")
})
cat("Reference loaded. Cell types available:\n")
print(table(ref$label.main))
cat("\n")

query_data <- GetAssayData(seurat, assay = "SCT", layer = "data")
cat("Query data dims:", paste(dim(query_data), collapse=" x "), "\n\n")

singler_results <- SingleR(test = query_data, ref = ref,
                           labels = ref$label.fine,
                           clusters = Idents(seurat))

cat("SingleR results per cluster:\n")
print(data.frame(cluster = rownames(singler_results),
                 label = singler_results$labels,
                 pruned_label = singler_results$pruned.labels))

annotation_df <- data.frame(
  cluster = rownames(singler_results),
  auto_annotation = ifelse(is.na(singler_results$pruned.labels),
                           paste0(singler_results$labels, " (low confidence)"),
                           singler_results$labels)
)
write.csv(annotation_df, file.path(DATADIR, "P18excl_cluster_annotations_SingleR_fine.csv"), row.names = FALSE)
cat("\nSaved P18excl_cluster_annotations_SingleR_fine.csv\n\n")

cluster_labels <- setNames(annotation_df$auto_annotation, annotation_df$cluster)
seurat@meta.data$auto_annotation <- cluster_labels[as.character(Idents(seurat))]

cat("############################################################\n")
cat("STEP 6: Cross-validation against author's major_cluster label\n")
cat("############################################################\n")
myeloid_hits <- grep("Macrophage|Monocyte", cluster_labels, value = TRUE, ignore.case = TRUE)
cat("Cluster label(s) matched to myeloid identity:", paste(unique(myeloid_hits), collapse=", "), "\n")
if (length(myeloid_hits) > 0) {
  our_cluster_cells <- sum(seurat@meta.data$auto_annotation %in% myeloid_hits)
  overlap <- sum(seurat@meta.data$auto_annotation %in% myeloid_hits &
                seurat@meta.data$major_cluster == "Myeloid")
  author_total <- sum(seurat@meta.data$major_cluster == "Myeloid")
  cat(sprintf("Our Myeloid-labeled cells: %d\n", our_cluster_cells))
  cat(sprintf("Author's Myeloid label (this subset): %d cells\n", author_total))
  cat(sprintf("Overlap: %d cells\n", overlap))
  cat(sprintf("Precision: %.1f%%\n", 100 * overlap / our_cluster_cells))
  cat(sprintf("Recall: %.1f%%\n\n", 100 * overlap / author_total))
} else {
  cat("WARNING: no cluster matched to Macrophage/Monocyte identity.\n\n")
}

saveRDS(seurat, file.path(DATADIR, "seurat_P18excl_annotated_fine.rds"))
cat("Saved seurat_P18excl_annotated_fine.rds (final object)\n\n")

cat("############################################################\n")
cat("STEP 7: Reference comparison against ORIGINAL (P18-included,\n")
cat("MSigDB C8-based) annotation -- READ-ONLY, context only.\n")
cat("############################################################\n")
original_annot <- tryCatch(
  read.csv(file.path(ORIGDIR, "cluster_annotations.csv")),
  error = function(e) NULL)
if (!is.null(original_annot)) {
  cat("Original (P18-included, C8) labels:\n"); print(original_annot)
  cat("\nNew (P18-excluded, SingleR) labels:\n"); print(annotation_df)
}
cat("\n")

cat("############################################################\n")
cat("STEP 8: Four-panel UMAP -- PR-Pre, PR-Post, PD-Pre, PD-Post\n")
cat("############################################################\n")
seurat@meta.data$panel_group <- paste(seurat@meta.data$recist, seurat@meta.data$timepoint_cell)
panels <- list(
  "PR - Pre-treatment"  = "PR Pre-treatment",
  "PR - Post-treatment" = "PR Post-treatment",
  "PD - Pre-treatment"  = "PD Pre-treatment",
  "PD - Post-treatment" = "PD Post-treatment"
)

umap_df <- as.data.frame(Embeddings(seurat, "umap"))
umap_df$panel_group <- seurat@meta.data$panel_group
umap_df$auto_annotation <- seurat@meta.data$auto_annotation

plots <- list()
for (title in names(panels)) {
  grp <- panels[[title]]
  sub_df <- umap_df %>% filter(panel_group == grp)
  cat(title, "-- cells:", nrow(sub_df), "\n")
  p <- ggplot(sub_df, aes(x = umap_1, y = umap_2, color = auto_annotation)) +
    geom_point(size = 0.3, alpha = 0.6) +
    theme_minimal() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"))
  plots[[title]] <- p
}
cat("\n")

combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "Whole-Tumor Clustering, P18 Excluded (SingleR-Annotated)",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(FIGDIR, "umap_P18excl_4panel_fine.png"),
       combined, width = 12, height = 10, dpi = 300)
cat("Saved umap_P18excl_4panel.png\n\n")

cat("############################################################\n")
cat("JOB 2b COMPLETE\n")
cat("############################################################\n")
cat("Total elapsed:", round(difftime(Sys.time(), t_start, units="mins"), 1), "min\n")
cat("Done:", as.character(Sys.time()), "\n")

EOF

echo "Job finished: $(date)"
