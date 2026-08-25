#!/bin/bash
#SBATCH --job-name=P18excl_final
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/final_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/final_%J.err
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
library(ggplot2)
library(patchwork)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
FIGDIR  <- "/ibex/user/alghamlk/HCC_P18_excluded/figures"

cat("############################################################\n")
cat("FINAL ANNOTATION -- pure CellTypist majority_voting, as-is.\n")
cat("Labels + save + 4-panel UMAP in one job. Single source of\n")
cat("truth downstream: seurat_P18excl_final.rds, group by cell_type.\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Cells:", ncol(seurat), " Clusters:", length(unique(Idents(seurat))), "\n\n")

# ── stored UMAP, do not recompute ───────────────────────────────
if (!("umap" %in% Reductions(seurat))) {
  cat("NOTE: no stored UMAP -- computing once (becomes committed).\n")
  seurat <- RunUMAP(seurat, dims = 1:30, seed.use = 42, verbose = FALSE)
} else {
  cat("Using STORED UMAP embedding from clustering checkpoint.\n")
}
cat("\n")

# ── CellTypist majority_voting -> per-cluster label ─────────────
preds <- read.csv(file.path(DATADIR, "celltypist_predictions.csv"))
cluster_map <- data.frame(barcode = colnames(seurat),
                          cluster  = as.character(Idents(seurat)),
                          stringsAsFactors = FALSE)
merged <- merge(cluster_map, preds, by = "barcode")
cat("Merged:", nrow(merged), "of", ncol(seurat),
    ifelse(nrow(merged) == ncol(seurat), "(FULL MATCH)\n", "(MISMATCH!)\n"))

cluster_majority <- merged %>%
  group_by(cluster) %>%
  count(majority_voting, sort = TRUE) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.integer(cluster))

cat("\nPer-cluster committed CellTypist label:\n")
print(as.data.frame(cluster_majority), row.names = FALSE)

lab <- setNames(cluster_majority$majority_voting, cluster_majority$cluster)
cell_lab <- lab[as.character(Idents(seurat))]
cell_lab[is.na(cell_lab)] <- "Unassigned"
cell_lab <- unname(cell_lab)                       # FIX: strip cluster-id names
seurat$cell_type <- factor(cell_lab, levels = sort(unique(cell_lab)))

cat("\nCommitted cell_type distribution:\n")
print(sort(table(seurat$cell_type), decreasing = TRUE))
cat("\n")

# ── guarded cross-check vs author major_cluster ─────────────────
if ("major_cluster" %in% colnames(seurat@meta.data)) {
  ct <- tryCatch(
    table(cell_type = seurat$cell_type,
          author    = as.character(seurat@meta.data$major_cluster)),
    error = function(e) { cat("cross-tab skipped:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ct)) { cat("cell_type vs author major_cluster:\n"); print(ct); cat("\n") }
} else {
  cat("No major_cluster column -- cross-tab skipped.\n\n")
}

# ── save committed object (BEFORE plotting, so a plot error can't
# cost you the object) ──────────────────────────────────────────
saveRDS(seurat, file.path(DATADIR, "seurat_P18excl_final.rds"))
cat("Saved seurat_P18excl_final.rds  <-- downstream loads THIS\n")
write.csv(data.frame(barcode = colnames(seurat),
                     cluster = as.character(Idents(seurat)),
                     cell_type = as.character(seurat$cell_type)),
          file.path(DATADIR, "P18excl_final_cell_type.csv"), row.names = FALSE)
cat("Saved P18excl_final_cell_type.csv\n\n")

# ════════════════════════════════════════════════════════════════
# 4-panel UMAP by RECIST x timepoint, colored by cell_type.
# Clean layout: 2x2 grid, one shared legend collected on the right.
# ════════════════════════════════════════════════════════════════
cat("############################################################\n")
cat("4-panel UMAP colored by committed cell_type\n")
cat("############################################################\n")
seurat$panel_group <- paste(seurat$recist, seurat$timepoint_cell)
panels <- c("PR - Pre-treatment"  = "PR Pre-treatment",
            "PR - Post-treatment" = "PR Post-treatment",
            "PD - Pre-treatment"  = "PD Pre-treatment",
            "PD - Post-treatment" = "PD Post-treatment")

umap_df <- as.data.frame(Embeddings(seurat, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$panel_group <- seurat$panel_group
umap_df$cell_type   <- seurat$cell_type

all_types <- levels(seurat$cell_type)
palette <- setNames(scales::hue_pal()(length(all_types)), all_types)

make_panel <- function(title, show_legend = FALSE) {
  sub_df <- umap_df %>% filter(panel_group == panels[[title]])
  cat(sprintf("%-22s cells: %d\n", title, nrow(sub_df)))
  p <- ggplot(sub_df, aes(UMAP_1, UMAP_2, color = cell_type)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_manual(values = palette, drop = FALSE, name = "Cell type") +
    theme_minimal() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme(plot.title = element_text(size = 12, face = "bold"))
  if (show_legend) p + guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
  else p + theme(legend.position = "none")
}

p1 <- make_panel("PR - Pre-treatment")
p2 <- make_panel("PR - Post-treatment", show_legend = TRUE)  # this one carries the legend
p3 <- make_panel("PD - Pre-treatment")
p4 <- make_panel("PD - Post-treatment")
cat("\n")

combined <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Whole-Tumor UMAP by Response and Timepoint (P18 Excluded, CellTypist majority_voting)",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  legend.position = "right"))

ggsave(file.path(FIGDIR, "umap_P18excl_4panel_FINAL.png"),
       combined, width = 15, height = 11, dpi = 300)
cat("Saved umap_P18excl_4panel_FINAL.png\n\n")

cat("############################################################\n")
cat("DONE. Annotation committed; figure written.\n")
cat("############################################################\n")
cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
