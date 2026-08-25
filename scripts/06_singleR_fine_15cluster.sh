#!/bin/bash
#SBATCH --job-name=P18excl_sr15
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/sr15_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/sr15_%J.err
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
cat("COMMITTED FIGURE: 15 Leiden clusters, SingleR-fine (HPCA)\n")
cat("names, labeled 'C{cluster}: {fine label}' so repeated names\n")
cat("stay as DISTINCT colors. 15-class publication palette.\n")
cat("NOTE: several clusters share the fine label 'CD4+ effector\n")
cat("memory' (incl. the Treg cluster) -- this is HPCA-fine's known\n")
cat("collapse and is retained here by the user's explicit choice.\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
cat("Cells:", ncol(seurat), " Clusters:", length(unique(Idents(seurat))), "\n\n")

# stored UMAP, do not recompute
if (!("umap" %in% Reductions(seurat))) {
  seurat <- RunUMAP(seurat, dims = 1:30, seed.use = 42, verbose = FALSE)
  cat("Computed UMAP (none stored).\n")
} else cat("Using stored UMAP embedding.\n")
cat("\n")

# ── Read SingleR-fine per-cluster labels from disk ──────────────
sr <- read.csv(file.path(DATADIR, "P18excl_cluster_annotations_SingleR_fine.csv"),
               stringsAsFactors = FALSE)
cat("SingleR-fine CSV columns:", paste(colnames(sr), collapse=", "), "\n")
print(sr)
cat("\n")

# Robustly identify the cluster column and the label column.
# (script adapts to common column namings; edit here if needed.)
clust_col <- intersect(c("cluster","Cluster","clusters"), colnames(sr))[1]
lab_col   <- intersect(c("auto_annotation","label","labels","pruned_label",
                         "pruned.labels","annotation"), colnames(sr))[1]
if (is.na(clust_col) || is.na(lab_col))
  stop("Could not auto-detect cluster/label columns. Columns present: ",
       paste(colnames(sr), collapse=", "),
       " -- edit clust_col/lab_col in the script.")
cat("Using cluster column:", clust_col, "| label column:", lab_col, "\n\n")

sr[[clust_col]] <- as.character(sr[[clust_col]])

# ── Build 'C{n}: {label}' display name per cluster ──────────────
# order by numeric cluster id so legend reads C1..C15
sr <- sr[order(as.integer(sr[[clust_col]])), ]
sr$display <- sprintf("C%s: %s", sr[[clust_col]], sr[[lab_col]])
label_map <- setNames(sr$display, sr[[clust_col]])

# apply to cells (positional; strip names to avoid barcode-match error)
cell_lab <- label_map[as.character(Idents(seurat))]
cell_lab[is.na(cell_lab)] <- "unlabeled"
cell_lab <- unname(cell_lab)
# factor levels ordered C1..C15
lev <- sr$display[order(as.integer(sr[[clust_col]]))]
seurat$cluster_fine <- factor(cell_lab, levels = lev)

cat("cluster_fine levels (legend order):\n"); print(levels(seurat$cluster_fine))
cat("\nDistribution:\n"); print(table(seurat$cluster_fine)); cat("\n")

saveRDS(seurat, file.path(DATADIR, "seurat_P18excl_final_SRfine15.rds"))
cat("Saved seurat_P18excl_final_SRfine15.rds (committed: 15-cluster SR-fine)\n\n")

# ════════════════════════════════════════════════════════════════
# Publication palette for up to 15 categories.
# Uses a curated, maximally-distinct qualitative set (Okabe-Ito
# extended + Polychrome-style additions), not hue_pal (which makes
# 15 similar hues). Colorblind-friendlier than default ggplot.
# ════════════════════════════════════════════════════════════════
pub15 <- c(
  "#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00",
  "#CC79A7","#000000","#999999","#AD7700","#1C91D4","#007756",
  "#D5C711","#005685","#A04700"
)
n_lev <- nlevels(seurat$cluster_fine)
palette <- setNames(pub15[seq_len(n_lev)], levels(seurat$cluster_fine))

# ── 4 panels by RECIST x timepoint ──────────────────────────────
seurat$panel_group <- paste(seurat$recist, seurat$timepoint_cell)
panels <- c("PR - Pre-treatment"  = "PR Pre-treatment",
            "PR - Post-treatment" = "PR Post-treatment",
            "PD - Pre-treatment"  = "PD Pre-treatment",
            "PD - Post-treatment" = "PD Post-treatment")

umap_df <- as.data.frame(Embeddings(seurat, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1","UMAP_2")
umap_df$panel_group  <- seurat$panel_group
umap_df$cluster_fine <- seurat$cluster_fine

make_panel <- function(title, legend = FALSE) {
  sub_df <- umap_df %>% filter(panel_group == panels[[title]])
  cat(sprintf("%-22s cells: %d\n", title, nrow(sub_df)))
  p <- ggplot(sub_df, aes(UMAP_1, UMAP_2, color = cluster_fine)) +
    geom_point(size = 0.35, alpha = 0.75) +
    scale_color_manual(values = palette, drop = FALSE, name = "Cluster (SingleR-fine)") +
    theme_bw(base_size = 12) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme(plot.title = element_text(size = 13, face = "bold"),
          panel.grid.minor = element_blank())
  if (legend) p + guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
  else p + theme(legend.position = "none")
}

p1 <- make_panel("PR - Pre-treatment")
p2 <- make_panel("PR - Post-treatment", legend = TRUE)
p3 <- make_panel("PD - Pre-treatment")
p4 <- make_panel("PD - Post-treatment")
cat("\n")

combined <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Whole-Tumor UMAP by Response and Timepoint (P18 Excluded, 15 clusters, SingleR-fine)",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  legend.position = "right"))

ggsave(file.path(FIGDIR, "umap_P18excl_4panel_SRfine15_FINAL.png"),
       combined, width = 16, height = 11, dpi = 300)
cat("Saved umap_P18excl_4panel_SRfine15_FINAL.png\n\n")

cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
