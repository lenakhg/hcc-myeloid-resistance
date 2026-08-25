#!/bin/bash
#SBATCH --job-name=P18excl_breakdown
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/breakdown_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/breakdown_%J.err
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
cat("PART 1: FULL per-cluster label breakdown (not just the\n")
cat("plurality winner) -- showing the real internal composition\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_clustered.rds"))
preds <- read.csv(file.path(DATADIR, "celltypist_predictions.csv"))
cluster_map <- data.frame(barcode = colnames(seurat), cluster = as.character(Idents(seurat)))
merged <- merge(cluster_map, preds, by = "barcode")
cat("Merged rows:", nrow(merged), "of", ncol(seurat), "\n\n")

for (cl in sort(as.numeric(unique(merged$cluster)))) {
  sub <- merged %>% filter(cluster == as.character(cl))
  cat(sprintf("=== Cluster %d (n=%d total cells) ===\n", cl, nrow(sub)))
  breakdown <- sub %>% count(majority_voting, sort = TRUE) %>%
    mutate(pct = round(100*n/sum(n), 1))
  print(as.data.frame(breakdown))
  cat("\n")
}

write.csv(merged %>% count(cluster, majority_voting, sort = TRUE),
    file.path(DATADIR, "P18excl_full_label_breakdown_per_cluster.csv"), row.names = FALSE)
cat("Saved P18excl_full_label_breakdown_per_cluster.csv (complete,\n")
cat("not collapsed to a single winner per cluster)\n\n")

cat("############################################################\n")
cat("PART 2: 4-panel UMAP (PR-Pre, PR-Post, PD-Pre, PD-Post),\n")
cat("colored by TRUE PER-CELL CellTypist label (majority_voting),\n")
cat("not the cluster-collapsed version -- this correctly shows\n")
cat("heterogeneity within a single Leiden cluster if it exists\n")
cat("############################################################\n")
seurat <- RunUMAP(seurat, dims = 1:30, verbose = FALSE)

label_map <- setNames(merged$majority_voting, merged$barcode)
seurat@meta.data$celltype_percell <- label_map[colnames(seurat)]

seurat@meta.data$panel_group <- paste(seurat@meta.data$recist, seurat@meta.data$timepoint_cell)
panels <- list(
  "PR - Pre-treatment"  = "PR Pre-treatment",
  "PR - Post-treatment" = "PR Post-treatment",
  "PD - Pre-treatment"  = "PD Pre-treatment",
  "PD - Post-treatment" = "PD Post-treatment"
)

umap_df <- as.data.frame(Embeddings(seurat, "umap"))
umap_df$panel_group <- seurat@meta.data$panel_group
umap_df$celltype <- seurat@meta.data$celltype_percell

# Locked, consistent color palette across all 4 panels
all_types <- sort(unique(na.omit(umap_df$celltype)))
palette <- setNames(scales::hue_pal()(length(all_types)), all_types)

plots <- list()
for (title in names(panels)) {
  grp <- panels[[title]]
  sub_df <- umap_df %>% filter(panel_group == grp)
  cat(title, "-- cells:", nrow(sub_df), "\n")
  p <- ggplot(sub_df, aes(x = umap_1, y = umap_2, color = celltype)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_manual(values = palette, drop = FALSE) +
    theme_minimal() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    theme(plot.title = element_text(size = 12, face = "bold"),
         legend.position = "none")
  plots[[title]] <- p
}
cat("\n")

# Shared legend panel
legend_plot <- ggplot(data.frame(x=1,y=1,celltype=all_types), aes(x,y,color=celltype)) +
  geom_point() + scale_color_manual(values = palette) +
  theme_void() + guides(color = guide_legend(override.aes = list(size=4), ncol=1)) +
  theme(legend.text = element_text(size = 10))

combined <- (wrap_plots(plots, ncol = 2) | 
            (legend_plot + plot_layout(guides = "collect"))) +
  plot_layout(widths = c(4, 1)) +
  plot_annotation(title = "Whole-Tumor UMAP by Response and Timepoint (P18 Excluded, CellTypist-Annotated)",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))

ggsave(file.path(FIGDIR, "umap_P18excl_4panel_CellTypist_percell.png"),
       combined, width = 16, height = 11, dpi = 300)
cat("Saved umap_P18excl_4panel_CellTypist_percell.png\n\n")

cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
