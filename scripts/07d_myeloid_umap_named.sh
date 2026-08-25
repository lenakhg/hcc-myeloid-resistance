#!/bin/bash
#SBATCH --job-name=P18excl_myeUMAP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=00:45:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/myeUMAP_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/myeUMAP_%J.err
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
cat("Myeloid UMAP (res 0.4), descriptive names. Single + 4-panel.\n")
cat("Contamination (subcluster 6) shown in grey. Subcluster 7\n")
cat("merged into cDC2. Publication palette.\n")
cat("############################################################\n")
mye <- readRDS(file.path(DATADIR, "seurat_P18excl_myeloid_res04.rds"))
Idents(mye) <- mye@meta.data[["SCT_snn_res.0.4"]]
cat("Subclusters:", length(unique(Idents(mye))), " cells:", ncol(mye), "\n\n")

# ── Descriptive names (analyst interpretation of markers).
# 7 -> cDC2 (merged). 6 -> contamination (shown, greyed). ────────
name_map <- c(
  "1"  = "C1Q/APOE-TAM",
  "2"  = "cDC2",
  "3"  = "FCN1-monocyte",
  "4"  = "MRC1/CD9-TAM",
  "5"  = "S100A8-monocyte",
  "6"  = "Contamination (T/hepatocyte)",
  "7"  = "cDC2",                       # merged into 2
  "8"  = "Inflammatory-TAM",
  "9"  = "cDC1",
  "10" = "Perivascular-TAM",
  "11" = "mregDC"
)
lab <- unname(name_map[as.character(Idents(mye))])
# order: biological groups first, contamination last
bio_order <- c("C1Q/APOE-TAM","Perivascular-TAM","MRC1/CD9-TAM","Inflammatory-TAM",
               "FCN1-monocyte","S100A8-monocyte","cDC1","cDC2","mregDC",
               "Contamination (T/hepatocyte)")
mye$myeloid_label <- factor(lab, levels = bio_order)

cat("Named myeloid subcluster distribution:\n")
print(table(mye$myeloid_label))
cat("\n")

# ── Publication palette: Okabe-Ito-based distinct colors for the 9
# biological groups; contamination forced to grey. ──────────────
bio_cols <- c(
  "C1Q/APOE-TAM"     = "#E69F00",
  "Perivascular-TAM" = "#D55E00",
  "MRC1/CD9-TAM"     = "#CC79A7",
  "Inflammatory-TAM" = "#B22222",
  "FCN1-monocyte"    = "#0072B2",
  "S100A8-monocyte"  = "#56B4E9",
  "cDC1"             = "#009E73",
  "cDC2"             = "#66C2A5",
  "mregDC"           = "#9370DB",
  "Contamination (T/hepatocyte)" = "grey75"
)
palette <- bio_cols[levels(mye$myeloid_label)]

umap_df <- as.data.frame(Embeddings(mye, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1","UMAP_2")
umap_df$label       <- mye$myeloid_label
umap_df$panel_group <- paste(mye$recist, mye$timepoint_cell)

# centroids for on-plot labels (single UMAP)
cent <- umap_df %>% group_by(label) %>%
  summarise(x = median(UMAP_1), y = median(UMAP_2), .groups="drop")

# ═══════════════ SINGLE UMAP ═══════════════
p_single <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = label)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = palette, name = "Myeloid subset") +
  ggrepel::geom_text_repel(data = cent, aes(x, y, label = label),
                           color = "black", size = 4, fontface = "bold",
                           bg.color = "white", bg.r = 0.15, seed = 42) +
  theme_bw(base_size = 13) +
  labs(title = "Myeloid Compartment, Sub-clustered (P18 Excluded, res 0.4)",
       x = "UMAP 1", y = "UMAP 2") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggsave(file.path(FIGDIR, "umap_myeloid_named_single.png"),
       p_single, width = 11, height = 8, dpi = 300)
cat("Saved umap_myeloid_named_single.png\n")

# ═══════════════ 4-PANEL by RECIST x timepoint ═══════════════
panels <- c("PR - Pre-treatment"  = "PR Pre-treatment",
            "PR - Post-treatment" = "PR Post-treatment",
            "PD - Pre-treatment"  = "PD Pre-treatment",
            "PD - Post-treatment" = "PD Post-treatment")
make_panel <- function(title) {
  sub_df <- umap_df %>% filter(panel_group == panels[[title]])
  cat(sprintf("%-22s cells: %d\n", title, nrow(sub_df)))
  ggplot(sub_df, aes(UMAP_1, UMAP_2, color = label)) +
    geom_point(size = 0.45, alpha = 0.8) +
    scale_color_manual(values = palette, drop = FALSE, name = "Myeloid subset") +
    theme_bw(base_size = 12) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "none")
}
q1 <- make_panel("PR - Pre-treatment")
q2 <- make_panel("PR - Post-treatment")
q3 <- make_panel("PD - Pre-treatment")
q4 <- make_panel("PD - Post-treatment")

# legend carrier
leg <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = label)) +
  geom_point(size = 0.45) +
  scale_color_manual(values = palette, drop = FALSE, name = "Myeloid subset") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
  theme(legend.title = element_text(face = "bold"))

combined <- (q1 | q2) / (q3 | q4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Myeloid Subsets by Response and Timepoint (P18 Excluded, res 0.4)",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  legend.position = "right"))
combined <- combined & theme(legend.position = "right")
# collect a single legend from the panels
combined <- (q1 | q2) / (q3 | q4) +
  plot_layout(guides = "collect") &
  scale_color_manual(values = palette, drop = FALSE, name = "Myeloid subset") &
  theme(legend.position = "right")
combined <- combined + plot_annotation(
  title = "Myeloid Subsets by Response and Timepoint (P18 Excluded, res 0.4)",
  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))

ggsave(file.path(FIGDIR, "umap_myeloid_named_4panel.png"),
       combined, width = 15, height = 11, dpi = 300)
cat("Saved umap_myeloid_named_4panel.png\n\n")

saveRDS(mye, file.path(DATADIR, "seurat_P18excl_myeloid_res04_named.rds"))
cat("Saved seurat_P18excl_myeloid_res04_named.rds (myeloid_label added)\n")
cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
