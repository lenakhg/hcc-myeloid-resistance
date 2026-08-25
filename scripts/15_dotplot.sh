#!/bin/bash
#SBATCH --job-name=P18_dotplot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48GB
#SBATCH --time=00:25:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/dotplot_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/dotplot_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)
suppressMessages({library(Seurat); library(ggplot2); library(dplyr)})

DATA <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
FIG  <- "/ibex/user/alghamlk/HCC_P18_excluded/figures"

mye <- readRDS(file.path(DATA,"seurat_P18excl_myeloid_res04_named.rds"))
DefaultAssay(mye) <- "RNA"
mye <- NormalizeData(mye, verbose=FALSE)

# identity column: detect the named label column
labcol <- NULL
for (cn in c("myeloid_label","cell_type","subset","ident")) {
  if (cn %in% colnames(mye@meta.data)) { labcol <- cn; break }
}
cat("Label column:", labcol, "\n")
cat("Levels:\n"); print(table(mye@meta.data[[labcol]]))
Idents(mye) <- mye@meta.data[[labcol]]

# drop contamination for the figure
keep_id <- setdiff(levels(Idents(mye)), grep("ontam", levels(Idents(mye)), value=TRUE))
mye2 <- subset(mye, idents=keep_id)
Idents(mye2) <- droplevels(Idents(mye2))

# curated canonical markers per state (ordered to give a clean diagonal)
markers <- c(
  # C1Q/APOE-TAM
  "C1QA","C1QC","APOE","TREM2",
  # Perivascular-TAM (angiogenesis)
  "NRP1","STAB1","VCAM1","LYVE1","FOLR2",
  # MRC1/CD9-TAM
  "MRC1","CD9","SELENOP",
  # Inflammatory-TAM
  "CXCL9","CXCL10","IL1B","CCL3",
  # FCN1-monocyte
  "FCN1","VCAN","S100A12",
  # S100A8-monocyte
  "S100A8","S100A9","LYZ",
  # cDC1
  "CLEC9A","XCR1","BATF3",
  # cDC2
  "CD1C","FCER1A","CLEC10A",
  # mregDC
  "LAMP3","CCR7","IDO1","FSCN1"
)
markers <- unique(markers[markers %in% rownames(mye2)])
cat("Markers present:", length(markers), "\n")

p <- DotPlot(mye2, features=markers, dot.scale=7) +
  scale_color_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0) +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=11, face="italic"),
        axis.text.y = element_text(size=12),
        axis.title = element_blank(),
        legend.position = "right",
        panel.grid.major = element_line(color="grey92")) +
  labs(title="Myeloid subset marker expression")

ggsave(file.path(FIG,"dotplot_myeloid_markers.png"), p,
       width=13, height=6.5, dpi=300, bg="white")
cat("Saved dotplot_myeloid_markers.png\n")
cat("Done:", as.character(Sys.time()), "\n")
EOF
echo "Job finished: $(date)"
