#!/bin/bash
#SBATCH --job-name=volcano65
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --time=00:15:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/volcano65_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/volcano65_%J.err

cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({library(ggplot2); library(ggrepel)})

res <- read.csv("/ibex/user/alghamlk/HCC_P18_excluded/data/pseudobulk_deseq2_PRvsPD_P18excl.csv",
                row.names = 1)
res$gene <- rownames(res)
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
cat("Total genes plotted:", nrow(res), "\n")

res$sig <- "Not significant"
res$sig[res$padj < 0.05 & res$log2FoldChange < 0] <- "Up in PD (non-responders)"
res$sig[res$padj < 0.05 & res$log2FoldChange > 0] <- "Up in PR (responders)"
cat("\nClassification:\n"); print(table(res$sig))

label_genes <- c("CD14","FCGR1A","FCGR3A","SYK","MERTK","ADORA3","FLT3","GATD3A")
res$lab <- ifelse(res$gene %in% label_genes & res$padj < 0.05, res$gene, NA)

cols <- c("Not significant"            = "#C8CDD2",
          "Up in PD (non-responders)"  = "#A01414",
          "Up in PR (responders)"      = "#0064A0")

res$sig <- factor(res$sig, levels = c("Not significant",
                                      "Up in PD (non-responders)",
                                      "Up in PR (responders)"))
bg  <- res[res$sig == "Not significant", ]
sig <- res[res$sig != "Not significant", ]

p <- ggplot() +
  # background: faint grey dots
  geom_point(data = bg, aes(x = log2FoldChange, y = -log10(padj)),
             color = "#B8BEC4", size = 1.0, alpha = 0.28, stroke = 0) +
  # significant: shape 21 = colored outer ring + faint same-color fill
  geom_point(data = sig, aes(x = log2FoldChange, y = -log10(padj),
                             fill = sig, color = sig),
             shape = 21, size = 2.9, stroke = 0.7, alpha = 0.45) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_color_manual(values = cols, name = NULL, guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = 0.35, color = "#555555") +
  geom_vline(xintercept = 0, linetype = "dotted",
             linewidth = 0.3, color = "#999999") +
  geom_text_repel(data = sig, aes(x = log2FoldChange, y = -log10(padj), label = lab),
                  size = 3.5, fontface = "italic", color = "#222222",
                  max.overlaps = Inf, box.padding = 0.7, point.padding = 0.4,
                  min.segment.length = 0, segment.color = "#888888",
                  segment.size = 0.3) +
  labs(x = expression(log[2]~fold~change~(PR/PD)),
       y = expression(-log[10]~adjusted~italic(p)),
       title = "Baseline differential expression (pre-treatment PR vs PD)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(size = 13, hjust = 0, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.text = element_text(size = 11),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.7, size = 3.2)))

ggsave("/ibex/user/alghamlk/HCC_P18_excluded/figures/volcano_P18excl_65.png",
       p, width = 7.6, height = 6.6, dpi = 300, bg = "white")
ggsave("/ibex/user/alghamlk/HCC_P18_excluded/figures/volcano_P18excl_65.pdf",
       p, width = 7.6, height = 6.6, bg = "white")
cat("\nSaved volcano_P18excl_65.png and .pdf\n")
EOF

echo "Job finished: $(date)"
