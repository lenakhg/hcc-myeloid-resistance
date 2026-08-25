#!/bin/bash
#SBATCH --job-name=gsea_P18excl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --time=01:30:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/gsea_P18excl_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/gsea_P18excl_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({library(msigdbr); library(fgsea); library(dplyr)})

cat("Loading P18-excluded DESeq2 results:", as.character(Sys.time()), "\n")
deseq_res <- read.csv(
    "/ibex/user/alghamlk/HCC_P18_excluded/data/pseudobulk_deseq2_PRvsPD_P18excl.csv",
    row.names = 1)
cat("Genes:", nrow(deseq_res), "\n\n")

# Rank by DESeq2 Wald statistic (identical method to original)
deseq_res <- deseq_res %>% filter(!is.na(stat))
ranks <- setNames(deseq_res$stat, rownames(deseq_res))
ranks <- sort(ranks, decreasing = TRUE)
cat("Ranked genes:", length(ranks), "\n\n")

categories <- c("H", "C2", "C6", "C7", "C8")
category_labels <- c(H="Hallmark", C2="Curated_Pathways", C6="Oncogenic_Signatures",
                     C7="Immunologic_Signatures", C8="Cell_Type_Signatures")

cat("=== GSEA on P18-EXCLUDED ranking, per-category ===\n")
summary_table <- data.frame()

for (cat_code in categories) {
  label <- category_labels[[cat_code]]
  cat(sprintf("\n=== %s (%s) ===\n", cat_code, label))
  cat_sets <- msigdbr(species = "Homo sapiens", category = cat_code)
  pathway_list <- split(cat_sets$gene_symbol, cat_sets$gs_name)
  cat("Gene sets:", length(pathway_list), "\n")

  set.seed(42)
  fgsea_res <- fgsea(pathways = pathway_list, stats = ranks,
                     minSize = 10, maxSize = 500, eps = 0) %>% arrange(padj)

  cat("Significant (padj<0.05):", sum(fgsea_res$padj<0.05, na.rm=TRUE), "\n")
  cat("Top 10 pathways:\n")
  print(as.data.frame(fgsea_res %>% select(pathway, NES, padj, size) %>% head(10)))

  write.csv(fgsea_res %>% select(-leadingEdge),
      sprintf("/ibex/user/alghamlk/HCC_P18_excluded/data/gsea_PRvsPD_P18excl_%s.csv", cat_code),
      row.names = FALSE)

  summary_table <- rbind(summary_table, data.frame(
    category=cat_code, label=label, n_tested=length(pathway_list),
    n_sig_p05=sum(fgsea_res$padj<0.05, na.rm=TRUE),
    top_hit=fgsea_res$pathway[1], top_NES=round(fgsea_res$NES[1],3),
    top_padj=signif(fgsea_res$padj[1],3)))
}

cat("\n=== SUMMARY (P18-excluded) ===\n")
print(summary_table)
write.csv(summary_table,
    "/ibex/user/alghamlk/HCC_P18_excluded/data/gsea_PRvsPD_P18excl_summary.csv",
    row.names = FALSE)
cat("\nDone:", as.character(Sys.time()), "\n")
EOF
echo "Job finished: $(date)"
