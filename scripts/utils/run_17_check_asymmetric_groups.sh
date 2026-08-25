#!/bin/bash
#SBATCH --job-name=P18_cc_asymcheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=00:20:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/ccasym_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/ccasym_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
library(CellChat)
library(dplyr)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("For each of the 27 REAL/GRADED/SIGNIFICANT pathways, check\n")
cat("whether its significant L-R edges involve mregDC or\n")
cat("S100A8-monocyte (the two groups asymmetric between arms)\n")
cat("############################################################\n")

cellchat <- readRDS(file.path(DATADIR, "cellchat_preonly_merged_PRvsPD.rds"))
classified <- read.csv(file.path(DATADIR, "cellchat_preonly_rankNet_CLASSIFIED.csv"))
real_pathways <- classified %>% filter(status == "REAL, GRADED, SIGNIFICANT") %>% pull(name)
cat("Real pathways to check:", length(real_pathways), "\n")
print(real_pathways)
cat("\n")

df_list <- subsetCommunication(cellchat)
cat("Datasets in merged object:", length(df_list), "\n\n")

asymmetric_groups <- c("mregDC", "S100A8-monocyte")

for (i in seq_along(df_list)) {
  cat("=== Dataset", i, "===\n")
  df <- df_list[[i]]
  flagged <- df %>%
    filter(pathway_name %in% real_pathways,
          source %in% asymmetric_groups | target %in% asymmetric_groups) %>%
    select(pathway_name, source, target, prob, pval) %>%
    distinct()
  if (nrow(flagged) > 0) {
    cat("Real pathways WITH edges involving mregDC/S100A8-monocyte:\n")
    print(flagged)
  } else {
    cat("NONE of the real pathways involve mregDC or S100A8-monocyte\n")
    cat("as sender or receiver in this dataset.\n")
  }
  cat("\n")
}

cat("Done:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
