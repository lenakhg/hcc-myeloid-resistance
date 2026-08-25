#!/bin/bash
#SBATCH --job-name=P18_cc_rigor
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=00:30:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/ccrigor_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/ccrigor_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
library(CellChat)
library(dplyr)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"

cat("############################################################\n")
cat("Rigorous pathway significance check on the pre-only PR vs PD\n")
cat("CellChat comparison. rankNet(do.stat=TRUE); separate real,\n")
cat("graded signal from degenerate one-arm-only artifacts.\n")
cat("############################################################\n")

cellchat <- readRDS(file.path(DATADIR, "cellchat_preonly_merged_PRvsPD.rds"))

gg <- rankNet(cellchat, mode = "comparison", measure = "weight",
              stacked = FALSE, do.stat = TRUE)
rankNet_data <- as.data.frame(gg$data)

# ── SAVE FIRST, before any printing (so a display bug can't lose it)
write.csv(rankNet_data,
    file.path(DATADIR, "cellchat_preonly_rankNet_raw.csv"), row.names = FALSE)
cat("Saved cellchat_preonly_rankNet_raw.csv (full, unfiltered)\n\n")

cat("Columns returned:\n"); print(colnames(rankNet_data)); cat("\n")
cat("Full raw table:\n")
print(rankNet_data)                 # plain data.frame print, no n= arg
cat("\n")

# ── which groups existed in BOTH arms vs one arm (cell-count flag) ──
# PR dropped these for <10 cells: mregDC, Perivascular-TAM, S100A8-monocyte
pd_only_groups <- c("mregDC","S100A8-monocyte","Perivascular-TAM")
cat("Groups present in PD arm but dropped from PR (cell-count):",
    paste(pd_only_groups, collapse=", "), "\n")
cat("=> any pathway routing mainly through these is cell-count-suspect,\n")
cat("   even if statistically 'significant'.\n\n")

cat("############################################################\n")
cat("Classify each pathway: real/graded vs degenerate artifact\n")
cat("############################################################\n")
per_pathway <- rankNet_data %>%
  group_by(name) %>% slice(1) %>% ungroup() %>%
  mutate(status = case_when(
    !is.finite(contribution.relative.1)      ~ "DEGENERATE (Inf ratio)",
    contribution.relative.1 == 0             ~ "DEGENERATE (zero ratio)",
    contribution.relative.1 == 1             ~ "NEAR-NULL (ratio=1)",
    pvalues >= 0.05                          ~ "NOT SIGNIFICANT",
    TRUE                                     ~ "REAL, GRADED, SIGNIFICANT"
  )) %>%
  arrange(status, pvalues)

cat("Pathway count by classification:\n")
print(as.data.frame(table(status = per_pathway$status)))
cat("\n")

real_sig <- per_pathway %>% filter(status == "REAL, GRADED, SIGNIFICANT")
cat("============================================================\n")
cat("REAL, GRADED, SIGNIFICANT pathways (p<0.05, finite ratio):\n")
cat("============================================================\n")
print(as.data.frame(real_sig[, c("name","group","contribution",
                                 "contribution.relative.1","pvalues")]))
cat("\n")

deg <- per_pathway %>% filter(grepl("DEGENERATE|NEAR-NULL", status))
cat("============================================================\n")
cat("DEGENERATE / artifact pathways (do NOT report as findings):\n")
cat("============================================================\n")
print(as.data.frame(deg[, c("name","status","contribution.relative.1","pvalues")]))
cat("\n")

write.csv(as.data.frame(per_pathway),
    file.path(DATADIR, "cellchat_preonly_rankNet_CLASSIFIED.csv"), row.names = FALSE)
cat("Saved cellchat_preonly_rankNet_CLASSIFIED.csv\n\n")

cat("############################################################\n")
cat("Cross-check vs the naive setdiff() lists from the first run\n")
cat("############################################################\n")
pd_only_raw <- c("VCAM","CD23","CD70","CD80","XCR","CD46","MK","GDF",
                 "CD137","GRN","PVR","NOTCH","SN","DESMOSOME","NEGR")
pr_only_raw <- c("LCK","ANXA1","VEGI","SEMA6","CDH")
cat("Original 'PD-only' names -- how they classify now:\n")
print(as.data.frame(per_pathway %>% filter(name %in% pd_only_raw) %>%
        select(name, status, pvalues)))
cat("\nOriginal 'PR-only' names -- how they classify now:\n")
print(as.data.frame(per_pathway %>% filter(name %in% pr_only_raw) %>%
        select(name, status, pvalues)))

cat("\nDone:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
