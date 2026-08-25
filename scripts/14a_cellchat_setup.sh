#!/bin/bash
#SBATCH --job-name=P18_cc_prePRPD
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --time=04:00:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/ccprePRPD_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/ccprePRPD_%J.err
#SBATCH --mail-user=alghamlk@kaust.edu.sa
#SBATCH --mail-type=ALL

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'

.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
set.seed(42)
library(Seurat)
library(CellChat)
library(dplyr)
library(patchwork)
options(stringsAsFactors = FALSE)

DATADIR <- "/ibex/user/alghamlk/HCC_P18_excluded/data"
FIGDIR  <- "/ibex/user/alghamlk/HCC_P18_excluded/figures"

cat("############################################################\n")
cat("PRE-TREATMENT ONLY, PR vs PD comparison.\n")
cat("CAVEATS (user-acknowledged): PR=3 patients, PD=2 patients\n")
cat("(patient-confounded); small groups may drop per arm;\n")
cat("merge may be partial. Individual objects saved before merge.\n")
cat("############################################################\n")
seurat <- readRDS(file.path(DATADIR, "seurat_P18excl_combined.rds"))

# pre-treatment only, exclude contamination and SD
seurat <- seurat[, seurat$timepoint_cell == "Pre-treatment" &
                   seurat$cell_type_combined != "Contamination (T/hepatocyte)" &
                   seurat$recist %in% c("PR","PD")]
seurat$cell_type_combined <- droplevels(factor(seurat$cell_type_combined))
cat("Pre-treatment PR+PD cells:", ncol(seurat), "\n")
print(table(recist = seurat$recist))
cat("\nCells per group per arm:\n")
print(table(seurat$cell_type_combined, seurat$recist))
cat("\n")

DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat, verbose = FALSE)

run_arm <- function(grp) {
  cat("\n===== ", grp, " =====\n")
  sub <- seurat[, seurat$recist == grp]
  tb <- table(droplevels(factor(sub$cell_type_combined)))
  small <- names(tb)[tb < 10]
  if (length(small) > 0) cat("DROPPED (<10 cells) in", grp, ":", paste(small, collapse=", "), "\n")
  keep <- sub$cell_type_combined %in% names(tb)[tb >= 10]
  sub <- sub[, keep]
  sub$cell_type_combined <- droplevels(factor(sub$cell_type_combined))
  cat(grp, "groups kept:", nlevels(sub$cell_type_combined), "| cells:", ncol(sub), "\n")

  data.input <- GetAssayData(sub, assay = "RNA", layer = "data")
  meta <- data.frame(labels = sub$cell_type_combined, row.names = colnames(sub))
  cc <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  cat(grp, "pathways:", length(cc@netP$pathways), "\n")
  cc
}

cc.PR <- run_arm("PR")
saveRDS(cc.PR, file.path(DATADIR, "cellchat_preonly_PR.rds"))
cc.PD <- run_arm("PD")
saveRDS(cc.PD, file.path(DATADIR, "cellchat_preonly_PD.rds"))
cat("\nIndividual objects saved (safe even if merge fails).\n")

# ── liftCellChat to a common group set, then merge ──────────────
tryCatch({
  all_groups <- union(levels(cc.PR@idents), levels(cc.PD@idents))
  cat("\nUnion of groups across arms:", length(all_groups), "\n")
  cc.PR2 <- liftCellChat(cc.PR, all_groups)
  cc.PD2 <- liftCellChat(cc.PD, all_groups)
  object.list <- list(PR = cc.PR2, PD = cc.PD2)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  saveRDS(cellchat, file.path(DATADIR, "cellchat_preonly_merged_PRvsPD.rds"))
  saveRDS(object.list, file.path(DATADIR, "cellchat_preonly_objlist_PRvsPD.rds"))
  cat("MERGE SUCCEEDED. Generating comparison figures...\n")

  p1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2))
  p2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "weight")
  ggsave(file.path(FIGDIR, "cc_preonly_total_PRvsPD.png"), p1+p2, width=10, height=5, dpi=300)

  png(file.path(FIGDIR, "cc_preonly_diffheatmap_PRvsPD.png"), width=2400, height=2000, res=300)
  print(netVisual_heatmap(cellchat)); dev.off()

  gg <- lapply(seq_along(object.list), function(i)
    netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i]))
  ggsave(file.path(FIGDIR, "cc_preonly_roles_PRvsPD.png"),
         wrap_plots(gg, ncol=2), width=12, height=6, dpi=300)

  gg1 <- rankNet(cellchat, mode="comparison", stacked=TRUE, do.stat=TRUE)
  ggsave(file.path(FIGDIR, "cc_preonly_rankNet_PRvsPD.png"), gg1, width=8, height=10, dpi=300)

  cat("PR-only pathways:", paste(setdiff(cc.PR@netP$pathways, cc.PD@netP$pathways), collapse=", "), "\n")
  cat("PD-only pathways:", paste(setdiff(cc.PD@netP$pathways, cc.PR@netP$pathways), collapse=", "), "\n")
}, error = function(e) {
  cat("\n!!! MERGE/COMPARISON FAILED:", conditionMessage(e), "\n")
  cat("Individual PR and PD CellChat objects are still saved and usable.\n")
})

cat("\nDone:", as.character(Sys.time()), "\n")
EOF

echo "Job finished: $(date)"
