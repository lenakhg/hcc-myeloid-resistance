#!/bin/bash
#SBATCH --job-name=myeloid_abund_activ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=00:30:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/abund_activ_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/abund_activ_%J.err

cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({library(Seurat); library(dplyr); library(Matrix)})

obj <- readRDS("/ibex/user/alghamlk/HCC_P18_excluded/data/seurat_P18excl_myeloid_res04_named.rds")

# ---- Filter: pre-treatment, PR/PD only, drop contamination ----
obj <- subset(obj, subset = timepoint_cell == "Pre-treatment" &
                            recist %in% c("PR","PD") &
                            myeloid_label != "Contamination (T/hepatocyte)")
cat("Pre-treatment PR/PD myeloid cells:", ncol(obj), "\n")
cat("\nCells per patient x response:\n")
print(table(obj$patient, obj$recist))

# =========================================================
# ANALYSIS 1: DIFFERENTIAL ABUNDANCE (proportions per patient)
# =========================================================
cat("\n\n===== ANALYSIS 1: MYELOID STATE ABUNDANCE (PR vs PD) =====\n")
md <- obj@meta.data %>% select(patient, recist, myeloid_label)

# per-patient proportion of each state
prop_tbl <- md %>%
  group_by(patient, recist, myeloid_label) %>%
  summarise(n = n(), .groups="drop") %>%
  group_by(patient) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# mean proportion by response group
summary_tbl <- prop_tbl %>%
  group_by(recist, myeloid_label) %>%
  summarise(mean_prop = mean(prop), .groups="drop") %>%
  tidyr::pivot_wider(names_from = recist, values_from = mean_prop, values_fill = 0) %>%
  mutate(PD_minus_PR = PD - PR) %>%
  arrange(desc(PD_minus_PR))

cat("\nMean proportion of each myeloid state, by response group:\n")
cat("(positive PD_minus_PR = more abundant in non-responders)\n\n")
print(as.data.frame(summary_tbl))

write.csv(as.data.frame(summary_tbl),
  "/ibex/user/alghamlk/HCC_P18_excluded/data/myeloid_abundance_PRvsPD.csv", row.names=FALSE)

# =========================================================
# ANALYSIS 2: MACROPHAGE-RESTRICTED PSEUDOBULK DE
# (does the signature survive when we compare macrophage-to-macrophage?)
# =========================================================
cat("\n\n===== ANALYSIS 2: MACROPHAGE-RESTRICTED PSEUDOBULK DE =====\n")
suppressMessages(library(DESeq2))

# define macrophage/TAM states (the 4 TAM + could include monocytes; here: TAM states)
macro_states <- c("C1Q/APOE-TAM","Perivascular-TAM","MRC1/CD9-TAM","Inflammatory-TAM")
mac <- subset(obj, subset = myeloid_label %in% macro_states)
cat("Macrophage-only cells:", ncol(mac), "\n")
cat("Macrophage cells per patient x response:\n")
print(table(mac$patient, mac$recist))

# pseudobulk: sum raw counts per patient
raw <- GetAssayData(mac, assay="RNA", layer="counts")
pats <- unique(mac$patient)
pb <- sapply(pats, function(p){
  cells <- rownames(mac@meta.data)[mac@meta.data$patient == p]
  Matrix::rowSums(raw[, cells, drop=FALSE])
})
colnames(pb) <- pats

coldata <- mac@meta.data %>% distinct(patient, recist) %>%
  arrange(match(patient, colnames(pb)))
coldata <- data.frame(row.names=coldata$patient,
                      recist=factor(coldata$recist, levels=c("PD","PR")))
cat("\nMacrophage pseudobulk design:\n"); print(coldata)

dds <- DESeqDataSetFromMatrix(round(pb), coldata, design = ~ recist)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("recist","PR","PD"))
res <- res[order(res$pvalue),]

cat("\nMacrophage-restricted DE: top 20 genes\n")
print(head(as.data.frame(res), 20))

# check YOUR nominated genes specifically
cat("\n=== Do nominated genes stay elevated in PD within macrophages? ===\n")
cat("(negative log2FC = still up in PD when comparing macrophage-to-macrophage)\n\n")
for (g in c("CD14","FCGR1A","FCGR3A","SYK","MERTK","ADORA3","FLT3")) {
  if (g %in% rownames(res)) {
    r <- res[g,]
    cat(sprintf("%-8s log2FC=%+.2f  padj=%.2e\n", g, r$log2FoldChange, ifelse(is.na(r$padj),NA,r$padj)))
  } else cat(sprintf("%-8s not tested (filtered)\n", g))
}

write.csv(as.data.frame(res),
  "/ibex/user/alghamlk/HCC_P18_excluded/data/macrophage_restricted_DE_PRvsPD.csv", row.names=TRUE)

cat("\nDone.\n")
EOF
