# A Baseline Myeloid Immunosuppression Signature Associated with Resistance to Anti–PD-1 Plus Lenvatinib in Hepatocellular Carcinoma

Single-cell RNA-sequencing re-analysis nominating druggable myeloid targets of immunotherapy resistance in HBV-positive HCC.

**Author:** Lena Alghamdi¹²³ · **Supervisor:** Dr. Hugo Tavares³
¹College of Pharmacy, King Khalid University · ²KAUST Academy · ³Department of Genetics, University of Cambridge

---

## Overview

Most patients with advanced hepatocellular carcinoma (HCC) derive no durable benefit from anti–PD-1 plus lenvatinib, and no validated pre-treatment biomarker of response exists. This project re-analyses a public single-cell RNA-sequencing dataset (GSE235863; nine HBV-positive HCC patients) to test whether the **baseline myeloid landscape** distinguishes responders from non-responders, and whether the resulting signature is **therapeutically actionable**.

**Key result:** A patient-level pseudobulk comparison of pre-treatment tumours (3 partial responders vs 2 progressors) identifies a **65-gene signature** (46 elevated in progressive disease, 19 in partial response). The progression-associated genes form a coherent myeloid/macrophage programme (CD14, FCGR1A, FCGR3A, SYK, MERTK, ADORA3), enriched for Kupffer-cell and TNFα/NF-κB inflammatory biology, and nominate **three druggable, repurposable targets**: SYK (fostamatinib), MERTK/FLT3 (gilteritinib) and ADORA3 (namodenoson).

> ⚠️ These findings are **hypothesis-generating**, limited by a small cohort (n = 5 in the primary comparison), and require prospective validation.

---

## Data availability

| Resource | Location |
|---|---|
| **Raw data** (scRNA-seq, 9 patients) | [GEO: GSE235863](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235863) (Guo et al., *Cancer Cell* 2025) |
| **Processed Seurat objects** (`.rds`) | Zenodo: [(https://doi.org/10.5281/zenodo.22131800)] |
| **Result tables & figures** | This repository (`results/`, `figures/`) |

The raw count matrices are downloaded from GEO. The processed Seurat objects (post-normalisation, annotation, and sub-clustering) are large and are archived on Zenodo; the scripts below regenerate them from the raw data.

---

## The complete-responder exclusion

The dataset contains one complete responder (patient P18) contributing ~18% of cells with an atypical plasma-cell population. To avoid biasing shared normalisation and the PCA embedding, **P18 is excluded** from the analyses in this repository (hence the `P18excl` naming).

---

## Pipeline

Scripts are numbered in run order and are SLURM batch scripts for the KAUST Ibex HPC cluster (R 4.5.0). Paths inside the scripts are absolute (`/ibex/user/alghamlk/HCC_P18_excluded/...`); adjust these to your environment before running.

| # | Script | Purpose |
|---|---|---|
| 00 | `00_preprocess_sct_pca.sh` | SCTransform normalisation + PCA (shared embedding) |
| 01 | `01a_cluster_markers.sh`, `01b_singleR_annotate.sh` | Leiden clustering + initial SingleR annotation |
| 02 | `02_silhouette_sweep.sh` | Silhouette sweep (106 k/resolution combinations) → continuous structure |
| 03 | `03a–03e_annotate_*.sh` | Cell-type annotation: CellTypist, Monaco, Azimuth, scType, SCINA |
| 04 | `04_annotation_breakdown_umap.sh` | Cross-method annotation concordance + UMAP |
| 05 | `05_final_annotation.sh` | Consensus whole-tumour annotation (7 cell types) |
| 06 | `06_singleR_fine_15cluster.sh` | Fine 15-cluster annotation |
| 07 | `07a–07d_myeloid_*.sh` | Myeloid sub-clustering (9 states) + named UMAP |
| 08 | `08_build_combined_object.sh` | Build combined labelled object (`seurat_P18excl_combined.rds`) |
| **09** | **`09_pseudobulk_deseq2.sh`** | **Primary: pseudobulk DESeq2, PR vs PD → 65-gene signature** |
| 10 | `10_volcano.sh` | Volcano plot |
| 11 | `11_gsea.sh` | GSEA across MSigDB H, C2, C6, C7, C8 |
| 12 | `12_ontreatment_dynamics.sh` | Paired pre/post pseudobulk in responders |
| 13 | `13_dgidb_druggability.sh` | Unbiased DGIdb query of all 65 genes |
| 14 | `14a–14b_cellchat_*.sh` | Myeloid–T-cell communication (CellChat, preliminary) |
| 15 | `15_dotplot.sh` | Myeloid marker dotplot |
| 16 | `16_heatmaps.sh` | Signature/marker heatmaps |

`scripts/utils/` contains sanity-check scripts (cell counts, timepoints, metadata) used during development, kept for transparency.

---

## Reproducibility map

Every figure and number in the manuscript maps to a script and an output file:

| Manuscript element | Script | Output file |
|---|---|---|
| 65-gene signature; Table 2 | `09_pseudobulk_deseq2.sh` | `results/deseq2_PRvsPD_P18excl_SIG.csv` |
| Full DE (all 15,805 genes) | `09_pseudobulk_deseq2.sh` | `results/pseudobulk_deseq2_PRvsPD_P18excl.csv` |
| Volcano figure | `10_volcano.sh` | `figures/volcano_P18excl_65.png` |
| GSEA; Table 3 | `11_gsea.sh` | `results/gsea_PRvsPD_P18excl_{H,C2,C6,C7,C8}.csv` |
| On-treatment dynamics; §3.5 | `12_ontreatment_dynamics.sh` | `results/pseudobulk_deseq2_PR_PrePost_P18excl.csv` |
| Druggability; Table 4 | `13_dgidb_druggability.sh` | `results/dgidb_druggability_P18excl_65.csv` |
| Myeloid UMAP | `07d_myeloid_umap_named.sh` | `figures/umap_myeloid_named_single.png` |
| Myeloid marker dotplot | `15_dotplot.sh` | `figures/dotplot_myeloid_markers.png` |
| CellChat (supplementary) | `14a–14b` | `results/cellchat_preonly_rankNet_CLASSIFIED.csv` |
| Abundance/activation analysis (Discussion) | `17_abundance_activation.sh` | `macrophage_restricted_DE_PRvsPD.csv`, `myeloid_abundance_PRvsPD.csv` |
| Silhouette (Methods 2.3) | `02_silhouette_sweep.sh` | `results/silhouette_sweep_P18excl.csv` |

---

## Key methods notes

- **Pseudobulk, not per-cell:** raw counts are summed per patient (not per cell), giving n = 5 patient-level samples, which avoids pseudoreplication.
- **DESeq2 contrast:** `~ recist` with PD as the reference level; **negative log2FC = elevated in PD** (non-responders).
- **Filter:** genes retained with ≥10 counts in ≥2 of the 5 patients (15,805 genes tested).
- **GSEA ranking:** genes ranked by the DESeq2 Wald statistic; GSEA is reported as a finding characterising the enriched biology (it is not independent of the DESeq2 ranking).
- **Druggability:** all 65 significant genes are queried against DGIdb; of the 29 druggable hits, three leads are prioritised by pre-specified criteria (PD-elevated + approved/phase III agent + myeloid-immunosuppressive mechanism). DGIdb queried 2026-08-25.
- **DepMap co-dependency** is pan-cancer and presented as supporting, not primary, evidence.

---

## Environment

- **R** 4.5.0 (GNU 12.2.0)
- Key packages: Seurat 5.5.1, SCTransform, DESeq2, fgsea, msigdbr, CellChat, SingleR, CellTypist, scType, SCINA, httr2
- Full session details: `environment/sessionInfo.txt`
- Random seed: 42 (fixed throughout)
- Compute: KAUST Ibex HPC cluster (SLURM)

---

## Citation

If you use this code or the derived results, please cite:

> Alghamdi L, Tavares H. *A Baseline Myeloid Immunosuppression Signature Associated with Resistance to Anti–PD-1 Plus Lenvatinib in Hepatocellular Carcinoma.* 2026. GitHub: https://github.com/lenakhg/hcc-myeloid-resistance

And the original data:

> Guo X, Nie H, Zhang W, et al. Contrasting cytotoxic and regulatory T cell responses underlying distinct clinical outcomes to anti-PD-1 plus lenvatinib therapy in cancer. *Cancer Cell.* 2025;43(2):248-268.e9. doi:10.1016/j.ccell.2025.01.001

See `CITATION.cff` for machine-readable citation metadata.

---

## License

Code released under the MIT License (see `LICENSE`). The underlying data (GSE235863) is subject to its original terms of use at GEO.

---

## Acknowledgements

The author thanks the KAUST Supercomputing Core Laboratory for access to the Ibex HPC cluster, and the KAUST Academy for programme support.
