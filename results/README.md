# Results

Derived result tables from the P18-excluded analysis. Each maps to a script in `../scripts/`.

| File | Description | Script |
|---|---|---|
| `pseudobulk_deseq2_PRvsPD_P18excl.csv` | Full DESeq2 result, all 15,805 tested genes (PR vs PD, pre-treatment) | `09_pseudobulk_deseq2.sh` |
| `deseq2_PRvsPD_P18excl_SIG.csv` | The 65 significant genes (padj < 0.05) | `09_pseudobulk_deseq2.sh` |
| `deseq2_PRvsPD_P18excl_ranked.csv` | Genes ranked by Wald statistic (input for GSEA) | `09_pseudobulk_deseq2.sh` |
| `gsea_PRvsPD_P18excl_H.csv` | GSEA, Hallmark collection | `11_gsea.sh` |
| `gsea_PRvsPD_P18excl_C2.csv` | GSEA, curated pathways | `11_gsea.sh` |
| `gsea_PRvsPD_P18excl_C6.csv` | GSEA, oncogenic signatures | `11_gsea.sh` |
| `gsea_PRvsPD_P18excl_C7.csv` | GSEA, immunologic signatures | `11_gsea.sh` |
| `gsea_PRvsPD_P18excl_C8.csv` | GSEA, cell-type signatures | `11_gsea.sh` |
| `gsea_PRvsPD_P18excl_summary.csv` | GSEA top-hit summary across collections | `11_gsea.sh` |
| `pseudobulk_deseq2_PR_PrePost_P18excl.csv` | On-treatment dynamics (paired pre/post in responders) | `12_ontreatment_dynamics.sh` |
| `dgidb_druggability_P18excl_65.csv` | All drug–gene interactions for the 65 genes | `13_dgidb_druggability.sh` |
| `dgidb_druggable_genes_P18excl_65.csv` | Per-gene druggability summary (29 druggable genes) | `13_dgidb_druggability.sh` |
| `silhouette_sweep_P18excl.csv` | Silhouette width across 106 k/resolution combinations | `02_silhouette_sweep.sh` |
| `cellchat_preonly_rankNet_CLASSIFIED.csv` | CellChat pathway comparison (preliminary) | `14a–14b` |
| `P18excl_myeloid_subcluster_markers_res0.4_final.csv` | Myeloid sub-cluster markers (9 states) | `07c_myeloid_res04_finalize.sh` |

**Note on DESeq2 direction:** the contrast is PR vs PD with PD as reference, so **negative log2FoldChange = elevated in Progressive Disease (non-responders)**.
