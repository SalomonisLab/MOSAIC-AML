**Supplementary Table 3 | Performance of each model on each cohort.**

| Model | Cohort | Assay | Categories evaluated | Sensitivity | Specificity | AUROC |
|---|---|---|---|---|---|---|
| Bulk RNA-only (BeatAML-trained) | BeatAML (5-fold CV) | bulk | 50 | 0.532 | 0.956 | 0.830 |
| Bulk RNA-only (BeatAML-trained) | Leucegene (external) | bulk | 18 | 0.546 | 0.956 | 0.855 |
| Bulk RNA-only (BeatAML-trained) | Held-out single-cell | single-cell | 3 | 0.605 | 0.929 | 0.858 |
| Bulk RNA-only (BeatAML-trained) | All single-cell | single-cell | 14 | 0.303 | 0.946 | 0.695 |
| MOSAIC-AML multimodal | All single-cell (CV) | single-cell | 28 | 0.505 | 0.954 | 0.908 |
| MOSAIC-AML multimodal | BeatAML / Leucegene | bulk | n/a | falls back to the bulk caller | — | — |

The bulk model performs well within bulk (AUROC 0.830–0.855) but degrades on single-cell input (0.695). The multimodal model requires single-cell modalities and therefore has no independent bulk-cohort performance; on bulk input the platform operates as its bulk caller.
