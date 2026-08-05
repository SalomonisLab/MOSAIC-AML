**Table 1 | Pooled held-out performance across three single-cell cohorts.**

| Cohort | Specimens | Scored calls | Drivers with ≥1 positive | TP | FP | FN | TN | Sensitivity | Specificity | Precision | F1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Sealed held-out (internal) | 29 | 728 | 22 | 52 | 43 | 51 | 582 | 0.505 | 0.931 | 0.547 | 0.525 |
| Trumpp/Waclawiczek (external) | 16 | 928 | 17 | 32 | 91 | 28 | 777 | 0.533 | 0.895 | 0.260 | 0.350 |
| GSE281087 (external) | 15 | 390 | 17 | 11 | 54 | 19 | 306 | 0.367 | 0.850 | 0.169 | 0.232 |
| **All cohorts pooled** | **60** | **2046** | — | **95** | **188** | **98** | **1665** | **0.492** | **0.898** | **0.336** | **0.399** |

Every (specimen × driver) call with a known label is counted; no minimum-positive filter is applied, so a driver observed once still contributes. GSE281087 is scored panel-honestly: genes its targeted panel never assayed are left unlabelled rather than counted as wild-type. Sensitivity and specificity are micro-averaged over calls.
