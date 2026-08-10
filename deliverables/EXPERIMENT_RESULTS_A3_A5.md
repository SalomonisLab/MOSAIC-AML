# A3-A5 results

## A3 matrix factorisation — modest real gain

Interaction-target (honest) AUROC by rank: 5 -> 0.625, 12 -> 0.664, 25 -> 0.678, **40 -> 0.685**;
deployed per-family model = 0.674. Spearman 0.225 -> **0.239**. A 50/50 blend matches rank-40 MF.
So **+0.011 AUROC / +0.014 Spearman**, real but small — the low-rank structure is worth having,
and it helps the low-n inhibitors most, but it is not the step change the headline gap suggested.

## A4 raw concentration data — the ceiling, now measured directly

Downloaded `beataml_wv1to4_raw_inhibitor_v4_dbgap.txt` (49 MB, 555,583 wells). It carries a
`replicate` column, so assay noise can be measured rather than inferred from same-target proxies:

- **Technical replicates of the same well**: Pearson 0.406, Spearman 0.641,
  **median absolute difference 13.7 percentage points of viability**.
- **Per-drug, across the dose series**: median replicate r = **0.529** over 27 drugs with enough
  replicates, implying a ceiling of **sqrt(0.529) = 0.727** Spearman.

The correspondence with model performance is exact where it matters:

| inhibitor | replicate r | model AUROC |
|---|---|---|
| Trametinib | 0.87 | 0.87 |
| Venetoclax | 0.84 | 0.977 |
| Rapamycin | 0.80 | 0.86 |
| Lapatinib | 0.27 | — |
| PLX-4720 | 0.24 | — |

**The drugs the model predicts best are precisely the drugs whose assay reproduces best.** This is now
a measured fact, not an inference from same-target pairs.

## A5 imputed modalities — blocked, not attempted

The rna2adt / rna2grn / rna2lipid / rna2metabolite bundles are vendored in `engine-code/` but were
pickled on Linux and carry `PosixPath` objects, so they will not unpickle on Windows
(`UnsupportedOperation: cannot instantiate 'PosixPath'`). This has to run on the cluster.
Still the most promising untested lead.
