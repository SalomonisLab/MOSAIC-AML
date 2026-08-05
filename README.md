# CIPHER-AML

**C**ell-state **I**nference of **P**athogenic **H**its from **E**xpression and **R**egulation —
the **mutation predictor** of [MOSAIC-AML](https://github.com/SalomonisLab/MOSAIC-AML).

> This is a component branch. It carries the mutation-prediction layer and the shared platform it runs
> on; the drug layer (COMPASS-AML) lives on the `compass-aml` branch, and the complete system is on
> `main`.

## What it does

Infers driver lesions from expression, cell state and regulon activity — no sequencing required.
Two complementary callers:

| caller | trained on | covers |
|---|---|---|
| **bulk variant-level** (deployed primary) | BeatAML2, n=707 with full WES | 50 variant-level categories, mean CV AUROC 0.829 |
| **single-cell multimodal** | the marrow atlas, 8 modality blocks | 26 drivers + 8 cytogenetic events, fused OOF AUROC 0.875 |

One classifier per (driver × modality), differential feature selection, percentile calibration, and
non-negative least-squares fusion floored to the best single block. Donor-grouped cross-validation
throughout — specimens from one patient never span folds — with the operating threshold chosen on
held-out folds, not the fold it is applied to.

## Run it

```
python pipeline/train_bulk_predictor.py      # -> bulk_mutation_predictor.pkl + model card
python pipeline/train_predictor.py           # -> the multimodal fused predictor
python pipeline/ingest_patient.py --sample <10x dir | .h5ad>  --name PATIENT
python pipeline/ingest_patient.py --bulk <expression.tsv>     --name PATIENT
```

Validation pages: `gui/validation.html`, `gui/calibration.html`, `gui/evidence.html`.
Methods and the full validation battery: `deliverables/METHODS_CIPHER-AML.md`.

**Standing caveat.** A call marked *predicted* is inferred from expression, not sequenced. The system
recommends a confirmatory assay for a predicted lesion; it never converts one into a therapy.
