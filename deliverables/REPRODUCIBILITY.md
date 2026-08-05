# Reproducibility manifest — MOSAIC-AML mutation caller

Everything needed to re-derive every number in the validation deliverables. Paths are relative to the
project root (`AML-multimodal/`), mirrored on the cluster at
`/data/salomonis-archive/LabFiles/Nicholas/AML-multimodal/`.

## Models (frozen)

| Artifact | What | Trained on |
|---|---|---|
| `pipeline/mutation_predictor.pkl` | Multimodal caller (8 modalities) | Single-cell atlas (n≈357 labeled) |
| `pipeline/model_card.json` | Multimodal per-mutation card (AUROC, threshold, weights) | — |
| `pipeline/bulk_mutation_predictor.pkl` | Bulk RNA-only caller (50 categories) | BeatAML2 (n=707) |
| `pipeline/bulk_model_card.json` | Bulk per-category card (CV AUROC, threshold) | — |

## Data splits

- `pipeline/holdout_samples.txt` — the **29 sealed held-out** single-cell samples (masked during training, never in any CV fold).
- Bundle `../aml-bakeoff/bundle_data.npz` — `ba_*` (BeatAML 707), `lg_*` (Leucegene 367), `sc_*` (atlas 387), `drivers` (58 categories), `holdout` (29). BeatAML RNA↔DNA crosswalk: strip/replace trailing `R`→`D`.

## Recipe (deterministic)

- **Multimodal**: per-modality `StandardScaler` → differential top-500 select → `LinearSVC(C=0.02, class_weight="balanced")` → percentile → **donor-grouped 3-fold CV-OOF** → ridge-NNLS fusion → F1-max threshold.
- **Bulk**: `log2(x+1)` → z-score(cohort) → top-2500-variance → `LogisticRegression(C=1.0, class_weight="balanced")` → **5-fold StratifiedKFold CV-OOF** → F1-max threshold; cohort-matched score references.
- **Seeds**: `random_state=0` (bootstrap, StratifiedKFold), `random_state=1` (permutation). `GroupKFold` is deterministic (donor order). Bootstrap B=1000, permutation P=1000.

## Analysis scripts → outputs

| Script | Produces |
|---|---|
| `pipeline/train_predictor.py` | multimodal model + per-held-out board reports |
| `pipeline/eval_oof_metrics_v2.py` | `scratchpad/oof_metrics_v2.json` (nested, CIs, perm p, ladder) + `oof_scores_v2.json` |
| `pipeline/bench_matrix.py` | `deliverables/bulk_matrix.{json,tsv}` (bulk × 4 cohorts) |
| `pipeline/build_full_matrix.py` | `deliverables/MOSAIC-AML_model_x_cohort_matrix.pdf` + `full_matrix.tsv` |
| `pipeline/build_provability_figs.py` | `deliverables/provability_A/B/C.*` |
| `pipeline/vaf_stratified.py` | `deliverables/vaf_stratified.*` |
| `pipeline/feature_attribution.py` | `deliverables/feature_attribution.{json,tsv}` |
| `pipeline/calibration_dca.py` | `deliverables/calibration_reliability.*`, `decision_curves.*` |

## Metrics tables (machine-readable)

`scratchpad/oof_metrics_v2.json`, `deliverables/{bulk_matrix,full_matrix,vaf_stratified,feature_attribution}.*`,
`gui/validation_stats.json` (per-mutation AUROC + CI + perm p, served to the board).

## Environment

- Cluster `anaconda3-2020` **Python 3.12**, **numpy 2.4**, scikit-learn / scipy / pandas (compute-node env — importable only under LSF, not the login node's old glibc).
- Compute: LSF `bsub -q test`; every analysis runs on a compute node (login node cannot import numpy).

## To reproduce end-to-end

```bash
cd pipeline
bsub -q test -W 150 -M 32000 -R "rusage[mem=32000]" -o ev2.log python eval_oof_metrics_v2.py
bsub -q test -W 90  -M 24000 -R "rusage[mem=24000]" -o bm.log  python bench_matrix.py
bsub -q test -W 60  -M 24000 -R "rusage[mem=24000]" -o vaf.log python vaf_stratified.py
bsub -q test -W 20  -M 12000 -R "rusage[mem=12000]" -o fa.log  python feature_attribution.py
# figures render locally (matplotlib): build_full_matrix.py, build_provability_figs.py, calibration_dca.py
```

## Still to deposit (open-science)

- ☐ Public code repository + the two `.pkl` models + `holdout_samples.txt`.
- ☐ Per-center panel gene lists (see `genes_assayed_scaffold.md`).
