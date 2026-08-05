# TRIPOD+AI reporting checklist — MOSAIC-AML mutation caller

Maps each TRIPOD+AI (2024) item to where it is addressed in this project. Status: ✅ done ·
�repartial · ☐ to-do. This is the scaffold for a manuscript methods section.

| # | TRIPOD+AI item | Where / status |
|---|---|---|
| 1–2 | Title / abstract identify it as a multivariable prediction model study | ☐ (manuscript) |
| 3a | Background & rationale | AML driver detection from expression; ✅ project docs |
| 3b | Objectives | Predict 58 variant-level drivers from RNA (bulk) or single-cell (multimodal) |
| 4a | Source of data (cohorts) | ✅ `truth_provenance.md` — BeatAML2, Leucegene, single-cell atlas, Trumpp, MDS |
| 4b | Study dates / setting | ◑ per contributing center; to compile |
| 5a | Eligibility criteria | Mutations with ≥8 pos / ≥8 neg (well-powered); ✅ in `eval_oof_metrics*.py` |
| 5b | Treatment/care details | n/a (retrospective genomic) |
| 6a | Outcome definition | Known driver genotype (variant- or gene-level); ✅ `truth_provenance.md` |
| 6b | Outcome assessment blinding | Labels from independent sequencing, not model-derived; ✅ |
| 7a | Predictors | 8 single-cell modalities / bulk RNA; ✅ `feature_attribution.*` lists top drivers |
| 7b | Predictor assessment blinding | Features built without touching labels (donor-grouped CV) ✅ |
| 8 | Sample size | ✅ n per mutation reported (n⁺ in every figure); CIs reflect it |
| 9 | Missing data | Modality blocks fillna(0); imputed modalities fidelity-filtered ✅ |
| 10a | Predictors handling | StandardScaler → differential top-500 → LinearSVC → percentile ✅ |
| 10b | Model type / building | Per-modality linSVM → ridge-NNLS fusion; ✅ `train_predictor.py` |
| 10c | **Internal validation** | **Donor-grouped CV-OOF + nested-CV threshold** ✅ `eval_oof_metrics_v2.py` |
| 10d | **Class imbalance / calibration** | Isotonic calibration + reliability/Brier/ECE ◑ `calibration_reliability.*` |
| 11 | Risk groups | Calls at F1-max (nested) threshold ✅ |
| 12 | Development vs validation | Bulk: train BeatAML / test Leucegene+sc ✅ `bulk_matrix.*` |
| **AI-1** | **Discrimination + uncertainty** | **AUROC with 95% bootstrap CI** ✅ `provability_B` |
| **AI-2** | **Signal vs chance** | **Per-mutation permutation null (p)** ✅ `provability_B`, board badges |
| **AI-3** | **Independence of inputs** | **Modality-ablation ladder (imputation-independent gain)** ✅ `provability_A` |
| **AI-4** | **Clinical utility** | Decision-curve / net benefit ◑ `decision_curves.*` |
| AI-5 | Failure-mode analysis | VAF-stratified misses (subclonal) ✅ `vaf_stratified.*` |
| AI-6 | Fairness / subgroups | By-cohort/center performance ◑ (matrix by cohort ✅; per-center to add) |
| AI-7 | Data provenance & quality | ✅ `truth_provenance.md`; label-gap (not-assayed) flagged |
| AI-8 | Code availability | ✅ pipeline scripts versioned; ☐ public deposit |
| AI-9 | Compute / reproducibility | ✅ `REPRODUCIBILITY.md` (seeds, splits, model cards, versions) |
| 13–17 | Results / model presentation | ✅ figures + TSVs in `deliverables/`; GUI Validation page |
| 18 | Limitations | ✅ threshold optimism (corrected), label not-assayed gap, bulk→sc shift, small held-out |
| 19 | Interpretation / generalizability | ✅ "match the model to the assay"; external cohorts |
| 20 | Supplementary / data access | ◑ deposit models + splits |

## The honest gaps (what a reviewer will still flag)

- ◑ **Calibration & decision curves** — computing now (`calibration_reliability.*`, `decision_curves.*`).
- ◑ **Per-center subgroup performance** (fairness) — the cohort matrix exists; a per-contributing-center breakdown is the next cut.
- ☐ **Prospective / orthogonal validation** — the one thing retrospective CV cannot substitute (re-sequence model-positive / label-absent samples).
- ☐ **Public code + model deposit** — assemble from `REPRODUCIBILITY.md`.
