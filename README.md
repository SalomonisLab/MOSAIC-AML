# COMPASS-AML

**C**ell-state **O**riented **M**odelling of **P**harmacologic **A**ssay **S**ensitivity —
the **drug predictor** of [MOSAIC-AML](https://github.com/SalomonisLab/MOSAIC-AML).

> This is a component branch. It carries the drug-response layer and the shared platform it runs on;
> the mutation predictor (CIPHER-AML) lives on the `cipher-aml` branch, and the complete system is on
> `main`.

## What it does

Predicts which small-molecule inhibitors a patient's leukaemia is sensitive to, learned from the
BeatAML2 **ex-vivo** functional screen: 48,998 dose–response curves over 520 specimens from 479
patients and 118 inhibitors, after curve QC and inclusion filtering.

Three deliberately **separate** models:

- **A — patient-level response.** Hierarchical: a ridge over every (specimen × inhibitor) row in a
  target-pathway family, plus a shrunk per-drug residual. Four feature blocks (RNA / differentiation
  state / mutations / clinical) late-fused by non-negative least squares on inner donor-grouped
  out-of-fold predictions, floored to the best single block.
- **B — state-resolved.** The same model re-applied to each cell state's pseudobulk, giving blast and
  LSC-like coverage, the most resistant state above 1% abundance, and bulk-versus-single-cell
  disagreement.
- **C — mechanistic.** Target expression, pathway output readouts, BCL2-family dependency, genetic
  activation and measurable resistance proxies — kept *out* of A so agreement and disagreement between
  an empirical and a mechanistic line of evidence stay informative.

Then a treatment-utility score with every penalty itemised (uncertainty, resistance burden,
infeasibility, out-of-distribution), rankings produced **per clinical tier**, and eight expert agents —
all non-voting with respect to the anchored subtype call.

## Validation

| | |
|---|---|
| per-inhibitor mean AUROC (118 drugs, donor-grouped CV) | **0.774** |
| approved-agent subset (n=42) | **0.809** |
| per-patient top-1 retrieval | **34%** vs 10% matched chance |
| calibration | ECE **0.012**, Brier 0.185 (baseline 0.249) |
| leave-wave-out / leave-centre-out | 0.72 / 0.73–0.89 |
| permutation null (100 shuffles) | observed 0.367 vs 0.001 ± 0.019 — **19 null SDs**, p = 0.0099 |
| sealed hold-out (15% of *patients*) | 0.784 |

Model B passed a pre-registered test: fitted only on BeatAML **bulk** and never shown a cell-state
label, it predicts higher venetoclax sensitivity in primitive than in monocytic states in **93%** of
387 atlas samples (p = 8×10⁻⁵¹) — **rank 1 of 118 inhibitors**, so not a generic artefact.

## Run it

Needs the BeatAML2 inhibitor screen in `data/external/beataml/`
(fetch from https://github.com/biodev/beataml2.0_data):

```
python pipeline/build_state_signatures.py    # lineage signatures from the atlas
python pipeline/train_drug_model.py          # Model A -> drug_response_model.pkl
python pipeline/build_drug_score_refs.py     # cohort-matched score references (REQUIRED)
python pipeline/eval_drug_model.py           # the validation battery
python pipeline/predict_drugs.py --atlas-sample "<Dataset::Sample>"
```

Per-patient view: `gui/therapy.html`. Validation page: `gui/rx_validation.html`.
Methods: `deliverables/METHODS_COMPASS-AML.md`.

**Standing caveat.** Ex-vivo sensitivity is an experimentally grounded prioritisation signal, **not** an
estimate of clinical benefit: culture conditions, pharmacokinetics, the marrow microenvironment,
combination therapy and toxicity are not represented. Nothing here is validated against patient
outcomes.
