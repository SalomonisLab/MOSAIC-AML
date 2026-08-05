# MOSAIC-AML — consolidated answers

**Read this first:** several numbers below **supersede** what was sent earlier. Two things changed the
results after those first answers were written:

1. **The NYU-2 metadata correction** — 149 mutation positives had been recorded as wild-type.
2. **BeatAML augmentation** — now built, validated, and deployed (your suggestion; it works).

A "supersedes" note marks every figure that has moved. §6 lists them in one place.

---

## 0. The architecture, in one picture

There are **two models**, one per input type. This resolves most of the confusion.

```
BULK RNA input                        SINGLE-CELL input
(BeatAML / Leucegene / any bulk)      (the atlas, or an uploaded scRNA sample)
        │                                     │
        ▼                                     ▼
  Bulk RNA-only caller              MOSAIC-AML multimodal caller
  trained on BeatAML2 (n=707)       trained on the single-cell atlas
                                    + (NEW) BeatAML-augmented shared blocks
```

- **BeatAML's role** = the training corpus for the bulk caller — and, as of this round, an
  **augmentation source** for the multimodal caller too.
- **BeatAML is never "imputed onto" itself as a source of new information.** The imputers matter
  because they give bulk and single-cell a **shared feature space** (see §5.2).

---

## 1. The five primary questions

**Q1. You mentioned 92.86% accuracy. Which result is that — single-cell held-out, BeatAML, or Leucegene?**

**Single-cell held-out.** It is **39/42 = 92.86%**, the *gene-level* score of one individual held-out
single-cell sample (8 of the 29 sealed held-out samples scored exactly 39/42). Not BeatAML, not
Leucegene. It reads high because gene-level accuracy counts every correctly-called-*absent* gene.

**Q2. Is bulk RNA-seq alone far more accurate than modality imputation?**

**No — the opposite.** On single-cell data the multimodal model is better on every metric except
sensitivity, which ties:

| | Bulk RNA alone | Multimodal |
|---|---|---|
| Accuracy | 83.1% | **94.1%** |
| AUROC | 0.71 | **0.88** |
| F1 | 0.26 | **0.44** |
| Specificity | 0.85 | **0.96** |
| Sensitivity | 0.51 | 0.51 (tied) |
| False positives | 1,272 | **325** |

**Q3. Does MOSAIC-AML multimodal use BeatAML + single-cell, or just BeatAML?**

**At the time you asked: neither — it used only the single-cell atlas, zero BeatAML.** All 8 modality
blocks were single-cell-derived. **This has now changed at your suggestion** — see §5.3. Note also that
the "bulk RNA-seq alone" bar in the comparison PDFs is the *single-cell atlas collapsed to pseudobulk*,
**not** BeatAML.

**Q4. Are the statistics for the held-out 50 single-cell samples or all 350? Is 50 correct?**

Those figures were **all ~350 (n = 357 labeled)**, by donor-grouped cross-validated out-of-fold
prediction — not a held-out subset. **50 was not correct**; the sealed held-out set is **29**.
*Superseded:* we now also report a **pooled held-out of 60 samples** (§3).

**Q5. inv(16) CBFB-MYH11 — bulk or single-cell? It says TP=2, FN=2, but BeatAML has 50 such samples.**

**Single-cell.** Three different numbers were being compared:

| Number | What it is |
|---|---|
| **11** | inv(16) samples in the single-cell atlas — the delivered file: **TP=11, FP=4, FN=0, sensitivity 1.000** |
| **4** | inv(16) samples in the 29-sample sealed held-out *board* — model caught 2, missed 2 (**your TP=2/FN=2**) |
| **~50** | inv(16) samples in **BeatAML** (bulk, n=707) — a dataset that figure doesn't use |

inv(16) is one of our **best**-captured events (CV AUROC 0.997).

---

## 2. Model × cohort performance (the eight secondary questions)

### Bulk RNA-only model (BeatAML2-trained)

| Cohort | Overall sensitivity | Overall specificity | AUROC |
|---|---|---|---|
| **BeatAML** (5-fold CV, n=707, 50 categories) | 0.53 | 0.96 | 0.83 |
| **Leucegene** (external bulk, n=367, 18 categories) | 0.55 | 0.96 | 0.855 |
| **All scRNA** (n=387, 14 categories) | **0.30** | 0.95 | **0.695** |
| Held-out scRNA | *superseded — see §3* | | |

Per-mutation highlights — BeatAML: NPM1 0.96/0.96, DNMT3A 0.88/0.90, RUNX1 0.87/0.94, KMT2A 0.84,
IDH1 0.82; weak KRAS 0.53/0.79, PTPN11 0.35, FLT3-ITD 0.31. Leucegene: NPM1 0.87/0.99, CEBPA 0.73,
IDH2 0.67; FLT3-ITD 0.00 (a truth-definition mismatch, not purely a model failure), TET2 0.23.

**The headline here:** the bulk model **collapses on single-cell** (AUROC 0.83 → 0.695). It cannot be
substituted for a single-cell-native model.

### MOSAIC-AML multimodal model

| Cohort | Overall sensitivity | Overall specificity | AUROC |
|---|---|---|---|
| **All scRNA** (donor-grouped CV, 28 mutations) | 0.505 | 0.954 | **0.908** |
| **Held-out scRNA** (29 samples, 22 mutations) | 0.505 | 0.931 | — |
| **BeatAML / Leucegene** | falls back to the bulk caller → **0.53 / 0.55 sens, 0.96 spec** | | |

On a bulk cohort there are no single-cell modalities to add, so the multimodal *system* operates as its
bulk caller — its BeatAML and Leucegene performance **is** the bulk model's row above.

---

## 3. "Why only 3 mutations? The held-out should be ~50 samples."

**You were right, and this was a reporting flaw on our side.** That analysis required **≥3 positives AND
≥3 negatives per mutation** before scoring it — reasonable for a stable *per-mutation* estimate, wrong
for an *overall* number, since a mutation seen once still contributes a true positive or a false negative.

Removing the filter and pooling every scored call:

| Cohort | Samples | Mutations with ≥1 positive | Sensitivity | Specificity |
|---|---|---|---|---|
| Sealed held-out scRNA | 29 | **22** (not 3) | 0.505 | 0.931 |
| Trumpp / Waclawiczek | 16 | 17 | 0.533 | 0.895 |
| GSE281087 (panel-honest) | 15 | 17 | 0.367 | 0.850 |
| **All pooled** | **60** | — | **0.492** | **0.899** |

2,046 scored calls; TP 95, FP 188, FN 98, TN 1,665 (precision 0.336, F1 0.399, accuracy 0.860).
**This is the number to quote for held-out performance.** GSE281087 is the weakest, which is expected
for a fully external cohort.

*On GSE281087 labels:* its panel only assayed 24 of 67 genes, so genes it **never assayed** are left
**unlabelled** rather than counted as wild-type. Counting them as wild-type would have inflated
specificity from 0.775 to 0.881 on that cohort.

---

## 4. "For pseudobulk classification you don't use BeatAML?"

**Correct — it didn't, and that was worth changing.** Three experiments now define exactly what BeatAML
can and cannot contribute.

---

## 5. The three BeatAML experiments

### 5.1 Does imputation add information *within* BeatAML? — **No**

Imputed all four modalities from BeatAML bulk RNA, then predicted BeatAML's own mutations (47 mutations,
5-fold CV):

| Arm | AUROC | Accuracy |
|---|---|---|
| RNA only | 0.802 | 0.921 |
| RNA + imputed | 0.860 | 0.939 |
| **RNA + random nonlinear (control)** | **0.865** | **0.937** |
| imputed − random | **−0.004 (p = 0.42, n.s.)** | +0.002 (n.s.) |

Imputation *looks* like it helps (+0.058) — but a **random nonlinear transform of the same RNA helps
equally**. It is nonlinear feature-expansion, not new biological information. The same lift is available
from a nonlinear RNA model with no imputation.

*(Also scored by accuracy at your request. Note accuracy is a poor metric here: at 5.6% prevalence the
trivial "call everything absent" classifier scores 0.944 — higher than every model.)*

### 5.2 Does BeatAML-with-imputation predict *single-cell* mutations better than RNA alone? — **Modestly**

Train on BeatAML (707 bulk), test on single-cell pseudobulk (34 mutations):

| Arm | AUROC | Sensitivity | Specificity | F1 |
|---|---|---|---|---|
| RNA only | 0.663 | 0.524 | 0.772 | 0.210 |
| **+ imputed modalities** | **0.672** | 0.513 | **0.814** | **0.277** |

**Yes, imputation helps the transfer — but modestly, and via specificity/precision, not sensitivity.**
Both arms sit far below a single-cell-trained model (0.908), so BeatAML cannot replace single-cell training.

### 5.3 Can BeatAML *augment* single-cell training? — **Yes, and this is the win**

This is the version of your idea that pays off. The imputers give bulk and single-cell the **same
feature names**, so 707 bulk samples can be pooled into the training side of the shared blocks
(BeatAML joins training only, never the test fold; each cohort z-scored on its own statistics).

| Model | AUROC | Sensitivity | Specificity | F1 |
|---|---|---|---|---|
| Deployed (8 atlas modalities) | 0.889 | 0.493 | 0.950 | 0.428 |
| **+ BeatAML augmentation** | **0.908** | 0.505 | 0.954 | **0.456** |

**22/28 mutations improved; Wilcoxon p = 0.016.** Gains concentrate exactly where single-cell positives
are scarce: **WT1 +0.132**, KRAS +0.087, SRSF2 +0.062, CEBPA +0.058, IDH2 +0.056, NRAS +0.053,
IDH1 +0.044, DNMT3A +0.034. Small residual losses: GATA2 −0.071 (only 8 atlas positives), RUNX1 −0.033.

**This is now deployed.** Poolable blocks: BulkRNA (14,237 features), ADT (129), Metabolite (995),
GRN (7,486). Lipid could not be pooled — the atlas used the `newnormelastic` imputer bundle, which is
missing from the cluster, so its lipid species don't correspond to the available bundle. **If that
bundle can be located, Lipid becomes poolable too.**

**On including VAF:** VAF can't be an *input feature* — it comes from the sequencing the model is meant
to predict, so the model would need the answer to produce the answer. Implemented instead as
**VAF-weighted training** (BeatAML positives weighted by clonality), which is legitimate and gave the
**best sensitivity (0.536)**. Consistent with our VAF-stratified finding: sensitivity 0.29 for subclonal
(<10% VAF) rising to 0.69 for clonal (≥40%).

---

## 6. What changed since the first answers

| Metric | Originally sent | Now | Why |
|---|---|---|---|
| Multimodal AUROC (all scRNA) | 0.885 | **0.908** | NYU-2 label fix (+0.004), BeatAML augmentation (+0.019) |
| Held-out evaluation | 29 samples, 3 mutations | **60 samples, 22/17/17 mutations** | Removed the wrong ≥3-positive filter; added Trumpp + GSE281087 |
| Held-out sensitivity / specificity | 0.61 / 0.93 (3 mutations) | **0.492 / 0.899 (pooled, 60 samples)** | Properly pooled, many more mutations counted |
| Trainable mutations | 26 | **28** | SF3B1 and SRSF2 became trainable |
| Does BeatAML feed the multimodal model? | No | **Yes** | Your suggestion, built and validated |

**The NYU-2 correction:** all 28 NYU-2 samples had been recorded as wild-type for every gene, while the
harmonized v2 metadata supplies **149 positives** — which were being used as *true negatives*. 25 of the
28 have atlas expression. Effect: **IDH2 +0.120 AUROC** (its positives went 14 → 35), SF3B1 and SRSF2
became trainable, and F1 rose across every arm. Leakage check came back clean — the six samples
appearing under both NYU-1 and NYU-2 share a Donor_ID and already land in the same CV fold.

---

## 7. Figures

| File | What it shows |
|---|---|
| `sens_spec_panels.pdf` | **Your requested graphic** — one dot per mutation, sensitivity × specificity, panels per model × dataset |
| `precision_recall_panels.pdf` | Same in the style of the journal figure you sent; precision discriminates far better on rare mutations |
| `MOSAIC-AML_model_x_cohort_matrix.pdf` | Two models × four cohorts, overall + per-mutation |
| `beataml_impute_experiment.pdf` | §5.1 — imputation vs the random control |
| `provability_A_independence_ladder.pdf` | 69% of the multimodal gain needs no imputed modalities |
| `provability_B_auroc_forest.pdf` | Every mutation's AUROC + 95% CI vs chance |
| `provability_C_nested_operating_point.pdf` | Honest (nested-CV) vs optimistic operating point |
| `vaf_stratified.pdf` | Sensitivity rises with clonality — misses are subclonal |
| `calibration_reliability.pdf` | Calibration: ECE 0.454 → 0.014 |
| `decision_curves.pdf` | Clinical net benefit vs treat-all / treat-none |

A note on the sensitivity-vs-specificity graphic: **specificity saturates near 1** for rare mutations
(every model correctly calls most samples absent), so those panels compress into a thin band. The
precision-recall version spreads the data out and is the better choice for a manuscript figure.

---

## 8. Open items

- **The LLM gateway API key must be rotated.** It was committed to the public repo; it has been removed
  from the code, but removal does not invalidate an already-exposed key.
- **The `newnormelastic` lipid imputer bundle** is missing from the cluster — locating it would add a
  fifth poolable modality.
- **EBI (62 samples) and NIH (4)** still have no mutation calls; and the harmonized file's own dictionary
  notes `0 = absent OR not reported`, so other cohorts may carry the same hidden defect NYU-2 revealed.
- **Prospective validation** — re-sequencing a handful of model-positive / label-absent samples is the
  one thing retrospective cross-validation cannot substitute for.
