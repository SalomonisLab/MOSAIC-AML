# Cohort adapters

A cohort adapter lets each mutation category use whichever model was better on one cohort. Routed
categories are called by a model fitted on that cohort. Every other category keeps the shipped
BeatAML-trained call, unchanged. Code: `amlmm/cohort_adapter.py`. Training and evaluation:
`train_cohort_adapter.py`.

```bash
python ingest_patient.py --bulk sample.tsv --bulk-ref sc --name X --cohort-adapter cohort_adapters/32266596.json.gz
```

It is opt-in, and it refuses samples that don't look like its cohort; those keep every shipped call.
The refusal fires on gene coverage below 90% or on an expression profile that correlates poorly with
the cohort's. It rejected all 64 10x atlas MDS samples.

## 32266596 (BD Rhapsody WTA, 184 AML)

| file | trained on | how it is validated |
|---|---|---|
| `32266596.json.gz` (default) | 129 patients | 55 sealed patients it never saw (IDs stored in the file) + nested CV |
| `32266596_alldata.json.gz` | all 184 | nested CV only |

**Routing rule.** Within training folds only, the cohort model and the shipped model are scored on the
same held-out patients. A category is routed if the cohort model gains ≥0.03 AUROC on average, is better
in ≥70% of folds, and has F1 no worse. A new category (cytogenetics the shipped caller lacks) is added
if its cross-validated AUROC is ≥0.75 and stable across folds.

**Measured on patients the adapter never trained on** (pooled over the routed-eligible units):

| | shipped | routed |
|---|---|---|
| sealed 55: F1 / sensitivity / precision | 0.273 / 0.49 / 0.19 | **0.332 / 0.59 / 0.23** |
| sealed 55: new cytogenetic calls | — | 11 correct, 5 false, 3 missed |
| nested CV, all 184: F1 | 0.299 | 0.368 |

**Known limits.**
- In nested CV, routing *every* category to the cohort model scored higher still (F1 0.398 vs
  0.368). At 130–184 patients, which categories get worse is not stable from sample to sample, so the
  conservative rule gives up some gain to avoid regressions.
- Categories with no labeled truth in the cohort (KRAS, SF3B1, CREBBP, …) cannot be routed. They keep
  the shipped over-calling, so healthy controls still receive ~12 "present" calls.
- `32266596_alldata.json.gz` routes `+8` (trisomy 8), which scored 0.51 in nested CV. Treat that category
  as unvalidated.
- The files contain derived statistics from dataset 32266596: per-gene means/SDs, coefficients,
  out-of-fold score distributions, and the sealed patient IDs. Check that dataset's data-use terms
  before committing them to a public repository.

## Custom per-mutation models (format v2)

| file | fitted on | use |
|---|---|---|
| `32266596_custom.json.gz` | all 184 AML patients | the one to use |
| `32266596_custom_dev.json.gz` | 129 development patients | the blend recipe as scored once on the 55 sealed patients |

Two layers, both applied by `CohortAdapter.apply`:

1. Every shipped category is re-referenced to the cohort (z-scored against the cohort's own per-gene
   mean/sd, percentiled against the cohort's score distribution). No labels are used.
2. One custom model per mutation: an equal-weight ensemble of linear SVM, PLS and logistic models
   (chromosome-arm dosage for -7 / +8 / complex karyotype, cohort RNA plus the 50 re-referenced
   BeatAML scores for everything else), blended 1:1 with the shipped model's own re-referenced score
   where the shipped caller has that category. Calls use a prevalence-matched cut.

### No mutation on a worse model than before

The all-184 file serves the custom recipe by default and switches a mutation to a different model
generation only where that generation beat the custom model by >= 0.03 AUROC in >= 80% of patient
bootstraps, on patients neither model was fitted on. Four mutations switch:

- TET2, IDH2: the re-referenced BeatAML score alone (training on this cohort adds nothing there and
  dragged the blend down: 0.68 vs 0.74 and 0.68 vs 0.76).
- TP53: the cohort-trained model alone (0.90 vs 0.87 for the blend).
- SRSF2: the BeatAML SRSF2 model fine-tuned on this cohort (0.92 vs 0.87; see below).

U2AF1 and KMT2A-r (5 positives each) stay on the re-referenced shipped call. All 7 lesions the
shipped caller cannot call are served (held-out AUROC 0.81-0.99). Each unit records which model
serves it (`served`), why (`reason`) and every generation's held-out score (`held_out_auroc_mcc`).

Held-out results on all 184 patients (5-fold x 2 nested CV, whole runs held out):

| | shipped | last round's adapter | re-referencing alone | custom recipe |
|---|---|---|---|---|
| mean AUROC, 22 shipped-covered mutations | 0.669 | 0.723 | 0.785 | 0.799 |
| pooled F1 of the calls | 0.280 | 0.362 | 0.468 | 0.526 |
| pooled precision | 0.187 | 0.257 | 0.428 | 0.507 |

- Served model vs the shipped caller: higher in 22 of 22 mutations.
- Served model vs the better of shipped / last round's cohort model: higher in 25 of 29. NPM1, IDH1,
  t(15;17) and inv(16) are 0.002-0.013 lower, which is within chance.
- Why not simply take the top-scoring model per mutation: tested, and it is worse on held-out
  patients (0.774 vs 0.799). A strict "keep the old model unless the new one wins on the training
  patients" rule left more mutations below their predecessor on held-out patients (7 of 29, up to
  -0.03) than not applying it (5 of 29, up to -0.02). With 5-50 mutated patients per mutation, a
  0.02 difference between two models is a handful of patient pairs.

### Searching for the best model these data supports

The principled near-optimal method was also tried for every mutation the shipped caller covers:
treat the BeatAML model as a prior and let this cohort update it only as far as its data supports
(maximum a posteriori logistic model, penalty centered on the BeatAML model, prior strength chosen by
cross-validation per mutation or shared across mutations). On held-out patients it averaged 0.76
against 0.79 for the served models, and was lower for 17 of 19 mutations. Even choosing the best
prior strength per mutation in hindsight, after seeing the held-out answers, reaches only 0.80. The
served models are therefore within about 0.01 of what these model families can get from 184
patients; the limit is the number of patients, not the search. The learning curve gains about 0.03
AUROC per doubling of patients and has not flattened. SRSF2 was the one clear exception and now uses
the fine-tuned model.

Sealed 55 patients, development-fit blend, 14 shipped-covered mutations with >= 3 positives
(mean AUROC / pooled F1 / precision): shipped 0.729 / 0.251 / 0.164; v1 routed adapter
0.750 / 0.305 / 0.202; re-referencing alone 0.795 / 0.427 / 0.353; custom 0.803 / 0.483 / 0.406.
A sealed number and a cross-validation number for the same mutation are the same model measured on
different patients: RUNX1 reads 0.88 on the development patients and 0.71 on the sealed ones, and the
untrained re-referenced BeatAML score shows the same gap (0.87 vs 0.59).

Caveats: no untouched patients remain (the sealed set is spent and the three switches were chosen on
the all-184 held-out scores, so their numbers are slightly optimistic); calls are prevalence-matched
to this cohort; the platform guard refuses 7 of 184 in-cohort samples and all 64 real 10x samples in
`labels/mds_bulk_equiv.npz` (their best profile correlation is 0.838 against a limit of 0.85, a thin
margin); the files embed statistics derived from dataset 32266596, so check its data-use terms before
publishing them.

### 2026-10-14: highest-scoring model per mutation

At the user's request every mutation now uses the model with the highest held-out AUROC among all
generations tried, even where the margin is within chance. Changed: NPM1, t(15;17), inv(16) -> last
round's cohort model; BCOR, IDH1, WT1, PTPN11, FLT3-TKD -> BeatAML model fine-tuned on this cohort;
STAG2, NRAS, KRAS, CEBPA -> re-referenced shipped call; U2AF1 -> re-referenced score with a cohort cut.
Gains are 0.001-0.043 AUROC. These were chosen on the same held-out scores they are reported with,
so the small ones may not hold on new patients. Each unit's `served` and `reason` fields say which
model it uses and why.

### 2026-10-14: BeatAML2 added as training data

BeatAML2 (Bottomly et al. 2022, github.com/biodev/beataml2.0_data; 596 AML patients, one RNA
specimen each, labels from clinical NPM1/FLT3-ITD calls, exome mutation calls, fusions and
karyotype) was trained together with this cohort, each cohort standardized against itself. Scored
on held-out patients of this cohort (same folds as above), it beat the served model for 12 of 29
mutations, which now use it: TET2, ASXL1, WT1, PTPN11, FLT3-ITD, CEBPA-bZIP, KMT2A-r, KRAS, DNMT3A,
t(8;21), inv(3), +8 (gains 0.005-0.057). For the other 17 the earlier model stayed better.
Mean held-out AUROC over the 22 shipped-covered mutations is now 0.826 (start: 0.669). The joint
models use an extra block `rna_joint` (this cohort's 2,500 most variable genes that BeatAML2 also
measures). The raw BeatAML2 files are in Downloads\beataml2 and are not part of this repository;
check BeatAML2's data-use terms before publishing models trained on it.

### 2026-10-14: choosing BeatAML2 patients one by one

Instead of the whole BeatAML2 dataset, each mutation's training set was also built from individually
weighted or selected BeatAML2 patients: (a) weights from how closely each patient's expression
resembles this cohort (label-free density ratio), (b) KNN-Shapley data valuation of each patient's
contribution to predicting the mutation in this cohort's training patients, dropping the least
useful quarter, and (c) both. Selection was always done inside the training folds. On average it was
not better than using the whole dataset (0.80-0.82 vs 0.82), but for 14 mutations the best selected
version scored higher on held-out patients and is now used: KIT (+0.031), complex karyotype (+0.031),
-7 (+0.011), t(8;21), DNMT3A, PTPN11, FLT3-ITD, KRAS, WT1, TP53, t(15;17), and ASXL1, CEBPA-bZIP, +8
by less than 0.001. Most of these margins are within chance.

### Predicted response to intensive induction (`response_models`)

`32266596_custom.json.gz` also carries a treatment-response model, applied by
`CohortAdapter.response()`. `ingest_patient.py --cohort-adapter` writes its output to the report's
`treatment_response` field, and the decision board shows it as "Predicted response to treatment".

- **Endpoint:** complete remission (CR or CRi) after first-line intensive induction chemotherapy.
- **Training:** the 134 intensive-induction patients of 32266596 with a recorded first-line response,
  plus 320 BeatAML2 patients sampled at initial diagnosis and treated with standard induction
  chemotherapy. Model: linear SVM + logistic regression on 2,415 genes, each cohort standardized against
  itself, Platt-calibrated on out-of-fold scores.
- **Held-out AUROC** (32266596, 5-fold x 2, whole sequencing runs held out): 0.842. On the same patients:
  - COMPASS predicted cytarabine sensitivity (what the drug layer offered before): 0.552
  - ELN 2022 from true genetics: 0.783
  - LSC17: 0.696
  - the same model trained on BeatAML2 alone, never shown a label from this cohort: 0.787
  - with age and true ELN added: 0.864
- **Calibration:** by tertile, predicted vs observed CR was 0.28 vs 0.24-0.29, 0.67 vs 0.61-0.66 and 0.88
  vs 0.93; Brier score 0.15-0.16 against 0.24 for the base rate.
- **Not applicable to:** APL (the board says so when t(15;17) is called), venetoclax or HMA regimens,
  and specimens the control gate calls healthy. The prediction does not see performance status,
  comorbidities or the exact regimen.

Fix: earlier versions of the two v2 files stored 51 NaN values in their evaluation notes. Reports built
with them were not valid JSON, so the board could not open those patients. The values are now null, and
the adapter replaces any non-finite value before it reaches a report.

### Overall survival calibrated to the cohort (`survival_models`)

`CohortAdapter.survival()` returns overall survival at 1, 2 and 3 years on this cohort's own scale.
The report field is `cohort_survival` and the board card is "Prognosis". Two ridge Cox models (2,327
genes, block `rna_joint3`) were trained on three cohorts: 184 patients of 32266596 (96 deaths), 440
BeatAML2 patients sampled at initial diagnosis (243 deaths), and 250 GSE106291 patients (German AML Cooperative Group (AMLCG) trial patients;
intensively treated, 150 deaths). The partial likelihood is stratified by cohort, and the baseline
survival curve is estimated on this cohort's out-of-fold risk scores. The genes are those of
`rna_joint` that GSE106291 also measures.

| Held-out C-index (32266596, 5-fold x 2, whole runs held out) | C |
|---|---|
| MOSAIC survival layer, RNA upload only | 0.608 |
| Two-cohort model (BeatAML2 + 32266596), expression only | 0.665 |
| Three-cohort model, expression only | 0.682 |
| Two-cohort model, expression + age | 0.700 |
| MOSAIC survival layer, age + ELN + labs | 0.703 |
| Three-cohort model, expression + age | 0.718 |

- **Three cohorts against two:** with age, +0.018 (95% CI +0.005 to +0.030, paired bootstrap over
  patients, 300 resamples, seed 0); expression only, +0.018 (-0.000 to +0.038). The three-cohort models
  replaced the two-cohort models on 2026-11-04. GSE106291 plus this cohort without BeatAML2 scored 0.724
  with age, with a wider interval (-0.010 to +0.056), so it was not chosen.
- **Calibration:** MOSAIC's general layer predicted about 5% 2-year survival on this cohort against 55%
  observed. The three-cohort model predicts 0.54-0.55 overall, and 0.31/0.57/0.75 against observed
  0.33/0.58/0.73 by risk tertile (expression only, repeat 1 of 2).
- **Against the general layer:** with age alone, the three-cohort model matches the general layer that
  also needs ELN risk and blood counts (0.718 against 0.703).
- **GSE106291 for mutations:** as a third training cohort for RUNX1 (46 mutated) and t(8;21) (7 fused),
  it lowered held-out AUROC (RUNX1 0.819-0.822 against 0.839; t(8;21) 0.932-0.933 against 0.945), so
  both keep their models.
- **Treatment-suggestion cut:** choosing each mutation's call cut for F1 (cross-fitted over 5 groups
  of sequencing runs) instead of matching prevalence lowered pooled suggestion F1 from 0.685 to 0.662,
  and a cut per drug class to 0.646, so the prevalence-matched cut stays.
- **Remission model not raised:** predicted genetics, age and a 5,000-gene variant were tested on top of
  the remission model; none improved it (0.831-0.843 against 0.842).

Survival-layer fixes in `survival_layer.py` and `amlmm/survival/data.py`:

- **Age-only input returned NaN.** Clinical input with age but no labs produced a NaN risk score: the
  one-patient clinical block filled a missing lab with the median of that one patient. Missing values are
  now filled with the training-cohort medians. `train_survival_model.py` stores these as `clin_medians`,
  with constants as a fallback for the current bundle.
- **`arm_note` misdescribed the model.** It said the expression-only model was used whenever the arm was
  `deployed`. It now names the arm used, the inputs supplied and the ones filled.

### Independent test on GSE106291

GSE106291 is the AMLCG cohort: 250 intensively treated AML patients with RNA-seq, 235 of them with an
induction response. The first remission recipe (linear SVM + logistic on BeatAML2 + 32266596) scored AUROC
0.563 on it, below age alone there (0.716).

### Remission model that transfers between cohorts (2026-11-04)

Every recipe was scored leave-one-cohort-out over the three cohorts with an induction response
(32266596, 134 patients; BeatAML2, 320; GSE106291, 235): trained on two, scored on the third. 88 recipes
were compared (gene sets, ridge strength, one model per cohort, age stacked on top); scripts r7, r7b, r8 and
r8b in `Downloads/32266596_experiments/scripts`.

- **Where the drop comes from:** 32266596 and BeatAML2 transfer to each other (AUROC 0.73-0.83 in both
  directions). Any model scored on GSE106291 reached 0.59-0.63, and GSE106291's own model reached
  0.62-0.68 on the other two. Within GSE106291, cross-validated, RNA reached 0.789. The cohorts have
  similar ages (median 56, 57 and 58 years) and CR rates (61%, 69%, 70%), so the difference is in how
  expression relates to response in GSE106291 (AMLCG regimens, resistance as its endpoint), not in the
  patient mix.
- **Served now:** one ridge logistic model (C=0.005) per cohort, each standardized on its own training
  scores, averaged, trained on all three cohorts (2,327 genes, block `rna_joint3`).

| Held-out AUROC, complete remission | Previous model | Served now |
|---|---|---|
| 32266596, 5-fold x 2, whole runs held out | 0.842 | 0.842 (difference +0.001, 95% CI -0.045 to +0.048) |
| 32266596 left out of training | 0.785 | 0.828 |
| BeatAML2 left out of training | 0.726 | 0.721 |
| GSE106291 left out of training | 0.563 | 0.641 |
| BeatAML2, its own patients held out | 0.684 | 0.728 |
| GSE106291, its own patients held out | 0.774 | 0.762 |

- **Still open:** on GSE106291 as an unseen cohort, no RNA model beat age alone (0.716). Adding age to the
  model raised GSE106291 to 0.659 left out and 0.786 held out, but lowered 32266596 to 0.828 (paired
  difference -0.013), where age does not predict response (0.487), so it is not served.
- **Calibration on 32266596 (out-of-fold):** predicted vs observed CR by tertile 0.26/0.29, 0.69/0.61,
  0.89/0.93; Brier 0.159 against 0.237 for the base rate.
