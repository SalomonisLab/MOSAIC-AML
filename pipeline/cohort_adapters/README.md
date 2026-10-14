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
- Categories with no labelled truth in the cohort (KRAS, SF3B1, CREBBP, …) cannot be routed. They keep
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

### Searching for the best model this data supports

The principled near-optimal method was also tried for every mutation the shipped caller covers:
treat the BeatAML model as a prior and let this cohort update it only as far as its data supports
(maximum a posteriori logistic model, penalty centred on the BeatAML model, prior strength chosen by
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
karyotype) was trained together with this cohort, each cohort standardised against itself. Scored
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
