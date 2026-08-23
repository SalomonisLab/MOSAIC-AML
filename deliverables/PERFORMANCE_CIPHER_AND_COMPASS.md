# CIPHER-AML and COMPASS-AML — complete performance reference

Every number here is pulled from a saved artefact, named at the end of each section. Where two figures
exist for the same quantity, both are given with the reason they differ.

**One-line summary.** Both layers rank well and decide badly. CIPHER-AML reaches mean AUROC 0.908
(single-cell) / 0.830 (bulk) but pooled precision 0.336 at its deployed thresholds. COMPASS-AML reaches
0.775 apparent AUROC, 0.671 once the patient main effect is removed — against a *measured* assay ceiling
of 0.727, so ~92% of the recoverable signal is already taken.

---

# Part 1 — CIPHER-AML (mutation / driver-lesion prediction)

## 1.1 Two models, different jobs

| | single-cell multimodal | bulk variant-level (**primary caller**) |
|---|---|---|
| unit | 28 gene-level drivers | 50 variant-level categories |
| input | 89 cell states × 8 modality blocks | bulk RNA (or single-cell bulk-equivalent) |
| mean AUROC | **0.908** | **0.830** |
| trained on | atlas single-cell + 707 BeatAML bulk | BeatAML2 (n=707, WES) |
| role now | cytogenetics and driver context | every mutation call in a patient report |

### Single-cell multimodal, per-driver (28 drivers, donor-grouped CV)

| | AUROC | sensitivity | specificity | F1 |
|---|---|---|---|---|
| deployed | 0.889 | 0.494 | 0.950 | 0.428 |
| fused (all modalities) | **0.908** | 0.505 | 0.954 | 0.456 |
| fused (selective) | 0.902 | 0.497 | 0.955 | 0.453 |

Distribution: median 0.918, range 0.769–0.997. **17 of 28 drivers ≥ 0.90**, 6 at 0.85–0.90, 4 at
0.80–0.85, 1 below 0.80.

- Best: inv(16)/*CBFB-MYH11* 0.997, *KMT2A* 0.974, del(7) 0.970, *KMT2A*-rearrangement 0.965, *SF3B1* 0.959
- Worst: trisomy 8 **0.769**, *KRAS* 0.804, *FLT3* 0.826, *GATA2* 0.846, *PTPN11* 0.849

### Bulk variant-level caller, per-category (50 categories)

Mean CV AUROC **0.830**, median 0.858. **32 categories ≥ 0.80**, 8 at 0.70–0.80, 9 at 0.60–0.70, 1 below.

- Best: *U2AF1*_S34 1.000, *SF3B1*_K666 0.997, *SRSF2* 0.996, *SF3B1*_K700 0.983, NPM1 exon12 frameshift 0.981
- Worst: *CBL* 0.578, *DNMT3A*_nonR882 0.605, *RAD21* 0.652, *SUZ12* 0.660, *CREBBP* 0.665

## 1.2 Held-out performance — the numbers that matter

Pooled across three external single-cell cohorts, **60 specimens, 2,046 scored calls**:

| cohort | n | calls | sensitivity | specificity | precision | F1 | accuracy |
|---|---|---|---|---|---|---|---|
| sealed single-cell hold-out | 29 | 728 | 0.505 | 0.931 | 0.547 | 0.525 | 0.871 |
| Trumpp/Waclawiczek | 16 | 928 | 0.533 | 0.895 | 0.260 | 0.350 | 0.872 |
| GSE281087 (panel-honest) | 15 | 390 | 0.367 | 0.850 | **0.169** | 0.232 | 0.813 |
| **ALL held-out** | **60** | **2,046** | **0.492** | **0.898** | **0.336** | **0.399** | 0.860 |

**Performance degrades monotonically with externality.** That gradient is the honest expectation for
transfer, and the most external cohort is the most informative estimate of real-world behaviour.

**Drivers actually recovered** (gene-level, the informative half of the old "N/M correct" pill):
Leucegene **210/364 = 0.577**, sealed single-cell **47/97 = 0.485**, Trumpp **30/44 = 0.682**.

## 1.3 Detection is clonality-limited

Sensitivity by variant allele fraction (bulk caller, 1,312 positive calls):

| VAF stratum | positives | called | sensitivity |
|---|---|---|---|
| subclonal < 0.10 | 31 | 9 | **0.290** |
| 0.10 – 0.25 | 212 | 94 | 0.443 |
| 0.25 – 0.40 | 399 | 258 | 0.647 |
| clonal ≥ 0.40 | 670 | 460 | **0.687** |

Per-gene sensitivity tracks median VAF: *NPM1* 0.957 (VAF 0.295), *SRSF2* 0.917, *SF3B1* 0.857,
*U2AF1* 0.818, *IDH1* 0.816, *DNMT3A* 0.793 — versus *FLT3* 0.304, *NRAS* 0.379, *GATA2* 0.414.

## 1.4 Threshold placement, not ranking, is the binding constraint

Mean AUROC 0.908 coexists with pooled precision 0.336. Fitting per-category thresholds on half a cohort
and transferring them to the other half:

| cohort | deployed F1 | cohort-percentile | prevalence-matched | **split-half oracle** | in-sample oracle |
|---|---|---|---|---|---|
| sealed single-cell | 0.525 | 0.455 | 0.500 | **0.569** | 0.755 |
| Trumpp | 0.418 | 0.416 | 0.380 | **0.520** | 0.723 |
| GSE281087 | 0.250 | 0.296 | 0.319 | **0.407** | 0.694 |
| Leucegene | 0.558 | 0.457 | 0.432 | **0.586** | 0.667 |

Two things follow. Threshold recalibration is worth **+0.03 to +0.16 F1**, biggest where the caller is
worst. And **label-free rules do not capture it** — both are worse than the shipped thresholds on three
of four cohorts. The in-sample oracle (0.67–0.76) is inflated: up to 29 free parameters on as few as
390 calls.

## 1.5 Where the modality gain comes from

- Bulk-trained → single-cell transfer collapses: **0.664 → 0.673** (RNA-only vs multimodal, both far
  below the 0.908 achieved when trained natively). Assay-matched modelling is required.
- Bulk **augmentation** of single-cell training helps: atlas-only 0.830 → atlas + 707 BeatAML **0.852**.
- Imputed modalities add **nothing beyond a matched random nonlinear control**: rna_only 0.802,
  rna + imputed 0.860, rna + *random* 0.865. The imputed blocks supply capacity, not information.

## 1.6 Reference-mismatch failure mode

Scores are percentiled against a cohort-matched reference. A specimen unlike any available reference
inherits that reference's offset. On GSE281087 (CD34-sorted), **54 of 65 positive calls were wrong**.

The z-scoring supplies its own diagnostic: a specimen the reference describes has **mean |z| ≈ 0.8** by
construction. Measured — BeatAML vs `beataml` 0.739, atlas single-cell vs `sc` 0.770, **GSE281087 vs
`sc` 21.25**, atlas single-cell vs the wrong (`beataml`) reference 5.35. A guard now fires above 2.0.
Within-cohort ranking on GSE281087 still reaches mean AUROC 0.688, so the information survives; the
threshold does not.

*Sources: `production_fused_model.json`, `pooled_heldout_eval.json`, `vaf_stratified.json`,
`full_matrix.tsv`, `exp_cipher_thresholds.json`, `beataml_{transfer,augment,impute}_experiment.json`,
`exp_ood_guard.json`, `bulk_mutation_predictor.pkl`.*

---

# Part 2 — COMPASS-AML (ex-vivo drug response)

## 2.1 Cohort and headline

**520 specimens × 118 inhibitors, 479 subjects, 48,998 measurements, 79.9% of the matrix observed.**
Donor-grouped 5-fold CV, tail classes at the 20th/80th percentile.

| metric | out-of-fold | sealed hold-out |
|---|---|---|
| mean per-drug Spearman | 0.369 | 0.372 |
| mean per-drug AUROC | **0.775** | 0.788 |
| mean per-drug AUPRC | 0.750 | 0.824 |

## 2.2 The headline is inflated — and by exactly how much

Two-way decomposition of the response matrix:

| component | share of variance |
|---|---|
| patient main effect (some specimens die readily in culture) | **15.4%** |
| drug main effect | 45.9% |
| **patient × drug interaction** (the only part that is a recommendation) | 46.8% |

The model's mean prediction per specimen correlates **0.540** with the patient main effect. Re-scored
against the interaction term alone:

| metric | as reported | **interaction only** |
|---|---|---|
| mean per-drug Spearman | 0.365 | **0.223** |
| mean per-drug AUROC | 0.775 | **0.671** |

**0.671 is the figure to quote for drug-specific skill.**

## 2.3 The assay is the ceiling, not the model

Estimated from technical replicates: median per-drug reliability **0.529** over 27 drugs with enough
replicates → ceiling √0.529 = **0.727**. Replicates of the same well differ by a median **13.7
percentage points of viability**. COMPASS at 0.671 is **~92% of the recoverable signal**.

Two independent confirmations:

- **Predictability tracks reproducibility**: Spearman **0.288, P = 0.0017** across 116 inhibitors, with
  a monotone gradient by reliability tier — reproducible 0.315 → moderate 0.231 → weak 0.200 →
  unreliable 0.184.
- **Predictability does NOT track training-set size**: Spearman 0.112, **P = 0.23**. The weak inhibitors
  are *hard*, not under-trained.

Reliability tiers across 118 inhibitors: **17 reproducible, 32 moderate, 29 weak, 38 unreliable,
2 unmeasurable** (median reliability 0.274, median implied ceiling 0.520).

## 2.4 Deployed matrix-factorisation blend

Rank-40 soft-impute factorisation with a features → patient-factors regression, blended 50/50:

| arm | interaction AUROC | interaction Spearman |
|---|---|---|
| deployed (per-family only) | 0.671 | 0.221 |
| MF rank 40 alone | 0.680 | 0.232 |
| **blend 50/50, rank 40 (deployed)** | **0.682** | **0.234** |
| blend, rank 60 | 0.682 | 0.234 |

Blend-weight sweep is flat between 0.5 and 0.75 (0.682–0.683), so 0.5 is used rather than chasing 0.001
on the sweep that selected it. The raw-target headline is unchanged (0.772), as expected.

## 2.5 Block ablations

| blocks | OOF AUROC | Spearman |
|---|---|---|
| clinical only | 0.579 | 0.106 |
| mutations only | 0.641 | 0.185 |
| cell state only | 0.658 | 0.211 |
| state + mut + clin (no RNA) | 0.707 | 0.273 |
| **RNA only** | **0.766** | 0.354 |
| RNA + state | 0.763 | 0.350 |
| RNA + mut | 0.773 | 0.364 |
| **all four fused** | **0.772** | 0.363 |

**RNA carries almost all of it.** Fusion adds ~0.006 over RNA alone; every non-RNA block alone is far
weaker. This is why the fusion is floored to the best single block.

## 2.6 Per-patient ranking — the task the product actually performs

468 specimens: mean per-patient Spearman 0.309.

| | hit@k | chance | precision@k |
|---|---|---|---|
| top-1 | **0.336** | 0.100 | 0.336 |
| top-3 | 0.618 | 0.274 | 0.294 |
| top-5 | 0.763 | 0.416 | 0.277 |
| top-10 | 0.925 | 0.669 | 0.246 |

**3.4× chance at rank 1.** This is the more trustworthy figure, and it is what the deployed system uses.

## 2.7 Calibration, abstention and clinically actionable subsets

**Calibration** (n = 18,919): ECE **0.0117**, Brier 0.185 against a 0.250 baseline.

**Abstention** — error falls monotonically as low-confidence calls are withheld:

| coverage | error rate | AUROC |
|---|---|---|
| 100% | 0.283 | 0.794 |
| 75% | 0.228 | 0.835 |
| 50% | 0.167 | 0.879 |
| 25% | 0.097 | 0.927 |

**Clinically actionable subset** (42 approved agents): mean Spearman **0.418**, mean AUROC **0.809** —
materially better than the all-118 average. Venetoclax: Spearman 0.766, **AUROC 0.977** (n = 363).
Quizartinib 0.869.

**By cell-state stratum**: primitive 0.762 AUROC, intermediate 0.752, monocytic/mature 0.721.

## 2.8 The clinical translation is null

| drug received | n scored | complete-response AUROC | P | OS C-index age+ELN → +model |
|---|---|---|---|---|
| **Cytarabine** | 149 | **0.435** (n=131, 93 CR) | 0.249 | 0.739 → 0.743 |
| Azacitidine | 50 | 0.550 (n=40) | 0.595 | 0.669 → 0.694 |
| Decitabine / sorafenib / midostaurin / venetoclax | 16–25 | underpowered | — | — |

**Predicted ex-vivo sensitivity does not translate into clinical response in this cohort.** Cytarabine
is on the wrong side of chance. Caveats that do not rescue it: cytarabine is wave-unstable, single-agent
ex-vivo cytarabine at 72 h is a poor model of 7+3 combination induction, and n = 131 is limited power.

*Sources: `drug_model_card.json`, `drug_model_validation.json`, `exp_compass_deep.json`,
`drug_assay_reliability.tsv`, `drug_model_card_abl_*.json`, `EXPERIMENT_RESULTS.md`.*

---

# Part 3 — What can be done to improve them

Ordered by measured expected value, not by appeal.

## 3.1 CIPHER-AML

**A. Per-cohort threshold calibration — largest measured gain in the platform.**
Measured: +0.03 to +0.16 F1, precision 0.182 → 0.458 on the worst cohort. Requires labelled specimens
from the target cohort; `calibrate_to_cohort.py` already does it (on GSE281087, fitting on 7 specimens
raised held-out F1 0.240 → 0.333). **This is a deployment protocol decision, not a modelling problem** —
it means asking each new site for 10–20 sequenced specimens. *Expected: large. Cost: organisational.*

**B. Raise sensitivity on the low-VAF tail.** Sensitivity is 0.290 below 10% VAF against 0.687 above
40%. Options that have not been tried: a VAF-aware target (train to predict clonal lesions only, and
report subclonal detection as out of scope), or explicitly modelling clonal fraction as a covariate so
the caller can abstain rather than miss. *Expected: moderate, and it converts a silent failure into a
declared one. Cost: 2–3 days.*

**C. Fix the categories that are structurally weak, or retire them.** *CBL* 0.578, *DNMT3A*_nonR882
0.605, *RAD21* 0.652 — nine categories sit below 0.70 CV AUROC. A caller that reports a 0.58-AUROC
category alongside a 0.99 one invites equal weight being given to both. Either improve them or move them
to a declared "research-only, do not act" tier. *Expected: no metric change, real reduction in
misreading. Cost: 1 day.*

**D. Build score references for new assay types.** The GSE281087 failure was a reference that did not
describe the specimen. Adding a sorted/enriched-population reference would close that specific case, and
the mean|z| guard now detects the general one. *Expected: moderate on sorted cohorts. Cost: 2 days.*

**E. Do NOT bother with label-free recalibration.** Measured and rejected: within-cohort percentile
matching and prevalence matching are worse than the shipped thresholds on three of four cohorts.

## 3.2 COMPASS-AML

**A. Nothing computational moves the headline much — the ceiling is real.** 0.671 against 0.727 leaves
~8% of recoverable signal. The rank-40 blend already took +0.011 of it. Anyone proposing a new model
here should first be shown §2.3.

**B. The one thing that would actually raise the ceiling: more replicate wells.** Only 15.4% of
concentration cells currently have more than one well. Replicate-averaging failed as a post-hoc fix
(reliability 0.192 → 0.168) precisely because the replicates mostly do not exist. This is a wet-lab
request, not an analysis. *Expected: raises the ceiling itself. Cost: assay time.*

**C. Report the actionable subset separately, and lead with it.** 42 approved agents reach mean AUROC
0.809 and Spearman 0.418, versus 0.775/0.369 across all 118. Venetoclax alone is 0.977. The product
recommends from the approved tier anyway, so the all-118 average understates the deployed use case.
*Expected: no model change, a fairer headline. Cost: hours.*

**D. Use abstention as a product feature.** At 50% coverage the error rate halves (0.283 -> 0.167) and
AUROC rises to 0.879. A report that ranks 10 drugs and declines to score the other 108 is both more
accurate and more honest than one that scores all 118. *Expected: large practical gain. Cost: 1 day.*

**E. The clinical question is the real gap, and it needs data we do not have.** The cytarabine null
(0.435) bounds every claim this layer can make. Closing it requires a cohort with drug-matched outcomes
and adequate n per agent — ideally venetoclax-treated patients, where the ex-vivo assay is strongest
(AUROC 0.977) and where BeatAML has only 21. *Expected: transformative if positive, decisive if null.
Cost: a new cohort.*

**F. Drop or down-weight the 38 unreliable inhibitors.** Their mean interaction Spearman is 0.184 and
their assay does not reproduce. Scoring them at all implies a precision the measurement cannot support;
the utility score already penalises low reliability, but they still appear in rankings.

## 3.3 Cross-cutting

The same sentence describes both layers: **ranking is strong, thresholding is weak, and the external
gradient is real.** The highest-value work in both cases is not a better model — it is (i) per-cohort
calibration with a small labelled set, (ii) abstention instead of scoring everything, and (iii) honest
tiering so a 0.58-AUROC category or a non-reproducing inhibitor cannot be read like a 0.98 one.
