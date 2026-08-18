# MOSAIC-AML — survival layer: complete results

Everything measured, including the results that do not favour the model. Reproduce any number with the
commands in §9.

---

## 1. Bottom line

| question | answer |
|---|---|
| Can it rank who dies sooner? | **Yes.** C-index **0.756 ± 0.029** across 60 hold-out re-draws (the single sealed split gave 0.787); **0.706** on TCGA-LAML with every coefficient frozen. |
| Does it beat age + ELN 2017? | **Yes, but only in combination.** +0.059 [+0.030, +0.088], P = 0.001. |
| Does molecular data beat age + ELN on its own? | **No.** −0.003 in BeatAML *and* −0.003 in TCGA, independently. |
| Calibrated probability of surviving to a horizon? | **Yes.** Mean absolute calibration gap 0.043–0.075. |
| How long will a *group* of similar patients live? | **Yes.** Median per risk tertile accurate to **6–7 weeks**. |
| How long will **one person** live? | **No.** MAE 0.87 y; honest 80% interval spans **22.6 years**. |
| Does it work in non-intensively-treated patients? | **No.** C-index **0.554** (n=71). Do not use it there. |

**Deployed output is a risk group, horizon probabilities and restricted mean survival — never a single
predicted number of months for one patient.**

> **Stability of the headline.** Re-drawing the sealed hold-out 60 times gives **0.756 ± 0.029** (5th–95th percentile 0.707–0.799); the single sealed split's 0.787 sits at the 86th percentile, so 0.756 is the honest expectation. The **gain** is the firmer claim: **+0.062 over age + ELN, positive on 100% of 60 draws**.

---

## 2. Cohort

BeatAML2, restricted three ways, each of which changes the answer:

| step | n |
|---|---|
| clinical rows with expression | 671 |
| drop unknown vital status (follow-up but no outcome — neither event nor censoring) | −22 → 649 |
| initial-diagnosis specimens only (a relapse specimen leaks prognosis: the patient already survived to relapse) | → 446 |
| one specimen per patient (earliest) | **444** |

**444 patients · 245 deaths (55%) · 199 censored · median follow-up of the living 2.26 y.**
Sealed hold-out: **89 patients, 51 deaths**, drawn once and never touched by any fitting, selection or
calibration step. Patients never span folds (5-fold, patient-grouped). The feature space — PCA,
z-reference, imputation medians, age-spline knots — is refit inside every fold.

Censoring is handled properly throughout: a patient last known alive contributes as censored at their
follow-up time, never as a death.

---

## 3. What each block is worth

Age and ELN 2017 are free and available on day one, so the number that matters is the **gain over
age + ELN**, not the C-index.

| arm | CV (out-of-fold) | sealed hold-out | AUC 2 y | Brier 2 y | Δ vs age+ELN (95% CI) | P(no gain) |
|---|---|---|---|---|---|---|
| ELN 2017 alone | 0.580 | 0.624 | 0.710 | 0.202 | — | — |
| cell state | 0.576 | 0.569 | 0.568 | 0.239 | −0.116 [−0.168, −0.069] | — |
| mutations | 0.626 | 0.579 | 0.672 | 0.217 | −0.066 [−0.112, −0.017] | — |
| age alone | 0.686 | 0.736 | 0.797 | 0.188 | — | — |
| bulk RNA | 0.690 | 0.715 | 0.816 | 0.181 | −0.002 [−0.041, +0.041] | — |
| molecular (RNA + state + mut) | 0.689 | 0.716 | 0.829 | 0.172 | **−0.003** [−0.042, +0.036] | 0.538 |
| **age + ELN 2017** | **0.692** | **0.725** | 0.834 | 0.164 | *baseline* | — |
| clinical block | 0.709 | 0.745 | 0.822 | 0.176 | +0.018 [−0.006, +0.043] | — |
| full (molecular + clinical) | 0.726 | 0.752 | 0.852 | 0.162 | +0.034 [+0.005, +0.059] | 0.010 |
| **deployed** (+ age spline, baseline induction) | **0.751** | **0.787** | **0.872** | **0.149** | **+0.059 [+0.030, +0.088]** | **0.001** |

**Read this honestly: age carries most of the signal, ELN adds a little, and molecular data adds a real
but modest amount only in combination.** Cell state and mutations *on their own* are materially worse
than the clinical baseline.

Hold-out risk tertiles separate at log-rank χ² = 38.8, **p = 4.7 × 10⁻¹⁰**.

### The treatment covariate is the single largest gain anywhere in the platform
Adding **baseline induction type** takes CV C-index 0.726 → 0.751 — more than all molecular data
combined. Only **diagnosis-time** treatment is used: `n_regimens` correlates 0.571 with follow-up time
(a patient must survive to accumulate regimens) and including it would read the outcome. Transplant is
excluded for the same reason (immortal-time bias).

---

## 4. External validation — TCGA-LAML (the credibility result)

Patient-grouped CV and a sealed hold-out protect against learning a *patient*. Neither protects against
learning a *cohort*. This is the test that separates them.

**Frozen:** Cox coefficients, PCA rotation, variable-gene selection, NNLS fusion weights — all loaded
unchanged, no refitting, no recalibration against TCGA outcomes.
**Refit on TCGA:** the per-gene z-reference (median/MAD) only, because RSEM log2 and BeatAML units are
not the same scale. Same cohort-matched-reference mechanism the drug layer uses for single-cell input.

**149 patients · 92 deaths · 57 censored · median follow-up 1.00 y · 89.4% gene overlap**
(UCSC Xena `TCGA.LAML.sampleMap/HiSeqV2`, Illumina HiSeq RSEM)

| arm | C-index | 95% CI | AUC 1 y | AUC 2 y | Δ vs age + cytogenetics | P(Δ ≤ 0) |
|---|---|---|---|---|---|---|
| **deployed** | **0.706** | 0.654 – 0.758 | 0.741 | 0.806 | **+0.035** | **0.030** |
| full | 0.689 | 0.636 – 0.741 | 0.720 | 0.787 | +0.018 | 0.182 |
| molecular | 0.669 | 0.607 – 0.728 | 0.716 | 0.756 | −0.003 | 0.522 |
| age + cytogenetics | 0.670 | 0.613 – 0.729 | 0.696 | 0.752 | *baseline* | — |
| clinical block | 0.664 | 0.608 – 0.723 | 0.681 | 0.735 | — | — |

**Risk tertiles** (cut at TCGA's own tertiles of predicted risk):

| group | n | deaths | 2-y survival | median |
|---|---|---|---|---|
| low | 50 | 17 | **71.7%** | not reached |
| intermediate | 49 | 31 | 41.1% | 1.42 y |
| high | 50 | 44 | **13.7%** | 0.75 y |

Log-rank low vs high: χ² = 38.0, **p = 7.0 × 10⁻¹⁰**. Curves ordered, non-crossing (Fig. Sv4).

**Three caveats, stated not buried:**
1. TCGA mutation calls are absent from this download, so the `mut` block enters as unreported — the same
   degraded path a real upload without a mutation caller takes. The fused arms run one block short.
2. TCGA-LAML is **treatment-homogeneous**, so the +0.035 that induction type contributes in BeatAML
   cannot be validated there at all. The deployed arm's edge over `full` here is the age spline alone.
3. CALGB cytogenetic risk is mapped onto ELN 2017 (Favorable → Favorable, Intermediate/Normal →
   Intermediate, Poor → Adverse). Close, not identical, so the bar cleared is approximate.
   Median follow-up is only 1.0 y, so the 5-y AUCs rest on a thin tail.

**The molecular-alone result reproduced exactly: −0.003 in BeatAML, −0.003 in TCGA, independently.**

---

## 5. Calibration (sealed hold-out, n = 89)

| horizon | predicted | observed (KM) | mean abs. gap across predicted quartiles |
|---|---|---|---|
| 1 y | 0.540 | 0.486 | 0.064 |
| 2 y | 0.395 | 0.387 | **0.043** |
| 5 y | 0.301 | 0.370 | 0.075 |

Quartile detail at 2 y — predicted vs observed: 0.018/0.000, 0.173/0.243, 0.531/0.496, 0.837/0.787.

---

## 6. Group timing works; individual timing does not

**By risk tertile** (predicted vs observed median survival):

| group | n | deaths | predicted median | observed KM median | error |
|---|---|---|---|---|---|
| low | 30 | 8 | 6.93 y | not reached | — |
| intermediate | 29 | 16 | 1.29 y | 1.15 y | **0.14 y** |
| high | 30 | 27 | 0.44 y | 0.33 y | **0.11 y** |

**For one patient** (evaluated only on the 47 who died inside follow-up — censored patients have no
known survival time, and including them would flatter the error):

| metric | value |
|---|---|
| median absolute error | 0.39 y |
| MAE | 0.87 y |
| within 6 months | 59.6% |
| within 1 year | 80.9% |

But observed survival **inside a single predicted-risk band** spans, p10–p90: low 0.38–1.42 y,
intermediate 0.04–1.43 y, high 0.04–0.86 y. And conformal prediction intervals — which achieve exact
empirical coverage — are **22.6 years wide at 80%** (52.8 y at 90%).

That is irreducible uncertainty in the underlying data, not a tuning failure. It is why the deployed
report leads with **restricted mean survival** (average of the next N years survived, defined even when
the median is never reached) rather than the median: the low-risk median of 6.93 y is extrapolated past
the end of follow-up and should not be quoted.

---

## 7. Limits, measured rather than asserted

- **Close to useless without intensive induction.** C-index 0.724 in the 284 intensively treated,
  **0.554** in the 71 who were not — where age + ELN is itself *below chance* at 0.481. Small n, but the
  model should not be trusted in that group.
- **Not causal about treatment.** Induction type removes treatment as a *confounder*; it does not license
  counterfactuals. 529 of 535 patients received the same standard induction, and induction intensity
  partly encodes clinician judgement about fitness rather than a randomised choice.
- **Single-cell survival is a bridge, not a fit.** The atlas has vital status for **0** samples, so the
  layer is BeatAML-bulk-trained and reaches single-cell patients only through the bulk-equivalent
  bridge — a transfer not itself validated against observed survival.
- **Refuses thin input.** Below 80% gene coverage the risk *group* starts to flip (percentile drift
  0.06 → 0.23 at 75% coverage), so it refuses rather than returning a confident number; 80–95% warns.
- **Assumes AML on bulk input.** The healthy-vs-diseased gate needs single cells. A healthy pooled-CD34
  control once scored at the 97th risk percentile; the layer now refuses when the gate says healthy, and
  states the assumption when no gate could run.

---

## 8. Two defects found and fixed while running the external validation

1. **The bundle never stored the age-spline knots**, so inference fell back to a hardcoded `[45, 60, 71]`
   while the coefficients had been fitted against the training quartiles. Every upload rebuilt the spline
   basis in the wrong place. Knots are now computed and carried in the model: **[47.0, 61.0, 70.0]**.
2. **`build_blocks` accepted a `train_idx` argument that neither call site passed**, so knots were placed
   using quantiles computed over the sealed hold-out patients as well as the training ones.

Both fixed; model retrained. **Every metric is unchanged to four decimal places** (CV 0.751, hold-out
0.787, AUC 2 y 0.8715, Brier 0.1494) — the leak was immaterial (three quantiles of one covariate across
444 patients). Verified the argument now actually takes effect: train-only knots give 70.0 where
whole-cohort quantiles would give 70.5.

Separately, the deployed arm was **missing from the incremental-gain bootstrap**, so any report that
selected it showed a blank where its value over age + ELN goes. Now computed: **+0.059, P = 0.001.**

---

## 9. Files and reproduction

| file | what it is |
|---|---|
| `METHODS_survival.md` | full methods |
| `VALIDATION_TCGA_LAML.md` | the external validation write-up |
| `survival_model_card.json` | every arm, CV + hold-out + incremental, cohort provenance |
| `validation_tcga_laml.json` | all TCGA numbers |
| `survival_time_validation.json` | calibration, tertile medians, per-patient error |
| `Sv1_survival_discrimination.png` | discrimination |
| `Sv2_risk_groups_and_calibration.png` | risk groups + calibration |
| `Sv3_group_versus_individual_timing.png` | why group timing works and individual does not |
| `Sv4_tcga_external_validation.png` | TCGA-LAML KM by predicted risk tertile |

```bash
python pipeline/train_survival_model.py      # cohort, all arms, sealed hold-out, model card
```

```bash
python pipeline/validate_tcga_laml.py        # fetches TCGA-LAML and runs the frozen transfer
```

`lifelines` / `scikit-survival` are not installed, so Cox partial likelihood, Breslow baseline hazard,
Harrell's C-index, Kaplan-Meier, log-rank, time-dependent AUC and IPCW Brier are implemented directly in
`pipeline/amlmm/survival/coxph.py`, which runs **six known-answer self-checks** (`python coxph.py`): it
recovers a simulated β = 1.0 as 0.987, returns exactly 1.0/0.0/0.5 for perfect/reversed/uninformative
concordance, and reproduces a hand-worked Kaplan-Meier to four decimals. Nothing downstream is trusted
until those pass.
