# MOSAIC-AML — survival layer

Can the platform predict whether a patient survives, and for how long? **Partly.** This document says
exactly how far it goes and where it stops.

---

## 1. The short answer

| question | answer |
|---|---|
| Can it rank who dies sooner? | **Yes.** C-index **0.752** on a sealed hold-out of patients (0.726 cross-validated). |
| Can it give a calibrated probability of surviving to a horizon? | **Yes.** Predicted vs observed 0.54/0.49 at 1 y, 0.40/0.39 at 2 y, 0.30/0.37 at 5 y; mean absolute gap across predicted quartiles 0.043–0.075. Two-year AUC **0.852**. |
| Can it say how long a *group* of similar patients will live? | **Yes.** Median survival per predicted risk tertile is accurate to ~6 weeks (1.29 y predicted vs 1.15 observed for intermediate risk; 0.44 vs 0.33 for high). |
| Can it say how long **one person** will live? | **Not reliably.** Median absolute error 0.39 y, MAE 0.87 y; 60% within six months, 81% within a year — but actual survival inside a single predicted-risk band spans ~1–1.4 years between the 10th and 90th percentile. |
| Does adding **treatment** help? | **Yes, most of all.** Baseline induction type takes C-index 0.726 → **0.750** (+0.059 over age+ELN vs +0.034 without it). Survives the reverse-causation check. |
| Does the molecular data beat age + ELN 2017? | **Not on its own.** Expression alone adds **−0.002** C-index over age+ELN. Only the *combination* adds: **+0.034, 95% CI [+0.006, +0.059], P(no gain) = 0.008.** |

**So the deployed output is a survival curve, horizon probabilities and a risk group — never a single
predicted number of months for one person.**

---

## 2. Cohort

BeatAML2, restricted three ways, each of which changes the answer:

1. **Known vital status** — 22 specimens with expression are `Unknown` and are dropped (they have a
   follow-up time but no outcome, so they are neither events nor censorings).
2. **Initial-diagnosis specimens only** — `overallSurvival` runs from diagnosis. A relapse specimen
   answers a different question *and* leaks prognosis, because such a patient has by definition already
   survived to relapse. This takes 649 usable specimens down to 446.
3. **One specimen per patient** — earliest kept.

**Final: 444 patients, 245 deaths (55%), 199 censored, median follow-up of the living 2.3 years.**
Sealed hold-out: **89 patients (20%), 51 deaths**, drawn once and never touched by any fitting,
selection or calibration step.

Censoring is respected everywhere: a patient last known alive contributes as censored at their
follow-up time, never as a death.

---

## 3. Model

`lifelines` / `scikit-survival` are not installed, so the primitives are implemented directly in
`amlmm/survival/coxph.py` — ridge-penalised Cox partial likelihood, Breslow baseline hazard, Harrell's
C-index, Kaplan-Meier, log-rank, time-dependent AUC, IPCW Brier. The module runs **six known-answer
self-checks** (`python coxph.py`): it recovers a simulated β = 1.0 as 0.987, returns exactly 1.0/0.0/0.5
for perfect/reversed/uninformative concordance, reproduces a hand-worked Kaplan-Meier to four decimals,
and so on. Nothing downstream is trusted until those pass.

Structure mirrors COMPASS-AML, for the same reason (blocks differ by orders of magnitude in width):

    risk = stack( cox_rna(x), cox_state(x), cox_mut(x), cox_clin(x) )

Per-block Cox models, stacked by a small Cox fitted on **inner out-of-fold** block predictions, floored
to the best single block if stacking does not beat it. The feature space (PCA, z-reference) is refit
inside every fold; patients never span folds.

---

## 4. What each arm is worth

| arm | CV C-index | hold-out | Δ vs age+ELN (95% CI) |
|---|---|---|---|
| ELN 2017 alone | 0.580 | 0.624 | — |
| cell state | 0.576 | 0.569 | −0.116 [−0.168, −0.069] |
| mutations | 0.626 | 0.579 | −0.066 [−0.112, −0.017] |
| age alone | 0.686 | 0.736 | — |
| RNA | 0.690 | 0.715 | −0.002 [−0.041, +0.041] |
| molecular (RNA+state+mut) | 0.689 | 0.716 | −0.003 [−0.043, +0.037] |
| **age + ELN 2017** | **0.692** | **0.725** | *baseline* |
| clinical block | 0.709 | 0.745 | +0.018 [−0.006, +0.043] |
| **full (molecular + clinical)** | **0.726** | **0.752** | **+0.034 [+0.006, +0.059]** |

Risk tertiles separate on the hold-out with log-rank **p = 2.8 × 10⁻⁹**.

The honest reading: **age carries most of the signal, ELN adds a little, and molecular data adds a
real but modest amount only in combination.** Cell state and mutations *on their own* are worse than
the clinical baseline.

---

## 5. Deployment behaviour, and what it refuses to do

- **It picks the arm it can actually feed.** `full` needs age and ELN; most uploads have neither. When
  they are missing it uses `molecular` (0.716) and says so, rather than imputing a median age and
  quoting the stronger model's accuracy.
- **Cohort-matched references.** A single-cell bulk-equivalent z-scored against BeatAML lands far
  outside that distribution and collapses the curve — the first smoke test told a patient they had two
  weeks to live. The model now carries a single-cell expression reference and the single-cell risk
  distribution; a single-cell patient is ranked among single-cell samples, then mapped onto the BeatAML
  risk scale, whose baseline hazard is the only one tied to observed survival.
- **It refuses below 80% gene coverage.** Withholding genes drifts one sample's risk percentile
  0.06 → 0.08 (90%) → 0.23 (75%) → 0.58 (50%): below ~80% the risk *group* flips. Between 80% and 95%
  it warns. An all-zero vector from a gene-ID mismatch previously returned a confident 97th-percentile
  prognosis; that is now impossible.
- **It refuses for a specimen the control gate called healthy.** A healthy pooled-CD34 control was
  scoring at the 97th risk percentile with 0.4% one-year survival.
- **For bulk input it says the number assumes AML.** The healthy-vs-diseased gate needs single cells,
  so it cannot fire on a bulk upload.

---

## 6. Limitations

- **Prognosis, not prophecy — and no longer fully treatment-blind, but not causal either.** Baseline
  induction type is now a covariate (C-index 0.726 → 0.750), which removes treatment as a *confounder*.
  It does **not** license counterfactuals: 529 of 535 patients received the same standard induction, and
  induction intensity partly encodes clinician judgement about fitness rather than a randomised choice.
  The model cannot say what would have happened under a different therapy.
- **It is close to useless in non-intensively-treated patients.** Stratified: C-index 0.724 in the 284
  who received intensive induction, but **0.554** in the 71 who did not (where age+ELN is itself below
  chance at 0.481). Small n, but the model should not be trusted in that group.
- **An individual lifespan interval is decades wide.** Conformal prediction achieves exact coverage, and
  the honest 80% interval spans **22.6 years** (90%: 52.8 y). This is irreducible uncertainty in the
  baseline data, not a tuning failure.
- **Individual timing is weak** (§1). A single number of months for one person is not supportable from
  this data, which is why the deployed output does not produce one as its headline.
- **The low-risk median is extrapolated.** More than half the low-risk group is still alive, so the
  observed median is genuinely not reached; the predicted 6.9 years sits past the follow-up window.
- **Single-cell survival cannot be learned here.** The atlas has essentially no survival metadata
  (vital status non-null for 0 samples), so the layer is BeatAML-bulk-trained and reaches single-cell
  patients only through the bulk-equivalent bridge — a transfer that is not itself validated against
  observed survival in single-cell-profiled patients.
- ~~**No external cohort.**~~ **Resolved — see §7.** The model has now been transferred, frozen, to
  TCGA-LAML. It holds up (C-index 0.706), with the caveats recorded there.

---

## 7. External validation — TCGA-LAML

Full report: [`VALIDATION_TCGA_LAML.md`](VALIDATION_TCGA_LAML.md) · reproduce with
`python pipeline/validate_tcga_laml.py`

The Cox coefficients, PCA rotation, variable-gene selection and NNLS fusion weights were loaded
unchanged and applied to **149 TCGA-LAML patients (92 deaths)** profiled on a different platform
(Illumina HiSeq RSEM), at different institutions, in an earlier treatment era. Only the per-gene
z-reference was cohort-matched, as it is for single-cell input.

| | BeatAML CV | BeatAML sealed hold-out | **TCGA-LAML (frozen transfer)** |
|---|---|---|---|
| C-index, deployed arm | 0.751 | 0.787 | **0.706** (95 % CI 0.655 – 0.758) |
| AUC at 2 y | — | 0.872 | 0.806 |
| gain over age + ELN/cytogenetics | **+0.059** [+0.030, +0.089], P = 0.001 | — | +0.035, P = 0.030 |

(The deployed arm's gain over age + ELN was previously absent from the model card — the bootstrap ran
over six arms and `deployed` was not one of them, so every report that selected it showed a blank where
its added value goes. It is now computed: **+0.059, 95 % CI [+0.030, +0.089], P(no gain) = 0.001**,
against an age + ELN out-of-fold C-index of 0.692.)

Risk tertiles separate **71.7 % vs 13.7 %** two-year survival (log-rank p = 7.0 × 10⁻¹⁰), curves
ordered and non-crossing.

Two findings that cut against the model are on the record there rather than omitted: the **molecular
blocks alone do not beat age + cytogenetics** in TCGA (ΔC = −0.003), and the treatment covariate that
produces the largest gain in BeatAML **cannot be validated** in a cohort where everyone received the
same induction.

---

## Reproducibility

| step | script |
|---|---|
| survival primitives + self-checks | `pipeline/amlmm/survival/coxph.py` |
| cohort assembly and provenance | `pipeline/amlmm/survival/data.py` |
| stacked per-block Cox | `pipeline/amlmm/survival/model.py` |
| train + all arms + incremental bootstrap | `pipeline/train_survival_model.py` → `deliverables/survival_model_card.json` |
| calibration, group vs individual timing | `pipeline/eval_survival_time.py` → `deliverables/survival_time_validation.json` |
| figures | `pipeline/build_survival_figures.py` → `deliverables/figures/Sv1–Sv3` |
| per-patient inference hook | `pipeline/survival_layer.py` (called from `ingest_patient.py`) |
