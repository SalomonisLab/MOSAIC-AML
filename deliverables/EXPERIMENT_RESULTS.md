# What the experiments actually showed

Results of running the optimisation roadmap. Every arm is patient-grouped CV on the same split as the
deployed models, so the numbers are directly comparable. **Negative results are reported at the same
weight as positive ones** — three of the most useful findings here are negatives.

---

## 1. ⚠ The COMPASS headline number was inflated

**Experiment A2.** The response matrix decomposes as:

| component | share of variance |
|---|---|
| patient main effect (some specimens just die readily in culture) | 15.4% |
| drug main effect | 45.9% |
| **patient × drug interaction** (the only part that is a *recommendation*) | 46.8% |

The model's mean prediction per specimen correlates **0.540** with the patient main effect — it is
partly learning which specimens are globally sensitive. Re-scoring against the interaction term alone:

| metric | as reported | interaction only |
|---|---|---|
| mean per-drug Spearman | 0.365 | **0.223** |
| mean per-drug AUROC | 0.774 | **0.672** |

**61% of the signal survives.** The drug-specific component is real and well above chance, but the
headline was inflated. **All documentation should quote 0.672 as the drug-specific figure**, with 0.774
labelled as including the patient main effect. The per-patient ranking result (top-1 34% vs 10% chance)
is unaffected, because ranking within a patient removes the patient effect by construction — and that
is the number the deployed system actually uses.

## 2. ⚠⚠ Ex-vivo sensitivity did NOT predict clinical response

**Experiment A1**, the most important test in the roadmap. Out-of-fold predicted sensitivity to the drug
the patient *actually received*, against their real outcome:

| drug received | n scored | complete response AUROC | *p* | OS C-index: age+ELN → +model |
|---|---|---|---|---|
| **Cytarabine** | 149 | **0.435** (n=131, 93 CR) | 0.249 | 0.739 → 0.743 |
| Azacitidine | 50 | 0.550 (n=40) | 0.595 | 0.669 → 0.694 |
| Decitabine / Sorafenib / Midostaurin / Venetoclax | 16–25 | underpowered | — | — |

**In this cohort, predicted ex-vivo sensitivity does not translate into clinical response.** Cytarabine
is on the wrong side of chance; azacitidine is null.

Fair caveats, none of which rescue the result: cytarabine is already flagged wave-unstable; single-agent
ex-vivo cytarabine at 72 h is a poor model of 7+3 induction, which is a combination; and n=131 is
limited power. But 0.435 is not a trend in the right direction.

**Consequence for how the tool is described.** COMPASS-AML predicts an *ex-vivo assay*, and the one
chance we had to show that this tracks patient benefit came back null. The existing caveat — "a
prioritisation signal for trial matching or laboratory validation, not an estimate of clinical benefit"
— is now not merely cautious wording but a *tested and confirmed* limitation.

---

## 3. ★ Treatment information is the single biggest survival gain — and it is legitimate

**Experiments L1/L2.** This directly answers the treatment-blindness question.

| arm | C-index | Δ vs age+ELN (95% CI) |
|---|---|---|
| age + ELN 2017 | 0.692 | *baseline* |
| full (molecular + clinical) | 0.726 | +0.034 [+0.006, +0.060] |
| **full + baseline induction type** | **0.750** | **+0.059 [+0.029, +0.092]** |
| full + baseline + post-baseline treatment | 0.759 | +0.067 [+0.035, +0.098] |

**The reverse-causation check matters and it passed.** `n_regimens` correlates **0.571** with follow-up
time — you must survive to accumulate regimens — so a naive treatment block would be partly reading the
outcome. Restricting to **induction type, which is decided at diagnosis**, keeps almost the entire gain
(+0.059 of the +0.067). Post-baseline exposure adds only +0.008.

*Honest qualification:* induction intensity partly encodes clinician judgement about fitness, so this is
legitimate prognostic information available at diagnosis rather than a causal treatment effect. It does
not license counterfactual statements about what would happen under a different therapy.

**Treatment-stratified (L2)** exposes a real limitation:

| stratum | n | deaths | C (full) | C (age+ELN) |
|---|---|---|---|---|
| intensive induction | 284 | 139 | 0.724 | 0.691 |
| non-intensive | 71 | 55 | **0.554** ¹ | 0.481 |

¹ **Superseded — this was a pooled-fitting artefact.** The 0.554 came from a model trained on the whole cohort and evaluated inside this subgroup. Fitting *within* the stratum gives **0.681**, and adding the published TP53/FLT3-ITD/NRAS/KRAS rule gives **0.701** (`ELN2022_RISK_BENCHMARK.md` §4.2). What survives is the finding about the *guideline*: ELN 2017 scores 0.496 and ELN 2022 0.462 in this stratum — at or below chance.


The model works in fit patients receiving standard therapy and **is close to useless in the
non-intensive group**, where age+ELN is actually *below* chance. Small n, but it should be stated.

---

## 4. Four clean negatives that save future effort

| experiment | result | reading |
|---|---|---|
| **B1 LSC17** | 0.585 alone; **no gain** in combination (0.723 vs 0.726) | As implemented — 13/17 genes mapped, cohort-wide z-scoring rather than the published normalisation — it did not help. A faithful implementation might; this one did not. |
| **B3 variant detail** (VAF, karyotype, FLT3-ITD allelic ratio, ELN 2022) | **0.717, slightly worse than 0.726** | Forty-plus noisy columns on 194 events overfits. The information may be real but needs strong selection, not concatenation. |
| **B7 gradient boosting** | 0.591 | Nonlinearity does not help at n=355. Consistent with the mutation layer's earlier bake-off finding linear models optimal. |
| **C8/C3 proportional hazards** | Schoenfeld min *p* = **0.57**; log-normal AFT gives an identical **0.726** | PH holds. **The accelerated-failure-time and discrete-time work in the roadmap is unnecessary** — that is two or three days saved. |

Also null: univariate RNA screening (0.687 vs 0.690 — 60 PCs are not the bottleneck).

Small positive: **age as a spline**, 0.726 → **0.733**. Cheap, keep it.

---

## 5. ★ The lifespan question is now settled

**C1 — conformal prediction intervals.** Coverage is exact (80% target → 80% empirical; 90% → 90%), so
the method is working correctly. The intervals it produces:

| target coverage | median interval width |
|---|---|
| 80% | **22.6 years** |
| 90% | **52.8 years** |

**A statistically valid individual survival interval is too wide to be clinically useful.** This is not
a modelling failure to be optimised away — it is a measurement of how much irreducible uncertainty
there is in an individual AML survival time given baseline data (the log-normal AFT scale parameter is
σ = 2.03). It definitively answers "can it say how long *this person* will live": **no, and no amount of
model tuning will change that.**

**C2 — RMST works and should replace the median.** Defined for every risk group, including the low-risk
group whose median is never reached:

| risk group | RMST at 2 y | RMST at 5 y |
|---|---|---|
| low | 1.73 y | **3.72 y** |
| intermediate | 1.26 y | 2.38 y |
| high | 0.61 y | 0.76 y |

Clean monotone separation, no extrapolation past follow-up. **Recommend switching the reported quantity
from median survival to RMST.**

---

## 6. What this changes

**Deploy now:** baseline treatment covariates (+0.024 C-index over the current model), age spline
(+0.007), RMST in place of median survival.

**Change in the documentation:** the drug-specific COMPASS figure is 0.672, not 0.774; the ex-vivo →
clinical link was tested and was null; the survival model is near-useless in non-intensively-treated
patients; individual lifespan intervals are decades wide.

**Do not bother with:** AFT, discrete-time survival, gradient boosting, RNA feature screening, LSC17 as
currently implemented, or bulk concatenation of variant detail.

**Still worth running (not yet done):** COMPASS multi-task matrix factorisation (A3), the imputed
modalities for drug response (A5), the raw concentration-level data for a clinically-anchored AUC (A4),
and TCGA-LAML as an external survival cohort (B2). A3 and A5 are the two most likely remaining sources
of real gain.

---

*Reproduce:* `pipeline/exp_compass_diagnostics.py`, `pipeline/exp_survival_sweep.py` →
`deliverables/exp_compass_diagnostics.json`, `deliverables/exp_survival_sweep.json`.
