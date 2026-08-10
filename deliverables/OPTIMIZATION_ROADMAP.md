# MOSAIC-AML — what would actually make it better

An audit of every improvement worth trying to COMPASS-AML (drug response), the survival model, and the
lifespan estimate. Each entry says what to do, why it should help, roughly how much, what it costs, and
**whether the data to run it already exists** — because a catalogue of experiments we cannot run is
worthless.

Current state, for reference:

| | metric |
|---|---|
| COMPASS-AML | per-inhibitor AUROC 0.774 · per-patient top-1 retrieval 34% vs 10% chance |
| survival discrimination | C-index 0.752 sealed hold-out · +0.034 over age+ELN |
| lifespan (group) | median survival per risk tertile accurate to ~6 weeks |
| lifespan (individual) | MAE 0.87 y; ~1–1.4 y spread inside one risk band |

---

# 0. "Treatment-blind" — what it means, and yes, it can largely be fixed

## What it means now

The survival model learns *P(death | tumour biology at diagnosis)*. Two patients with identical biology
who received completely different therapy are, to the model, the same observation. So it silently
conflates two things:

- **prognosis** — how aggressive this disease is, and
- **what happened to be done about it** — which is not random.

Three consequences, in increasing order of how much they matter:

1. It can only say *"patients who look like this tend to do X"*, never *"this patient will do X **if
   treated with Y**"*.
2. Part of what looks like molecular prognostic signal may be a **treatment-assignment** signal. Older
   and less fit patients get less intensive therapy and die sooner; the model reads whatever biology
   correlates with that assignment and is rewarded for it. This inflates the apparent value of the
   biology.
3. It cannot support the question a clinician actually has, which is comparative: *which* treatment.

## The data exists to fix most of this

Checked, in the 671 BeatAML specimens with expression:

| field | non-null |
|---|---|
| `typeInductionTx` | 535 |
| `responseToInductionTx` | 533 |
| `responseDurationToInductionTx` | 535 |
| `cumulativeTreatmentTypes` (includes transplant) | 611 |
| `cumulativeTreatmentRegimens` (free-text drug names) | 611 |

Named-drug exposure recoverable from the regimen text, among initial-diagnosis patients: **cytarabine
303, azacitidine 93, decitabine 64, sorafenib 33, midostaurin 29, venetoclax 21**.

## Three levels of fix, honestly graded

**L1 — treatment-adjusted prognosis. Easy, do it. Modest effect.**
Add treatment covariates to the Cox: intensive chemotherapy vs hypomethylating vs palliative,
transplant, number of regimens. This removes treatment as a *confounder* so the reported molecular
effect is cleaner.
*Caveat that limits the gain:* **529 of 535 patients received "Standard Chemotherapy" induction.** At
induction, treatment is very nearly a constant, so there is little to adjust away. The exception is
transplant, which does vary and matters enormously.

> **⚠ Transplant must be time-dependent.** A patient has to survive long enough to be transplanted.
> Coding transplant as a baseline covariate makes it look magically protective — classic immortal-time
> bias. Use a time-dependent covariate or a landmark analysis at 3–6 months. Getting this wrong would
> produce a dramatic, completely spurious result.

**L2 — treatment-stratified prognosis. Easy, genuinely useful.**
Fit separately within strata (intensive vs non-intensive; transplanted vs not). Answers *"given this
treatment path, what is the prognosis"*, which is a real clinical question and needs no causal
assumptions.

**L3 — counterfactual "survival if given A instead of B". Possible to attempt, but I would not trust it here.**
IPTW, T-learners and causal forests all exist, but they require no-unmeasured-confounding, and
treatment assignment in AML is driven overwhelmingly by fitness and physician judgement — much of which
is not in this dataset. With n=444 and near-constant induction therapy, an estimated treatment effect
would mostly be a restatement of the confounding. **Recommendation: do L1 and L2, and be explicit that
L3 needs a randomised or much larger, richer observational cohort.**

---

# 1. COMPASS-AML — treatment-recommendation accuracy

## ★ The three that matter most

### A1. Validate against real clinical outcomes — the biggest gap in the whole platform
Not an accuracy optimisation; a credibility transformation. COMPASS currently predicts an **assay**.
Nobody has asked whether it predicts a **patient**.

The cohort exists: **254 patients with an inhibitor screen, a clean complete-response/refractory label,
and survival.** Of the labelled patients, **272 received cytarabine** and **93 received an HMA**.

Three tests, all runnable now:
- Does predicted ex-vivo sensitivity to **cytarabine** separate CR from refractory in the 272 who
  actually received it?
- Does it predict overall survival in that group, adjusting for age and ELN?
- Same for the HMA group (n=93) and, underpowered but worth reporting, venetoclax (n=21).

*Why it matters:* a positive result is the difference between "an interesting ex-vivo model" and
"evidence that ex-vivo prioritisation tracks clinical benefit". A null result is equally important and
should be published as such — it would set a hard ceiling on how the tool may be described.
**Cost: ~1 day. Data: exists.**

### A2. ⚠ Decompose the response matrix — the headline AUROC may be partly an artefact
The current per-drug AUROC (0.774, "given this drug, which patients respond") may substantially reflect
a **patient main effect**: some patients' cells simply die readily in culture, for reasons including
blast fraction, sample handling and viability. A model that learns "this specimen is globally
sensitive" scores well on every drug without knowing anything drug-specific.

**The experiment:** decompose the 520 × 118 AUC matrix into `patient effect + drug effect +
interaction`, then re-evaluate the model against the **interaction term alone**. That number is the
honest measure of drug-*specific* prediction.

*Why it matters:* if most of the signal is the patient main effect, the headline is inflated and the
per-patient ranking task (top-1 34% vs 10%) is the more trustworthy figure — which happens to be the
one the deployed system uses anyway. Either way we would then know.
**Cost: half a day. Data: in hand. Run this before optimising anything else.**

### A3. Multi-task matrix factorisation across drugs
Current design is per-family ridge + per-drug residual. The response matrix is 520 × 118 and **~87%
observed** — close to ideal for collaborative filtering with side information (patient features + drug
target descriptors). Low-rank sharing across drugs is the standard winning approach on GDSC/CTRP-style
data and should particularly help the low-n inhibitors that the per-drug residual currently barely
trains on.
*Expected: the largest single modelling gain, plausibly +0.02–0.04 AUROC.* **Cost: 2–3 days.**

## Data and target improvements

### A4. Download and use the raw concentration-level inhibitor data
`beataml_wv1to4_raw_inhibitor_v4_dbgap.txt` is **not currently downloaded**. It enables:
- **A clinically-anchored AUC** — integrate the viability curve only over concentrations plausibly
  reachable in a patient, instead of over the whole tested range including implausibly high doses. This
  is arguably a better target than the released AUC, and directly attacks the "sensitivity only appears
  at the top decade" problem the pharmacology agent currently flags after the fact.
- Alternative endpoints: IC50, E_max, Hill slope, area over the curve.
- **Per-observation curve-quality weights** from actual fit residuals, replacing the current
  binary QC drop with a weighted fit (keeps information instead of discarding 603 rows).
*Expected: modest but foundational — a better target raises everyone's ceiling.* **Cost: 1–2 days.**

### A5. Test the imputed modalities on drug response — a promising untested lead
COMPASS currently uses RNA / cell-state / mutations / clinical. The platform also has **ADT,
metabolite, lipid and GRN imputers** (`rna2*`) that were never tried for drug response. In the mutation
layer these were the *best single block* for many drivers — metabolite and lipid especially. Drug
sensitivity is a metabolic and apoptotic phenotype, so this is a well-motivated guess rather than a
fishing trip.
*Expected: unknown, plausibly the second-largest gain. Also a clean negative if it fails.*
**Cost: 1–2 days. Data: imputers exist.**

### A6. Two-way normalisation of the response target
Currently the AUC is standardised **within drug**. Also standardising **within patient** would make the
model predict drug-specific deviation rather than general sensitivity — which is exactly what a
recommendation needs. Directly related to A2; run them together.

### A7. Wave-aware modelling instead of wave-exclusion
Three clinically important inhibitors — **Cytarabine, Nutlin-3a, GDC-0941** — are currently demoted for
acquisition-wave instability. Adding wave as a covariate, or a ComBat-style correction on the AUC,
could recover them into the primary panel. Cytarabine is the induction backbone; having it properly
modelled matters more than any single AUROC point.

## Model and calibration refinements

| # | experiment | expected | cost |
|---|---|---|---|
| A8 | Drug-similarity kernel (target-overlap Jaccard + learned embedding) replacing the 16 coarse family groups — helps low-n drugs | +0.005–0.015 | 1 d |
| A9 | Nonlinear learners (gradient-boosted trees, small MLP with drug embedding). *Prior:* the mutation-layer bake-off found linear already optimal, so this may well be a null — worth one controlled test, not a campaign | 0 to +0.01 | 1 d |
| A10 | Pathway/gene-set scores instead of raw PCs — often better-powered and interpretable | +0.005–0.015 | 1 d |
| A11 | Per-drug-family hyperparameter tuning (currently one global alpha per block) | +0.005 | 0.5 d |
| A12 | Seed/fold ensembling — cheap variance reduction | +0.005–0.01 | 0.5 d |
| A13 | Ordinal or censored-tail modelling of the indeterminate middle 60%, instead of excluding it from the binary task | +0.005–0.01 | 1 d |
| A14 | Isotonic vs Platt calibration head-to-head (Platt won for mutations; untested here) | calibration only | 0.5 d |
| A15 | Self-supervised pretraining on all 707 BeatAML + 387 single-cell samples, including the **187 specimens with expression but no drug screen** that are currently unused | speculative | 2–3 d |
| A16 | External ex-vivo cohort (e.g. a functional-precision-medicine cohort with matched RNA) — the true generalisation test | validation only | depends on access |

*Note on the 187 unused specimens:* they were checked — every one has `analysisDrug = n`, i.e. no
inhibitor screen was ever run. They cannot be recovered as labelled training data, only used
unsupervised (A15).

---

# 2. Survival discrimination (who dies sooner)

Realistic ceiling first: published AML survival models from baseline data plateau around **C 0.75–0.80**.
We are at 0.752. So expect **+0.02–0.05 from everything below combined**, not +0.15. The large gains in
this field come from *post-treatment* information (response, MRD), which answers a later question.

## ★ Highest value

### B1. LSC17 and other published prognostic signatures
LSC17 is a 17-gene AML stemness score and one of the strongest published transcriptomic prognostic
markers. It is directly computable from the expression we already hold, and the platform already has
LSC infrastructure. Use it as (a) a benchmark the model must beat and (b) a feature.
*Why:* 60 RNA PCs on 194 events is badly over-parameterised; a compact, externally-validated signature
is exactly the right prior. **Cost: 0.5 day. Expected: +0.01–0.03, plus a much stronger paper.**

### B2. External validation on TCGA-LAML
~200 patients, public, with expression and survival — a genuinely independent cohort. Current validation
is a sealed hold-out *within* BeatAML, which does not test cohort transfer at all. Also pool the two for
training to increase the event count.
**Cost: 2 days. Expected: the credibility gain is larger than the accuracy gain.**

### B3. Variant-level detail the model currently throws away
Three specific, well-evidenced prognostic features that exist in the data and are not used:
- **FLT3-ITD allelic ratio** — the `allelic_ratio` column is right there in the clinical file, and high
  vs low ratio is prognostically distinct.
- **Mutation VAF and clonality** — `t_vaf` is in `mutations.txt`. TP53 VAF and allelic state in
  particular are strongly prognostic; presence/absence discards that.
- **Karyotype-derived features** — complex and monosomal karyotype, plus `consensusAMLFusions`. The
  free-text `karyotype` field is currently unparsed.

**Cost: 1–2 days. Expected: +0.01–0.03**, and these are the features a haematologist will immediately
ask why you left out.

### B4. ELN 2022 instead of ELN 2017
ELN 2022 is the current standard and reclassifies a meaningful fraction of patients. We already derived
it (`labels/eln2022_beataml.tsv`) — **available for 339 of 671**, which is the catch: it is less
complete than ELN 2017. Test both; possibly use ELN2022 where available with ELN2017 as fallback.

## Modelling

| # | experiment | expected | cost |
|---|---|---|---|
| B5 | Age as a spline rather than linear — age is the single dominant predictor and is currently modelled linearly | +0.005–0.015 | 0.5 d |
| B6 | Elastic net / univariate pre-screening on the RNA block instead of 60 unscreened PCs | +0.01 | 0.5 d |
| B7 | Random survival forest and gradient-boosted Cox — nonlinearity and interactions | +0.01–0.02 | 1 d |
| B8 | Competing-risks model (death from AML vs other causes) — `causeOfDeath` exists | correctness | 1 d |
| B9 | Deep survival (DeepSurv / Cox-nnet). *Prior: n=444 is small; expect a null* | 0 | 1 d |
| B10 | Joint multi-task training with the mutation and drug layers (shared representation) | speculative | 3 d |
| B11 | Landmark analysis at 3 and 6 months — handles transplant and early deaths properly | correctness | 1 d |
| B12 | Dynamic prognosis conditioning on induction response (533 labelled) — **large C-index jump expected, but it answers a later question and must not be compared to the baseline model** | large but different question | 2 d |

---

# 3. Remaining-lifespan prediction (how long)

This is the weakest part and, in my view, the most interesting to work on — the current honest position
is that group timing is good and individual timing is not.

### ★ C1. Conformal survival prediction — the right tool for this exact question
Gives **finite-sample-valid prediction intervals** with almost no distributional assumptions. Instead of
"about 14 months", the output becomes "with 80% confidence, between 5 and 26 months" — which is both
honest and genuinely usable. This is the single best answer to "how long, and how sure are you".
**Cost: 2 days. Expected: does not improve accuracy; it makes the output *correct*, which matters more.**

### ★ C2. Restricted mean survival time (RMST) instead of median
RMST is well-behaved under censoring and is **defined even when the median is not reached** — which
fixes the current low-risk group returning "not reached" and the model extrapolating 6.9 years past the
follow-up window. RMST at 2 and 5 years is a cleaner, more comparable quantity.
**Cost: 0.5 day. Expected: removes an existing embarrassment.**

### C3. Accelerated failure time (AFT) models
Cox models the *hazard*; AFT models log(time) directly. For a "how long" question AFT is the more
natural parameterisation and often better calibrated in the tails. Test log-normal and Weibull AFT
against the Cox.
**Cost: 1 day. Expected: +0.05–0.15 y on MAE.**

### C4. Discrete-time survival / multi-horizon classifiers
Instead of one proportional-hazards model, train a set of "alive at 6/12/24/60 months" classifiers with
shared features (or a discrete-time neural hazard model). This drops the proportional-hazards
assumption, which is very likely violated here — early chemotherapy deaths and late relapse deaths have
different hazard shapes.
**Cost: 1–2 days. Expected: better-calibrated horizon probabilities, especially at 5 years where the
current model is pessimistic (0.30 predicted vs 0.37 observed).**

### C5. Censoring-aware quantile regression
Predict the 10th, 50th and 90th percentiles of survival time directly. Gives an interval natively.
Complementary to C1.

| # | experiment | expected | cost |
|---|---|---|---|
| C6 | Per-horizon recalibration — the 5-year prediction is systematically pessimistic | calibration | 0.5 d |
| C7 | Smoothed/spline baseline hazard instead of the Breslow step function | cosmetic + slight | 0.5 d |
| C8 | Explicitly test the proportional-hazards assumption (Schoenfeld residuals) — if violated, C3/C4 are not optional | diagnostic | 0.5 d |
| C9 | Cause-specific time-to-event modelling | correctness | 1 d |

---

# 4. If I could only run ten things, in order

| # | experiment | why it is first | cost |
|---|---|---|---|
| 1 | **A2** decompose the response matrix (patient vs drug vs interaction) | may show the headline AUROC is partly a patient main effect — changes what everything else is optimising | 0.5 d |
| 2 | **A1** clinical-outcome validation of COMPASS (254 patients) | converts an assay model into evidence about patients; null result equally valuable | 1 d |
| 3 | **C2** RMST + **C1** conformal intervals | makes the lifespan output honest and usable rather than merely caveated | 2.5 d |
| 4 | **B1** LSC17 as benchmark and feature | strongest published prior; exposes whether 60 PCs are earning their keep | 0.5 d |
| 5 | **B3** FLT3-ITD allelic ratio, VAF, karyotype features | already in the files; the omissions a haematologist notices immediately | 1.5 d |
| 6 | **L1+L2** treatment adjustment and stratification (time-dependent transplant) | answers the treatment-blindness objection directly | 1.5 d |
| 7 | **A5** imputed modalities for drug response | well-motivated, untested, potentially the second-biggest COMPASS gain | 1.5 d |
| 8 | **A4** raw concentration data → clinically-anchored AUC | better target raises every downstream ceiling | 1.5 d |
| 9 | **A3** multi-task matrix factorisation | likely the largest pure-modelling gain for COMPASS | 3 d |
| 10 | **B2** TCGA-LAML external validation | the only true test of generalisation | 2 d |

Roughly three weeks of work. Items 1–5 are about half of it and carry most of the credibility gain;
items 7–9 carry most of the raw accuracy gain.

---

# 5. What will *not* fix this

Worth stating so effort is not wasted:

- **Bigger neural networks.** n=444 for survival and 520 specimens for drug response. The binding
  constraint is sample size and label quality, not model capacity.
- **More RNA features.** 14,237 genes on 194 events is already far past the point where adding features
  helps; the wins are in better *targets* and better *priors*, not wider inputs.
- **Chasing the ex-vivo AUROC upward.** Past roughly 0.80 the remaining variance is assay noise —
  culture conditions, viability, handling. A1 (does it predict patients?) is worth more than another
  0.02 of assay AUROC.
- **Predicting individual survival time to the month.** The ~1-year spread inside a risk band is
  substantially irreducible from baseline data: it reflects treatment, complications and chance. The
  fix is honest intervals (C1), not a better point estimate.
