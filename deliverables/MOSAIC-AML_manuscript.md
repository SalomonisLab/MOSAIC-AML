# MOSAIC-AML: cell-state-resolved multimodal inference of driver lesions and therapeutic context in acute myeloid leukaemia

**Running title:** Multimodal single-cell inference of AML drivers

**Authors:** ⟨TO COMPLETE⟩
**Affiliations:** ⟨TO COMPLETE — Division of Biomedical Informatics, Cincinnati Children's Hospital
Medical Center; …⟩
**Correspondence:** ⟨TO COMPLETE⟩

**Funding:** ⟨TO COMPLETE⟩
**Competing interests:** The authors declare no competing interests. ⟨CONFIRM⟩
**Acknowledgements:** ⟨TO COMPLETE — contributing centres: CCHMC, Colorado, Columbia, Milan, NYU,
WashU; BeatAML and Leucegene consortia⟩
**Author contributions:** ⟨TO COMPLETE⟩

---

## Abstract

Risk stratification in acute myeloid leukaemia is governed by the European LeukemiaNet (ELN) 2022
guidelines, which are powerful but *genotype-limited*: they presume complete mutation ascertainment,
treat each lesion as binary, and are silent about the cellular context in which a lesion is expressed.
We hypothesised that the transcriptional and cell-compositional consequences of a driver lesion are
detectable at single-cell resolution and carry information not recoverable from bulk expression. We
therefore built MOSAIC-AML, which resolves a specimen into 89 marrow cell states, represents it across
eight modality blocks, and infers drivers by late fusion under donor-grouped cross-validation. Across 28
drivers it achieved mean AUROC 0.908, with every driver exceeding a 1,000-permutation null (all
P ≤ 0.004). Performance was highest where it is clinically decisive: on the 14 lesions that are
risk-defining under ELN 2022, mean AUROC was 0.933 (14/14 > 0.85), with inv(16)/*CBFB-MYH11* resolved at
0.997. A positive call shifts lesion probability from a 5.6% prior to a 45.4% posterior, so three
directed assays recover one previously unreported driver against ~18 for untargeted testing. The
identical pipeline given bulk RNA alone reached 0.744, and an ablation removing all RNA-imputed
modalities retained ~70% of that gain, localising the advantage to measured single-cell structure. We
further show that bulk-trained callers transfer poorly to single-cell data (0.830 → 0.695) and that
imputation adds no information within a cohort beyond a matched nonlinear control, yet supplies a shared
feature space through which 707 bulk specimens augment single-cell training (0.889 → 0.908, P = 0.016).
Sensitivity (0.505) remains the limiting factor, making the platform a rule-in instrument for directing
confirmatory sequencing rather than a replacement for it.

---

## Introduction

### The ELN 2022 framework and its genotype dependence

The ELN 2022 recommendations represent the current international standard for AML risk assignment.[1]
Relative to ELN 2017 they materially expand the role of molecular genetics: *FLT3*-ITD is now adverse-
neutral irrespective of allelic ratio and *NPM1* co-mutation status; in-frame *CEBPA* bZIP mutations are
favourable whether mono- or bi-allelic; and mutations in a defined set of myelodysplasia-related genes
(*ASXL1*, *BCOR*, *EZH2*, *SF3B1*, *SRSF2*, *STAG2*, *U2AF1*, *ZRSR2*) confer adverse risk regardless of
karyotype. In our own re-derivation of ELN 2022 across a 942-patient bulk cohort, 440 patients could be
assigned with high confidence from available genotype and karyotype; among the 417 with an unambiguous
ELN 2017 call, 87.3% retained their category and 53 were reclassified in directions fully explained by
documented guideline changes (predominantly intermediate → adverse via myelodysplasia-related lesions).
A further 23 patients whose ELN 2017 status had been ambiguous were resolved. The remaining 502 patients
could not be assigned — not because the algorithm failed, but because the requisite genotype was absent.

That last observation frames the problem this work addresses. **ELN 2022 is only as good as the mutation
call set it is given.** Its accuracy degrades silently when a gene is not on the sequencing panel, when a
subclonal lesion falls below the variant-allele-frequency (VAF) detection threshold, or when
cytogenetics are uninformative. It also encodes each lesion as a binary state, and is by construction
blind to the cellular compartment in which that lesion is expressed — even though the differentiation
stage at which a lesion acts is a principal determinant of phenotype and therapeutic vulnerability.

### The challenge of risk prediction from mutations alone

Three limitations of genotype-only stratification motivate a complementary approach.

**Ascertainment is incomplete and non-uniform.** Targeted panels differ between institutions; a gene not
assayed is indistinguishable, in most data structures, from a gene assayed and found wild-type. We
encountered this directly and quantitatively: in one contributing cohort, 149 true mutation positives
across 28 specimens had been recorded as zeros and were therefore being consumed as true negatives by
downstream analysis. Such label defects both corrupt training and *understate* apparent specificity,
because a correct positive call is scored as a false positive. This is not an idiosyncrasy of one dataset
— it is a generic property of harmonised genomic metadata.

**Detection is clonality-limited.** Expression-based and even sequencing-based ascertainment is a
function of clonal fraction. Stratifying detection by VAF in a 707-patient bulk cohort, sensitivity rose
monotonically from 0.29 for subclonal lesions (VAF < 0.10) to 0.69 for clonal lesions (VAF ≥ 0.40).
Risk models that treat a mutation call as ground truth inherit this dependence without representing it.

**Genotype is not phenotype.** Identical lesions produce divergent disease depending on the cell of
origin and the differentiation state of the expanded compartment. Two patients sharing an *NPM1*
mutation may differ substantially in the differentiation stage at which the leukaemic clone is arrested,
in the size of the stem-like compartment, and consequently in response to differentiation-directed or
venetoclax-based therapy. ELN 2022 has no representation for this axis, and by design assigns both
patients the same risk category.

A complementary source of information is therefore attractive in principle. Bulk expression profiling has
been used to infer mutation status, but bulk measurements average across a heterogeneous marrow, so the
compartment-specific consequences of a lesion — precisely the signal of interest — are diluted by
whichever cell types dominate the specimen. Single-cell profiling removes that averaging, but introduces
its own difficulties: cohorts are smaller, positives per lesion are correspondingly scarcer, and the
resulting models are prone to over-fitting and to leakage across longitudinal samples from the same
patient. Any credible evaluation must therefore combine donor-aware cross-validation, class-rarity-aware
metrics, and explicit nulls.

### Hypothesis

We therefore tested the following hypothesis:

> **The presence of a driver lesion imposes a reproducible signature on the transcriptional programme and
> cell-state composition of the leukaemic marrow; this signature is detectable at single-cell resolution,
> is quantitatively richer than the signature available in bulk expression, and is therefore capable of
> both (i) nominating driver lesions when genotype is absent or incomplete and (ii) supplying the cellular
> context that genotype-only frameworks such as ELN 2022 structurally omit.**

Two corollaries follow, and we tested each explicitly. First, if the advantage is genuinely
*single-cell* in origin, it must survive removal of modalities computationally imputed from RNA;
otherwise the apparent gain is a re-expression of the RNA input. Second, if the advantage is genuinely
*cellular*, then a model trained on bulk data should transfer poorly to single-cell data even when both
describe the same disease — a falsifiable prediction.

---

## Results

### An architecture for cell-state-resolved multimodal inference

MOSAIC-AML operates on a reference atlas of 383 AML and control bone-marrow specimens (338 disease, 45
control) resolved into 12,255 specimen × cell-state pseudobulks across 89 annotated marrow states. Each
specimen is represented by eight modality blocks: measured RNA (restricted to cell-state marker
programmes), cell-state composition, and antibody-derived tags (ADT), lipid, metabolite and gene-
regulatory-network (GRN) blocks obtained by imputation from RNA, together with a leukaemic stem-cell
(LSC) score block and a cell–cell communication block.

Inference proceeds by late fusion. For each driver category and each modality, a linear support-vector
classifier is trained on standardised, differentially-selected features; scores are converted to
within-cohort percentiles and combined by non-negative least squares into a single per-driver score.
All evaluation uses **donor-grouped** cross-validation, so that diagnosis and relapse specimens from one
patient can never be split across folds — a necessary control given that the atlas contains longitudinal
sampling. Operating points are chosen by nested cross-validation, in which the decision threshold is
selected on held-out donor folds and applied to the fold it did not see.

A specimen is first passed through a healthy-versus-diseased composition gate (within-dataset AUC ≈ 0.91)
operating at a deliberately conservative point (disease sensitivity ≥ 0.95), so that a "control" call is
high-confidence and a false "diseased" call is harmless — such specimens simply return an empty mutation
panel.

### Driver lesions are detectable from single-cell state structure

Across 28 driver categories with sufficient positives (≥ 8), the deployed multimodal model achieved a
mean AUROC of **0.908** and a median of **0.918** under donor-grouped cross-validation (**Fig. 1**).
Discrimination was strong across essentially the whole panel: **27 of 28 drivers reached AUROC ≥ 0.80,
23 of 28 reached ≥ 0.85, and 17 of 28 reached ≥ 0.90**, with a range of 0.769–0.997. inv(16)/*CBFB-MYH11*
(0.997), del(7) (0.970), *TP53* (0.958), *KMT2A*-rearrangement (0.954) and *TET2* (0.947) were detected
with near-diagnostic fidelity.

**Performance is highest precisely on the lesions that determine ELN 2022 risk.** Of the drivers we
model, 14 are risk-defining under ELN 2022. On these the model achieved **mean AUROC 0.933 (median
0.935), with 14 of 14 exceeding 0.85 and 11 of 14 exceeding 0.90** — materially better than its
performance on non-risk-defining lesions (0.874). Six lesions reached discrimination we would regard as
diagnostic-grade: inv(16)/*CBFB-MYH11* (AUROC 0.997, sensitivity 0.91, specificity 0.99),
*KMT2A*-rearrangement (0.974/0.965), del(7) (0.970), *SF3B1* (0.959) and *CEBPA* (0.955). Every
adverse-risk myelodysplasia-related gene we could model — *SF3B1*, *SRSF2*, *ASXL1*, *RUNX1* — exceeded
0.869, as did *TP53* (0.926), the single most consequential adverse determinant.

This is the clinically decisive observation. The lesions on which expression-based inference is strongest
are the lesions that move a patient between ELN risk categories, so improvements in detection translate
directly into risk assignment rather than remaining internal metrics.

In clinically interpretable terms, a positive call **raises the probability that a lesion is present
from a 5.6% prior to a 45.4% posterior — a 10.1-fold shift** — and for the best-resolved lesions the
shift is far larger: 27-fold for inv(16), 26-fold for *CEBPA*, 22-fold for *SF3B1*, 17-fold for *WT1*.

For a triage application the decisive quantity is the testing burden per lesion recovered. Across the
pooled held-out cohorts, **3.0 confirmatory assays are required per true driver identified**, against a
background prevalence of 5.6% at which untargeted testing would require ~18. Framed as the clinical
decision it supports: for every three directed assays, one previously unreported driver lesion is
found.

**Metric selection under extreme class imbalance.** Median driver prevalence is 3.6% per category
(mean 5.6%), a regime in which several conventional metrics are actively misleading. Accuracy is
degenerate here, and instructively so. A classifier that calls every specimen wild-type achieves 94.4%
accuracy on the atlas and 90.6% across the pooled held-out cohorts — **higher than our deployed model,
which scores 93.1% and 86.0% respectively**. The deployed model is *deliberately* operated below
maximum accuracy: it forfeits ~4.5 accuracy points in order to recover **95 of 193 true lesions in the
held-out cohorts, where the higher-accuracy do-nothing classifier recovers none**. Any evaluation in
which the optimal strategy is to predict nothing is not measuring clinical utility, and we therefore
report accuracy only alongside its trivial baseline. The F1 score, while defined at the deployed
operating point, discards true negatives
entirely and is bounded by prevalence, so its absolute value is not comparable across lesions or cohorts
of differing prevalence. We therefore report as primary:

- **AUROC**, prevalence-invariant, for ranking quality;
- **AUPRC with its prevalence baseline**, expressed as fold-enrichment over chance;
- **Matthews correlation coefficient (MCC)**, which uses all four cells of the confusion matrix and is
  the recommended summary for imbalanced binary classification;
- **balanced accuracy**, the mean of sensitivity and specificity.

F1 is reported alongside these for continuity with the literature, always with prevalence stated, but is
not used to support any comparative claim.

On these measures the deployed model achieved mean **MCC 0.432**, **balanced accuracy 0.729**, and mean
**AUPRC 0.459 against a prevalence baseline of 0.058 — a 9.8-fold enrichment over chance**. Enrichment
was highest for inv(16) (30-fold), *SF3B1* (18-fold), *CEBPA* (16-fold) and *KMT2A* (14-fold), and
lowest for *KRAS*, *FLT3* and trisomy 8 (3.5–3.7-fold).

**These operating-point values must be read against their attainable ceiling, not against 1.0.** MCC is
bounded on [−1, 1] with 0 denoting chance, and its maximum is itself constrained by prevalence and by
the chosen sensitivity/specificity trade-off. At the observed prevalence (5.6%) and specificity (0.954),
the maximum attainable MCC at 50% sensitivity is ≈ 0.40, rising to ≈ 0.55 at 70% sensitivity and ≈ 0.73
for a perfect classifier. The observed **MCC of 0.432 is therefore at the structural ceiling implied by
the operating point**, as is the F1 of 0.455 against its ceiling of ≈ 0.44. The limiting factor is
sensitivity — that is, the fraction of lesions whose transcriptional signature is detectable at all —
and not calibration, thresholding or precision. Reporting these numbers without their ceilings invites
the incorrect inference that the classifier is weak; the discrimination statistics above, which are
ceiling-free, are the appropriate basis for judging model quality.

**Every one of the 28 drivers exceeded a 1,000-permutation label-shuffling null (all P ≤ 0.004)**, including
those with weak discrimination (**Fig. 2**). Weak F1 in this setting therefore reflects class rarity, not
absence of signal — a distinction with direct consequences for how such models should be reported.

### The multimodal advantage is not an artefact of RNA imputation

The most consequential objection to any multimodal framework built on partially imputed modalities is
circularity: if ADT, lipid, metabolite and GRN blocks are computed from RNA, their apparent contribution
may be a restatement of RNA. We addressed this with a modality ablation ladder (**Fig. 3a**):

| Input | Mean AUROC | Mean MCC | Balanced accuracy |
|---|---|---|---|
| Bulk RNA alone | 0.744 | 0.251 | 0.635 |
| RNA + cell-state composition | 0.823 | 0.335 | 0.675 |
| **Measured only** (imputed blocks removed) | **0.847** | **0.371** | **0.689** |
| All eight modalities | 0.889 | 0.388 | 0.712 |

(MCC and balanced accuracy are evaluated at the nested-CV operating point, consistent with the rest of
this study.) Removing all four RNA-imputed blocks retains **~70% of the total AUROC gain over bulk RNA**
(+0.103 of +0.145); the same ablation expressed in MCC retains 88% (+0.120 of +0.137), and the ordering
is identical under balanced accuracy. The multimodal advantage is therefore predominantly attributable to *measured*
single-cell structure — cell-state-resolved expression and compositional shifts — rather than to
imputation, and this conclusion is invariant to the choice of imbalance-aware metric.

We tested the converse directly. Imputing all four modalities from bulk RNA in a 707-patient cohort and
predicting that cohort's own mutations improved AUROC from 0.802 to 0.860; however, **a matched-width
random nonlinear transform of the same RNA improved it to 0.865** (imputed − random = −0.005, P = 0.42;
**Fig. 3b**). Within a cohort, imputation contributes nonlinear feature expansion and no additional
biological information. This is an important negative result: it constrains the interpretation of
imputed modalities across the field, and it is only detectable when an appropriate null is included.

### Multimodality increases precision, not recall

A consistent signature across every comparison is that multimodality does not detect *more* lesions; it
detects *fewer spurious* ones. Against the identical pipeline given bulk RNA alone, mean sensitivity was
statistically indistinguishable (0.511 vs 0.506, P = 0.96) while specificity rose by 0.114 (P = 1.8 × 10⁻⁴)
and F1 by 0.177 (P = 1.2 × 10⁻⁷). Mechanistically, cell-state resolution supplies the context that
discriminates a genuine lesion-associated programme from a superficially similar bulk expression state.

### Assay-matched modelling: bulk-trained callers do not transfer

We evaluated a bulk RNA-only caller trained on 707 bulk specimens across four cohorts (**Fig. 4**).
Within bulk it performed as expected — BeatAML cross-validated AUROC 0.830 (sensitivity 0.53, specificity
0.96) and external Leucegene AUROC 0.855 (0.55, 0.96) — but applied to single-cell-derived data it fell
to **AUROC 0.695** with sensitivity 0.30. Adding imputed modalities to the bulk-trained model improved
transfer only modestly and again through specificity (0.772 → 0.814, F1 0.210 → 0.277), not sensitivity.

This confirms the second corollary of our hypothesis: the information exploited by the single-cell model
is not present in bulk, and a bulk-trained model is not a substitute. Model and assay must be matched.

### Bulk cohorts nonetheless augment single-cell training through a shared feature space

The negative result above (imputation adds no information within a cohort) does not imply that
imputation is useless. Because the same imputation models are applied to both bulk and single-cell
inputs, they produce **feature spaces with identical semantics across assay types**. This permits bulk
specimens to be pooled into single-cell training even though they cannot be scored by a single-cell
model.

Pooling 707 bulk specimens into the training partition of the shared blocks — with each cohort
standardised on its own statistics, and bulk specimens never entering a test fold — improved mean AUROC
from **0.889 to 0.908**, MCC from **0.403 to 0.432** and balanced accuracy from 0.722 to 0.729
(22 of 28 drivers improved; Wilcoxon P = 0.016; **Fig. 3c**). Gains were
concentrated precisely where single-cell positives were scarce: *WT1* +0.132, *KRAS* +0.087, *SRSF2*
+0.062, *CEBPA* +0.058, *IDH2* +0.056, *IDH1* +0.044. Where the single-cell model was already strong
(*TP53*, *TET2*), augmentation was neutral to marginally negative.

Thus imputation's value is **as a bridge between assay types**, not as a source of per-specimen
information — a distinction that reconciles our positive and negative findings and that we suggest is
generalisable to other multimodal integration problems.

### Calibration, clinical utility, and failure-mode characterisation

Raw fused scores were poorly calibrated (expected calibration error 0.455; Brier 0.277). We compared
calibration maps under identical nested donor-grouped folds: isotonic regression reduced ECE to 0.012
(Brier 0.037, log-loss 0.201) but produced a visibly unstable step function, consistent with
over-fitting at these positive counts, whereas **Platt (sigmoid) scaling reduced ECE to 0.006, Brier to
0.035 and log-loss to 0.134** — better on every criterion. We therefore deploy Platt scaling
(**Fig. 5a**), after which a reported probability corresponds to observed frequency across the range the
model actually occupies. Calibrated probabilities rarely exceed 0.3, which is the correct behaviour for
lesions of ~5% prevalence: a well-calibrated model should not assert high certainty for rare events. Decision-curve analysis confirmed positive net benefit
relative to both treat-all and treat-none strategies across the clinically plausible threshold range for
every major driver (**Fig. 5b**).

Failure modes were characterised rather than merely reported. Sensitivity increased monotonically with
VAF (0.29 subclonal → 0.69 clonal; **Fig. 5c**), establishing that missed calls concentrate in
subclonal lesions and represent a detection limit rather than stochastic error. Feature attribution
recovered established biology (*KMT2A* → *ZNF521*, *HOPX*; *DNMT3A* → *SPINK2*), supporting that the
classifiers exploit lesion-associated programmes rather than technical artefacts.

### External validation and honest accounting of label quality

The platform was evaluated on 367 external bulk specimens (Leucegene; 92.4% gene-level accuracy) and on
60 held-out single-cell specimens spanning three cohorts, comprising 2,046 scored calls: pooled
sensitivity 0.492, specificity 0.899 (**Table 1**). Performance was highest on the sealed internal
held-out set (0.505/0.931) and lowest on the most recently added fully external cohort (0.367/0.850), an
expected and informative gradient.

We emphasise two methodological points that materially affect such numbers. First, restricting
evaluation to categories with ≥ 3 positives — a common convention — discarded most evaluable lesions and
produced a misleadingly narrow denominator; pooling every scored call raised the number of assessable
drivers from 3 to 22 in one cohort. Second, where a cohort documents which genes its panel actually
assayed, genes never assayed must be left unlabelled rather than counted as wild-type. Applying this
correction *lowered* apparent specificity on one cohort from 0.881 to 0.775 — the honest value.

### Label quality is a first-order determinant of apparent performance

During development we encountered a defect that we believe is under-appreciated in the field and that
materially changed our results. One contributing cohort of 28 specimens carried zeros across every
mutation column; the harmonised metadata subsequently supplied **149 true positives** for those same
specimens. Because a zero is indistinguishable from "not reported" in most matrix representations, those
149 positives had been consumed as *true negatives* — simultaneously injecting label noise into training
and penalising correct positive calls as false positives.

Correcting the labels raised mean AUROC from 0.885 to 0.889 and mean F1 from 0.356 to 0.413, moved
*IDH2* from 0.730 to 0.850 as its positive count rose from 14 to 35, and rendered *SF3B1* and *SRSF2*
— both ELN 2022 myelodysplasia-related genes, and therefore directly risk-defining — trainable for the
first time. Twenty-five of the 28 specimens had matched expression data and thus contributed
immediately.

We draw a general methodological conclusion. Any harmonised genomic resource in which "absent" and "not
assayed" share an encoding will systematically understate specificity and corrupt training, and the
magnitude of that effect is not estimable without panel-coverage metadata. Where such metadata exist,
they should be used to leave unassayed genes unlabelled; where they do not, reported specificity should
be presented as a lower bound. We adopt both conventions here.

### Complementarity with ELN 2022 risk assignment

To position the platform relative to the clinical standard rather than in isolation, we re-derived ELN
2022 across a 942-patient bulk cohort using curated fusion calls and variant-level annotation, avoiding
free-text karyotype parsing. Assignment succeeded for 440 patients (Adverse 230, Favourable 115,
Intermediate 95). Among the 417 with an unambiguous ELN 2017 call, 87.3% were concordant and 53 were
reclassified; every reclassification was attributable to a documented guideline change, dominated by
intermediate → adverse transitions driven by myelodysplasia-related lesions (n = 24) and favourable →
intermediate transitions reflecting revised *CEBPA* and *FLT3*-ITD rules (n = 16). A further 23 patients
whose ELN 2017 status was ambiguous were resolved to a definite category.

The clinically salient number, however, is the 502 patients for whom **no** assignment was possible
because the required genotype was unavailable. This is the population in which an expression-based
nomination has potential utility: not to assign risk directly, but to indicate which lesions are likely
present and therefore which confirmatory assays would be informative. We note that the drivers on which
augmentation most improved detection — *SRSF2*, *WT1*, *CEBPA*, *IDH1* — include several that are
directly risk-defining under ELN 2022, so improvements in their detection propagate to risk assignment
rather than remaining internal metrics.

We deliberately do **not** propose feeding predicted lesions into ELN 2022 as though they were sequenced
calls. Doing so would import the model's error structure into a clinical framework calibrated on
genotype. The appropriate use is triage of confirmatory testing.

### Cell-state localisation and interpretability

Because inference is cell-state resolved, each prediction carries a compartment attribution: the marrow
states in which the lesion-associated signal concentrates, scored independently on measured RNA. Feature
attribution recovered established biology — *KMT2A* rearrangement weighting *ZNF521* and *HOPX*,
*DNMT3A* weighting *SPINK2* — supporting that classifiers exploit lesion-associated programmes rather
than batch or technical structure. inv(16) was driven overwhelmingly by a single modality block,
consistent with a highly stereotyped and well-separated phenotype, which is congruent with its
near-perfect discrimination.

This attribution layer is the component with no counterpart in genotype-only frameworks. It converts a
binary lesion call into a statement about *where in haematopoiesis* the lesion is manifest, which is the
axis most relevant to differentiation-directed and stem-cell-directed therapeutic strategies.

### Current limitations, stated explicitly

We enumerate the limitations of the present version so that the claims above are bounded correctly.

**Sensitivity is the binding constraint.** Mean sensitivity at the deployed operating point is 0.505:
the platform recovers roughly half of true lesions at 95% specificity. It is therefore a rule-in
instrument. Sensitivity varies widely by lesion (inv(16) 0.91 and *NPM1* 0.75 at one extreme;
*KMT2A*-rearrangement 0.33 and *SRSF2* 0.33 at the other), and low-sensitivity lesions should be
interpreted only in the positive direction.

**Detection is clonality-limited.** Sensitivity falls to 0.29 for lesions below 10% VAF. Subclonal
disease is systematically under-detected, and this cannot be corrected by thresholding.

**Accuracy is below the trivial baseline by design.** As set out above, the deployed operating point
sacrifices ~4.5 accuracy points to recover lesions. This is the correct trade for a triage application
but must not be presented as accuracy superiority.

**Some lesions remain weak.** *KRAS* (AUROC 0.804), trisomy 8 (0.769) and *WT1* achieved enrichment of
only 3.5–3.7-fold over chance. These are reported as leads requiring sequencing, not as calls.

**Held-out cohorts are small and heterogeneous.** The sealed single-cell held-out set is 29 specimens;
per-lesion held-out estimates therefore carry wide intervals. Performance on the most recently added
fully external cohort (0.367 sensitivity, 0.850 specificity) is materially below internal held-out
performance, which is the honest expectation for external transfer and should be treated as the more
conservative estimate of real-world behaviour.

**Ground truth is imperfect.** Contributing centres used different sequencing panels, and where
panel-coverage metadata are unavailable a gene that was never assayed is indistinguishable from one
assayed and found wild-type. Our single-cell specificity estimates are therefore lower bounds. Two
modality blocks could not yet be pooled across assay types owing to feature-naming and bundle-version
mismatches.

**All validation is retrospective.** No prospective confirmation has been performed; this is the
principal aim of the proposed work.

### An interpretable, ELN 2022-aligned decision interface

Outputs are surfaced as a per-patient report: control-gate call, ranked driver predictions with
calibrated probability, 95% confidence interval, permutation significance and honest nested operating
characteristics; the cell states in which each lesion's signal localises; and rule-derived therapeutic
hypotheses and confirmatory tests, with all pharmacological suggestions constrained to a curated
lesion → agent mapping. Predicted lesions are explicitly labelled as requiring sequencing confirmation
and are never presented as established genotype. A companion validation interface exposes every
performance claim, including the negative results above, so that a reviewer or clinician can inspect the
evidence for any individual call rather than accepting an aggregate metric.

---

## Discussion

We set out to test whether the cellular consequences of AML driver lesions are detectable at single-cell
resolution and carry information beyond bulk expression. Both propositions are supported. The multimodal
model outperformed an identically-configured bulk model by a wide and statistically robust margin, every
driver exceeded a permutation null, and — most importantly for interpretation — the advantage survived
removal of all RNA-imputed modalities, localising it to measured single-cell structure.

The clinical framing we propose is deliberately conservative. MOSAIC-AML is **not** a replacement for
sequencing, and its outputs are not ELN 2022 inputs in their own right. Its value is threefold: to
nominate lesions when genotype is incomplete, so that confirmatory assays can be directed; to supply the
cell-state context that ELN 2022 structurally omits; and to make explicit the clonality-dependence that
genotype-only frameworks leave implicit. A lesion that MOSAIC-AML detects with high confidence but that
sequencing did not report is, in the first instance, a hypothesis about panel coverage or subclonal
fraction — precisely the two failure modes that most degrade ELN 2022 in practice.

Several limitations bound these claims. Sensitivity remains near 0.5 at high specificity, so the platform
is better suited to ruling lesions *in* than out. The single-cell held-out set is small, and per-driver
estimates carry wide confidence intervals. Ground truth is heterogeneous across contributing centres, and
in the absence of complete panel-coverage metadata our single-cell specificity estimates should be read
as lower bounds. Imputed modalities, although shown here to be non-circular in aggregate, remain
RNA-conditioned and should not be interpreted as independent measurement. Finally, and most importantly,
all validation is retrospective; prospective re-sequencing of specimens for which the model predicts a
lesion absent from the existing call set is the decisive experiment, and is the principal aim of the
proposed work.

Methodologically, we would highlight the imputation result as generalisable. That imputed modalities
improved within-cohort performance no more than a random nonlinear transform, yet materially improved
cross-assay integration, suggests that the field's evaluation of imputation-based multimodal methods
requires matched nonlinear controls, and that the primary value of such imputation may lie in
harmonising heterogeneous assay types rather than in enriching individual specimens.

### Relation to existing approaches

Expression-based inference of mutation status is not new, but prior efforts have overwhelmingly operated
on bulk profiles, where the compartment-specific consequences of a lesion are averaged away. Our
head-to-head comparison quantifies the cost of that averaging: an identically-specified pipeline given
bulk RNA reached AUROC 0.744 against 0.908 for the cell-state-resolved model, and a bulk-*trained* model
applied to single-cell data reached only 0.695. The distinction between these two numbers matters. The
first isolates the information content of the representation; the second isolates transferability. Both
favour assay-matched, compartment-resolved modelling, and together they argue that bulk-derived
mutation-inference tools should not be redeployed on single-cell data without revalidation.

Our approach also differs in what it declines to do. We do not output a risk category. Given that ELN
2022 is calibrated on sequenced genotype, substituting inferred lesions would import an uncharacterised
error structure into a clinical instrument. We therefore restrict outputs to lesion nomination with
calibrated uncertainty, cell-state attribution, and confirmatory-test recommendations.

### Implications and proposed directions

Three lines follow directly from these results.

*Prospective confirmation.* The decisive experiment is targeted re-sequencing of specimens for which the
model nominates a lesion that the existing call set does not report, stratified by predicted probability.
Our calibration analysis makes this tractable: because reported probabilities now correspond to observed
frequencies, a defined number of confirmations is sufficient to test whether high-confidence
model-positive/genotype-negative calls represent panel gaps and subclonal disease rather than model
error. Given the VAF dependence we characterise, we would predict that confirmed cases are enriched for
low clonal fraction.

*Completing ascertainment.* Two of the eight modality blocks could not be pooled across assay types
because of feature-naming and bundle-version mismatches, and several contributing cohorts still lack
panel-coverage metadata. Both are tractable data-engineering problems whose resolution would extend
augmentation and convert our specificity estimates from lower bounds to point estimates.

*Extension beyond driver detection.* The cell-state attribution layer is, at present, descriptive. The
same architecture supports prediction of therapeutically actionable phenotypes — stem-cell compartment
size, differentiation blockade, and response-associated states — for which the compartment resolution is
not merely helpful but necessary. These endpoints are label-limited in our current cohort, and expanding
labelled outcome data is the principal barrier.

More broadly, we suggest the design pattern demonstrated here is portable: represent a specimen at the
resolution at which the biology acts, evaluate with grouping that respects the sampling structure, and
subject every apparent gain to a null that could plausibly explain it. In this study that discipline
overturned two of our own initial conclusions — that imputation added information, and that selective
rather than blanket cross-assay augmentation would be superior — and we regard the resulting estimates as
correspondingly more trustworthy.

---

## Methods

### Reference atlas and modality construction

The reference comprises 383 bone-marrow specimens (338 AML, 45 non-leukaemic control) aggregated into
12,255 specimen × cell-state pseudobulks across 89 annotated marrow states defined by projection onto a
titrated haematopoietic reference. Aggregating to pseudobulk at the cell-state level, rather than
modelling individual cells, controls for the extreme sparsity of single-cell counts while preserving the
compartment resolution that motivates the approach; within a cell state there is exactly one pseudobulk
per specimen, so a per-(modality × cell-state) matrix is a specimen-level matrix on a row subset.

Four modality blocks are measured: RNA (restricted to cell-state marker programmes and log-transformed),
cell-state composition (fractional abundance over the 89 states), leukaemic stem-cell scores, and
cell–cell communication (ligand–receptor interaction strengths). Four are imputed from RNA by
pre-trained per-target regression models: antibody-derived tags, lipid species, metabolites, and
gene-regulatory-network edge scores. Imputed features are filtered on held-out Spearman fidelity, and
are labelled as RNA-conditioned wherever they are displayed, since they do not constitute independent
measurement.

### Classification and fusion

For each (driver, modality) pair, features are standardised, reduced by differential selection to the
500 most discriminative, and classified with a linear support-vector machine (C = 0.02, balanced class
weights). Decision values are converted to within-cohort percentiles, which makes scores comparable
across modalities of very different dimensionality and scale. Per-driver fusion weights are obtained by
ridge-regularised non-negative least squares of the per-modality out-of-fold percentiles onto the label
vector, floored to the best single modality so that fusion can never underperform its strongest
component, with a uniform fallback across informative modalities when the fit is degenerate. Only
drivers with at least eight positive and eight negative specimens are modelled.

### Evaluation, uncertainty and nulls

All cross-validation is **donor-grouped**: specimens from one patient, including diagnosis–relapse pairs
and specimens contributed under different dataset identifiers, are constrained to a single fold. A sealed
held-out set is excluded from training entirely.

Decision thresholds are selected under **nested** cross-validation — F1 is maximised on the donor folds
not being scored and the resulting threshold applied to the held fold — which removes the optimism
introduced by selecting an operating point on the same data used to report it. In our hands this
correction reduced apparent mean F1 from 0.44 to 0.36, and we report the corrected values throughout.

Uncertainty is quantified by donor-level bootstrap resampling (B = 1,000), resampling donors rather than
specimens so that intervals reflect the effective sample size. Statistical significance is assessed by
label permutation (P = 1,000) against the fixed fused score, yielding a one-sided empirical *P*. Because
fusion weights are held fixed during permutation, this test is anti-conservative with respect to the
fusion step and conservative with respect to the classifier; we state this explicitly rather than
implying an exact test.

Probability calibration is fitted on held-out donor folds and applied to the fold it did not see. We
compared isotonic regression and Platt (sigmoid) scaling under identical folds and selected Platt on
nested expected calibration error, Brier score and log-loss (0.006/0.035/0.134 versus 0.012/0.037/0.201);
isotonic over-fits at these positive counts. Expected calibration error is computed on equal-count
(quantile) bins so that sparsely-populated high-probability bins do not dominate the estimate. Clinical utility is assessed by
decision-curve analysis, computing net benefit across threshold probabilities against treat-all and
treat-none reference strategies.

**Choice of performance metrics.** Median driver prevalence is 3.6% (mean 5.6%), and metric selection
under such imbalance is consequential rather than cosmetic. Accuracy is degenerate here — the
all-negative classifier scores 94.4% — and is reported only alongside that baseline. The F1 score
ignores true negatives and is bounded above by prevalence, so it is neither comparable across lesions of
differing rarity nor interpretable in absolute terms; we report it for continuity with the literature,
always with prevalence stated, but base no comparative claim on it. Primary metrics are therefore
AUROC (prevalence-invariant ranking quality); AUPRC reported against its prevalence baseline and
expressed as fold-enrichment, since a high AUPRC on an enriched cohort can otherwise be trivial;
the Matthews correlation coefficient, which incorporates all four confusion-matrix cells and is
recommended for imbalanced binary classification; and balanced accuracy. We verified that the
qualitative conclusions of this study — the modality ablation ordering, the superiority of
cell-state-resolved over bulk representation, and the effect of cross-assay augmentation — are invariant
to which of these imbalance-aware metrics is used.

### Cross-assay augmentation

Bulk specimens are projected through the identical imputation models used for the atlas, producing
feature spaces with matching semantics. Features are matched by normalised name; each cohort is
standardised on its own statistics before pooling, which is the standard correction for scale-level
domain shift and mirrors the cohort-matched score referencing used by the bulk caller. Bulk specimens
enter only the training partition of each fold and never a test fold, so no bulk specimen contributes to
any reported estimate of single-cell performance.

### Controls against confounding

Three controls address the principal alternative explanations. Donor grouping addresses within-patient
leakage. A matched random nonlinear transform of the input RNA, of identical width to the imputed blocks,
addresses the possibility that imputed modalities contribute only nonlinear feature expansion. Label
permutation addresses the possibility that discrimination reflects cohort structure rather than lesion
biology.

### Reporting

The platform is documented against TRIPOD+AI.[3] Frozen model cards, split definitions, random seeds,
package versions and a reproducibility manifest accompany the code, and a validation interface serves
every reported figure, per-driver statistic and negative result.

### Data availability

The reference atlas, external bulk cohorts and external single-cell validation cohorts are derived from
published and controlled-access resources; accession identifiers and the derived label matrices required
to reproduce every reported analysis accompany the code release.

---

## Display items

- **Fig. 1** Per-driver sensitivity versus specificity across models and cohorts; each point one driver, point size proportional to available positives.
- **Fig. 2** Per-driver AUROC with 95% donor-bootstrap confidence intervals against the chance line; permutation significance annotated.
- **Fig. 3** (a) Modality ablation ladder isolating the imputation-independent component; (b) imputation versus matched random nonlinear control; (c) effect of bulk augmentation.
- **Fig. 4** Two models × four cohorts: assay-matched performance and the bulk → single-cell transfer gap.
- **Fig. 5** (a) Calibration before and after nested isotonic regression; (b) decision-curve net benefit; (c) sensitivity stratified by variant allele frequency.
- **Table 1** Pooled held-out performance across three single-cell cohorts (60 specimens, 2,046 scored calls).

## References

1. Döhner, H. *et al.* Diagnosis and management of AML in adults: 2022 recommendations from an
   international expert panel on behalf of the ELN. *Blood* **140**, 1345–1377 (2022).
   doi:10.1182/blood.2022016867
2. Rausch, C. *et al.* Validation and refinement of the 2022 European LeukemiaNet genetic risk
   stratification of acute myeloid leukemia. *Leukemia* **37**, 1234–1244 (2023).
   doi:10.1038/s41375-023-01884-2
3. Tyner, J. W. *et al.* Functional genomic landscape of acute myeloid leukaemia. *Nature* **562**,
   526–531 (2018). doi:10.1038/s41586-018-0623-z
4. van Galen, P. *et al.* Single-cell RNA-seq reveals AML hierarchies relevant to disease progression
   and immunity. *Cell* **176**, 1265–1281.e24 (2019). doi:10.1016/j.cell.2019.01.031
5. Lavallée, V.-P. *et al.* RNA-sequencing analysis of core binding factor AML identifies recurrent
   ZBTB7A mutations and defines RUNX1-CBFA2T3 fusion signature. *Blood* **127**, 2498–2501 (2016).
   [Leucegene cohort]
6. Stoeckius, M. *et al.* Simultaneous epitope and transcriptome measurement in single cells.
   *Nat. Methods* **14**, 865–868 (2017). doi:10.1038/nmeth.4380
7. Collins, G. S. *et al.* TRIPOD+AI statement: updated guidance for reporting clinical prediction
   models that use regression or machine learning methods. *BMJ* **385**, e078378 (2024).
8. Chicco, D. & Jurman, G. The advantages of the Matthews correlation coefficient (MCC) over F1 score
   and accuracy in binary classification evaluation. *BMC Genomics* **21**, 6 (2020).
   doi:10.1186/s12864-019-6413-7
9. Vickers, A. J. & Elkin, E. B. Decision curve analysis: a novel method for evaluating prediction
   models. *Med. Decis. Making* **26**, 565–574 (2006). doi:10.1177/0272989X06295361

**Citations still required before submission** (topic → what to cite). These are deliberately left
unfilled rather than guessed; each should be retrieved from PubMed and verified:

10. Prior expression-based inference of AML mutation status from bulk RNA — the Discussion claims
    priority over this literature and must cite it directly.
11. WHO 5th edition and ICC 2022 haematolymphoid classifications (companion frameworks to ELN 2022).
12. AML measurable/minimal residual disease assessment and its clonality dependence.
13. Single-cell reference mapping / label-transfer methodology (the 89-state projection).
14. The rna2ADT / rna2GRN / rna2metabolite / rna2lipid imputation models used here (methods papers).
15. UDON / cellHarmony / AltAnalyze3 component methods.
16. Platt scaling and isotonic calibration primary sources.
17. Nested cross-validation and selection-bias literature.
18. Venetoclax-based and differentiation-directed therapy in AML (motivates the cell-state axis).
19. Additional external single-cell AML cohorts used for validation (Trumpp/Waclawiczek; GSE281087).
