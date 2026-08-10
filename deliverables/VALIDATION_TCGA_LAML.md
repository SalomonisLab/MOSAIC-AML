# External validation — TCGA-LAML

**Roadmap item B2.** Every survival number reported before this came from BeatAML2: one cohort, one
sequencing pipeline, one set of clinical conventions. Patient-grouped CV and a sealed hold-out protect
against learning a *patient*; neither protects against learning a *cohort*. This is the test that
separates the two.

    reproduce with:  python pipeline/validate_tcga_laml.py

| | |
|---|---|
| cohort | TCGA-LAML, UCSC Xena hub (`HiSeqV2`, Illumina HiSeq RSEM, gene symbols) |
| patients scored | **149** (173 with expression ∩ 186 with follow-up, one specimen each, follow-up > 0) |
| events | 92 deaths, 57 censored |
| follow-up | median 1.00 y, maximum 7.8 y |
| gene overlap | 12,733 / 14,237 of the model's space (**89.4 %**) |
| cytogenetic risk | 90 intermediate, 31 favorable, 26 adverse, 2 unknown |

## What was frozen and what was not

Loaded unchanged from `survival_model.pkl`: **the Cox coefficients, the PCA rotation, the variable-gene
selection, and the NNLS fusion weights.** No refitting, no recalibration against TCGA outcomes.

Refit on TCGA: **the per-gene z-reference (median/MAD) only.** RSEM log2 counts and BeatAML's units are
not the same scale, and a model applied across platforms without cohort-matched normalisation measures
the platform rather than the patient. This is the same cohort-matched-reference mechanism the drug
layer already uses for single-cell input — a normalisation, not a predictor.

## Result

| arm | blocks | C-index | 95 % CI | AUC 1 y | AUC 2 y | AUC 5 y |
|---|---|---|---|---|---|---|
| **deployed** | rna+state+mut+clin+age-spline+tx | **0.706** | 0.655 – 0.758 | 0.742 | 0.806 | 0.865 |
| full | rna+state+mut+clin | 0.689 | 0.636 – 0.741 | 0.720 | 0.787 | 0.859 |
| molecular | rna+state+mut | 0.669 | 0.607 – 0.728 | 0.716 | 0.756 | 0.727 |
| age + cytogenetics | — | 0.670 | 0.613 – 0.729 | 0.696 | 0.752 | 0.886 |
| clinical block | — | 0.664 | 0.608 – 0.723 | 0.681 | 0.735 | 0.889 |

**Gain over age + cytogenetic risk**, bootstrapped over patients (both arms scored on the same people):

| arm | ΔC | 95 % CI | P(ΔC ≤ 0) |
|---|---|---|---|
| deployed | **+0.035** | −0.001 – +0.073 | **0.029** |
| full | +0.018 | −0.023 – +0.058 | 0.182 |
| molecular | −0.003 | −0.070 – +0.063 | 0.522 |

**Risk tertiles** (deployed arm, cut at the TCGA tertiles of predicted risk):

| group | n | deaths | 2-y survival | median survival |
|---|---|---|---|---|
| low | 50 | 17 | **71.7 %** | not reached |
| intermediate | 49 | 31 | 40.9 % | 1.42 y |
| high | 50 | 44 | **13.7 %** | 0.75 y |

Log-rank low vs high: χ² = 38.5, **p = 5.4 × 10⁻¹⁰**. Curves are ordered and do not cross
(`fig_Sv4_tcga_km.png`).

## Reading it honestly

**It generalises.** 0.706 in an independent cohort against 0.751 CV / 0.787 sealed hold-out in BeatAML.
A drop of that size across a different platform, different institutions and a decade's difference in
treatment era is the expected cost of transfer, not evidence of overfitting — the ranking survives
intact and the tertile separation (71.7 % vs 13.7 % at two years) is larger than most of what the
clinical variables achieve alone.

**The molecular blocks alone do not beat age + cytogenetics here** (ΔC = −0.003). That is a real
negative and it is not hidden: the added value in TCGA comes from combining molecular and clinical
information, not from expression on its own. The same pattern holds in BeatAML, where `molecular` runs
~0.04 below `full`.

**The one significant gain is the deployed arm** (+0.035, P = 0.029) — nominally, one-sided, on 149
patients with 92 events. It is a genuine result, and it is also a single comparison at the edge of
significance in a modest cohort. It should be read as *consistent with* the BeatAML finding, not as
independent confirmation of the effect size.

### Limitations, stated rather than buried

1. **Mutation calls are absent** from this download, so the `mut` block enters as unreported (all-zero
   with the trailing fraction-missing column set to 1). This is exactly the degraded path a real upload
   without a mutation caller takes, so the numbers are honest for that input — but the fused arms are
   running one block short of what BeatAML gave them.
2. **TCGA-LAML is treatment-homogeneous** (essentially uniform standard induction). The +0.035 C-index
   that baseline induction type contributes in BeatAML therefore *cannot* be validated here; the
   `deployed` arm's treatment block is constant across all 149 patients and contributes nothing. Its
   advantage over `full` in this table comes from the age spline alone.
3. **CALGB cytogenetic risk is mapped onto ELN 2017** (Favorable → Favorable, Intermediate/Normal →
   Intermediate, Poor → Adverse). Close, not identical, so the bar being cleared is approximate.
4. **Follow-up is short** — median 1.00 y. The 5-year AUCs rest on few patients still at risk, and the
   clinical arms' apparently strong 5-y AUC (0.886) sits on that thin tail.
5. Xena's `HiSeqV2` is log2(RSEM+1); it is exponentiated back to linear before the reference is built,
   because `FeatureSpace.add_reference` expects linear space as it does for BeatAML.

## Related correction found while running this

The bundle did not store the age-spline knots, so `survival_layer.py` fell back to a hardcoded
`[45, 60, 71]` while the model had been fitted against the training quartiles — every upload rebuilt
the spline basis in slightly the wrong place. Separately, `build_blocks` accepted a `train_idx`
argument that **neither call site passed**, so knots were placed using quantiles computed over the
hold-out patients as well as the training ones.

Both are fixed. The knots are now computed on training patients only and carried in the bundle
(`age_knots = [47.0, 61.0, 70.0]`; the whole-cohort quantiles would have been `[47.0, 61.0, 70.5]`,
which is how we know the argument is now actually in effect — the last time an ignored parameter
produced byte-identical numbers in this project it was a real bug).

The model was retrained with the leak closed. **Every metric is unchanged to four decimal places**
(CV 0.751, sealed hold-out 0.787, AUC 2 y 0.8715, Brier 0.1494), which says the leak was immaterial —
three quantiles of a single covariate across 444 patients. The fix is still correct, and the
inference-time miscalibration it removes was the part that actually mattered.
