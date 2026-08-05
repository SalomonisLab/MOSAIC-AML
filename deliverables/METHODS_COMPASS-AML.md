# COMPASS-AML — predicting ex-vivo drug response from the transcriptome

The drug layer of **MOSAIC-AML** (*Multimodal Omics and State-Aware Inference of Cancer Drivers in
Acute Myeloid Leukemia*). It answers a different question from the mutation caller and is deliberately
built as a **parallel** layer, not a downstream one: mutation calls enter as evidence, never as a drug
look-up table.

---

## 1. What the BeatAML2 inhibitor data are, and are not

`beataml_probit_curve_fits_v4_dbgap.txt` holds one fitted dose–response curve per (specimen ×
inhibitor) from the Beat AML **ex-vivo** screen: fresh mononuclear cells exposed to seven graded
concentrations plus an untreated control, fitted with a probit model. These are **not** treatments
given to patients. The response variable is the area under the fitted viability curve (`auc`, on a
0–300 scale here) — **lower AUC = more sensitive**.

That number is not comparable across inhibitors: concentration ranges, potency and cohort
distributions all differ. Every use below is therefore of a **within-drug standardised** value, as the
original Beat AML heat maps did it.

**Cohort actually used.** 63,395 single-agent curves → 53,571 have an RNA-seq specimen present in the
14,237-gene shared feature space → after curve QC and inhibitor inclusion, **48,998 measurements over
520 specimens from 479 patients and 118 inhibitors**.

---

## 2. Curve QC, normalisation and response classes

**QC flags** are attached per row and the conjunction is `qc_pass`; nothing is dropped silently.
Excluded: non-convergent fits (10), an off-panel concentration window (10), and a within-drug extreme
deviance above that inhibitor's own 99th percentile (603). *Increasing* curves — viability rising with
dose, 7,092 rows — are **kept** (they are genuine non-responders) but flagged, and can never be
labelled sensitive.

**Within-drug normalisation** is robust (median / 1.4826·MAD) rather than mean/SD, because the AUC
distribution for a potent inhibitor is strongly left-skewed and a handful of exquisitely sensitive
specimens would otherwise set the scale. The modelled quantity is `sens = −z(AUC)`, so larger always
means more sensitive.

**Response classes** are the tails, not a global cutoff: **sensitive** = bottom 20% of within-drug AUC,
**resistant** = top 20%, and the middle 60% is `indeterminate` — retained for regression, excluded from
the binary task. One universal AUC cutoff across 118 inhibitors would mostly encode drug potency.

**Inhibitor inclusion.** Of the 165 screened compounds, **115 enter the primary panel**
(≥120 specimens, ≥110 patients, IQR ≥ 20 AUC units, ≤60% increasing curves, ≥25 specimens in each
tail, |wave shift| ≤ 0.6 SD). Three more — **Cytarabine, Nutlin 3a, GDC-0941** — fail *only* the
acquisition-wave stability test; rather than delete the induction backbone we place them in a
`wave_conditional` tier, model them, and attach the caveat everywhere they are reported. 47 are
excluded outright, with the reason recorded per compound.

---

## 3. Model A — the patient-level response model

    sens_ij  =  f_group(X_i, D_j)  +  w_j · f_j(X_i)

**Hierarchical**, because 118 inhibitors with ~400 trainable specimens each is exactly the regime where
118 independent models overfit and one pooled model underfits. `f_group` is a ridge fitted over every
(specimen × inhibitor) row in a **target-pathway family group** (16 groups: RTK, FLT3, MAPK,
PI3K/AKT/mTOR, JAK-STAT, cell cycle, epigenetic, apoptosis, …), with a **drug descriptor** — the
z-expression of that inhibitor's own annotated target genes — supplied as interaction features. It
therefore learns "what kind of leukaemia is sensitive to RTK inhibition" and can say something about an
inhibitor it has little data for. `f_j` is a per-drug ridge on the residual, shrunk by
`w_j = n_j/(n_j + 150)`.

**Four feature blocks, late-fused.**

| block | contents |
|---|---|
| `rna` | 100 PCs of the cohort-z-scored transcriptome (top-4,000-variance genes) + the drug's target descriptor |
| `state` | 10 atlas-derived lineage-signature scores + an explicit primitive→mature axis |
| `mut` | 58 variant-level driver categories (observed WES for BeatAML; the mutation layer's calibrated calls at inference) |
| `clin` | age, sex, BM/PB blast %, WBC, specimen type, disease stage, relapse/de-novo/transformed, prior MDS, ELN 2017 |

A single ridge over the concatenated blocks **loses to RNA alone** (0.742 vs 0.766 mean AUROC): the
blocks differ by orders of magnitude in dimensionality and scale, so one shared penalty over-penalises
the small ones. Late fusion is the remedy the mutation layer already uses — per-block ridges combined
by **non-negative least squares on inner donor-grouped OOF predictions**, floored to the best single
block so fusion can never underperform its best part. That recovers **0.772**.

**The `state` block.** BeatAML2's own analysis found differentiation state to be a broad determinant of
ex-vivo response, so the model is given that axis explicitly rather than left to hope the RNA PCs
encode it. Full 89-state deconvolution of bulk was previously shown to be too collinear to trust, so
the 89 atlas states are collapsed into 10 lineage groups and a marker signature is learned per group
from the atlas itself (top 60 genes by specificity over the best competing lineage). The signatures
recover textbook markers unprompted — AVP/CRHBP/MECOM for HSC/MPP, FLT3 for LMPP/CLP,
PRTN3/ELANE/AZU1/MPO for granulocytic, S100A8/9/FCN1 for monocytic, HBB/HBA/AHSP for erythroid,
PPBP/PF4/ITGA2B for MEP/Mk, CXCL12/VCAM1/FN1 for stroma.

**Calibration** is per-drug Platt scaling fitted on the OOF **percentile**, not the raw score. This is
the cross-assay lesson the mutation caller learned the hard way: a single-cell bulk-equivalent occupies
a different band of the raw decision scale than a BeatAML bulk specimen, so a Platt curve fitted on raw
BeatAML scores returns ≈0.99 for essentially every inhibitor on a single-cell sample — a model that
enthusiastically recommends all 118 drugs. Percentiles are cohort-invariant provided the reference is
matched, so three references are stored: `beataml` (bulk specimens), `sc_sample` (387 single-cell
bulk-equivalents) and `sc_state` (6,000 individual cell-state pseudobulks). A cell state's profile is
systematically more extreme than any whole sample; scored against the sample-level reference every
state looks like an outlier.

---

## 4. Model B — the state-resolved layer

Because Model A consumes one expression vector in the shared gene space, it can be handed a single
cell state's pseudobulk exactly as it would be handed a whole sample. Per (patient × inhibitor) that
yields: the prediction for each state; the abundance-weighted mean; **coverage** of the blast and
LSC-like compartments; the **most resistant state above 1% abundance** (the escape candidate);
dispersion across states; and `bulk_vs_sc`, the disagreement between the bulk view and the
single-cell-weighted one.

**"Blast compartment" is a compartment, not a malignant-cell call.** No per-cell somatic genotype is
used. The blast compartment is the set of states AML blasts occupy (HSC/MPP, LMPP/CLP, GMP,
granulocytic, monocytic, MEP/Mk, erythroid, DC); lymphoid and stromal states are treated as
presumed-normal bystanders and excluded from coverage. The LSC-like proxy is the primitive compartment.

---

## 5. Model C — the mechanistic model, kept separate

Model A learns "patients who look like this responded ex vivo"; Model C asks "is there a mechanistic
reason this compound should work here?" — target expression percentile, a short curated
**transcriptional output** readout of the target pathway (output genes, not pathway members, since a
kinase's mRNA says little about its activity), the anti-apoptotic balance
BCL2/(BCL2+MCL1+BCL2L1+BCL2A1), genetic activation of the target pathway tagged observed-vs-predicted,
and measurable proxies for the compound's curated bypass routes.

Fusing the two into one score would destroy exactly the signal that makes the pair useful: agreement
between an empirical and a mechanistic line of evidence is much stronger than either alone, and
disagreement is flagged as a conflict rather than averaged away.

---

## 6. The treatment-utility score

    S_ij = w1·P(sensitive) + w2·coverage_blast + w3·coverage_LSC + w4·mechanistic + w5·clinical
           − w6·uncertainty − w7·resistance − w8·infeasibility − w9·OOD

Each negative term exists because some real failure mode would otherwise rank first. **Uncertainty**
combines whether the inhibitor can be predicted at all (its own held-out AUROC), how decisive this
particular call is, how much data backs it, and the spread across cell states. **Resistance** combines
flagged bypass mechanisms with an explicit escape term — a bulk-sensitive prediction over a resistant
primitive subclone is a relapse, not a response. **Infeasibility** returns 0 and says "not assessed"
when no clinical covariates were supplied, rather than inventing a penalty. **OOD** is measured against
the *matched* cohort, so a single-cell patient is asked "are you unusual among single-cell AML
samples?" — the assay-level distance from BeatAML is reported separately as a fixed caveat, since
penalising every single-cell patient for it would just subtract a constant.

The positive part is renormalised by the weights that were **evaluable**, so a bulk patient with no
cell-state data is not silently penalised for missing coverage terms.

**Rankings are produced per clinical tier** — approved in AML / approved for another indication /
clinical-trial agent / research-only compound — never as one merged list. A research compound with a
clinically available analogue (AGI-6780 → enasidenib) is scored up and the analogue named.

---

## 7. The eight agents

All eight are `therapeutic`-domain and therefore **non-voting** with respect to the anchored subtype
call: a drug recommendation must never be able to reach back and alter the genetic anchor.

| agent | what it contributes |
|---|---|
| `drug_response` | predicted response, calibrated probability, percentile, uncertainty, and the **measured** AUCs of the 20 nearest BeatAML specimens — an observation a reader can check without the model |
| `cell_state_coverage` | which states respond, which escape, and the largest bulk-vs-single-cell disagreements |
| `molecular_mechanism` | Model C's itemised evidence |
| `pharmacology` | exposure plausibility, multi-kinase promiscuity, and — computed from the curve fits, not curated — whether sensitivity only appears in the top decade of the tested concentration range |
| `clinical_evidence` | tier separation and eligibility filtering |
| `combination` | complementary-coverage hypotheses across **different** pathways only; it explicitly refuses to add single-agent scores, because BeatAML2 contains no combination measurements |
| `skeptic` | nine specific challenges: differentiation-state confound, small training set, assay batch, modality conflict, curve quality, out-of-distribution, no clinical route, weak model, and predicted activity in normal haematopoiesis |
| `reporting` | the structured recommendation; no free narrative — every field traces to an evidence item |

---

## 8. Validation

**Split discipline.** Patients, never specimens (34 BeatAML subjects contributed >1 drug-tested
specimen). A sealed hold-out of ≥15% of *patients*, drawn once and never touched by any fitting,
selection, threshold or calibration step. Everything refits inside every fold — the cohort z-reference,
the PCA, the clinical imputation medians, the within-drug normalisation constants and the tail
cut-points. Normalising AUC against the whole cohort before splitting is the subtle leak that makes
ex-vivo response models look better than they are.

**Two tasks are evaluated and they are not interchangeable.** Per-drug across patients ("who responds
to this drug?") is what a biomarker paper reports. Per-patient across drugs ("which drug for this
patient?") is what the deployed system is actually asked, and it is harder — it requires predictions to
be comparable *across* inhibitors.

| | result |
|---|---|
| per-inhibitor mean AUROC (donor-grouped CV, 118 drugs) | **0.774** |
| per-inhibitor mean Spearman | 0.365 (median 0.362); 100% reach *p* < 0.05 |
| per-inhibitor mean AUPRC | 0.748 (prevalence baseline 0.475) |
| approved-agent subset (n = 42) | mean AUROC **0.809** |
| sealed hold-out (72 patients, 74 specimens, 6,648 measurements) | mean AUROC 0.784, Spearman 0.371 |
| per-patient top-1 retrieval | **34%** vs 10% matched chance (**3.4×**) |
| per-patient top-5 retrieval | 76% vs 42% chance |
| per-patient ranking concordance | mean Spearman 0.309 |
| calibration | ECE **0.012**, Brier 0.185 (baseline 0.249) |
| abstention | error 28% at full coverage → **4.7%** at 10% coverage; AUROC 0.795 → 0.963 |
| leave-wave-out | AUROC 0.722 / 0.733 |
| leave-centre-out (4 centres) | AUROC 0.731 – 0.890 |
| differentiation-state strata | AUROC 0.721 – 0.762 in every stratum |
| permutation null (specimen↔expression re-pointing, 100 shuffles) | observed 0.367 vs null 0.001 ± 0.019 — **19 null SDs**, *p* = 0.0099 |

**Best-predicted inhibitors** include Venetoclax (AUROC 0.977, AUPRC 0.970 vs 0.482 baseline),
Dasatinib 0.936, Tivozanib 0.927, Rapamycin 0.926, Trametinib 0.916, Ponatinib 0.916.

**Model B's falsifiable test.** Venetoclax response in AML is repeatedly reported to track
differentiation state — monocytic AML relatively resistant, primitive AML relatively sensitive. Model B
was fitted entirely on BeatAML **bulk** and never saw a cell-state label. Handed the 89 atlas cell
states across all 387 atlas samples, it predicts higher venetoclax sensitivity in primitive than in
monocytic states in **93.2%** of samples (mean +0.662 z, Wilcoxon *p* = 8.4 × 10⁻⁵¹) — and this is
**rank 1 of 118 inhibitors**, so it is not a generic "primitive states look sensitive to everything"
artefact. ABT-737, the other BH3 mimetic in the panel, ranks 8th. Internal consistency between the
bulk prediction and the abundance-weighted state mean is Spearman 0.970.

---

## 9. Limitations, stated plainly

- **Ex vivo is not clinical benefit.** Culture conditions, pharmacokinetics, the marrow
  microenvironment, combination therapy and toxicity are not represented. Nothing here is validated
  against patient outcomes.
- **The state-resolved layer has no ground truth.** No assay measured per-cell-state response, so
  Model B is checked against a pre-registered biological expectation and internal consistency, not
  accuracy.
- **Single-cell input is an assay transfer.** 100% of single-cell bulk-equivalents lie beyond the 95th
  percentile of BeatAML's own distance distribution. Cohort-matched references make the transfer
  well-behaved, but the transfer itself is not validated against measured ex-vivo response in
  single-cell-profiled patients — that would require a functional-precision-medicine cohort with paired
  scRNA.
- **Combinations are coverage hypotheses**, not synergy estimates.
- **47 of 165 inhibitors are excluded**, and 3 more carry a wave-instability caveat.
- **The mutation block is observed WES in training and predicted probabilities at inference.** Its
  average contribution is small (+0.006 AUROC over RNA alone), so this substitution is unlikely to
  matter much, but headline numbers are also reported without it.

---

## Reproducibility

| step | script |
|---|---|
| BeatAML2 curve QC, normalisation, inclusion | `pipeline/amlmm/drug/data.py` |
| inhibitor target / mechanism / clinical-tier curation | `pipeline/amlmm/knowledge/drug_annotation.tsv`, `amlmm/drug/targets.py` |
| lineage signatures from the atlas | `pipeline/build_state_signatures.py` → `labels/cellstate_signatures.json` |
| feature blocks | `pipeline/amlmm/drug/features.py` |
| Model A | `pipeline/amlmm/drug/model.py`, `pipeline/train_drug_model.py` |
| cohort-matched score references | `pipeline/build_drug_score_refs.py` |
| Model B | `pipeline/amlmm/drug/statemodel.py`, `pipeline/validate_state_response.py` |
| Model C | `pipeline/amlmm/drug/mechanism.py` |
| utility score | `pipeline/amlmm/drug/utility.py` |
| agents | `pipeline/amlmm/drug/agents.py` |
| validation battery | `pipeline/eval_drug_model.py` → `deliverables/drug_model_validation.json` |
| figures | `pipeline/build_drug_figures.py` → `deliverables/figures/Rx1–Rx6` |
| per-patient entry point | `pipeline/predict_drugs.py`, `pipeline/drug_layer.py` (pipeline hook) |
| GUI | `gui/therapy.html`, `gui/rx_validation.html` |

Data: `beataml_probit_curve_fits_v4_dbgap.txt`, `beataml_drug_families.xlsx`,
`beataml_wv1to4_clinical.xlsx`, `beataml_waves1to4_sample_mapping.xlsx` — all from the public BeatAML2
data repository (`github.com/biodev/beataml2.0_data`), the same source as the expression and mutation
tables already in `data/external/beataml/`.
