# CIPHER-AML — how the mutation-prediction models are actually used

Supplementary methodological detail, keyed to the four-layer architecture diagram (Layers A–D).

> **Two corrections to the diagram before anything else**, both worth fixing before it goes in a grant:
> 1. **Layer C, Calibration** reads *"OOF Youden's J Threshold"*. The deployed code no longer uses
>    Youden's J. It was replaced by **F1-maximisation**, because Youden (TPR − FPR) ignores base rate:
>    for drivers at ~3% prevalence it selects a low threshold that keeps FPR nominally small while
>    admitting many false positives. The threshold is additionally **nested** — chosen on held-out
>    donor folds, applied to the fold it did not see. (`train_predictor.py`, "was Youden's J".)
> 2. **Layer A** reads *"Bulk RNA NOT Supported"*. True of the single-cell ingest path, but the
>    platform now has two bulk relationships worth showing: a separate **BeatAML-trained bulk caller**
>    handles bulk specimens, and **707 BeatAML specimens are pooled into the training partition** of
>    the four blocks whose feature space is shared across assay types (bulk-equivalent RNA, ADT,
>    metabolite, GRN). Bulk is unsupported as *input to Layer A*, not absent from the system.

---

## 1. The fusion framework and modality weighting (Layer C)

**One classifier per (driver × modality), never a single joint model.** For each driver and each of the
eight modality blocks independently: `StandardScaler` → differential feature selection (top 500) →
`LinearSVC(C = 0.02, class-balanced)` → decision score → **cohort percentile** (rank against the
cohort's sorted scores). The percentile step is what makes fusion possible at all: blocks differ by
three orders of magnitude in dimensionality (4 LSC probabilities versus 14,237 genes) and in scale, and
percentiling puts every block on a common distribution-free footing before combination.

**Weights are fitted, not assigned.** Per-driver, the per-modality **donor-grouped 3-fold OOF**
percentiles are stacked and regressed onto truth by **ridge-regularised non-negative least squares**.
Non-negativity is deliberate — a modality may contribute nothing but never negatively, which would be
fitting noise. The ridge term stops one block absorbing all weight when blocks are correlated, which
they are, since four are imputed from RNA.

**Two guards prevent fusion from underperforming its parts.** The fused blend is compared on the same
OOF data against the best single modality; if it does not match or exceed it, weights collapse to that
single block ("floored"). If the NNLS fit is degenerate, the model falls back to a uniform mix over
modalities with individual CV AUROC > 0.52 ("uniform fallback"). No driver trains empty or ends worse
than its strongest single view.

**Weights are lesion-specific and biologically interpretable** — not a global modality ranking. GATA2 is
carried by surface protein (ADT weight 0.56; ADT alone AUROC 0.919 versus 0.75 for cell-state RNA),
FLT3-TKD by cell-state RNA (0.66), RUNX1 by metabolite (0.48), TET2 by GRN, ASXL1 by cell–cell
communication, SRSF2 by lipid. **Composition and LSC are never the best single modality for any driver**,
though both contribute inside the fusion — a negative result we report explicitly.

**Standalone strength and fusion weight need not agree.** Weights are fitted jointly, so a strong but
redundant block is down-weighted in favour of one adding independent information (for SRSF2, lipid wins
standalone at 0.846 while the fusion leans on ADT at 0.39). Supplementary Table 5 gives both quantities
for every driver so this is auditable rather than surprising.

**Calibration.** The fused percentile is mapped to probability by **Platt scaling** fitted on held-out
donor folds, selected over isotonic regression on nested ECE / Brier / log-loss
(0.006 / 0.035 / 0.134 versus 0.012 / 0.037 / 0.201); isotonic over-fits at these positive counts.

**Cell-type-resolved layer.** In parallel, 7,005 per-(mutation × modality × cell-state) OOF classifiers
localise each driver's signal to specific marrow compartments. This is an interpretation layer — it does
not alter the call, it explains where the evidence sits.

---

## 2. How the expert witnesses contribute (Layer D)

This is the part the diagram gets right and my earlier description got wrong, so it is worth stating
precisely: **the descriptive witnesses do not vote.**

**Voting is determined by DOMAIN, not by witness name.** The arbiter admits a weighted subtype signal
from exactly three domains — `genetic` (which also sets the anchor), `predictive`, and `cell_state`.
Every other entry contributes only to the *additive harvest*: biomarkers, validation claims and
descriptive context that feed therapies and recommendations but **can never change the anchored call**.

This distinction matters when reading the architecture diagram, because two witnesses are commonly
mislabelled:

| Witness | Domain | Votes? |
|---|---|---|
| `genetic` | `genetic` | **yes** — and sets the anchor |
| `composition` | `predictive` | **yes** (`ingest_patient.py`: `AgentResult("composition", "predictive", …)`) |
| `bulk_mutation` | `predictive` | **yes** |
| UDON / programme witnesses | `cell_state` | **yes** |
| LSC, surfaceome, metabolic, lipid, GRN-regulon, cell–cell communication | descriptive | **no** — additive harvest only |

So the **six** Phase B descriptive witnesses are the non-voting, anchor-invariant set. Composition and
the predictive/UDON witnesses do vote, subject to the weighting below and to the CV-gate (a
`predictive` witness that fails the permutation test is additionally down-weighted ×0.3 rather than
zeroed).

**Weighting of the voting witnesses.** Each carries a declared `grounding` (how the claim was obtained)
and `independence` (measurement versus inference), combined multiplicatively:

```
effective weight = reliability × GROUNDING_FACTOR × INDEPENDENCE_FACTOR × CV-gate
```

| GROUNDING | | INDEPENDENCE | |
|---|---|---|---|
| `deterministic_fact` (observed genotype) | 1.0 | `independent` | 1.0 |
| `honest_cv` (cross-validated) | 0.9 | `discovery` | 0.7 |
| `classifier_call` | 0.7 | `rna_derived` | 0.6 |
| `descriptive_aggregate` | 0.5 | `imputed_from_RNA` | 0.5 |

The **CV-gate** zeroes any witness whose permutation *P* ≥ 0.05 — a below-chance predictor cannot
influence the consensus however confidently it reports. The multiplicative form encodes the project's
central commitment: an imputed (0.5) classifier call (0.7) carries 0.35 of the weight of an observed
genotype before reliability is even considered. Imputed modalities can support a conclusion; they
cannot drive one.

**The genetic anchor is absolute.** If the genetic witness reports an *observed* driver, that lesion
sets the leading hypothesis and **cannot be outranked by any imputed prediction**. A disagreeing
prediction is recorded as an explicit conflict rather than averaged away. This is a hard constraint,
introduced after an imputed TP53 prediction outvoted an observed genotype in an earlier version.

**Therapies key to observed genetics only.** Recommendations come from the curated knowledge base
(`biomarker_drug.tsv`, `validation_rules.tsv`) keyed on *observed* drivers. A predicted-but-unconfirmed
lesion yields a **validation recommendation** — which confirmatory assay to order — never a direct
therapeutic suggestion. This is the mechanism preventing an expression-based inference from propagating
into a treatment decision.

---

## 3. The remaining agents

**Specimen / control gate (Layer A).** A RandomForest on 89-dimensional cell-state composition
classifies healthy versus diseased before any mutation calling. Two figures are quoted in different
places and both are real: **AUROC 0.965** on the full validated composition task, and **≈ 0.91
within-dataset**, the stricter estimate that controls for batch and is the one we use when arguing the
signal is biology rather than cohort. The operating point is deliberately conservative — disease
sensitivity ≥ 0.95 — so a "control" call is high-confidence while a false "diseased" call is harmless,
that specimen simply returning an empty mutation panel.

**Cell-state assignment (Layer A).** cellHarmony projection onto the `Hs-MarrowAtlas-L3M` reference
(CP10k + log1p, cosine matching) yields a per-cell state label and cosine score, aggregated to
per-(sample × cell-state) pseudobulks by summing raw counts.

**RAG knowledge graph (Layer D).** Supplies cited domain evidence to the supervisor. It is a curated
store, not free-web retrieval, which is what makes the therapy mapping auditable.

**Supervisor / arbiter (Layer D).** A **deterministic pre-pass** computes the weighted consensus,
applies the genetic anchor and fixes the conclusions; the LLM then narrates within them. It may lower
confidence but never raise it above the deterministic ceiling, and cannot add, drop or substitute a
therapy. On any LLM error the deterministic pre-pass *is* the returned consensus. Guarded
continuous-feedback loops defer deterministically to the anchor.

**Report and GUI decision board (Layer D).** Per-patient JSON + Markdown carrying the anchored call,
confidence, conflicts and deliberation trace; the web board renders these and accepts new patient
uploads, which re-enter the pipeline at Layer A.

---

## Reproducibility

Fusion and weighting: `pipeline/train_predictor.py` (deployed; `AMLMM_AUGMENT=1` enables cross-assay
augmentation), `pipeline/production_fused_model.py` (validation harness).
Arbiter and weighting factors: `pipeline/amlmm/arbiter.py`, `pipeline/amlmm/agent.py`.
Witness construction: `pipeline/amlmm/panel.py`. Control gate: `pipeline/control_gate.py`.
Knowledge base: `pipeline/amlmm/knowledge/`.
Per-modality standalone AUROC and fusion weights for every driver: Supplementary Table 5.
