# AML multimodal — local working copy

A curated local copy of the Salomonis-lab AML multimodal atlas, pulled for building
and testing an agent / modeling framework on top of it.

## The three names

| name | what it is |
|---|---|
| **MOSAIC-AML** | the platform — *Multimodal Omics and State-Aware Inference of Cancer Drivers in Acute Myeloid Leukemia*. The atlas, the ingest path, the witness panel, the arbiter, the decision board. |
| **CIPHER-AML** | the **mutation predictor** — *Cell-state Inference of Pathogenic Hits from Expression and Regulation*. Infers driver lesions from expression, cell state and regulon activity, without sequencing. |
| **COMPASS-AML** | the **drug predictor** — *Cell-state Oriented Modelling of Pharmacologic Assay Sensitivity*. Predicts ex-vivo inhibitor sensitivity from the BeatAML2 functional screen, resolved by cell state. |

CIPHER-AML and COMPASS-AML are deliberately **parallel** layers, not a chain: a mutation call enters the
drug layer as one piece of evidence among several, never as a look-up key for a therapy.


## Provenance
- **Data + scripts** copied from the cluster: `bmiclusterp-head:/data/salomonis-archive/LabFiles/Nicholas/AML-multimodal`
  (that folder is a *deposit* — outputs + deploy/eval scripts + reports, not the live engine).
- **Engine code** = `engine-code/altanalyze3/`, a fresh clone of
  https://github.com/SalomonisLab/altanalyze3 (current; UDON lives in `altanalyze3/components/udon/`).
- The deposit scripts have hardcoded macOS paths (`/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/...`)
  — the original working tree is on Nathan's machine, not here.

## Layout
```
AML-multimodal/
├─ engine-code/altanalyze3/   UDON normalize+cluster engine (components/udon/, ~37 modules)
├─ data/                      curated multimodal pseudobulk matrices (.h5ad), ~1.5 GB
│   ├─ RNA/pseudobulk_counts_hashed.h5ad            12,255 × 35,702   MEASURED
│   ├─ GRN/imputed_grn_all_pseudobulks.h5ad         12,255 × 7,486    imputed from RNA
│   ├─ Metabolite/pseudobulk_imputed_metabolite_aml.h5ad  12,255 × 2,486   imputed from RNA
│   ├─ Lipid/pseudobulk_imputed_lipid_aml.h5ad      12,255 × 1,009    imputed from RNA
│   ├─ ADT/pseudobulk_adt_imputed.h5ad              12,255 × 129      imputed from RNA
│   └─ cell-communication/combined_sample_by_interaction.h5ad   383 × 141,101  (fastComm)
├─ labels/                    everything small you cluster / predict on
│   ├─ AML_clinical_metadata_CLEAN.tsv              405 rows (SPARSE — see caveats)
│   ├─ sample_cellstate_counts.tsv                  397 × 90 cell-type composition (= cell frequency)
│   ├─ RNA_UDON_final_program_assignments.tsv       UDON program per pseudobulk
│   ├─ {ADT,GRN,Lipid,Metabolite}_udon_clusters.txt UDON clusters per modality
│   ├─ LSC_prediction_subtype_RF.tsv                deployed LSC subtype calls
│   ├─ cellstate_signatures.json                   10 lineage signatures learned from the atlas
│                                                  (used by the COMPASS-AML `state` feature block)
│   └─ cellcomm_sample_manifest.tsv
└─ scripts/                   deposit deploy/eval code + imputation reports
    ├─ LSC-prediction/algorithm/   honest_model_bakeoff.py (CV template), RF model .joblib
    ├─ cell-communication/algorithm/
    └─ {GRN,Lipid,Metabolite}/algorithm/   rna2*_report.docx (imputation validation)
```

## Key facts
- 383 specimens (45 healthy-marrow controls + 338 disease), 12,255 cell-state pseudobulks
  (median 30 states/sample, range 2–88). The 5 modality matrices are **row-aligned** on the
  same 12,255-pseudobulk index with identical `.obs`.
- The per-pseudobulk genetic label is `adata.obs['Annotation']` (28 classes: NPM1c, Inv16,
  TP53, FLT3-ITD, SF3B1, …) — richer than the clinical TSV, but only ~120–130 samples carry a
  *specific* driver; the rest are generic "AML"/"Pediatric-AML"/"Control".

## Caveats (carried over from the assessment — read before modeling)
1. **Imputed ≠ independent.** GRN/Metabolite/Lipid/ADT are deterministic regressions on RNA,
   trained on external bulk CPTAC (n≈85) then applied to scRNA pseudobulks. Honest held-out
   fidelity is modest (lipid median Spearman 0.36; metabolite 0.27, R² negative). Each `.h5ad`
   `var` carries `heldout_spearman` — filter on it. Treat them as interpretability layers, not
   extra evidence beyond RNA.
2. **Labels are the bottleneck.** In `AML_clinical_metadata_CLEAN.tsv`: vital_status 0/405,
   survival 20, WHO 42, FAB 49, ELN 75. Survival / drug-response are not trainable in-house yet.
3. **Three nested leakage traps:** 12,255 pseudobulks from 383 donors → group CV by donor;
   10 datasets → leave-one-cohort-out; don't CV imputed features alongside parent RNA.
   `scripts/LSC-prediction/algorithm/honest_model_bakeoff.py` is the reference harness
   (it caught a 0.85→0.59 selection-bias inflation; RF honest balanced acc ~0.59, 148/397 low-conf).

## Deliberately NOT downloaded (get from the cluster / Nathan if needed)
- `RNA/pseudobulk_scaled_log2_hashed.h5ad` (~1.0 GB) — redundant with counts; normalize yourself.
- Per-modality `clusters/udon_result.h5ad` intermediates (~1.1 GB total) — final labels are in
  the small `labels/*udon*` files instead.
- `cell-frequency/combined_annotations.txt` (224 MB raw per-cell) — already summarized in
  `labels/sample_cellstate_counts.tsv`.
- `scALABLE` imputation training code (rna2grn/lipid/metabolite/adt) — not on the cluster; on
  Nathan's machine. `cellHarmony-web` (fastComm engine) — likewise not on the cluster.

## Loading
Needs Python with `anndata`/`scanpy` (e.g. an AltAnalyze/scanpy conda env). Large matrices:
`ad.read_h5ad(path, backed='r')` to avoid loading fully into RAM.

## COMPASS-AML (the drug-response layer)

`pipeline/amlmm/drug/` predicts ex-vivo inhibitor sensitivity from the transcriptome, trained on the
BeatAML2 functional screen (`data/external/beataml/beataml_probit_curve_fits_v4_dbgap.txt`, fetched
from the public BeatAML2 repo). Three separate models — patient-level response (A), the same model
re-applied per cell state (B), and independent mechanistic target evidence (C) — combined into a
per-tier prioritisation. See `deliverables/METHODS_COMPASS-AML.md` for methods and validation, and
`gui/therapy.html` / `gui/rx_validation.html` for the two front ends.

**Ex-vivo sensitivity is a prioritisation signal for trial matching or laboratory validation, not an
estimate of clinical benefit.**
