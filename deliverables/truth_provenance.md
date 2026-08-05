# Ground-truth provenance — MOSAIC-AML mutation calling

Where every "known genotype" label comes from, its granularity, and its known gaps. This exists so a
reviewer can audit the labels rather than take them on faith. The single most important caveat is at the
bottom: an "absent" label sometimes means *"not assayed"*, not *"assayed and wild-type"* — those are
false negatives in the **labels**, not the model.

| Cohort | n | Genotype assay | Truth granularity | VAF? | Used as |
|---|---|---|---|---|---|
| **BeatAML2** | 707 | WES / targeted capture somatic calls (dbGAP phs001657, `mutations.txt`) | **Variant-level** (exact protein change, e.g. DNMT3A R882) | **Yes** (`t_vaf`) | Bulk model training + CV |
| **Leucegene** | 367 | dbGAP variant calls (RNA-seq + orthogonal validation) | Category / gene-level | No | Bulk model external test |
| **Single-cell atlas** | 387 (357 labeled) | Per-center: **karyotype** for cytogenetics; clinical/targeted panels or WES for SNVs | Mostly **gene-level** (0/1); 17 categories carry the exact variant | No | Multimodal training + CV; held-out |
| **Trumpp / Waclawiczek** | 16 | Published driver list (Table S4) | Gene-level, **positives only** (others inferred wild-type) | No | Single-cell external test |
| **MDS (Hs_MDS_UDON)** | 64 | Genotype embedded in sample `Group` string | Gene-level, 4 genes only (SRSF2 / U2AF1 / RUNX1 / TET2) | No | AML-precursor stress test |

## Granularity: variant-level vs gene-level

- The caller predicts **58 variant-level categories** (e.g. `DNMT3A_R882` vs `DNMT3A_nonR882`, `NRAS_G12` vs `G13` vs `Q61`).
- Only **BeatAML** has variant-level truth for all of them. **Leucegene / single-cell / Trumpp truth is mostly gene-level**, so single-cell performance is scored by aggregating a gene's variant categories to one gene-level call (present iff any variant category fires) — this is why the single-cell denominators are ~40 gene units, not 58 variant categories.
- Consequence: single-cell numbers should be read as **gene-level** ("did we detect a DNMT3A mutation") not variant-level ("did we detect R882 specifically").

## Per-cohort gaps (things a reviewer will ask about)

- **BeatAML** — capture panels differ (`capture_type`); a gene outside a sample's panel reads as wild-type. `t_vaf` is bulk clonality; a low-VAF subclonal driver may be under the model's expression-detection limit (see `vaf_stratified.*`).
- **Leucegene** — FLT3-ITD scored 0.00 sensitivity; the ITD calling/definition differs from BeatAML's, so this is very likely a **truth-definition mismatch**, not a pure model failure. Flagged, not hidden.
- **Single-cell atlas** — the highest-risk labels. Cytogenetics come from karyotype (reliable for the 6 recurrent events); SNVs come from whatever each of the 6 centers (CCHMC, Colorado, Columbia, Milan, NYU, WashU) sequenced. **A gene a center did not test is labeled "absent."** So single-cell specificity may be *understated* (some "false positives" are real mutations the center never looked for). A per-sample "genes assayed" panel from each center would close this — currently not available uniformly.
- **Trumpp** — only Table S4 drivers are positives; every other gene is assumed wild-type. Small n (16).
- **MDS** — only 4 genes have genotype; 28 samples with empty `Group` are **excluded as unlabeled** (not used as negatives), which is the honest choice.

## What would make the labels bulletproof

1. A per-sample **"genes assayed" matrix** for the single-cell cohorts, so "absent" can be split into *tested-wild-type* vs *not-tested* (removes the false-negative-label risk).
2. Variant-level truth (not just gene-level) for at least a subset of the single-cell samples, to validate the variant-resolution claims.
3. Orthogonal re-sequencing of a handful of model-positive / label-absent single-cell samples — if the mutation is really there, it converts a counted "false positive" into a corrected label **and** a novel true positive.
