# Risk prediction against ELN 2022 — inference, VAF sensitivity, and where the guideline fails

Prompted by Nathan's note: benchmark risk prediction against **ELN 2022**, which has to be *inferred*
in BeatAML, and test **10% vs 40% VAF** thresholds.

    python pipeline/eln2022.py --vaf 0.10   (and --vaf 0.40)
    python pipeline/exp_eln2022_benchmark.py

Papers used: Döhner *et al.* Blood 2022 (Table 6); Röllig/Bill Leukemia 2025 (MDS-related mutations);
Pollyea/Döhner Blood 2024 (ELN under venetoclax-azacitidine); Löwenberg/Döhner Blood 2024
(less-intensive therapies); Sekeres *et al.* ASH 2025 guidelines for older adults.

---

## 1. ELN 2022 had to be inferred — BeatAML ships only ELN 2017

Every benchmark in this platform had been scored against the shipped `ELN2017` column. The two standards
differ in exactly the patients we get wrong:

| | ELN 2017 | ELN 2022 |
|---|---|---|
| FLT3-ITD allelic ratio | used (AR ≥ 0.5 changed the category) | **dropped** |
| CEBPA | biallelic required | **in-frame bZIP, mono- or biallelic** |
| MDS-related genes | ASXL1, RUNX1 only | **9 genes**: ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2 |
| TP53 | not a criterion | **adverse at VAF ≥ 10%**, any allelic status |
| t(9;11) | adverse-leaning | **intermediate**, and takes precedence over concurrent adverse mutations |

`pipeline/eln2022.py` implements Table 6 in full, including the precedence footnotes (CBF unaffected by
KIT/FLT3; NPM1 + adverse cytogenetics → adverse; MR genes *not* counted adverse alongside a
favorable-risk subtype; KMT2A-PTD excluded from KMT2A-rearranged; hyperdiploidy excluded from complex
karyotype). Karyotypes are parsed from the free-text ISCN strings; fusions from `consensusAMLFusions`;
CEBPA bZIP from variant-level protein positions in `mutations.txt`.

**Agreement with the shipped ELN 2017 label: 0.776 (n = 548).** Every disagreement runs in the direction
the guideline change predicts:

| 2017 → 2022 | n | driven by |
|---|---|---|
| Adverse → Intermediate | 48 | FLT3-ITD allelic ratio dropped |
| Favorable → Intermediate | 37 | NPM1mut **with** FLT3-ITD no longer favorable |
| Intermediate → Adverse | 25 | the seven newly-added MDS-related genes |
| Favorable → Adverse | 11 | NPM1mut with adverse cytogenetics |
| → Favorable | 2 | monoallelic bZIP CEBPA newly favorable |

---

## 2. VAF threshold: 10% vs 40% reclassifies 4.7% of patients and changes nothing prognostically

ELN 2022 names 10% for TP53 and is silent for every other gene; 40% is the conventional proxy for a
clonal/biallelic event.

| | TP53-positive specimens | ELN 2022 distribution |
|---|---|---|
| VAF ≥ 10% | 80 | 193 Fav / 371 Int / 378 Adv |
| VAF ≥ 40% | 70 | 188 Fav / 401 Int / 353 Adv |

**30 of 638 specimens with mutation data (4.7%) change category**, all toward Intermediate — 25
Adverse → Intermediate (TP53 and MR-gene calls falling below 40%) and 5 Favorable → Intermediate
(CEBPA bZIP calls falling below 40%).

**But the discrimination is unmoved:**

| label | C-index (444 survival patients) |
|---|---|
| ELN 2017 (shipped) | 0.610 |
| ELN 2022 @ VAF ≥ 10% | **0.612** |
| ELN 2022 @ VAF ≥ 40% | 0.610 |

So the VAF choice is a real reclassification and a prognostic non-event at cohort level. It should be
stated as a fixed 10% (matching the guideline) rather than tuned, and the 4.7% who move should be
reported as an uncertainty band rather than treated as a modelling decision. One genuine improvement in
2022: the Intermediate group is cleaner — median survival 23.2 months versus 10.3 under 2017 — even
though the overall C-index does not move.

---

## 3. Where ELN fails, and what works instead

### 3.1 ELN does not stratify non-intensively-treated patients — reproduced

Pollyea/Döhner (Blood 2024) reported that ELN classifiers "did not provide clinically meaningful risk
stratification" under venetoclax-azacitidine, and that **TP53 / FLT3-ITD / NRAS / KRAS** status did
(26.5 / 12.1 / 5.5 months). That reproduces here:

| stratum | ELN 2017 | ELN 2022 | Pollyea 4-gene |
|---|---|---|---|
| intensive induction (n = 357, 177 deaths) | 0.631 | **0.642** | 0.574 |
| non-intensive / unknown (n = 87, 68 deaths) | 0.496 | **0.462** | **0.612** |

**ELN is at or below chance outside intensive induction**, and the published 4-gene rule works there
and only there — exactly its stated scope. Group medians in the non-intensive stratum: 5.9 / 0.9 / 0.5
months (lower-risk / ITD-RAS / TP53); in the intensive stratum 40.6 / 23.2 / 7.0.

*Caveat:* BeatAML records induction type, not venetoclax-azacitidine specifically, and this stratum
includes unknown induction. The very short medians suggest it contains patients who died before or
early into treatment. This is a proxy for Pollyea's population, not a replication of it.

### 3.2 Splitting ELN 2022 adverse by MDS-related genes — does NOT reproduce

Röllig/Bill (Leukemia 2025) reported that adverse patients whose adverse feature is an MR-gene mutation
outlive other adverse patients (14.7 vs 8.3 months, p < 0.001). Here:

| group | n | deaths | median survival |
|---|---|---|---|
| adverse, MR gene mutated | 126 | 86 | 9.5 months |
| adverse, no MR gene | 63 | 47 | 7.8 months |

Same direction, far smaller, **not significant (log-rank p = 0.38)**. On this cohort the ELN 2022
adverse category is not usefully heterogeneous along that axis, and splitting it is not free
discrimination. Reported as a negative.

---

## 4. What this changes for MOSAIC

### 4.1 Against the corrected bar

| | C-index (n = 444) |
|---|---|
| ELN 2022 @ 10% VAF | 0.612 |
| **MOSAIC deployed** | **0.771** |
| both together | 0.771 |

The model contains whatever ELN 2022 carries — adding the guideline label on top gains nothing.

### 4.2 A correction to our own published limitation

`METHODS_survival.md` has stated that the model is "close to useless in non-intensively-treated
patients (C-index 0.554)". That number came from a model trained on the **pooled** cohort and evaluated
inside the subgroup. Fitting **within** the stratum changes it:

| stratum | MOSAIC (stratum-fitted) | Pollyea 4-gene | MOSAIC + 4-gene |
|---|---|---|---|
| intensive (n = 357) | 0.729 | 0.574 | 0.731 |
| non-intensive (n = 87) | **0.681** | 0.612 | **0.701** |

So the failure was largely an artefact of pooled fitting, not an intrinsic limit: **0.554 → 0.681**
when the model is fitted within the stratum, and **0.701** when the published 4-gene rule is added.
The 4-gene rule adds nothing in intensively-treated patients (+0.003), which is the correct behaviour
for a classifier built for less-intensive therapy.

---

## 5. What to do

1. **Report ELN 2022 at a fixed 10% VAF** as the benchmark, with the 4.7% reclassification band noted.
   Do not tune the threshold — it does not change discrimination.
2. **Stratify by treatment intensity before fitting.** It is worth more than any relabeling: +0.127
   C-index in the group where the guideline is at chance.
3. **Carry the 4-gene rule as a component in the non-intensive arm only**, where it adds +0.020 over
   the stratified model and where ELN contributes nothing.
4. Do not split the adverse category by MR-gene status on this cohort — measured and null.
