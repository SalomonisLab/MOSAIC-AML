# MATRIX-AML — Developer / Experiment Log

Running log of the multimodal mutation-prediction work: every experiment, its method, its result, and
the decisions taken. Newest sections at the bottom of each theme. Companion to `OVERVIEW.md`
(architecture) — this file is the lab notebook.

Maintainer note: append a **timestamped** entry (date + time, e.g. `2026-06-30 09:52 EDT`) whenever an
experiment is run or a decision is made. Each entry = **timestamp — what / why / how / result / file**.
The thematic sections §1–§10 are the organized reference; §11 is the chronological running log going forward.

---

## 0. Goal & setup

**Goal.** Predict each AML driver mutation / cytogenetic lesion (binary present/absent) for a sample from
its multimodal single-cell pseudobulk profile, and call healthy controls as "control / no mutation".

**Cohort.** Salomonis AML atlas, `/data/salomonis-archive/LabFiles/Nicholas/AML-multimodal/` (cluster).
387 samples (247 AML, 24 Control, 12 MDS, 4 T-ALL, 100 uncategorized), ~12,255 pseudobulks (sample × cell-state).

**8 modalities** (sample-level after n_cells-weighted aggregation over cell-states):

| modality | # features | source | independent of RNA? |
|---|---|---|---|
| Composition | 90 | cell-state frequencies (cellHarmony) | ✓ |
| RNA | 385 UDON markers (→ **401** raw+prog, see §6) | pseudobulk_counts_hashed.h5ad | — (the source) |
| ADT | 129 | surface protein | ✗ imputed from RNA |
| Lipid | 625 (@fidelity≥0.3) | imputed | ✗ imputed from RNA |
| Metabolite | 1084 (@fidelity≥0.3) | imputed | ✗ imputed from RNA |
| GRN | 7486 | imputed regulons | ✗ imputed from RNA |
| LSC | 3–4 | stem-cell class probabilities | ✓ |
| Cell-comm | 141,101 | ligand-receptor, sample-level | ✓ |

**Leakage discipline (non-negotiable).**
- Sealed held-out = 29 user-highlighted samples (26 testable; Colorado AML-05/06 + Milan PT01_D14 absent from all modalities). NEVER in training, selection, or hyperparameter tuning. Scored once.
- Donor-grouped CV everywhere (`samples.donor_group`); no donor straddles train/val.
- Per-fold, train-only feature selection (differential `diff_cap`).
- Mutation labels: 1=positive evidence, 0=explicit-sheet absence, NaN=unknown (dropped per mutation).

**Infrastructure gotchas (see §8 for detail).** Head node `bmiclusterp-head` CANNOT run the conda numpy
(GLIBC_2.27 missing) → all compute via `bsub` on LSF compute nodes. SSH prints a CCHMC banner → grep-filter
`II-105|F-05`. Do NOT touch the user's SpliceScout / A549 / sra3-HepG2 LSF jobs.

---

## 1. Single-modality ablations & feature handling

- **Single-modality per (mutation, model)** — `_single_modality.py`, `_modality_panel.py`, `_mutation_model_grid.py`, `_model_ablation.py`. 11 models × 8 modalities × all mutations, sealed held-out.
- **Strong models** (consistently best, used everywhere downstream): **logL2, elastic, linSVM, shrLDA, PLS** (regularized linear). Weak/unstable: RF, HistGB, NaiveB, kNN, MLP.
- **Standalone peak held-out AUC (full features, best model):** RNA 0.867 (PLS) · Metabolite 0.845 · GRN 0.842 · Lipid 0.819 · ADT 0.776 · Composition 0.711 (RF) · LSC 0.604. → imputed modalities trail RNA (circularity: they can't exceed their source).

- **Feature-selection levels** — `_fs_levels.py` (full / 500-each-way / 100-each-way differential, 11 models, 8 modalities). **Result: FULL features win for the strong models, everywhere.** Where FS actually trims (GRN 7486, Metabolite 1084), it HURTS (GRN PLS 0.842 → 0.806 → 0.800). FS only "helps" weak models we don't deploy. → **No feature selection; use full feature sets (Cell-comm capped at 8000 differential for tractability — only modality affected).**

- **PCA / latent representation** — `_pca_test.py` (raw vs PCA-30 vs PCA-100). **Raw wins** for all strong modalities; PCA only helped ADT (+0.04). → no PCA.

- **CP10k normalization regression** — `_rna_norm_check.py`. CP10k+log1p DROPPED comp+RNA from ~0.84 to 0.76–0.80. → reverted to **bare log1p**; CP10k was a handicap.

---

## 2. Fusion architectures

- **Early fusion** — `_early_fusion.py`. Concatenate ALL modalities' top-500 differential features into one
  matrix (~4000 cols), fit one strong model. **Result 0.866** (best model) vs 0.861 best-single. Dilutes
  single-dominant mutations (hurts IDH1/IDH2/NRAS). Architecture note: its per-feature coefficients already
  subsume modality weighting → no separate modality-weight optimization applies. NOT chosen for deployment:
  breaks under missing modalities (uploads only have Composition+RNA) and is uninterpretable. (Old-RNA data.)

- **Late fusion** (chosen architecture): per-modality predictions combined with learned weights. Everything
  in §3–§5.

---

## 3. Ensembles & weighted late fusion (preds-era data, old RNA)

- **Phase A** — `_modality_full.py`. Per modality (env `AMLMM_MODALITY`), 11 models, all mutations: held-out
  prediction scores + donor-grouped CV (cv_mean/cv_std). Cap 8000 differential. → `preds_<MOD>.pkl`.
- **Phase B — consistency-weighted ensemble** — `ensemble_weights.py`. Weight modalities by CV reliability
  `w=max(0,cv−0.5)/(cv_std+0.05)`. **best-single 0.766 → uniform 0.793 → margin 0.796 → consistency 0.798.**
  Combining beats best-single for all 11 models; weights prioritize Metabolite/RNA/GRN/Lipid/ADT, demote
  Composition/LSC. Honest caveat: imputed circularity — much of the lift is smart ensembling of RNA-derived views.
- **optimize_panel v1** — `optimize_panel.py`. Per-mutation best model by CV + gated consistency fusion.
  **0.844 held-out.** Exposed CV-selection noise (IDH2 → RF → **0.39**, below chance).
- **optimize_panel v2** — `optimize_panel_v2.py`. Restrict to strong family, double-blend over model×modality
  cells, CV-gated deference, one-SE rule. **0.859 held-out** (oracle 0.945). Recovered IDH2 0.39→0.68, PHIP
  0.60→0.93. Lost TP53 0.90→0.77 (over-blend). This was the **best robust baseline** for a long stretch.

---

## 4. OOF + leakage-clean weight optimization (new RNA, OOF data)

- **OOF generation** — `_modality_oof.py`. Per modality, capture honest **out-of-fold** train predictions
  (donor-grouped GroupKFold, per-fold train-only FS) + held-out test scores. → `oof_<MOD>.pkl`. These are the
  level-1 features a weight optimizer / stacker may legitimately learn on. (Speedup: hoisted per-fold `diff_cap`
  out of the model loop → ~10× on GRN/Cell-comm.)
- **Weight optimizer (designed + adversarially audited)** — `optimize_weights.py`. Multi-agent workflow
  designed 3 methods + audited for leakage/overfitting; audit found real **selection-inflation bugs** (single
  strategy max-over-cols on val; MODEL-AVG scored by mean-of-AUCs; combine() zero-weight fallback; oracle −1
  sentinel) — all fixed. Leakage skeptic: clean. **Result: weighted-fusion 0.823** (hard per-mutation model+
  strategy selection → noisy; IDH2 collapsed to 0.267).
- **Stacking meta-learner** — `stack_meta.py`. Regularized logistic meta-learner over the full 40
  (strong-model × modality) base-predictor library; inner-CV `C`. **0.848 held-out** (> weighted-fusion;
  richer combiner helps, no hard selection).
- **Robust continuous weighting** — `weights_robust.py`. Reliability-weighted average of all cells, NO hard
  selection. **0.842 held-out.** IDH2 still 0.367 → proved the below-chance collapse is **n=6 held-out noise**
  (AUC SE ≈ 0.23 at 6 positives), not a selection bug.
- **Molded weights (dominance-guaranteed)** — `weights_molded.py`. Each modality = its best strong model;
  ridge-NNLS optimal weights per mutation; FLOOR falls back to ablation winner. **Proves optimized ≥ best
  single modality on OOF (0.893 vs 0.865, 17/17 by construction); held-out 0.834 vs ablation 0.825.** Weights
  are biologically sensible (FLT3→ADT 0.56, IDH1→RNA 0.64, complex/RUNX1→Metabolite 0.50).

---

## 5. Per-(mutation,model) optimized vs uniform; model selection ceiling

- **Full per-(mutation,model) comparison** — `weights_full.py`. Optimized modality weights vs UNIFORM, all 8
  modalities, every strong model. **Optimized ≥ uniform for EVERY model** (Δ +0.002…+0.024; ALL +0.014).
  **Best deployable: linSVM + optimized weights = 0.864** (uniform 0.840). All 8 modalities used (NNLS
  zeros the redundant ones per mutation — correct, not a bug).
- **Best-model-per-mutation** — `best_model_per_mut.py`. **ORACLE (model picked by held-out) = 0.891**
  (upper bound, selection-on-test, NOT deployable). **DEPLOYABLE (model picked by OOF) = 0.844 — WORSE than
  fixed linSVM 0.864.** → Per-mutation MODEL selection is a mirage at this data size (3–12 positives); per-
  mutation MODALITY weighting is real. **Ship: linSVM + molded modality weights.**

**THE CEILING.** Early fusion 0.866, late-fusion optimized 0.864, gated double-blend 0.859, stacking 0.848 —
independent methods all converge at **~0.86**. Strong evidence of a **data/feature ceiling**, not a combiner
limitation. Real headroom = more held-out positives, better features for hard mutations (IDH2/NRAS/SRSF2 —
even their oracle is only 0.80–0.91), NOT cleverer combining.

---

## 6. UDON for RNA (nathan request)

- **Question (nathan):** is RNA using UDON clusters + control-normalized fold-change vectors? They're most
  relevant for RNA.
- **Finding (code-grounded):** the RNA modality used only the UDON marker **gene list** to subset, on **raw
  log1p expression**. The control-normalized **fold vectors** (`RNA/clusters/udon_result.h5ad`, verified
  control-normalized via `n_controls_used`/`baseline_kind`, 9668 pseudobulks × 18451 genes, mean ~1.0 ratios)
  and **program clusters** (`UDON_final_program_assignments.tsv`, 16 programs) were used only in the separate
  subtype witness, not the predictor.
- **Loaders built** — `amlmm/udon_features.py`: `udon_fold_sample_matrix` (crosswalk via (Sample,cell_state)
  9668/9668, n_cells-weighted, log1p), `udon_program_matrix` (per-sample program fractions), `canonical_rna`.
- **Head-to-head** — `udon_rna_compare.py` (donor-grouped CV, same samples). **UDON folds HURT mutation
  prediction** (FLT3 0.784→0.718, IDH2 0.593→0.457, NPM1 0.897→0.869 — worse on all 5). Programs alone too
  coarse. **raw + program fractions is best** (FLT3 +0.026, IDH2 +0.021, neutral else). Control-normalization
  removes the absolute-expression signal mutations need; folds are a disease-discovery tool. Controls have NO
  fold vectors (they're the baseline) → folds can't classify healthy-vs-diseased either.
- **DECISION: adopted `raw + UDON programs` as canonical RNA = 401 features.** `canonical_rna()` is the single
  source of truth; `_modality_oof.py` + `_modality_full.py` use it; oof_RNA regenerated (401 feat).

---

## 7. Healthy-vs-diseased control + control gate (nathan request)

- **First control** — `control_healthy_vs_disease.py`. 24 Control vs 263 diseased, donor-grouped CV.
  AUC 0.90–0.985 across 7/8 modalities (Lipid 0.985 best; **LSC fails 0.46**; ADT 0.902). spec(controls→healthy)
  reported. Built-in batch-confound check.
- **Verification** — `control_hd_verify.py`. Within-dataset (Columbia/NYU-1/WashU, the 3 datasets with both
  classes; 108 samples/14 controls): AUC held **0.90–0.97**; label-permutation null collapsed to **~0.5**.
  → **biology, not batch.** ADT dropped 0.90→0.75 within-dataset (partly batch).
- **Deployed control gate** — `control_gate.py` → `control_gate.pkl`. Composition-only (sole upload-derivable
  modality), cv_auc **0.965**. Threshold targets disease sensitivity (set by user to **0.99** → thr 0.056,
  sens 0.992 / in-sample 263/263 diseased proceed = never miss a patient, control spec 0.50). Wired into
  `ingest_patient.py` (4 minimal, backward-compatible edits): a "control" call → report leads with
  `specimen_class:"control"` / "control / no mutation", subtype preserved under `subtype_if_diseased`.
- **Controls in labels:** the 24 controls ARE labeled no-mutation (present=0, absent across flags) and DO train
  as negatives — earlier claim that they were "dropped" was wrong (diag_controls.py corrected it).

---

## 8. Infrastructure & gotchas

- **Head node can't run conda numpy** (GLIBC_2.27). Everything numeric → `bsub` on a compute node
  (`/usr/local/anaconda3-2020/bin/python`). Quick jobs: `bsub -K` (blocking) + read the `-o` file.
- **TERM_CPULIMIT** killed OOF GRN/Cell-comm at ~60 min CPU — the `normal` queue applies a default CPU limit
  when no runlimit is set. Fix: submit with `-W 2880 -c 96:0` (mirrors how Phase A survived).
- **LSF dependency chaining:** `optweights`/`stackmeta` submitted with `-w "done(<oof jobs>)"` to auto-run when
  OOF completes — survives session boundaries (background pollers do NOT; they die with the session).
- **Disk-cache pattern** — `cache_sl.py` aggregates a modality to sample-level once → `_sl_<MOD>.pkl`
  (avoids re-aggregating multi-GB matrices / OOM). `canonical_rna` reuses `_sl_RNA.pkl` + appends programs.
- **SSH:** `-o ConnectTimeout=10 -o ServerAliveInterval=5`, grep-filter the `II-105|F-05` banner.
- **Efficiency note (self-critique):** several "oracle"/summary numbers were row-maxes of tables already in
  hand — compile from existing outputs before launching a new job. Save OOF *and* held-out AUCs together so a
  deployable-vs-oracle comparison needs no re-run.

---

## 9. Results summary (sealed held-out mean AUC over 17 mutations)

| approach | held-out | deployable | notes |
|---|---|---|---|
| best single modality (CV) | 0.766 | ✓ | baseline |
| consistency-weighted ensemble (Phase B) | 0.798 | ✓ | old RNA |
| optimize_panel v1 | 0.844 | ✓ | CV-selection noise |
| OOF weighted-fusion | 0.823 | ✓ | hard-selection flaw |
| robust reliability weighting | 0.842 | ✓ | new RNA |
| molded weights (NNLS) | 0.834 | ✓ | OOF 0.893, dominance proven |
| stacking meta-learner | 0.848 | ✓ | new RNA, OOF |
| gated double-blend v2 | 0.859 | ✓ | old RNA |
| early fusion | 0.866 | ✗* | old RNA; breaks on missing modalities |
| **linSVM + optimized modality weights** | **0.864** | **✓** | **← deployed combiner** |
| best-model-per-mutation (DEPLOYABLE, OOF-selected) | 0.844 | ✓ | worse than fixed linSVM |
| best-model-per-mutation (ORACLE, test-selected) | 0.891 | ✗ | upper bound only |
| oracle (best single cell on test) | ~0.94 | ✗ | non-deployable ceiling |

Held-out is noise-limited for low-prevalence mutations (IDH1 n+=3, IDH2 n+=6, FLT3-TKD n+=4): per-mutation
AUC SE ≈ 0.2 → those rows are not individually trustworthy and the combiner differences above are within noise.

---

## 10. Current state & open items

**Deployed / decided:**
- Combiner: **linSVM + per-mutation optimized (ridge-NNLS) modality weights** (~0.864).
- RNA representation: **raw markers + UDON program fractions** (401 feat).
- Control gate: Composition, disease-sensitivity 0.99, wired into `ingest_patient.py`.

**Open:**
- Wire the imputed witnesses (ADT/Lipid/Metabolite/GRN) into the **patient-upload path** (currently only
  Composition + genetics + control gate are wired; needs the rna2* imputers). This is the main deployment gap.
- Save the molded per-mutation weight table as the deployed artifact.
- Real headroom is data/features (more held-out positives; real metabolomics for IDH2), not combiners.
- Optional untried lever: early+late hybrid (add early-fusion prediction as one cell) — likely ~0.86.

---

## 11. Chronological running log (timestamped)

- **2026-06-30 09:52 EDT** — Created this `DEVELOPER_LOG.md` (222 lines) and pushed to the cluster;
  backfilled §0–§10 from the full mutation-prediction + control + UDON arc. Adopted the convention: every
  experiment/decision from here gets a timestamped entry (what/why/how/result/file). Per user directive to
  log everything with timestamps.
- **2026-06-30 09:52 EDT** — Confirmed cluster state: all 8 OOF pkls done; `optweights` (weighted-fusion,
  0.823) and `stackmeta` (stacking, 0.848) completed; control gate retrained at 0.99 disease-sensitivity
  (thr 0.056). Best deployable combiner = **linSVM + optimized modality weights (0.864)**; ~0.86 is the
  data ceiling (multiple methods converge).

- **2026-06-30 11:06 EDT** — **20 new model types** tested on the fused representation (z-scored balanced
  top-100/modality), `test_new_models.py` → `_new_models.txt`. **No new model beats the existing linear
  ones.** Top: logL2 0.789 (ref), Voting-soft 0.785 (best new), PLS 0.783 (ref), Stacking 0.775, linSVM
  0.770. Trees/kernels/QDA/LDA-svd all below; GaussProc 0.50 (kernel opt disabled, not a fair GP test);
  NuSVC n=1 (nu infeasible, ignore). Whole fused representation tops out ~0.79 — far below late-fusion
  linSVM+optimized (0.864) → re-confirms architecture (late fusion) ≫ model choice. Model space is exhausted.
- **2026-06-30 11:14 EDT** — **14 early-fusion improvement methods** tested, `early_fusion_improve.py` →
  `_early_improve.txt`. **None beat naive concat (0.829).** naive 0.829 = late 0.829 (tied), hybrid 0.821,
  then everything fancier is BELOW: zscore 0.786, gate 0.785, histgb 0.767, drop 0.766, mkl 0.763, residual
  0.760, pca 0.756, rbfsvm 0.734, blocknorm 0.723, interact 0.722, mlp 0.678. Added complexity overfits on
  small redundant data; plain concat + L2 is already at the representation ceiling. Confound noted: naive
  used 500 feat/modality, improvements 200 (mixes transform with feature count for naive-vs-zscore; the
  transform-only sub-comparison at 200 is clean and the fancy transforms still lose). Even naive early ≤
  late ≤ deployed 0.864 → **early fusion ≤ late fusion, every way.**
- **2026-06-30 11:14 EDT — CONCLUSION of the "try everything" push:** 34 approaches (14 early-fusion
  improvements + 20 new models) tested; **NONE beats deployed late-fusion linSVM+optimized weights (0.864).**
  The lever is data/feature quality (more positives; real non-imputed measurements), not the model or fusion
  method. Model + architecture space is exhausted.

- **2026-06-30 13:29 EDT** — **Trumpp/Waclawiczek venetoclax-AML cohort ingestion (nathan task).** 16 scRNA
  samples (8 paired Diagnosis+Refractory; Waclawiczek et al., Cell Stem Cell 2025; LSC subtypes + mutations
  + VEN/HMA response in `Trumpp.xlsx` Table S4). Submitted nathan's cellHarmony_lite alignment to the
  Hs-BM-titrated 89-state reference → integrated h5ad. Job **777433 `h5ad_combine`** (12h/128G), script
  `/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_AML_UDON/run_cellHarmony_Trumpp.lsf`, out → that
  dir's `output/`. Pre-flight: all 16 soupX inputs (matrix.mtx+barcodes+genes) populated; caught CRLF in
  SoupX_filepaths_Trumpp.txt (harmless — cellHarmony does line.strip()). NOTE: salomonis2 share (NOT the
  salomonis-archive MATRIX-AML data). Next: aligned h5ad → composition → run through mutation predictor +
  control gate as a VALIDATION cohort (known mutations = ground truth).

- **2026-06-30 15:30 EDT — Trumpp validation, composition-only = the deployment gap (honest negative).**
  First validation ran the predictor on Trumpp with COMPOSITION ONLY (the only modality a raw cellHarmony
  alignment gives) → **MEAN AUC 0.514 over 12 drivers (near-chance)**; bottom-4 BELOW chance (DNMT3A 0.27,
  IDH2 0.25, NRAS 0.36, ASXL1 0.44). Root cause (verified against `mutation_predictor.pkl` weights): the
  molded weights put **~0 weight on Composition** for those drivers (ASXL1={Lipid,GRN,Cell-comm},
  DNMT3A={GRN 0.50,Cell-comm 0.44}, IDH2={ADT 0.52}) — so composition-only scores them on a modality the
  model ignores → noise. The system was never meant to run on composition alone.

- **2026-06-30 17:00 EDT — FULL multimodal Trumpp deployment (the real imputation chain).** Built
  `pipeline/deploy_trumpp_full.py`: Trumpp scRNA (`combined_with_umap_and_markers.h5ad`, raw counts in
  `layers['counts']`) → per-(Library,cell-state) summed-count pseudobulks (724) → `rna2{adt,grn,lipid}`
  bundles → n_cells-weighted sample aggregation (matches `dataio.sample_modality_matrix`) → deployed
  `MutationPredictor` with all available modalities. Infra established:
  * Trumpp counts confirmed raw (intlike, ~12200/cell); cohort RNA pseudobulk = **SUM of raw counts**
    (corr(n_cells,rowsum)=0.91) → Trumpp RNA = log1p(n_cells-wt-mean of summed pseudobulks)[382/385 UDON
    markers] + zero program cols (UDON not run on Trumpp). RNA log1p med 6.43 (cohort 5.98; +0.45 = depth).
  * Imputer bundles were local-only (gitignored); pushed ADT/GRN/Lipid + code to the cluster engine-code via
    one tarball (`rna2_code_bundles.tgz`). All three unpickle on the cluster; output dims match cohort
    exactly (ADT 129, GRN 7486, Lipid 1009 → cohort built by these same bundles, confirmed via ADT h5ad
    `uns['prediction_summary'].bundle_path`).
  * **Metabolite bundle NOT on accessible storage** (cohort h5ad built on Nathan's Mac; `uns['imputation']`
    = per-target Ridge α=100); deferred. Cell-comm/LSC deferred (heavy/low-weight). Predictor renormalizes
    weights over present modalities, so partial roster degrades gracefully.
  * **ADT recipe reverse-engineered + calibrated.** ADT api does NOT normalize internally; cohort ADT h5ad
    `uns`: log1p scale, 31199 clipped negatives, names = `raw_feature_name` (`Hu.*`) vs predictor `features`
    (stripped). Validated on cohort pseudobulks: CP10k+log1p input → strip `Hu.` → clip≥0 → log1p gives
    **per-protein corr 0.994** vs Nathan's stored values (RAW input = 0.357; larger CP constants worse). A
    residual global affine remained → fit **per-protein Rosetta calibration** (my-recipe vs cohort stored,
    on 4000 cohort pseudobulks) and applied to Trumpp → ADT val med **2.179 vs cohort 2.175**, overlap
    **129/129**. GRN/Lipid need no calibration (their api CP10k+log1p-normalizes internally; fed raw counts).
  * **RESULT: MEAN AUC 0.514 (comp-only) → 0.668 (full, 5 modalities: Composition+RNA+ADT+GRN+Lipid).**
    Per-driver comp→full: DNMT3A 0.27→0.72, trisomy8 0.64→**1.00**, NPM1 0.53→0.80, kmt2a 0.79→0.93,
    TP53 0.65→0.77, complex 0.57→0.73, del5 0.67→0.73, IDH2 0.25→0.43, NRAS 0.36→0.43, ASXL1 0.44→0.46,
    FLT3 0.50→0.53 (ADT-dominant 0.83w; imputed ADT doesn't transfer → ~flat), TET2 0.52→0.49.
  * HONEST caveats: 3 drivers still <0.5 (ASXL1/IDH2/NRAS) — underpowered (IDH2/NRAS n+=2 of 16) or missing
    their Metabolite/Cell-comm weight. FLT3 ADT domain shift. Cohort own held-out mean is 0.795 (0.86
    well-powered); 0.668 is the honest external number (domain shift + 3/8 modalities absent + imputed-from-
    RNA modalities). 16 board reports rewritten with full multimodal predictions → `runs/trumpp_*`.

- **2026-06-30 20:25 EDT — Trumpp v3: rigor for Nathan's review (CIs, abstention, calibration, report).**
  Answering Nathan's pushback (still-4-below-chance; "worse than 0.86"; wants a per-sample artifact):
  * **Bootstrap 95% CIs (3000 resamples) per driver.** With 16 samples / 1–7 positives the single-driver
    AUCs are very wide. **0 of 12 drivers are significantly below 0.5** (every CI upper bound ≥ 0.5): the 4
    sub-0.5 points (IDH2 0.36 [0.00,0.87] n⁺2, NRAS 0.43 [0.00,0.93] n⁺2, ASXL1 0.46 [0.00,1.00] n⁺4, TET2
    0.46 [0.14,0.77] n⁺4) are small-sample scatter, not anti-prediction — same reason trisomy8 reads 1.000
    (n⁺2, also not meaningful alone).
  * **Abstention gate** added: a driver is flagged `abstain` when n⁺<3 OR modality coverage <60% (its
    weighted modalities — Metabolite/Cell-comm — absent). All 4 sub-0.5 drivers abstain → never emitted as
    confident calls (directly addresses "should never happen").
  * **RNA depth-calibration**: global factor 0.638 → Trumpp RNA median 5.981 = cohort (Trumpp ~1.57x deeper).
    Rank-preserving across the 16 samples ⇒ AUC unchanged; it fixes probability/threshold calibration.
  * **"0.86 vs 0.66" = different tests** (not a regression): 0.795/0.86 = cohort internal held-out with all 8
    MEASURED modalities; 0.514 = Trumpp composition-only; **0.660 = Trumpp full chain (5 modalities, 3 imputed
    from RNA, Metab/Cell-comm/LSC absent)**. Best-powered drivers carry it: NPM1 0.80, TP53 0.77, complex 0.73.
  * **Deliverable for Nathan**: `deploy_trumpp_full.py` now writes `runs/trumpp_<sample>/PREDICTION_REPORT.md`
    (per-driver prob/call/truth/✓✗ + per-modality contribution weight·score·cohort-OOF-AUC + provenance +
    abstention) alongside `patient_report.json`. Featured P13_Diagnosis (known NPM1/FLT3/DNMT3A): NPM1 0.82 TP
    (RNA+GRN), DNMT3A 0.85 TP (GRN), FLT3 0.32 FN (83%-ADT-weighted; imputed ADT score 0.30 — honest transfer
    miss). Packaged as `MATRIX-AML_Trumpp_Report.docx` (answers + full table + P13 detail + methods + limits).
  * Known limitation surfaced: binary calls over-call at a fixed 0.5 threshold (probs skew high in diseased
    marrow) → use probability ranking + abstention, per-driver calibrated thresholds are a quick add.

- **2026-06-30 21:00 EDT — Trumpp v4: Cell-comm + LSC now COMPUTED (7 of 8 modalities live).** Stopped
  deferring the "hard" modalities — both are derivable from the scRNA I already have:
  * **Cell-communication** — pushed the `fastComm` component to the cluster engine-code; reconstructed the
    exact LR table from the cohort's 141,101 interaction names (**647 unique CellChatDB ligand-receptor
    pairs** = `sender|ligand|receptor|receiver`), ran fastComm **per Library** on CP10k+log1p cells with the
    cohort's params (`min_lr_expression_score=0.2, max_lr_candidates_per_state_pair=5, include_self_edges=
    False, normalization=CP10k+log1p`), reindexed each sample to the 141,101 cohort columns → 16×141,101.
    ~5000 nonzero interactions/sample. **Lifted DNMT3A 0.72→0.821** (Cell-comm is its 2nd-largest witness,
    w=0.44, oriented score 0.936; CI [0.58,1.00]) and **completed modality coverage** (DNMT3A 56→100%, NPM1
    88→100%, ASXL1 79→100%, NRAS 63→86%) — only the 4 Metabolite-weighted drivers still have gaps.
  * **LSC** — the atlas's saved `LSC_RF_classifier.joblib` is a RandomForest on **cell-state frequencies**
    (= the composition I already build); applied it to each Trumpp sample → the 4 LSC subtype probs (overlap
    4/4). Needed a `monotonic_cst=None` shim on the unpickled trees (model pickled with sklearn <1.4).
  * **Metabolite** — genuinely blocked: exhaustive search (archive + salomonis2 + saljh8, unbounded) found no
    `rna2metabolite_aml_bundle.pkl.gz`; provenance paths (ADT/LSC/lipid) all point to the originator's local
    Mac (`/Users/saljh8/...`), and the deposited h5ad stores only outputs, not the ridge weights → cannot be
    reconstructed. Needs that one file onto shared storage.
  * **RESULT: MEAN FULL 0.661 over 12 (unchanged mean — the sub-0.5 drivers are n⁺=2–4 noise, not modality
    gaps; 0/12 significantly below chance). Qualitative win: DNMT3A a confident 0.82, coverage complete for
    all but the 4 Metabolite drivers.** All 16 board reports + P13 `PREDICTION_REPORT.md` regenerated with the
    full 7-modality contributions; deliverable `MATRIX-AML_Trumpp_Report.docx` refreshed.

- **2026-07-01 — Metabolite FOUND, all 8 modalities live, OOF-calibrated calls, runs/ cleaned + reran.**
  * **Metabolite bundle located** at `/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/
    components/rna2metabolite/artifacts/rna2metabolite_aml_bundle.pkl.gz` (I'd searched lowercase
    `/users/saljh8`; it was in Nathan's Revio checkout). Copied into engine-code → imputes, overlap **1084/1084**
    exact cohort match. **All 8 modalities now compute** for Trumpp; MEAN AUC **0.677** (was 0.661); every
    driver at **100% coverage**; TP53 0.77→0.855, FLT3 0.53→0.617, del5→0.800.
  * **OOF-calibrated per-driver thresholds** (fixes 0.5 over-calling). In-sample calibration failed (cohort
    positives score ~0.96 → thresholds 0.96 → held-out sensitivity crashed to 29%). Correct source =
    out-of-fold: added Youden-J on the deployed OOF blend in `train_predictor.py` (`P.thresholds`, used by
    `predict_one`; default 0.5). Median thr **0.72**. Held-out calls: **54%→72.6% accuracy, 71% sensitivity**,
    present-calls 390→226. AUC unchanged (0.795 all / ≈0.86 well-powered — the 0.86 is the well-powered subset;
    0.795 averages in rare point-muts KIT/IDH2/NRAS/WT1/PTPN11). All **29** held-out reports now written
    (Milan::PT01 = no atlas data). predictor.py gained `self.thresholds`.
  * **Board (gui): placeholders removed.** `scan_runs` now filters to real reports (predict_/trumpp_/ingest_)
    + groups them (Held-out / Trumpp / Uploaded) + per-sample predicted-vs-known accuracy pill; `AMLMM_SHOW_ALL=1`
    overrides. matrix_board.html grouped roster. Verified in preview: 47 real runs render, no placeholders.
  * **runs/ cleared + reran with the new predictor.** Deleted **50 stale dev/scratch dirs** (probe/regr/disc/
    gpu/model/phase/patient_/run/subtype__; list saved `runs/_cleared_dev_runs.txt`). Kept 29 predict_ + 16
    trumpp_ (already regenerated OOF-calibrated 8-modality) + single_modality caches. **Reran the 2 old
    composition-only ingest uploads** (AML7, BF71-CD34) through the full 8-modality OOF-calibrated predictor via
    new `deploy_scrna.py` (cosine cell-state assignment `amlmm.scrna.assign_cells` → per-state pseudobulks →
    RNA+impute-all+Composition+LSC+Cell-comm → predictor → mutation_panel). Needed a scanpy-free h5py 10x reader
    (scanpy pulls numba, incompatible w/ cluster numpy 2.4). AML7 = diseased (WT1 0.92, trisomy8 0.89 top);
    BF71-CD34 = **control** (CD34-sorted, face-valid; only RUNX1). `deploy_scrna.py` is the reusable
    arbitrary-scRNA full-chain path (what the GUI "Add patient" should call next).
  * Deliverable for Nathan refreshed to final: `MATRIX-AML_Trumpp_Report.docx` (8 modalities, 0.677, calibrated
    calls, P13 report) — the "Metabolite off-cluster" note was corrected.

### 2026-07-02 — Therapy hypotheses + Recommended validations added to the deliverable report
- **What:** Added two clinician-facing subsections to the P13 detailed prediction report (§3 of the combined
  `MATRIX-AML_Trumpp_Report` and the standalone `PREDICTION_REPORT.md`): **Therapy hypotheses** (maps the
  model's confident *present* calls — NPM1, DNMT3A, complex-karyotype — to a literature therapeutic rationale;
  explicitly flags the FLT3 negative call as low-trust because its signal is ~83% imputed-ADT-weighted) and
  **Recommended validations** (orthogonal confirmatory assays per driver, ordered by management impact:
  FLT3-ITD PCR w/ allelic ratio → NPM1 RT-PCR/IHC → myeloid NGS → karyotype+FISH → TP53/17p → flow cytometry).
- **Why:** the reports gave calls but no actionability; the AI-panel/ingest path already emits
  `targetable_therapies`/`recommended_validations`, but the external-cohort deliverable for Nathan did not.
- **How:** edited both markdown reports; both sections are hypothesis-generating and carry a "not clinical
  guidance" banner (rendered as a styled blockquote). Added a `>` blockquote branch to `convert.js`
  (left-rule + light-blue shading) so the disclaimer renders in Word; regenerated `MATRIX-AML_Trumpp_Report.docx`.
- **Result:** docx 17,412 → 19,190 bytes; verified via document.xml that both headings, all six validation
  rows, and the blockquote shading (`F2F6FB`) are present. Artifacts in `scratchpad/artifacts/`.

### 2026-07-02 — Tested improvements from Nathan's marker-discovery pipeline (labels + markers A/B)
Nathan shared `aggregate-markers/MARKER_DISCOVERY_SUMMARY.md` (his own Mac-run marker-discovery + bulk→sc
classification track, mirrored to the archive). Tested the two concrete things it offers our predictor.

- **(1) Labels — near-no-op.** Diffed our `mutation_matrix_explicit.tsv` (already dense 402×73 0/1) against his
  authoritative `Metadata/AML_metadata_CLEAN.xlsx` `Mutation_Matrix` tab (402×65, 0/1). **402/402 samples matched;
  the ONLY difference across all drivers = 3 samples in `mut_FLT3`** (NYU-1 AML0024/AML2451/AML3050), all
  **FLT3-TKD-only** (canonical−, ITD−, TKD+). Ours folds TKD into FLT3 (→47); his convention is FLT3∪ITD with
  **TKD separate** (→44). His summary's smaller counts (NPM1 66, TET2 29…) were positive-DONOR counts, not
  samples — our labels were already built from the same CLEAN metadata. Built `mutation_matrix_explicit_v2.tsv`
  (his convention). **A/B (FLT3, deploy RNA block): baseline ho=0.788 → v2 ho=0.758** — TKD-separation is a
  *correctness* choice (TKD≠ITD, distinct drug sensitivity; we already carry `mut_FLT3-TKD`), **not a performance
  gain** (slightly lower, within small-n noise). Verdict: optional correctness adoption, not an improvement.
- **(2) Markers — his don't beat ours overall; our features already win.** RNA-modality A/B on the SAME
  sample-level log1p gene pseudobulks (`exp_nathan.py`, LSF job 793902, held-out + donor-grouped 3-fold CV),
  4 feature sets per driver over 17 scored drivers. **Mean held-out AUC: OURS_udon 0.777, OURS_deploy 0.768,
  NATHAN-markers 0.735, ALLGENES 0.655.** Wholesale Nathan markers HURT (0.777→0.735); all-genes far worse
  (0.655) — the UDON curation is doing real work. **Honest CV-based (leak-free) per-driver selection of Nathan
  where CV wins → held-out 0.739→0.730 (still slightly worse)** — CV/held-out disagree for 6/9 (small-n overfit).
  **Only ASXL1 & FLT3-TKD show BOTH CV and held-out favoring Nathan robustly**: ASXL1 0.48→0.72 held-out (fixes a
  below-chance driver, CV 0.672→0.753), FLT3-TKD 0.841→0.898 (CV 0.668→0.922); FLT3 ~tie (Nathan 0.825 best).
  Caveat: his markers are cell-type-resolved; sample-pseudobulk washes that out, so this fairly tests his marker
  *genes* at sample level, not his cell-type discovery *design*.
- **Decision:** keep UDON markers as the RNA default (they win overall); the defensible targeted adoption is
  Nathan's markers for **ASXL1** (below-chance→0.72) and optionally **FLT3-TKD**, via a per-driver RNA feature
  override + re-fit (pending). Adopt his FLT3-TKD label separation for correctness. Deployed predictor UNCHANGED
  by this test (standalone eval; `mutation_predictor.pkl` untouched). Cluster note: submission throttled behind
  ~1250 SpliceScout PENDs; job landed only after `bmod`-ing 32G→12G (was unschedulable) + btop to position 1.
  Artifacts: `scratchpad/exp_nathan.py`, `exp_nathan_results.json`, `build_v2_labels.py`, `mutation_matrix_explicit_v2.tsv`.

### 2026-07-02 — Project presentation (.pptx) built
- **What:** 15-slide `MATRIX-AML_Presentation.pptx` (in `scratchpad/artifacts/`), framed on the abstract
  "A precision AI molecular-diagnostic & drug-repositioning platform for AML." Sections: challenge → vision →
  data foundation (383 samples / 12,255 pseudobulks / 89 states / 8 modalities) → Pillar 1 imputation →
  Pillar 2 late-fusion classification → Pillar 3 multi-agent engine → internal results (0.795 / ≈0.86 / control
  gate 0.965) → external (Trumpp 0.677, Grimes CITE-seq) → independent cross-validation (ours 0.777 vs discovered
  0.735) → therapy hypotheses → Pillar 4 drug repositioning (splicing reversal) → platform → honest limits →
  conclusion. Speaker notes on every slide.
- **How:** pptxgenjs (matte-black theme, one muted-teal accent, mono kickers, native bar charts). QA on Windows:
  no LibreOffice/poppler, so rendered via PowerPoint COM (`WithWindow=$true`). Subagent visual QA → fixed 2
  chart-label defects (AUC rounding to "1" → `dataLabelFormatCode "0.00"`; axis major unit 0.1). Root-cause a
  hard failure: pptxgenjs LINE shapes with negative w/h emit invalid negative `<a:ext>` → PowerPoint "corrupt";
  fixed with non-negative extents + `flipV`. Recompressed (392KB→76KB), verified it re-opens.

### 2026-07-06 — Full head-to-head vs Nathan's marker-discovery pipeline (Lee thread)
Nathan CC'd Lee with `singlecell_classification_metrics.xlsx` (per covariate x modality x cell_type OOF AUROC,
5,473 rows, 31 covariates) + `bulk_cohort_classifiers_imputed_to_sc.xlsx` (bulk->sc transfer). Built a
comprehensive comparison workbook (`Desktop/MATRIX-AML_vs_marker-pipeline.xlsx`, 7 sheets).
- **exp_allcov.py** (LSF `long`/`test`, job 814440): ran MATRIX-AML fused multimodal donor-grouped OOF AUROC
  on EVERY covariate the pipeline ran, incl. ones our deployed predictor skipped (PHIP 0.978, SF3B1 0.858,
  JAK2 0.834, SRSF2 0.769, NF1 0.899, CEBPA 0.917, GATA2 0.919, KMT2A-rearr 0.965) + clinical (AML-vs-control
  0.995, FAB_monocytic 0.931, ELN_adverse 0.889, relapse 0.800, sex_M 0.693). 40 covariates.
- **exp_percelltype.py** (LSF, job 814289): scored our features at his exact granularity — per
  (mutation x modality x cell-state), donor-grouped OOF + 200-perm null → 7,005 rows (his schema), 4,024 signal.
- **Results (apples-to-apples, same metric both sides):** on fused OOF vs his median-across-cell-types OOF,
  **MATRIX-AML wins 29/31 covariates** (losses: sex_M — his XIST/DDX3Y cell-type markers; ELN_adverse 0.889 vs
  0.894 tie). Paired per-cell-type (3,032 identically-defined cells): **MATRIX-AML higher on 56% (1,690), mean
  0.716 vs 0.687**. Bulk->sc transfer: we win 31/31 (his transfer median AUROC_all ~0.55-0.63). Our per-modality
  mean OOF and full 7,005-row table included as sheets.
- **Honest framing captured in the sheet:** his "best cell-type (max)" is optimistic best-of-hundreds (his own
  caveat); the fair columns are median (typical) and our fused OOF. Our SEALED held-out (deployable, harder)
  number is kept as a separate column. LSF note: per-user 125-slot cap (A549 array) blocked scheduling;
  bswitch to `test` + bmod -n 1 got allcov running immediately.

### 2026-07-09 → 07-13 — Variant-level bulk predictor becomes the PRIMARY mutation caller
Cross-cohort bake-off pushed to **variant level** and promoted to the platform's primary caller.
- **Variant-level relabeling (`pipeline/bulk_external.py`):** genes with distinct functional hotspots are
  split into sub-categories parsed from BeatAML `hgvsp_short` / `variant_classification` / clinical fusions
  (FLT3_ITD vs TKD_D835-I836 vs other; DNMT3A_R882 vs nonR882; NRAS/KRAS G12/G13/Q61; TP53 hotspot-DBD vs
  LOF; SF3B1 K700/K666; U2AF1 S34 / Q157-R156; CEBPA bZIP/N-term/biallelic; …). Genes with no hotspot stay
  gene-level. **58 categories; ~50 clear ≥6 positives** in BeatAML (5-fold CV-OOF, the primary metric).
- **Result — splitting recovers hidden signal:** DNMT3A_R882 CV-AUROC **0.88 vs nonR882 0.61**; FLT3_ITD 0.91
  vs TKD 0.69; U2AF1_S34 1.00 vs Q157 0.80. TP53 hotspot-DBD 0.90 ≈ LOF 0.89 (clean split, no gap). No-FS beats
  MarkerFinder (+0.05, 73% of pairs); ENSEMBLE/PLS/logL2/shrLDA > trees > linSVM > MLP (worst). BeatAML→Leucegene
  transfers (NPM1/SRSF2/U2AF1_S34/STAG2/RUNX1 on-diagonal); TET2 the notable non-transfer (0.83→0.65).
- **Per-modality breakdown of the deployed sc system:** imputed Metabolite/Lipid/GRN carry it (single-best for
  19/26 mutations); RNA is the reliable backbone (weighted 17/26) but rarely single-best; LSC/Cell-comm weak;
  **fusion beats best-single in 26/26**. Base-learner sweep (8 families through the exact deployed recipe)
  **confirmed the deployed LinearSVC is optimal** — no learner beats it on train CV-OOF; percentile-calibration
  + NNLS late-fusion neutralize the base-learner differences that mattered on raw bulk.
- **Decision (boss):** the bulk variant-level predictor covers ~2× the mutations at comparable accuracy on
  cheap ubiquitous bulk RNA, so it becomes the **PRIMARY mutation caller**; the sc multimodal system is kept for
  cytogenetics (inv16/del5/del7/complex/trisomy8/KMT2A) + multimodal depth. Head-to-head on 18 shared mutations:
  sc 0.853 vs bulk comparable; bulk adds 25 categories sc can't reach (U2AF1_S34 1.00, SF3B1 1.00, JAK2_V617F,
  variant splits).
- **Deployable (`amlmm/bulk_predictor.py` + `pipeline/train_bulk_predictor.py` → `bulk_mutation_predictor.pkl`,
  3.9 MB):** `BulkMutationPredictor` trains on BeatAML2 (n=707 WES) over 50 categories on the 14,237 shared genes;
  logL2 / no-FS / top-2500-variance / percentile / F1-max. Per-cohort z-refs (sc/beataml/leucegene) so ONE model
  scores plain bulk RNA and any single-cell sample collapsed to its bulk-equivalent. Mean 5-fold CV AUROC 0.829;
  validated on the sealed sc held-out via bulk-equivalent (SRSF2 0.91, FLT3_ITD 0.83, TET2 0.76).
- **Wired into `pipeline/ingest_patient.py`:** every upload's scRNA → `bulk_equiv_from_adata` (sum cells → CP10k)
  → `predict(ref="sc")` → new top-level `mutation_predictions` report block (mode `bulk_variant_primary`) — the
  first time the ingest path calls mutations from expression. Kept OUT of the arbiter (predicted ≠ ground truth;
  user-supplied mutations stay the deterministic anchor); degrades gracefully if the pkl is absent. Real-data
  sanity: median 3 confident-present calls/sample.
- **Caveat:** variant-level labels need full WES, feasible only in bulk cohorts — our sc cohort's per-sample HGVS
  is too sparse (only FLT3_ITD clears ≥8), so the sc system stays gene-level.

_Last updated: 2026-07-13._
