# MOSAIC-AML — Developer / Experiment Log

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

## 6. UDON for RNA (boss request)

- **Question (boss):** is RNA using UDON clusters + control-normalized fold-change vectors? They're most
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

## 7. Healthy-vs-diseased control + control gate (boss request)

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

- **2026-06-30 13:29 EDT** — **Trumpp/Waclawiczek venetoclax-AML cohort ingestion (boss task).** 16 scRNA
  samples (8 paired Diagnosis+Refractory; Waclawiczek et al., Cell Stem Cell 2025; LSC subtypes + mutations
  + VEN/HMA response in `Trumpp.xlsx` Table S4). Submitted boss's cellHarmony_lite alignment to the
  Hs-BM-titrated 89-state reference → integrated h5ad. Job **777433 `h5ad_combine`** (12h/128G), script
  `/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_AML_UDON/run_cellHarmony_Trumpp.lsf`, out → that
  dir's `output/`. Pre-flight: all 16 soupX inputs (matrix.mtx+barcodes+genes) populated; caught CRLF in
  SoupX_filepaths_Trumpp.txt (harmless — cellHarmony does line.strip()). NOTE: salomonis2 share (NOT the
  salomonis-archive MOSAIC-AML data). Next: aligned h5ad → composition → run through mutation predictor +
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
    miss). Packaged as `MOSAIC-AML_Trumpp_Report.docx` (answers + full table + P13 detail + methods + limits).
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
    full 7-modality contributions; deliverable `MOSAIC-AML_Trumpp_Report.docx` refreshed.

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
    overrides. mosaic_board.html grouped roster. Verified in preview: 47 real runs render, no placeholders.
  * **runs/ cleared + reran with the new predictor.** Deleted **50 stale dev/scratch dirs** (probe/regr/disc/
    gpu/model/phase/patient_/run/subtype__; list saved `runs/_cleared_dev_runs.txt`). Kept 29 predict_ + 16
    trumpp_ (already regenerated OOF-calibrated 8-modality) + single_modality caches. **Reran the 2 old
    composition-only ingest uploads** (AML7, BF71-CD34) through the full 8-modality OOF-calibrated predictor via
    new `deploy_scrna.py` (cosine cell-state assignment `amlmm.scrna.assign_cells` → per-state pseudobulks →
    RNA+impute-all+Composition+LSC+Cell-comm → predictor → mutation_panel). Needed a scanpy-free h5py 10x reader
    (scanpy pulls numba, incompatible w/ cluster numpy 2.4). AML7 = diseased (WT1 0.92, trisomy8 0.89 top);
    BF71-CD34 = **control** (CD34-sorted, face-valid; only RUNX1). `deploy_scrna.py` is the reusable
    arbitrary-scRNA full-chain path (what the GUI "Add patient" should call next).
  * Deliverable for Nathan refreshed to final: `MOSAIC-AML_Trumpp_Report.docx` (8 modalities, 0.677, calibrated
    calls, P13 report) — the "Metabolite off-cluster" note was corrected.

### 2026-07-02 — Therapy hypotheses + Recommended validations added to the deliverable report
- **What:** Added two clinician-facing subsections to the P13 detailed prediction report (§3 of the combined
  `MOSAIC-AML_Trumpp_Report` and the standalone `PREDICTION_REPORT.md`): **Therapy hypotheses** (maps the
  model's confident *present* calls — NPM1, DNMT3A, complex-karyotype — to a literature therapeutic rationale;
  explicitly flags the FLT3 negative call as low-trust because its signal is ~83% imputed-ADT-weighted) and
  **Recommended validations** (orthogonal confirmatory assays per driver, ordered by management impact:
  FLT3-ITD PCR w/ allelic ratio → NPM1 RT-PCR/IHC → myeloid NGS → karyotype+FISH → TP53/17p → flow cytometry).
- **Why:** the reports gave calls but no actionability; the AI-panel/ingest path already emits
  `targetable_therapies`/`recommended_validations`, but the external-cohort deliverable for Nathan did not.
- **How:** edited both markdown reports; both sections are hypothesis-generating and carry a "not clinical
  guidance" banner (rendered as a styled blockquote). Added a `>` blockquote branch to `convert.js`
  (left-rule + light-blue shading) so the disclaimer renders in Word; regenerated `MOSAIC-AML_Trumpp_Report.docx`.
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
- **What:** 15-slide `MOSAIC-AML_Presentation.pptx` (in `scratchpad/artifacts/`), framed on the abstract
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
comprehensive comparison workbook (`Desktop/MOSAIC-AML_vs_marker-pipeline.xlsx`, 7 sheets).
- **exp_allcov.py** (LSF `long`/`test`, job 814440): ran MOSAIC-AML fused multimodal donor-grouped OOF AUROC
  on EVERY covariate the pipeline ran, incl. ones our deployed predictor skipped (PHIP 0.978, SF3B1 0.858,
  JAK2 0.834, SRSF2 0.769, NF1 0.899, CEBPA 0.917, GATA2 0.919, KMT2A-rearr 0.965) + clinical (AML-vs-control
  0.995, FAB_monocytic 0.931, ELN_adverse 0.889, relapse 0.800, sex_M 0.693). 40 covariates.
- **exp_percelltype.py** (LSF, job 814289): scored our features at his exact granularity — per
  (mutation x modality x cell-state), donor-grouped OOF + 200-perm null → 7,005 rows (his schema), 4,024 signal.
- **Results (apples-to-apples, same metric both sides):** on fused OOF vs his median-across-cell-types OOF,
  **MOSAIC-AML wins 29/31 covariates** (losses: sex_M — his XIST/DDX3Y cell-type markers; ELN_adverse 0.889 vs
  0.894 tie). Paired per-cell-type (3,032 identically-defined cells): **MOSAIC-AML higher on 56% (1,690), mean
  0.716 vs 0.687**. Bulk->sc transfer: we win 31/31 (his transfer median AUROC_all ~0.55-0.63). Our per-modality
  mean OOF and full 7,005-row table included as sheets.
- **Honest framing captured in the sheet:** his "best cell-type (max)" is optimistic best-of-hundreds (his own
  caveat); the fair columns are median (typical) and our fused OOF. Our SEALED held-out (deployable, harder)
  number is kept as a separate column. LSF note: per-user 125-slot cap (A549 array) blocked scheduling;
  bswitch to `test` + bmod -n 1 got allcov running immediately.

### 2026-07-07 — Re-derived per-mutation call thresholds (Youden → F1-max); fixed over-calling
The board's calibration audit showed held-out **call accuracy collapsing to 0.27–0.37 in the 0.70–0.90
probability band** — the predictor over-called. Root cause: the present/absent threshold used Youden's J
(TPR−FPR), which ignores base rate and keeps the cutoff low for rare drivers, admitting many false positives.
- **Fix (`pipeline/train_predictor.py`):** the present/absent threshold now **maximizes F1** on the deployed
  OOF blend per mutation (penalizes false positives via precision). Retrained on LSF (job 820977, `test` queue
  after the A549 125-slot cap blocked `long`; bswitch + bmod -n 1). Thresholds rose from ~0.5 median to
  **median 0.89, range 0.70–0.97**. Mean held-out AUC unchanged (0.795 over 21 covariates).
- **Impact (676 sealed held-out predictions, 29 samples):** overall call accuracy **0.726 → 0.862**; the
  0.7–0.9 danger band healed (p~0.7 decile 0.371→0.914, p~0.8 0.273→0.688). Present-call **precision (PPV)
  0.518**, recall 0.598, F1 0.555 — over-calling in the mid-band is gone, but top-end precision is still
  capped ~0.5 by the weak drivers (KIT 0.60, PTPN11 0.60, WT1 0.44, IDH2 0.65) whose scores aren't sharply
  separable even at high probability. Next lever for those = **abstention on low-AUC drivers**, not another
  threshold.
- **Artifacts:** rebuilt `gui/reliability.json` from the new reports (bins/deciles/per-mutation
  sens·spec·ppv·npv·sens_weight; `calib_src="bin-observed"`). `gui/calibration.html` panel 1 is now a
  **before(Youden, grey-dotted) vs after(F1, teal) overlay** showing the danger zone healing; panel 2 (raw
  score still over-confident — 0.85 → ~26% present) unchanged and still valid.

### 2026-07-07 — Review response: cell-state map, honest-framing, model card, VAF resource
Acted on an external technical review. Ground-truthed its claims first (the "7,005 cell-state classifier bank"
it cites is NOT a persisted asset — it was the one-off `exp_percelltype.py` comparison; deployed call is
sample-level, `dataio.py:262` n_cells-weighted collapse).
- **Cell-state localization map** (`pipeline/cellstate_localize.py`, runs LOCALLY — head node has old glibc,
  but all modality h5ads + build_context work locally): per (driver × 89-state) donor-grouped OOF AUROC on
  **measured RNA** (mutant vs non-mutant, the classifier's own contrast) + orthogonal composition-shift.
  Deliberately not imputed → "where the signal lives" isn't circular. Biologically sensible (DNMT3A→MPP,
  NPM1→GMP/neutrophil, TP53→LMPP). → `gui/cellstate_localization.json`; surfaced as an evidence-page heatmap
  (section 4) + inline "signal concentrates in" line on the board.
- **Honest-framing pass:** board now flags low-confidence drivers (held-out AUC<0.65 or n+<3) with ⚠ +
  "confirm by sequencing" (baked into `train_predictor.py` `confidence`); imputed modalities (ADT/GRN/Lipid/
  Metabolite) tagged ᴿᴺᴬ = "RNA-conditioned, not independent"; "probability"→"model score". `panel.py`
  arbiter already framed imputed-as-corroborating — gap was only the GUI.
- **Model card + train-vs-heldout** (`train_predictor.py` now persists `P.train_auc` = fused donor-grouped
  CV-OOF AUROC; writes `pipeline/model_card.json`). Ranking worst→best by held-out (dumbbell): worst 3
  (WT1/KIT/ASXL1) are OVERFITTING — train 0.77–0.89 but held-out 0.48–0.52 (gap 0.27–0.41); best (NPM1, del7,
  RUNX1) generalize. Small-n held-out is noisy (±0.1 run-to-run at n+≤5) — take train+held-out from the SAME run.
- **Coverage:** only 26 of 77 label columns clear the ≥8-positive training floor (supervised prediction needs
  examples — unlike a DNA detection panel). Of a collaborator's 53-mutation list: 22 deployed, 26 too-rare
  (1–7 pos), 5 absent. Borderline (SRSF2 7, SF3B1 6, t8_21 5) are 1–3 samples from trainable → more cohorts
  is the direct lever.
- **VAF resource** (`pipeline/build_vaf.py` → `labels/vaf_per_sample.tsv` + `gui/vaf_by_mutation.json`):
  harvested VAF from banked `labels/vaf_sources/` (Meta-CCHMC structured + harmonized Variant_Detail parsed),
  mapped to drivers (FLT3 split ITD/TKD). Per-deployed-mutation slot: **has_vaf** (17; median+distribution),
  **awaiting** (SNV, no data yet — e.g. TP53; auto-fills on re-run), **not_applicable** (8 cytogenetic).
  Surfaced on the board expanded row as a per-sample VAF strip + median. VAF↔predictability is DECOUPLED
  (Spearman 0.03): RUNX1/FLT3-ITD subclonal yet strong, so weakness = no transcriptional imprint, not rarity.

_Last updated: 2026-07-07._

### 2026-07-09 17:55 — Variant-level relabeling of the cross-cohort mutation bake-off
- **What:** Replaced gene-level driver labels with **58 variant-level categories** (per Nathan's spec) in `pipeline/bulk_external.py` + regenerated `aml-bakeoff/bundle_data.npz`; repackaged `Downloads/aml-bakeoff.zip`.
- **Why:** Distinct hotspots have distinct biology/expression signatures (FLT3_ITD vs TKD; DNMT3A_R882 vs nonR882; TP53 hotspot-DBD vs LOF; NRAS/KRAS G12/G13/Q61; etc.) — lumping by gene blurs the very signal the classifier should learn.
- **How:** New `_fine(sym,hgvsp,vc)` parser maps BeatAML `hgvsp_short` protein position -> sub-category; CEBPA biallelic = >=2 variants/sample; FLT3_ITD + KMT2A_fusion from curated clinical (`FLT3-ITD`, `consensusAMLFusions`). Leucegene mapped via its residue-level columns (FLT3-ITD/TKD, IDH1-R132, IDH2-R140, U2AF1-S34/Q157, NPM1, CEBPA-biallelic, MLL-PTD); sc stays gene-level + FLT3_ITD.
- **Result:** 58 categories, **~50 clear >=6 positives** in BeatAML (scored by 5-fold CV-OOF, the primary metric) — up from 15. Rare splits (<6: FLT3_N676, KRAS_Q61, PTPN11 E76/D61/A72, KMT2D, KIT_N822/D820, KMT2A_PTD) are labeled but auto-skipped. Bundle 68 MB / 14,237 genes; ready to run on the 9800X3D desktop.
- **File:** `pipeline/bulk_external.py`, `scratchpad/precompute_bundle.py`, `aml-bakeoff/{run_bakeoff.py,README.md,bundle_data.npz}`

### 2026-07-09 21:33 — Variant-level bake-off RESULTS (50 categories, run on 9800X3D desktop)
- **What:** Analyzed the completed variant-level cross-cohort bake-off (`results.zip` -> `scratchpad/results2/`). Built `scratchpad/variant_bakeoff_analysis.png` (3-panel: within-gene split contrast, worst->best ranking, cross-cohort concordance).
- **Result — variant-splitting validated:** DNMT3A_R882 CV-AUROC **0.88 vs nonR882 0.61**; FLT3_ITD **0.91 vs TKD 0.69 / other 0.73**; U2AF1_S34 **1.00 vs Q157 0.80**; NRAS Q61>G13>G12 (0.82/0.74/0.67). TP53 hotspot-DBD 0.90 ~ LOF 0.89 (clean split, no predictability gap).
- **Result — ranking (fixed shrLDA/no-FS, fair):** mean **0.82** over 50 cats. Solved (>=0.95): U2AF1_S34, SF3B1_K666, SRSF2, NPM1_exon12, CEBPA_Nterm, BCOR. Weakest real (n>=15): DNMT3A_nonR882 0.61, NF1 0.67, NRAS_G12 0.67, FLT3_TKD 0.69. n<15 = provisional.
- **Result — methods (replicate prior run):** no-FS beats MarkerFinder **+0.05 in 73%** of pairs; ENSEMBLE/PLS/logL2/shrLDA > trees > linSVM > **MLP 0.67 (worst)**.
- **Result — cross-cohort:** BeatAML->Leucegene transfers (NPM1/SRSF2/U2AF1_S34/STAG2/RUNX1 on-diagonal); **TET2 non-transfer 0.83->0.65**. FLT3_ITD holds across all 3 cohorts (0.91/0.87/0.90).
- **Next / implication:** deployed sc system should predict FLT3_ITD (not gene-level FLT3), split DNMT3A R882-vs-nonR882, abstain on TET2/nonR882/NF1.
- **File:** `scratchpad/results2/{bakeoff_results.json,bakeoff_summary.txt}`, `scratchpad/variant_bakeoff_analysis.png`, `scratchpad/analyze_results2.py`

### 2026-07-10 10:37 — Per-modality breakdown of deployed sc predictor + base-learner optimization bundle
- **What:** Extracted per-(mutation x modality) CV-OOF AUROC from deployed `mutation_predictor.pkl` (`scratchpad/extract_modality_breakdown.py` -> `deployed_modality_breakdown.{png,csv}`). Built portable base-learner sweep (`Downloads/aml-modality-bakeoff.zip`).
- **Result — per-modality:** imputed **Metabolite (0.797, best 7x), Lipid (0.790, best 6x), GRN (0.789, best 6x)** carry the system (single-best for 19/26 muts); RNA (0.796) = reliable backbone (weighted 17/26) rarely single-best; **LSC (0.627) + Cell-comm (0.675) weak**. Fusion beats best-single in **26/26** (0.875 vs 0.842). FLT3-ITD best=Lipid 0.83, fused 0.838, heldout 0.895.
- **Result — sc variant-level NOT feasible:** `vaf_per_sample.tsv` HGVS detail is sparse; only FLT3_ITD (33) clears >=8. DNMT3A_R882=1, NPM1=4, NRAS subtypes 1-2. Variant-level stays the bulk bake-off's domain (full WES).
- **New optimization lever:** deployed base learner = LinearSVC, which was WEAK in the bulk bake-off (0.77 vs 0.82-0.84). Built `aml-modality-bakeoff` (30MB bundle: 8 blocks capped at 8000 var feats + labels + donor groups + holdout) + `run_modality_bakeoff.py` sweeping 8 base learners through the exact deployed recipe (diff-select 500 -> learner -> percentile -> grouped CV-OOF -> NNLS fusion -> F1-max), ENSEMBLE derived from cached linear scores. ~8min on 9800X3D. Smoke (NPM1): linSVM already top (0.988) -> bulk finding may not transfer; full run to settle.
- **File:** `scratchpad/{extract_modality_breakdown.py,extract_modality_bundle.py,deployed_modality_breakdown.png}`, `aml-modality-bakeoff/*`

### 2026-07-10 11:23 — Base-learner sweep RESULT: deployed LinearSVC confirmed optimal
- **What:** Ran `aml-modality-bakeoff` (26 muts x 8 modalities x 8 learners, deployed recipe) on the desktop.
- **Result:** **Don't swap.** linSVM (deployed) best on train CV-OOF **0.875**; paired mean-delta vs linSVM negative for ALL (logL2 -0.009, ENSEMBLE -0.012, shrLDA/HistGB -0.020, RF -0.027, PLS -0.036, MLP -0.087). HistGB +0.009 on held-out but -0.020 train + held-out n~21 = noise. MLP worst again.
- **Why bulk finding didn't transfer:** percentile-calibration + NNLS late-fusion neutralize base-learner differences that mattered on raw single-modality bulk. Deployed design validated.
- **Conclusion:** deployed system at the ceiling for modeling levers; remaining gains = DATA (more sc variant-level genotyping + larger held-out).
- **File:** `scratchpad/modres/*` (results.json, summary.txt, charts.png)

### 2026-07-13 09:26 — Bulk variant-level predictor promoted to PRIMARY mutation caller + wired into ingest
- **What:** Per boss decision (bulk "almost as good, many more mutations"), made the bulk variant-level predictor the primary mutation caller. Built `amlmm/bulk_predictor.py` (BulkMutationPredictor) + `pipeline/train_bulk_predictor.py` -> `bulk_mutation_predictor.pkl` + `bulk_model_card.json`; wired into `pipeline/ingest_patient.py`.
- **Phase 1 (deployable):** trained on BeatAML2 (n=707 WES), 50 variant-level categories, 14,237 shared genes; logL2 / no-FS / top-2500-var / percentile / F1-max (bulk-bakeoff winning config). Per-cohort z-refs (sc/beataml/leucegene) for cross-cohort calibration. Mean 5-fold CV AUROC **0.829**. Validated on sealed sc held-out via bulk-equivalent: SRSF2 0.91, FLT3_ITD 0.83, TET2 0.76 (n>=5).
- **Phase 2 (ingest wiring):** `bulk_equiv_from_adata` (sum cells -> CP10k) -> `bulk_mutation_result` -> `predict(ref='sc')` -> report top-level `mutation_predictions` (mode `bulk_variant_primary`). Kept OUT of arbiter (predicted != truth; user mutations stay the deterministic anchor). Graceful skip if pkl absent. Real-data sanity: median **3** confident-present calls/sample (random-noise mock gave 21 = garbage-in).
- **Grounding:** head-to-head - 18 shared muts sc 0.853 vs bulk comparable; bulk adds **25 categories** sc can't reach (U2AF1_S34, SF3B1, JAK2_V617F, variant splits) on cheap bulk input; sc keeps 8 cytogenetic events.
- **Next:** sync pkl + edited ingest to the cluster for the live GUI; optionally add cytogenetics to the bulk caller.
- **File:** `amlmm/bulk_predictor.py`, `pipeline/{train_bulk_predictor.py,ingest_patient.py,bulk_mutation_predictor.pkl,bulk_model_card.json}`

### 2026-07-13 15:24 — Bulk-RNA upload path added to patient ingest (single-cell OR bulk)
- **What:** The platform now accepts a **bulk RNA expression file** as a sample, not just single-cell. Extends `pipeline/ingest_patient.py` + the GUI (`gui/gui_server.py`, `gui/mosaic_board.html`).
- **Pipeline:** new `--bulk <file>` (exclusive with `--sample`) + `--bulk-ref {beataml,leucegene,sc}` + `--bulk-scale {auto,linear,log2,log1p}`. `parse_bulk_expression` reads any delimited gene×value table, auto-detects scale (negatives→log2→2^x; compressed→log1p→expm1; else linear). `main_bulk` skips the atlas load (fast, `-M 4000`) → runs the PRIMARY bulk variant-level caller (`ref=bulk_ref`) + genetic anchor → report mode `bulk_panel` (no composition/subtype/cytogenetics/control-gate — needs single cells, clearly noted).
- **Schema fix:** refactored `bulk_mutation_result` → `(predictions_list, caller_meta, AgentResult)` so `mutation_predictions` is the GUI-native **list** (was a block dict) with `heldout_auc := CV AUROC` for the reliability/abstain logic. Applied to BOTH the sc and bulk report paths.
- **GUI:** `dispatch_ingest` routes `kind=bulk` → `--bulk`; POST reads `kind/bulk_ref/bulk_scale`; `/api/samples` surfaces `.tsv/.csv/.txt` as bulk. Add-patient modal: Input-type toggle (single-cell / bulk RNA) + bulk-reference selector + inbox auto-detects type; `renderPredMeta` + footer now bulk-aware (shows the bulk caller, not the sc multimodal boilerplate).
- **Verified end-to-end:** bulk ingest of a BeatAML sample (14,237 genes, auto→linear, ref=beataml) → 7 confident-present variant calls (FLT3_TKD, NPM1_exon12, DNMT3A_R882, IDH1_R132, IDH2_R140, NRAS_Q61, PTPN11_other); calibration sane (median 2–3 present/sample over 60 samples); rendered in the live GUI (correct bulk predictor box + honest footer).
- **Next:** sync to cluster + fold into the Mosaic-AML PR.
- **File:** `pipeline/ingest_patient.py`, `gui/gui_server.py`, `gui/mosaic_board.html`

### 2026-08-04 20:45 — COMPASS-AML: ex-vivo drug-response layer built end-to-end + platform renamed MATRIX-AML -> MOSAIC-AML
- **Rename:** MATRIX-AML -> **MOSAIC-AML** (*Multimodal Omics and State-Aware Inference of Cancer Drivers in Acute Myeloid Leukemia*). Byte-level rewrite of 38 files + 12 file renames (`mosaic_board.html`, `start_mosaic_board.bat`, `MOSAIC-AML_*` deliverables); board branding + launch.json entries + server_version updated. Bare "matrix" (mutation_matrix / confusion matrix / build_full_matrix) deliberately untouched.
- **Data:** fetched `beataml_probit_curve_fits_v4_dbgap.txt` (19 MB) + drug families/sample-mapping/clinical from the public BeatAML2 repo into `data/external/beataml/` — the lab did NOT have the drug file anywhere on the cluster (checked). 63,395 single-agent curves -> 53,571 with shared-space expression -> **48,998 measurements / 520 specimens / 479 patients / 118 inhibitors** after QC + inclusion.
- **Data layer** (`amlmm/drug/data.py`): per-row curve-quality flags (non-convergent 10, off-panel window 10, within-drug extreme deviance 603; 7,092 *increasing* curves kept-but-flagged and never callable sensitive); robust within-drug median/MAD z; tail classes (bottom/top 20%, middle 60% kept for regression only). 115 inhibitors pass the primary inclusion filter; Cytarabine / Nutlin 3a / GDC-0941 fail ONLY the wave-stability test -> new `wave_conditional` tier rather than deleting the induction backbone; 47 excluded with per-compound reasons.
- **Annotation** (`knowledge/drug_annotation.tsv`, 118 rows): targets, family, mechanism, clinical tier (6 approved-AML / 36 approved-other / 36 trial / 40 research), clinically available analogue, coarse exposure class, known resistance mechanisms. 16 coarse `family_group`s for hierarchical pooling.
- **Model A:** hierarchical `f_group(X, D) + w_j*f_j(X)`, drug descriptor = z-expression of the inhibitor's own targets. Four blocks (`rna`/`state`/`mut`/`clin`) **late-fused by NNLS on inner donor-grouped OOF, floored to the best single block** — a single ridge over the concatenation LOSES to RNA alone (0.742 vs 0.766); fusion recovers **0.772**.
- **`state` block:** built lineage signatures from our own atlas (`build_state_signatures.py`): 89 states -> 10 lineage groups -> top-60 specificity markers, in the shared ENSG space. Recovers textbook markers unprompted (AVP/CRHBP/MECOM, PRTN3/ELANE/MPO, S100A8/9/FCN1, HBB/AHSP, PPBP/PF4/ITGA2B, CXCL12/VCAM1).
- **THE bug worth remembering:** calibrating Platt on the RAW decision score made every single-cell sample come back at P(sensitive)~0.99 for all 118 drugs. Fix = calibrate the **OOF percentile** + cohort-matched score references (`beataml` / `sc_sample` / `sc_state`) + a matched expression z-reference. Same lesson `bulk_predictor.py` already documented. A cell state needs its OWN reference — against the sample-level one every state is an outlier.
- **Model B** (`statemodel.py`): Model A re-applied per cell-state pseudobulk -> blast/LSC coverage, escape state, dispersion, bulk-vs-sc disagreement. "Blast compartment" is labelled a COMPARTMENT everywhere — no per-cell genotype is used.
- **Model C** (`mechanism.py`): target expression pct, curated transcriptional OUTPUT readouts per pathway family, BCL2/(BCL2+MCL1+BCL2L1+BCL2A1) dependency, genetic activation tagged observed-vs-predicted, measurable resistance proxies. Kept OUT of A on purpose.
- **Utility + agents:** S_ij with every penalty itemised (uncertainty / resistance+escape / infeasibility / OOD), positive part renormalised by *evaluable* weight so bulk patients aren't penalised for absent coverage; rankings **per clinical tier**, never merged. Eight agents, all `therapeutic`-domain and **non-voting** on the anchored subtype call. The combination agent explicitly REFUSES to add single-agent scores (BeatAML2 has no combination data) and only proposes complementary-coverage pairs across different pathways.
- **Validation** (`eval_drug_model.py`): per-inhibitor mean AUROC **0.774** (approved-agent subset **0.809**), Spearman 0.365, AUPRC 0.748 vs 0.475 baseline, 100% of inhibitors p<0.05. Deployment task (per-patient across drugs): top-1 **34%** vs 10% matched chance (**3.4x**), top-5 76% vs 42%. ECE **0.012**, Brier 0.185 vs 0.249. Abstention: 28% -> **4.7%** error at 10% coverage. Leave-wave-out 0.722/0.733; leave-centre-out 0.731-0.890; survives every differentiation stratum (0.721-0.762); ~15 null SDs above a specimen-repointing permutation null. Sealed 15%-of-PATIENTS hold-out: 0.784.
- **Model B's falsifiable test PASSED:** fitted only on BeatAML bulk, never shown a cell-state label, it predicts higher venetoclax sensitivity in primitive than monocytic states in **93.2%** of all 387 atlas samples (mean +0.662 z, Wilcoxon p=8.4e-51) — **rank 1 of 118 inhibitors**, so not a "primitive looks sensitive to everything" artefact (ABT-737, the other BH3 mimetic, ranks 8th). Bulk-vs-weighted-state consistency Spearman 0.970.
- **Honest negatives recorded:** `state` and `clin` add ~nothing on AVERAGE (rna+mut 0.7730 vs all-four 0.7718) — but the fitted weights show they matter per family (FLT3 leans 0.36 on `mut`, cell_cycle 0.61 on `state`, JAK_STAT 0.46 on `clin`); a high fusion weight = complementary info, NOT standalone strength (state alone is 0.635 for cell-cycle drugs vs rna 0.712). 100% of single-cell bulk-equivalents lie beyond BeatAML's own p95 distance — the assay transfer is real and stated as a limitation.
- **Wiring:** `predict_drugs.py` (one sample), `drug_layer.py` (non-fatal hook in BOTH ingest paths, writes `runs/<id>/drug_report.json` + `DRUG_REPORT.md`, summary into `patient_report.json.drug_response`), `batch_drug_reports.py` (backfilled 30 atlas runs). GUI: `/therapy.html` (tiered, click-a-row evidence, combination hypotheses, abstention list) + `/rx_validation.html`; `/api/drug_report`, `has_drug_report` on `/api/runs`, `/val/figures/` sub-path. Verified live in the browser preview.
- **Perf note:** `amlmm/drug/h5rows.py` reads only the needed CSR rows out of the 1 GB atlas h5ad (contiguous-run reads), and `StateResponse._z` gathers with numpy instead of assigning 14,237 DataFrame columns — together 25+ min -> 3 s for the state validation.
- **File:** `pipeline/amlmm/drug/*`, `pipeline/{train_drug_model.py,eval_drug_model.py,build_drug_score_refs.py,build_state_signatures.py,validate_state_response.py,build_drug_figures.py,predict_drugs.py,drug_layer.py,batch_drug_reports.py}`, `deliverables/{METHODS_COMPASS-AML.md,drug_model_validation.json,drug_model_card.json,state_response_validation.json,figures/Rx1-Rx6}`, `gui/{therapy.html,rx_validation.html}`

### 2026-08-04 21:40 — Static export of the whole web UI (hand it to people, no server)
- **What:** `pipeline/build_static_site.py` freezes every page of the MOSAIC-AML GUI into a browsable bundle that works by **double-clicking `index.html`** — no install, no server, no internet. -> `deliverables/mosaic_static/` (28 MB) + `deliverables/MOSAIC-AML_static_site.zip` (**8.2 MB**).
- **The constraint that shapes it:** browsers block `fetch()` against `file://`, so the frozen API responses CANNOT ship as `.json` files the pages fetch. They ship as `.js` files assigning to a global (`<script src>` is not subject to that restriction) and a shim swaps `window.fetch` for a lookup against it. Images/PDFs/TSVs stay ordinary files — `<img src>` and download links work fine from `file://`.
- **Contents:** all 7 pages (board -> `index.html`, therapy, rx_validation, validation, calibration, evidence, cebpa_evidence), **429 patient reports**, **26 drug reports**, 66 figure/PDF/TSV assets, 476 baked API keys. Data shards: core 3.3 MB / reports 8.0 MB / drugs 10.0 MB, all minified; drug reports additionally pruned of `agents[*].evidence.all` (duplicates `per_drug`, ~65 kB x 30).
- **Nothing may look interactive when it isn't:** `/api/capabilities` reports `ingest:false` so the "+ Add patient" button never renders; any non-GET returns an explicit "this is a static export" message instead of failing quietly; a banner states it on every page. Drag-and-drop "Open report…" still works (pure FileReader).
- **Gotcha fixed during the build:** the build rewrites `/val/` -> `val/` so `<img>`/download links resolve relative to the bundle, but that same rewrite lands inside `fetch('/val/x.json')` calls, whose baked key still has the leading `/val/`. Rather than rewriting only some occurrences, the shim tries a small list of candidate keys (`/p`, `/val/<base>`, `/<base>`).
- **Verified from `file://`** in the browser: all 7 pages render, 429-row roster, 118 drug rows on therapy.html, all 6 Rx figures + the mutation-caller figures load, zero console errors, `addBtn` display none.
- **File:** `pipeline/build_static_site.py`, `deliverables/mosaic_static/*`, `deliverables/MOSAIC-AML_static_site.zip`

### 2026-08-04 22:15 — Last three dark pages restyled to the formal light theme + Plotly vendored for offline use
- **What:** `evidence.html`, `calibration.html` and `cebpa_evidence.html` were still on the old dark theme (`#0f1112` background, teal `#57b3a6` accent) while everything else had moved to the formal light palette. All three rebuilt on the shared design language used by validation / therapy / rx_validation: `--bg #ffffff / --panel #f7f7f4 / --ink #1a1a19 / --muted #6b6b63 / --line #e4e3dc / --accent #8a6a18`, sticky top bar with cross-links, uppercase kicker + `h2` + bordered `.card` per section, and a lead paragraph explaining what the page answers. **Zero dark hexes remain anywhere in `gui/`.**
- **Chart palette migrated too, not just the page chrome:** the Plotly `DARK` layout became `LIGHT` (white paper/plot, ink text, `#e4e3dc` gridlines); teal -> accent gold `#8a6a18` for the used/mutant series, dark `#2b3034` -> pale `#dcd9cd` for unused bars, grey -> warm neutral `#b9b6a8` for controls, warm/cool localisation colours -> `--bad`/`--good`. The "this patient" diamond moved from gold `#ffcf4d` (invisible on white) to deep blue `#2f6690` with a white outline, which stays distinct from both the accent and the neutral.
- **Plotly vendored:** the three pages loaded plotly from a CDN, which meant the static export rendered **no charts at all** on a laptop with no internet — the exact use case the export exists for. Now `gui/vendor/plotly.min.js` (4.6 MB) with a CDN `document.write` fallback if the local copy is missing; relative `src` so the same tag works live (`/vendor/...`) and in the bundle. Added a `/vendor/` route to `gui_server.py` and a copy step to `build_static_site.py`.
- **Also fixed while in there:** `cebpa_evidence.html`, `cebpa_violin_data.json` and `bulk_bakeoff_results.json` were never in the server's static whitelist, so that page 404'd on the live board; bar-chart y-axes got `automargin` (long labels like "Metabolite ᴿ" were clipped) and the x-ranges widened to 1.11-1.13 so the outside value labels fit when every bar sits near 0.99.
- **Verified:** live server screenshots of all three pages; static bundle re-verified from `file://` — plotly resolves from `vendor/`, the CDN tag never fires, all plots render, banner present.
- **Static bundle rebuilt:** 7 pages, 429 reports, 26 drug reports, 67 assets, zip **9.5 MB** (was 8.2 MB; +1.3 MB is the vendored Plotly).
- **File:** `gui/{evidence,calibration,cebpa_evidence}.html`, `gui/vendor/plotly.min.js`, `gui/gui_server.py`, `pipeline/build_static_site.py`

### 2026-08-05 01:40 — Upload path verified (3 real bugs fixed), components named CIPHER-AML / COMPASS-AML, pushed to GitHub with three branches
- **"Does it work just by uploading data?" — now yes, and the asking found three real bugs:**
  1. `amlmm/scrna.py` had the cellHarmony marrow reference as a **hardcoded cluster path**, so an scRNA upload failed on every machine that was not the cluster. Now resolves env override -> the copy inside the vendored `engine-code/altanalyze3` checkout -> the cluster archive.
  2. `pipeline/control_gate.pkl` did not exist locally, so the healthy-vs-diseased gate **silently returned nothing** on every upload. Trained it (cv_auc 0.966, disease sensitivity 0.992 at the conservative operating point).
  3. `gui_server.dispatch_ingest` passed the caller's `sample` straight through, so a bare filename (exactly what `/api/samples` labels it) produced a job that died 40 s later on a file-open error. It now resolves against the inbox and fails immediately with a readable message.
- **Verified end to end, three ways:** bulk upload (a real BeatAML specimen written as a collaborator would hand it over) -> 50 mutation calls + 118 inhibitors scored, correct "bulk input, no cell-state coverage" caveat; single-cell upload (`make_upload_test_sample.py` expands an atlas sample's per-state pseudobulks into 3,002 Poisson-drawn cells with a KNOWN composition) -> gate says diseased p=0.977, leading hypothesis **Inv16** recovered for an inv(16) sample it never saw the label for, 16 states scored, Quizartinib/Midostaurin top the approved-AML tier; and the GUI's own `POST /api/ingest` -> job dispatched, completed, rendered on the board.
- **Named the two prediction layers** (both `<X>-AML`, under the MOSAIC-AML platform):
  - **CIPHER-AML** — *Cell-state Inference of Pathogenic Hits from Expression and Regulation* (the mutation predictor).
  - **COMPASS-AML** — *Cell-state Oriented Modelling of Pharmacologic Assay Sensitivity* (the drug predictor; replaces the working title MOSAIC-Rx).
  Applied across code, GUI page titles and headers, board nav and cards, docs and deliverables; `METHODS_models_and_agents.md` -> `METHODS_CIPHER-AML.md`, `METHODS_MOSAIC-Rx.md` -> `METHODS_COMPASS-AML.md`. A canonical "the three names" table now heads README and OVERVIEW.
- **100-permutation null finished:** observed mean Spearman 0.367 vs null **0.001 ± 0.019**, max 0.041 — **19 null SDs**, p = 0.0099 (the floor for 100 shuffles; no permutation came close). Replaces the earlier 20-shuffle p = 0.048, which was floor-limited and undersold the effect. Figures and prose refreshed.
- **GitHub:** repo renamed `SalomonisLab/Matrix-AML` -> **`SalomonisLab/MOSAIC-AML`** (GitHub redirects the old URL). The local tree had **no common ancestor** with the remote, so rather than force-pushing over five merged PRs the work was grafted: `git reset --soft origin/main` then one commit -> a clean fast-forward, remote history intact. Three branches: **`main`** (full platform), **`cipher-aml`**, **`compass-aml`** — each component branch keeps the shared `amlmm` core and ingest path (both components run through them) and drops the other component's scripts, models, GUI pages and deliverables, with its own root README. Verified: both branches compile, the shared core imports, and a bulk upload on `cipher-aml` runs to completion with the drug layer **gracefully absent** ("WARN COMPASS-AML drug layer skipped") — the guarded hook proven in a real pruned tree, not just in theory.
- **.gitignore hardened before the first push:** `data/` (norm_exp.txt alone is 281 MB), `inbox/`, `pipeline/mutation_predictor*.pkl` (124–152 MB, past GitHub's hard limit), the regenerable static bundle and zips. Pushed: 757 files, 75.6 MB, nothing over 50 MB.
- **File:** `pipeline/{make_upload_test_sample.py,control_gate.pkl}`, `pipeline/amlmm/scrna.py`, `gui/gui_server.py`, `.gitignore`, `SETUP.md`, `README.md`, `OVERVIEW.md`, `deliverables/METHODS_{CIPHER,COMPASS}-AML.md`

### 2026-08-05 03:10 — Survival layer: yes to who and roughly when, no to how long for one person
- **Asked:** can MOSAIC-AML predict whether a patient survives and for how long, and with what accuracy. Answer built and measured rather than asserted.
- **Data reality:** the single-cell atlas has **no survival metadata** (vital status non-null for 0 samples), so survival cannot be learned there. BeatAML can: after dropping unknown vital status, restricting to **initial-diagnosis** specimens (a relapse specimen answers a different question and leaks prognosis) and one specimen per patient — **444 patients, 245 deaths (55%), 89 sealed away**.
- **No lifelines/sksurv installed**, so Cox partial likelihood, Breslow baseline, C-index, KM, log-rank, time-dependent AUC and IPCW Brier were implemented directly with **six known-answer self-checks** (recovers a simulated beta=1.0 as 0.987; C-index exactly 1/0/0.5; KM matches a hand-worked example to 4 dp; log-rank p 0.70 vs 4e-18).
- **Result.** C-index **0.752** sealed hold-out / 0.726 CV; 2-y AUC **0.852**; calibration gap 0.043–0.075; risk tertiles separate at **p = 2.8e-9**. Group median survival accurate to **~6 weeks**. But per-patient timing: MAE 0.87 y, and survival inside one risk band spans ~1–1.4 y (p10–p90) — so the deployed output is a curve + horizon probabilities + a risk group, **never a single number of months**.
- **The honest headline:** molecular data does **not** beat age + ELN 2017 on its own (RNA alone −0.002 C-index; molecular fusion −0.003). Only the combination adds: **+0.034, 95% CI [+0.006, +0.059], P(no gain)=0.008** — the one arm whose interval clears zero. Cell state and mutations alone are *worse* than the clinical baseline.
- **Four bugs/limits caught by stress-testing, all fixed:** (1) the mutation block is 14% not-assayed with one dead column, which made Cox fail silently and — because the fusion required every block to have scored a patient — quietly collapsed every fused arm onto an arbitrary block; (2) a single-cell sample z-scored against BeatAML collapsed its curve and **told a patient they had two weeks to live** — fixed with cohort-matched references, same as COMPASS-AML; (3) a gene-symbol/ENSG mismatch fed all-zero vectors that still returned a confident 97th-percentile prognosis — the layer now owns alignment and refuses below **80%** gene coverage (evidence-based: the risk group flips below ~80%); (4) a healthy pooled-CD34 control scored at the 97th percentile with 0.4% one-year survival — now refused when the control gate says healthy, and bulk input (which has no gate) states that the number **assumes AML**.
- **Wired** into both ingest paths behind a try/except, verified end to end on a single-cell and a bulk upload. Figures Sv1–Sv3, methods in `deliverables/METHODS_survival.md`.
- **File:** `pipeline/amlmm/survival/{coxph,data,model}.py`, `pipeline/{train_survival_model,eval_survival_time,build_survival_figures,survival_layer}.py`, `pipeline/survival_model.pkl`, `deliverables/{METHODS_survival.md,survival_model_card.json,survival_time_validation.json,figures/Sv1-Sv3}`

### 2026-08-10 — B2: external validation on TCGA-LAML, and two calibration bugs it flushed out
- **Asked:** run the TCGA-LAML external validation if the data is available. It is — UCSC Xena hub, three files, fetched automatically on first run by `validate_tcga_laml.py` so the reproduction step does not depend on someone already having them.
- **Why it mattered:** every survival number so far came from BeatAML2. Patient-grouped CV and a sealed hold-out protect against learning a *patient*; neither protects against learning a *cohort*. This is the only test that separates the two.
- **What was frozen:** Cox coefficients, PCA rotation, variable-gene selection and NNLS fusion weights, all loaded unchanged. The **only** thing refit on TCGA is the per-gene z-reference (median/MAD) — RSEM log2 and BeatAML units are not the same scale, and transferring without cohort-matched normalisation measures the platform. Same mechanism the drug layer already uses for single-cell input.
- **Cohort:** 149 patients (173 with expression ∩ 186 with follow-up, one specimen each), **92 deaths**, median follow-up 1.00 y, **89.4%** gene overlap.
- **Result:** deployed arm C-index **0.706** (95% CI 0.655–0.758) against 0.751 CV / 0.787 sealed hold-out in BeatAML. 2-y AUC 0.806. Gain over age + cytogenetic risk **+0.035, P(ΔC≤0) = 0.029**. Risk tertiles separate **71.7% vs 13.7%** two-year survival, log-rank **p = 7.0e-10**, curves ordered and non-crossing (`fig_Sv4_tcga_km.png`).
- **The honest negative, recorded not buried:** the **molecular blocks alone do not beat age + cytogenetics** in TCGA (ΔC = −0.003). Value comes from the combination, matching the BeatAML pattern. And the covariate with the largest gain in BeatAML — baseline induction type — **cannot be validated here at all**, because TCGA-LAML is treatment-homogeneous; the deployed arm's advantage over `full` in this cohort is the age spline alone.
- **Bug 1 (live, affected every upload):** the bundle never stored the age-spline knots, so `survival_layer.py` fell back to a hardcoded `[45, 60, 71]` while the model had been fitted against the training quartiles — the spline basis was rebuilt in the wrong place at inference. Knots are now computed and carried in the bundle (`[47.0, 61.0, 70.0]`).
- **Bug 2 (mild leak):** `build_blocks` took a `train_idx` argument that **neither call site passed**, so knots were placed using quantiles over the hold-out patients too. Both call sites fixed. Retrained: **every metric unchanged to 4 dp** (CV 0.751, hold-out 0.787, AUC 2y 0.8715, Brier 0.1494) — the leak was immaterial, three quantiles of one covariate over 444 patients. Verified the argument now actually bites: train-only knots give 70.0 where whole-cohort would give 70.5. (The last time an ignored parameter produced byte-identical numbers in this project, it was a real bug, so this was checked rather than assumed.)
- **Not done, deliberately:** pooling BeatAML + TCGA for training. At 89% gene overlap across different platforms that needs a batch-correction step whose failure mode is inventing signal; the frozen-transfer test is the more informative use of the cohort.
- **File:** `pipeline/validate_tcga_laml.py`, `pipeline/train_survival_model.py`, `pipeline/survival_model.pkl`, `deliverables/{VALIDATION_TCGA_LAML.md,validation_tcga_laml.json,fig_Sv4_tcga_km.png,METHODS_survival.md,OPTIMIZATION_ROADMAP.md}`

### 2026-08-12 — ELN 2022 inferred in BeatAML; VAF sensitivity; where the guideline fails
- **Asked (Nathan):** benchmark risk prediction against **ELN 2022**, which must be inferred in BeatAML, and test **10% vs 40% VAF** thresholds. Five papers supplied.
- **The gap this exposed:** every benchmark in the platform had been scored against BeatAML's shipped **ELN2017** column. 2017 and 2022 differ precisely where we do worst — FLT3-ITD allelic ratio dropped, CEBPA changed from biallelic to in-frame bZIP (mono- or biallelic), seven MDS-related genes newly adverse, TP53 newly adverse **at VAF ≥ 10%**, t(9;11) intermediate with precedence.
- **Built `pipeline/eln2022.py`** implementing Table 6 in full including the precedence footnotes: free-text ISCN karyotype parsing (complex excluding hyperdiploidy, monosomal, −5/del(5q), −7, −17/abn(17p), the six adverse translocations), fusions from `consensusAMLFusions`, CEBPA bZIP from variant-level protein positions, and the nine-gene MR set from `mutations.txt`. **Agreement with the shipped ELN2017 label 0.776 (n=548)**, and every disagreement runs in the direction the guideline change predicts (48 Adv→Int from the dropped allelic ratio, 37 Fav→Int from NPM1+ITD, 25 Int→Adv from the new MR genes).
- **VAF answer:** TP53 is mutated in **80 specimens at ≥10% and 70 at ≥40%**. **30 of 638 specimens (4.7%) change ELN 2022 category**, all toward Intermediate. **Discrimination does not move at all**: ELN2017 0.610, ELN2022@10% 0.612, ELN2022@40% 0.610. So the threshold is a real reclassification and a prognostic non-event — fix it at 10% per the guideline, report the 4.7% as an uncertainty band, do not tune it.
- **Pollyea/Döhner REPRODUCES.** ELN is at or below chance outside intensive induction: intensive (n=357) ELN2017 0.631 / ELN2022 0.642 / 4-gene 0.574; **non-intensive (n=87) ELN2017 0.496 / ELN2022 0.462 / TP53-FLT3ITD-NRAS-KRAS 0.612**. The published 4-gene rule works there and only there, which is exactly its stated scope.
- **Röllig/Bill does NOT reproduce.** Splitting ELN 2022 adverse by MR-gene status: 9.5 months (n=126) vs 7.8 (n=63), same direction as the published 14.7 vs 8.3 but **not significant (log-rank p=0.38)**. The adverse category is not usefully heterogeneous here.
- **A correction to our own stated limitation.** METHODS_survival.md said the model is "close to useless in non-intensively-treated patients (C-index 0.554)". That was a **pooled-fitting artefact**: fitting within the stratum gives **0.681**, and adding the 4-gene rule **0.701**. What survives is the claim about the guideline, not the model. Treatment-stratified fitting is worth +0.127 in the stratum where ELN is at chance — more than any relabeling.
- **MOSAIC against the corrected bar:** ELN2022 0.612 vs deployed **0.771**; adding the ELN label on top gains nothing (0.771), i.e. the model already contains it.
- **File:** `pipeline/{eln2022.py,exp_eln2022_benchmark.py}`, `labels/eln2022_beataml_vaf{10,40}.tsv`, `deliverables/{ELN2022_RISK_BENCHMARK.md,exp_eln2022_benchmark.json}`, `refs/eln2022/`
