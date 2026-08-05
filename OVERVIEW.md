# MOSAIC-AML — Project Overview (master document)

***M**ultimodal **O**mics and **S**tate-**A**ware **I**nference of **C**ancer Drivers in **A**cute **M**yeloid **L**eukemia.*
(Formerly MATRIX-AML; renamed 2026-08-04 — same platform, same code, new name.)
Last updated 2026-06-15. This is the single source of truth: what the project is,
how everything works, where everything lives, and what is done vs. still ahead.

## The three names

| name | what it is |
|---|---|
| **MOSAIC-AML** | the platform — *Multimodal Omics and State-Aware Inference of Cancer Drivers in Acute Myeloid Leukemia*. The atlas, the ingest path, the witness panel, the arbiter, the decision board. |
| **CIPHER-AML** | the **mutation predictor** — *Cell-state Inference of Pathogenic Hits from Expression and Regulation*. Infers driver lesions from expression, cell state and regulon activity, without sequencing. |
| **COMPASS-AML** | the **drug predictor** — *Cell-state Oriented Modelling of Pharmacologic Assay Sensitivity*. Predicts ex-vivo inhibitor sensitivity from the BeatAML2 functional screen, resolved by cell state. |

CIPHER-AML and COMPASS-AML are deliberately **parallel** layers, not a chain: a mutation call enters the
drug layer as one piece of evidence among several, never as a look-up key for a therapy.


---

## 1. The goal

Build an analytical framework that **triages AML patients toward therapies** from a
single cheap input (**bulk RNA-seq or scRNA-seq**), by:

1. reconstructing a full multimodal single-cell profile (modalities imputed from RNA),
2. having a panel of **per-modality expert agents**, each an independent-ish
   "biological witness," reason over its modality's **honest, leakage-proof evidence**,
3. an **arbiter agent** reconciling them into a consensus that accounts for every
   modality, weighting genuinely independent evidence above imputed-from-RNA evidence,
4. producing a per-patient report: AML subtype, targetable driver programs
   (explainable), confidence-ranked **therapy hypotheses**, and the **validation assays**
   (flow markers from ADT, targeted metabolomics/lipidomics, ex-vivo drug tests) that
   would confirm them.

The scientific thesis: not "multimodal integration" but **mechanism-resolved triage** —
recommend action only where multiple independent modalities converge on a targetable state.

---

## 2. The data

**Atlas** (control-normalized pseudobulks; "up to 90 cell populations" per specimen):

| Modality | Shape (obs × var) | Status | Notes |
|---|---|---|---|
| RNA | 12,255 × 35,702 | **measured** | the source signal; `pseudobulk_counts_hashed.h5ad` |
| GRN (TF→target edge activity) | 12,255 × 7,486 | imputed from RNA | `var['heldout_spearman']` per feature |
| Metabolite | 12,255 × 2,486 | imputed from RNA | median held-out ρ≈0.27 |
| Lipid | 12,255 × 1,009 | imputed from RNA | median held-out ρ≈0.36 |
| ADT (surface) | 12,255 × 129 | imputed from RNA | → flow-panel nominations |
| cell-communication | 383 × 141,101 | derived (fastComm) | sample-level; not yet a witness |
| cell-frequency (composition) | 397 × 90 | **measured** | default predictive feature block |
| UDON RNA programs | 16 programs over 8,564 pseudobulks | computed | control-normalized disease programs |

- **Specimens loaded by the pipeline:** 387 (45 control-annotated, ~342 disease),
  11 datasets, 316 donor groups (63 donors have >1 specimen → grouped CV matters).
  Nathan's framing is 333 AML + 70 control (=403) — **reconcile the cohort definition** (open item):
  loaded breakdown = disease_category AML 247 + unlabeled 100 + Control 24 (metadata; 45 by obs
  annotation) + MDS 12 + T-ALL 4. So ~16 fewer specimens than the stated cohort and a control-count
  mismatch (45/24 vs 70) — align the canonical set with Nathan.
- **Genetic labels** live in the obs `Annotation` (28 driver classes) — richer than the
  clinical TSV. Canonicalized (NPM1c→NPM1, FLT3-ITD→FLT3, …): NPM1 43, Inv16 18, FLT3 14,
  TET2 12, TP53 10, then a long tail.
- **Mutation matrix** (derived by `genetics.py` from `Annotation` + free-text `karyotype`
  + ELN): 184/387 have any genetic data; prevalence NPM1 43, complex 28, trisomy8 23, del7
  19, inv16 19, del5 15, FLT3 14, TP53 10, KMT2A 9; targetable FLT3/IDH1/IDH2/NPM1/KIT.
- **Clinical metadata** is **sparse** and now an `.xlsx` (`Metadata/AML_metadata_CLEAN.xlsx`,
  402×22; upstream swapped it from `.tsv` on 2026-06-15 — the loader handles both):
  vital_status 0, overall_survival ~20, WHO ~40, FAB ~48, ELN ~70, age ~96, sex ~171,
  karyotype ~158, Donor_ID full.

**The central caveat (drives every design decision):** ADT/GRN/metabolite/lipid are all
**imputed from RNA**, so they are **not independent witnesses** — several agreeing is not
independent confirmation. Genuinely independent axes = **RNA, cell frequency, genetics,
and any externally-measured validation**. The framework treats imputed modalities as
interpretive lenses weighted by imputation fidelity, and the arbiter is told this explicitly.

---

## 3. Where everything lives

**Cluster** (`ssh bmiclusterp-head`, = bmiclusterp2.chmcres.cchmc.org, LSF 10.1, queue `normal`):
- Data deposit (original scattered layout): `/data/salomonis-archive/LabFiles/Nicholas/AML-multimodal/`
  (modality dirs `RNA/`, `GRN/`, … ; `Metadata/`; `LSC-prediction/algorithm/sample_cellstate_counts.tsv`;
  `RNA/clusters/UDON_final_program_assignments.tsv`).
- Pipeline code: `…/AML-multimodal/pipeline/` ; run outputs: `…/AML-multimodal/runs/<run_id>/`.
- Python: `/usr/local/anaconda3-2020/bin/python` (sklearn 1.5.2, anndata 0.10.8, numpy 1.26;
  **no openpyxl** → stdlib xlsx reader). **Canonical run env.**
- LLM gateway (LiteLLM): `http://bmiclusterp2.chmcres.cchmc.org:4000/v1`; models
  `nemotron-3-super` (used), `gpt-oss-120b`, `qwen3.6-35b-a3b`, `qwen3-embedding-4b`.
  Key embedded in `amlmm/llm.py` (override via `AMLMM_LLM_API_KEY`/`_BASE_URL`/`_MODEL`).
- Engine source (UDON etc.): repo **github.com/SalomonisLab/altanalyze3** (`components/udon/`);
  per-user checkouts exist on the cluster. **scALABLE imputation + cellHarmony-web are NOT on
  the cluster** — only Nathan's machine produced the imputed `.h5ad`.

**Local working copy** (Windows, dev): `C:\Users\krog5w\.gemini\antigravity\scratch\AML-multimodal\`
- `data/` (curated row-aligned matrices ~1.5 GB), `labels/` (metadata + composition + UDON +
  LSC TSVs), `engine-code/altanalyze3/` (GitHub clone), `pipeline/` (the code, mirrors cluster),
  `runs/` (local run outputs), `gui/` (zero-dep decision-board viewer over `runs/`; see `gui/README.md`),
  `README.md`, `OVERVIEW.md` (this file), `requirements-lock.txt`.
- Local Python has the full stack incl. anndata (sklearn 1.9 — newer than cluster; see §8).

**Memory** (persists across sessions): `…/.claude/projects/…/memory/aml-multimodal-project.md`.

---

## 4. Architecture — how it all works

Two layers, deliberately separated:

```
INPUT (bulk or scRNA-seq) → deconvolve frequency + impute modalities (AltAnalyze3)   [upstream]
        │
        ▼
DETERMINISTIC CORE  (amlmm/ — reproducible, leakage-proof, never agent-touched)
  dataio → feasibility → assemble_features → classify (grouped nested CV) → cluster_explore → report
        │                                   └── gate (DecisionHooks)
        ▼
AGENT LAYERS
  (a) decision-seam hooks: AgentHooks routes gate/feature/model/report decisions to Nemotron
  (b) per-modality panel: witnesses (predictive / genetic / UDON) → arbiter → consensus
        │
        ▼
OUTPUT  runs/<run_id>/ : feasibility.tsv, cv_result.json, final_model.joblib,
                         run_report.json + REPORT.md  (pipeline) ;  panel_report.json + PANEL.md (panel)
```

### 4a. The deterministic pipeline (one config = one LSF job)
`run.py` runs the step sequence for one configuration (a target × CV strategy × feature
blocks). `submit.py` fans out many configs as parallel `bsub` jobs (`bsub -L /bin/bash -n
<slots> -M <mb> -W <wall> -R span[hosts=1]`, default queue, thread-pinned, `bjobs`
fail-closed status). The headline metric is an **honest out-of-fold balanced accuracy**:

- **donor-grouped or leave-one-cohort-out** nested CV (`cv.py`); a hard assertion forbids a
  donor in both train and test; model+hyperparameters selected inside each fold; scaling
  fit inside each fold (sklearn Pipelines).
- a **group-level permutation p-value** is the chance reference (the trustworthy signal).
- the **gate** (`DecisionHooks.gate_result`) decides pass/rerun/abort on that p-value.

### 4b. Agent integration
- **Decision-seam hooks** (`hooks.py`): `DecisionHooks` = auditable defaults; `AgentHooks` =
  the same seams driven by Nemotron (gate/skeptic, feature-block choice, model choice, report
  synthesis), each validating + falling back to deterministic. `run.py --hooks agent`,
  `submit.py … --hooks agent`. The agent owns *which/whether/write-up*, never the CV math.
- **Per-modality panel** (`panel.py`): the MOSAIC-AML "tumor board." Witness kinds:
  - **predictive** (composition, or a modality block) — grounded in that block's honest CV.
  - **genetic** — mutation/cytogenetic context + targetable drivers (independent axis).
  - **UDON** — conserved program ↔ subtype associations (discovery; flags batch programs).
  - the **arbiter** reconciles all witnesses, weighting independent > imputed, distrusting
    batch-dominated programs, into a consensus + per-witness consistency + targetable summary
    + recommended validations. Every agent call capped at **500,000 tokens** (i.e. untruncated).

---

## 5. Component inventory (`pipeline/`)

| File | Role |
|---|---|
| `run.py` | orchestrator / one-config LSF job body (provenance, gate-abort, nonzero exit on error) |
| `submit.py` | fan-out LSF submitter + `status` (completed/errored/aborted/CRASHED?/CORRUPT) |
| `panel.py` | MOSAIC-AML panel CLI (witnesses + arbiter) |
| `amlmm/dataio.py` | layout-aware loader (local vs cluster), stdlib `.xlsx` reader, modality aggregation, cached cohort-matrix baseline + direct cell-communication reader |
| `amlmm/context.py` | `Config` + `Context` (tables, artifact store, hooks, atomic IO) |
| `amlmm/cv.py` | **leakage-proof grouped/nested CV** + group-level permutation baseline |
| `amlmm/models.py` | fold-safe model zoo (rf / logreg / gboost), thread-pinned |
| `amlmm/step.py` | `StepSpec`/`StepResult`/registry/`run_step` |
| `amlmm/targets.py` | prediction targets + label extraction |
| `amlmm/hooks.py` | `DecisionHooks` / `AgentHooks` / `GateDecision` — the decision seams |
| `amlmm/llm.py` | OpenAI-compatible client (stdlib) → Nemotron; 500k default cap |
| `amlmm/cluster.py` | LSF `bsub`/`bjobs` backend (lab conventions, fail-closed) |
| `amlmm/genetics.py` | mutation/cytogenetic matrix from annotation + karyotype + ELN |
| `amlmm/udon.py` | UDON program ↔ subtype/dataset associations (SATAY-UDON style) |
| `amlmm/panel.py` | the 9 witness gather/assess fns (predictive, genetic, UDON, LSC, surfaceome/ADT, metabolic, lipid, GRN-regulon, cell-communication) + `run_panel` (cohort) + `run_patient_panel` |
| `amlmm/agent.py` | `AgentSpec`/`AgentResult`/registry/`run_agent` (witness analogue of `step.py`); grounding + independence weight factors |
| `amlmm/agents/__init__.py` | witness roster: registers all 9 as AgentSpecs, `default_roster` (de-dups by name) |
| `amlmm/ledger.py` | shared **evidence ledger** (immutable evidence + revisable opinion, `deterministic_evidence_hash`) |
| `amlmm/arbiter.py` | **arbiter v2**: deterministic genetic anchor + concordance + additive harvest (surface therapies / validations / descriptive findings); LLM narrates only |
| `amlmm/feedback.py` | **Phase C** guarded deliberation loop: deterministic deference to the anchor over bounded rounds + round-0 baseline + convergence + groupthink-drift report |
| `amlmm/retrospective.py` + `validate.py` | **Phase D** retrospective clinical validation: do the engine's anchored-driver / LSC calls associate with ELN / response (Fisher + permutation, underpowered-aware, circularity-flagged) |
| `amlmm/knowledge.py` + `knowledge/*.tsv` | versioned curated KB (`kb-2026.07`): biomarker→drug, validation rules; LLM may only order within it |
| `amlmm/steps/*` | feasibility, assemble_features, classify, cluster_explore, report |
| `docs/AGENT_INTEGRATION.md` | how agents plug into the seams |
| `README.md`, `requirements-lock.txt` | usage + canonical versions |

---

## 6. How to run

```bash
# on the submit host (bmiclusterp-head), cd …/AML-multimodal/pipeline
PY=/usr/local/anaconda3-2020/bin/python

# deterministic pipeline, fanned out as LSF jobs
$PY submit.py configs                 # list the config matrix
$PY submit.py submit                  # fan out all configs (add --hooks agent for LLM seams)
$PY submit.py status                  # per-config state

# one config in-process (the job body)
$PY run.py --target subtype --strategy donor_kfold

# MOSAIC-AML per-witness panel (cohort)
$PY panel.py --target subtype                          # composition + genetic + UDON
$PY panel.py --target subtype --blocks composition,ADT,GRN   # add predictive modality witnesses
# MOSAIC-AML per-PATIENT tumor-board report
$PY panel.py --patient "WashU::10DD-1002__Diagnosis"             # subtype + therapies + validations
$PY panel.py --patient "WashU::10DD-1002__Diagnosis" --rounds 2  # + Phase C deliberation
```
Outputs land in `…/AML-multimodal/runs/<run_id>/`. Develop locally (same commands, layout
auto-detected); the cluster is the canonical run environment.

### GUI — decision board + uploader (`gui/`)
A zero-dependency board (stdlib `http.server` + one HTML file) that VIEWS per-patient reports and
UPLOADS new scRNA samples to process. Run it on the cluster login node and tunnel in:
```bash
python3 gui/gui_server.py 8765            # on the cluster login node (server is stdlib-only, 3.6-safe)
ssh -L 8765:127.0.0.1:8765 bmiclusterp-head   # from your laptop -> open :8765
```
The server runs under any python3; ingest jobs use the analysis python on the compute node
(`/usr/local/anaconda3-2020/bin/python`, override `AMLMM_PYTHON`).
Viewing shows the anchored leading hypothesis, ranked therapies/validations, the 9-witness evidence
panel with the grounding × independence weighting made visible, recorded conflicts (how the
observed-driver anchor overrules an imputed disagreement), and the Phase C deliberation. **+ Add
patient** uploads a sample (10x dir / `filtered_feature_bc_matrix.h5` / `.h5ad`; optional observed
mutations) → `POST /api/ingest` dispatches `ingest_patient.py` as an LSF job → live status in the
sidebar → the report joins the roster when done. Also opens standalone (drag a `patient_report.json`
onto the HTML, view-only). See `gui/README.md`.

### New-patient entry path (`ingest_patient.py`)
Turns a NEW external gene-level scRNA sample into a decision-board report:
gene-level scRNA → **composition** (cosine to the 89-state cellHarmony reference; `amlmm.scrna`,
sparse-aware) → **cohort-trained** subtype prediction (RF on the atlas composition→subtype, applied
out-of-cohort: `held_out=False`) + optional **genetic anchor** from supplied mutations → arbiter →
`runs/<id>/{patient_report.json, ledger.json, PATIENT.md, status.json}`. Validated end-to-end on the
Grimes **AML-CITE-Seq** 10x samples (AML-7 → FLT3/medium anchor-free; BF71-CD34 → NPM1/medium;
supplying a driver flips to anchored/high). Imputed descriptive witnesses (GRN/ADT/metabolic/lipid/
LSC/UDON/cell-comm) are not yet wired for uploads (listed under `ingest.deferred`); the arbiter
handles a partial roster. The LLM may lower confidence but can't raise it above the deterministic
ceiling for an unconfirmed (anchor-free) call.
```bash
python ingest_patient.py --sample <10x_dir|h5|h5ad> --name "Patient X" [--mutations "TP53,FLT3"]
```

---

## 7. Status — done vs. to-build

**Done & validated**
- Curated local data copy + organization; altanalyze3 engine cloned.
- Deterministic pipeline: dataio (incl. xlsx + layout auto-detect), leakage-proof grouped/nested
  CV + permutation baseline, model zoo, steps, registry, gate. Runs as **per-config LSF jobs**
  (validated: job 78149/78274 produced all artifacts).
- 5-dimension **adversarial audit** + fixes (permutation baseline, fail-closed status, atomic
  writes, single-donor-group class drop, LOCO guard, deterministic tie-break, provenance).
- **AgentHooks live** (Nemotron) — gate/skeptic overruled the deterministic rule on fold variance.
- **Mutation matrix** (`genetics.py`) and **UDON associations** (`udon.py`) — sensible output.
- **Per-modality panel v1** (`panel.py`): predictive + genetic + UDON witnesses, each
  Nemotron-assessed with sensible weights (composition ~0.85 measured, genetic ~0.48 down-weighted
  for sparse coverage, UDON ~0.82); **arbiter validated** — it produced a concordance-based consensus
  (Inv16 + NPM1 best-supported across composition + genetics + UDON), a targetable-driver list
  (FLT3/IDH1/IDH2/NPM1-menin/KIT) and recommended validations (flow/ADT panels + targeted seq + ex-vivo
  drug tests), with per-witness consistency. Per-agent cap = **500,000 tokens** (untruncated).
- **Per-patient mode** (`panel.py --patient <sample_key>`, atlas patients via honest out-of-fold
  predictions): each witness gives a per-patient read; the arbiter returns a subtype call + targetable
  therapies + recommended validations + conflicts (`PATIENT.md`). Validated on an NPM1 case (→ NPM1+Inv16
  signature, menin inhibitor + venetoclax/aza, NPM1/Inv16/flow validations) and a TP53 case (composition
  mispredicted FLT3 0.84, but the **genetic + UDON witnesses correctly caught TP53/del7** — the panel
  surfaced the conflict). **This TP53 mis-lead is now CLOSED in Phase A** (deterministic genetic anchor).
- **Phase A — blackboard backbone + arbiter v2 (done & validated 2026-06-16).** Agent registry
  (`agent.py`), shared **immutable evidence ledger** (`ledger.json` with a `deterministic_evidence_hash`),
  a versioned **knowledge layer** (`knowledge/*.tsv`, `kb-2026.06`; the LLM may only order within it,
  never invent), and **arbiter v2** with a **deterministic genetic anchor**: a *present* observed driver
  (mutation or cytogenetic — incl. APL/t(15;17)→ATRA+ATO, t(8;21), KMT2Ar→menin inhibitor) leads, and a
  single imputed prediction that disagrees is recorded as a conflict but **cannot outrank it**. Therapies
  key only on observed drivers; the LLM narrates rationale/confidence but cannot change the decision (the
  deterministic pre-pass IS the fallback). Verified on the cluster with real Nemotron on an LSF compute
  node: **the TP53 case now leads TP53** (composition still mispredicts FLT3 0.84, now a recorded
  conflict), NPM1 leads NPM1, both `confirmed_by_genetics`. Hardened by a 5-dimension **adversarial agent
  review** (16 confirmed findings fixed — incl. a *critical* shared-ctx DataFrame-truthiness bug that
  silently wiped the anchor on the 2nd patient onward); committed regression test
  `_phaseA_regression.py` (all pass, cluster + local).
- **Phase B — complete the 9-witness roster (done & validated 2026-06-16, `kb-2026.07`).** Added six
  specialized expert agents, each on its own modality: **LSC** (stemness/risk classifier call,
  `classifier_call`/`rna_derived`), **surfaceome/ADT** (imputed surface markers → flow/immunotherapy
  hypotheses; CD33→gemtuzumab, CD123→tagraxofusp, flow-pending), **metabolic** + **lipid** (imputed,
  fidelity-gated top features), **GRN-regulon** (imputed TF activity, fidelity UNKNOWN→flagged), and
  **cell-communication** (independent fastComm L-R signaling, direct sample-level reader). All six are
  **non-voting/descriptive**: they use new domain strings the arbiter's vote/anchor branches never read
  and contribute only through an **additive harvest** (therapy_biomarkers / validation_claims /
  descriptive_context) — surface therapies are kept in a SEPARATE flow-pending list, imputed witnesses
  never key a driver therapy, and harvested validations are allow-listed to descriptive claim types. The
  guarantee is mechanically verified: **leading hypothesis / driver therapies / concordance are
  byte-identical with vs without the 6 new witnesses** (anchor invariance). Verified on the cluster with
  real Nemotron (full 9-witness roster): NPM1→NPM1 (high, 0.829), TP53→TP53 (medium, 0.654, composition
  →FLT3 a recorded conflict), 6 descriptive findings each. Hardened by a 5-dimension **adversarial agent
  review** (12 confirmed findings fixed — incl. a *high-sev* GRN witness/​block name collision that
  crashed the run via the ledger guard, and a cluster-only cell-communication backed-CSR **segfault**
  fixed by a non-backed reader); committed regression test `_phaseB_regression.py` (anchor-invariance,
  surface-positive, no-anchor, missing-modality, cold-path determinism — all pass, cluster + local).
- **Phase C — guarded continuous feedback loops (done & validated 2026-06-16, opt-in `--rounds N`).**
  `amlmm/feedback.py` `deliberate()`: round 0 is the immutable independent **baseline**; in rounds
  1..N a *voting* witness whose round-0 vote conflicts with a **genetically-confirmed** leading
  hypothesis DETERMINISTICALLY **defers** (down-weights its opinion toward the observed driver), applied
  once → converges. Guards: evidence is immutable (the anchor is recomputed from evidence every round, so
  it has the final say); a `drift` report flags `groupthink_warning` only if the leading changed AND the
  final is not genetically anchored; `max_rounds` cap + convergence; `mode = continuous |
  conflict_triggered`. **What deliberation can do:** refine concordance/confidence (a conflicting imputed
  witness defers → agreement rises). **What it can NEVER do:** change the leading hypothesis, the driver
  therapies, or the anchor — proven by an adversarial test (a conflicting voter at weight 1e9 still can't
  flip the anchor). `run_patient_panel(..., max_rounds, mode)` (default 0 = single pass, fully backward
  compatible); CLI `panel.py --rounds N --mode ...`; `PATIENT.md` gains a deliberation section.
  Cluster-validated with real Nemotron: TP53→TP53 (0.655→0.792 over 1 deliberation round, converged),
  NPM1→NPM1 (0.851→0.919), `groupthink_warning=False`, evidence hash invariant to round count. Hardened
  by a 4-dimension **adversarial review** (12 findings, all low/nit; fixed the one real cosmetic
  inconsistency — `stop_reason` divergence between ledger and report — and added deliberation-round
  counting); regression `_phaseC_regression.py` (loop-safety incl. 1e9 adversarial, determinism,
  convergence/cap, immutability, drift watchdog, single-pass contract — all pass, cluster + local).
- **Phase D (partial) — retrospective clinical validation (done & cluster-validated 2026-06-16).**
  `amlmm/retrospective.py` + `validate.py`: tests whether the engine's deterministic cohort outputs
  (the genetic-anchored driver + its ELN-expected risk; the LSC stemness call) ASSOCIATE with the
  sparse clinical labels present (ELN_risk n=70, clinical_response n=52, survival n=19; `vital_status`
  is empty). Fisher-exact 2×2 + seeded label-permutation concordance, **underpowered-aware** (n + flag
  on every test), reproducible local==cluster. **Honest findings:** (1) the engine's anchored-driver
  risk reproduces the clinician's ELN bucket (concordance 0.744, perm p=5e-05; adverse-driver→ELN-adverse
  OR=144) — but this *shares logic* with how ELN is defined, so it validates the engine's driver
  **extraction + anchor correctness**, not novel prediction; (2) the genuinely independent cross-modality
  signal (RNA-derived p-LSC stemness ↔ genetic ELN-adverse) is **underpowered/null** (OR=2.6, p=0.45,
  n=44) — reported honestly, not as a discovery; (3) the favorable-driver→response result is flagged
  **suspect** (the `clinical_response` column is heterogeneously encoded across cohorts). A real
  environment-dependent bug (None-vs-NaN mask leak inflating an apparent p-LSC signal to a false p<1e-4)
  was caught + fixed so n's/ORs are identical local↔cluster. `_phaseD_regression.py` (stats helpers,
  binarization, explicit-support, coverage — all pass). Output: `runs/phaseD_validation/VALIDATION.md`.

**Partial**
- Genetic witness mutation matrix is still from **sparse free-text karyotype + annotation only**
  (184/387 have any genetic data); no structured variant feed.
- UDON witness uses RNA programs only; the per-modality UDON clusters (ADT/GRN/Lipid/Metabolite) are
  present locally but keyed by sample tokens that don't map cleanly to `sample_key` — deferred (would
  need a token→sample_key crosswalk + per-modality program→subtype association maps).
- Descriptive witnesses (metabolic/lipid/GRN/cell-comm) run **patient mode only**; cohort mode returns a
  skipped stub (the cohort path still uses `run_panel` + the legacy LLM reconcile).

**Not started / blocked (remaining Phase D stretch)**
- **New/external-patient mode (UNBLOCKED; integration build remaining)** — per-patient works for *atlas*
  patients (held-out OOF); a genuinely new bulk/scRNA upload needs the deconvolve→impute entry path. The
  assets are on the shared archive at `…/LabFiles/Nathan/Revio/altanalyze3/` (group-readable): all four
  `rna2{grn,adt,lipid,metabolite}` imputation bundles + `cellHarmony`. Confirmed loadable on the cluster
  Python — `from altanalyze3.components.rna2grn import load_bundle` (with `…/Revio/altanalyze3` on
  `sys.path`) → `.predict_from_adata/_h5ad/_dataframe/_10x_h5`. **Bulk keystone tested + found hard:**
  `amlmm/bulk.py` builds a cell-state×gene signature from the RNA atlas + NNLS deconvolution; the
  round-trip (`_bulk_validate.py`, compute node) shows it does NOT recover the 89 fine cell-states even
  in-sample (median r≈0.18–0.22, dominant-state match 6/40) — the states are too collinear for bulk
  deconvolution. So bulk needs **coarser composition or a dedicated/regularized deconvolver** (real
  research); the **scRNA path avoids it** (per-cell classification, not unmixing). **scRNA keystone
  built + VALIDATED END-TO-END (2026-06-17)** — `amlmm/scrna.py` cosine-assigns each query cell to the
  cellHarmony BM reference (`Hs-MarrowAtlas-L3M.txt`, 2870 markers × 89 populations = the atlas states
  *exactly*, verified) → composition. Proven on real gene-level scRNA: the **Grimes AML-CITE-Seq** 10x
  samples (`…/Nathan/Collaborators/Grimes/AML-CITE-Seq/{AML-7,BF71-CD34}/filtered_feature_bc_matrix.h5`),
  GEX-only via scanpy `read_10x_h5` (drops Antibody-Capture/ADT). AML-7 (8495 cells, 2833/2870 markers,
  cosine 0.57) → monocyte/myeloid-skewed; BF71-CD34 (10461 cells, cosine 0.62, 80/89 states) →
  stem/progenitor-dominated (HSC/MPP/GMP/MEP/LMPP) — matches the CD34 sort (face-valid) and uses the
  SAME 89-state names as the cohort reports (drop-in). `assign_cells` made sparse-aware → 300MB on a
  compute node (`_scrna_cite_validate.py`). [The old `Revio/MDS-AML-KINNEX-1/*.h5ad` are long-read
  JUNCTION data, not gene-level — use the 10x/gene-level path.] **Remaining:** per-state pseudobulk →
  `rna2*` impute → panel (a new patient has no CV OOF → composition witness needs a train-on-cohort
  predict-new path; no genetics → the Phase C no-anchor consensus applies).
- **Per-modality UDON crosswalk (DEFERRED, low value)** — the ADT/GRN/Lipid/Metabolite `*_udon_clusters`
  files key on lossy/ambiguous sample tokens (45–87% coverage, no Dataset prefix) and have no
  program→subtype association map; usable only as descriptive cluster-membership context.
- **Validation-panel generator** (ADT→flow panels, GRN→phospho-flow, metabolite/lipid→targeted assays)
  and **explainable counterfactual axis report** (per-patient driver-program decomposition).
- **Retrospective validation — done (Phase D above)** for ELN/response on the in-house cohort; still to
  do: **external independent test datasets** and a proper survival model (n=19 here is too small).
- Cohort-count reconciliation (387 loaded vs 333+70).

---

## 8. Caveats & scientific constraints (read before trusting any number)

1. **Imputed ≠ independent.** Down-weight imputed modalities; reserve "confirmation" for RNA,
   frequency, genetics, and measured validation. (The panel/arbiter already enforce this.)
2. **Labels are the bottleneck.** Survival/drug-response can't be *trained* (n≈20/24). Therapy
   outputs are **knowledge-grounded hypotheses + the validation to confirm them**, tested
   retrospectively — not learned predictors.
3. **Grouped CV is mandatory** (donor + cohort); enforced by assertion. The permutation p-value,
   not the raw accuracy, is the trustworthy signal.
4. **Small-n + version sensitivity.** Subtype n≈97; the point estimate shifts ~0.06 across sklearn
   versions (cluster 1.5.2 canonical; provenance recorded in every `run_report.json`).
5. **Honesty over flattery** is the design goal — the correct output on this dataset is often
   "significant but unstable; validate / get more data."

---

## 9. Roadmap

1. **Data foundation (in progress):** mutation matrix ✓, UDON-program grounding ✓; reconcile
   cohort counts; load per-modality UDON clusters + program marker genes.
2. **Full witness roster (done, Phase B):** all 9 witnesses wired — composition (predictive) + genetic
   + UDON + LSC + surfaceome/ADT + metabolic + lipid + GRN-regulon + cell-communication ✓, each with
   honest grounding/independence tags; the 6 new ones are non-voting/descriptive (anchor-invariant).
3. **Arbiter v2 (done, Phase A) + continuous feedback (done, Phase C):** deterministic concordance +
   genetic anchor + KB-grounded ranked therapies/validations + additive descriptive/surface harvest ✓;
   **guarded multi-round deliberation** (deterministic deference to the anchor; round-0 baseline +
   groupthink drift) ✓. Still to add: per-patient counterfactual driver-program decomposition.
4. **Action layer:** curated program→drug and program→flow-panel priors → ranked therapy
   hypotheses + tiered validation panels.
5. **Validation + bulk path (Phase D):** retrospective ELN/response association harness ✓ (done,
   `validate.py`; honest + underpowered-aware — it validates driver-extraction correctness and reports
   the independent cross-modality signal as a null at current n). Still to do: external independent test
   datasets, a survival model (n too small here), and the bulk/scRNA deconvolve→impute entry point —
   now UNBLOCKED: the `rna2*` bundles + cellHarmony are group-readable on the shared archive
   (`…/Nathan/Revio/altanalyze3/`) and load on the cluster Python. Bulk deconvolution keystone tested
   (`amlmm/bulk.py`) and found too weak at 89-state resolution (in-sample r≈0.2) → bulk needs coarser
   composition or a dedicated deconvolver; **scRNA route recommended** (cellHarmony → exact composition →
   `rna2*` impute → panel). Awaiting a direction call before wiring the front end.

---

## 10. Quick reference

- Run env: `/usr/local/anaconda3-2020/bin/python` on `bmiclusterp-head`.
- Data: `/data/salomonis-archive/LabFiles/Nicholas/AML-multimodal/` ; code: `…/pipeline/` ;
  outputs: `…/runs/`.
- LLM: `nemotron-3-super` @ `http://bmiclusterp2.chmcres.cchmc.org:4000/v1`.
- Engine: github.com/SalomonisLab/altanalyze3 (UDON in `components/udon/`).
- The honest signal = permutation p-value; the canonical metric env = cluster sklearn 1.5.2.

---

## 11. COMPASS-AML — the ex-vivo drug-response layer (added 2026-08-04)

A **parallel** layer to the mutation caller, not a downstream one: mutation calls enter as evidence,
never as a drug look-up. Trained on the BeatAML2 ex-vivo inhibitor screen
(`beataml_probit_curve_fits_v4_dbgap.txt`, fetched from the public BeatAML2 repo into
`data/external/beataml/`).

**Cohort.** 48,998 dose-response curves over **520 specimens / 479 patients / 118 inhibitors** (of 165
screened; 47 excluded for sample size, dynamic range, degenerate tails or acquisition-wave shift, and
3 — Cytarabine, Nutlin 3a, GDC-0941 — kept in a `wave_conditional` tier with the caveat attached).

**Three models, deliberately separate.**
- **A** patient-level response: hierarchical (target-pathway family + per-drug residual), four feature
  blocks (`rna`/`state`/`mut`/`clin`) late-fused by NNLS on inner donor-grouped OOF, floored to the
  best single block. Calibrated per drug on the OOF **percentile**, with cohort-matched score
  references (`beataml` / `sc_sample` / `sc_state`).
- **B** state-resolved: the same model applied per cell-state pseudobulk -> blast and LSC-like
  coverage, escape state, dispersion, bulk-vs-single-cell disagreement.
- **C** mechanistic: target expression, pathway output readout, BCL2-family dependency, genetic
  activation, measurable resistance proxies. Kept out of A so agreement/disagreement stays informative.

**Headline validation** (donor-grouped CV; 15% of *patients* sealed):
mean AUROC **0.774** over 118 inhibitors (**0.809** on the 42 approved agents); per-patient top-1
retrieval **34%** vs 10% matched chance; ECE **0.012**; leave-wave-out 0.722/0.733; leave-centre-out
0.731-0.890; survives inside every differentiation-state stratum; ~15 null SDs above a
specimen-repointing permutation null.

**Model B's falsifiable test.** Fitted only on BeatAML bulk and never shown a cell-state label, it
predicts higher venetoclax sensitivity in primitive than in monocytic states in **93.2%** of the 387
atlas samples (p = 8.4e-51) — **rank 1 of 118 inhibitors**, so not a generic artefact.

**Entry points.** `pipeline/predict_drugs.py` (one sample), `pipeline/drug_layer.py` (the hook
`ingest_patient.py` calls, non-fatal), `pipeline/batch_drug_reports.py` (backfill).
**GUI.** `/therapy.html` per-patient prioritisation (tiered, click a row for evidence),
`/rx_validation.html` the validation page. Full methods: `deliverables/METHODS_COMPASS-AML.md`.

**Standing caveat.** Ex-vivo sensitivity is an experimentally grounded prioritisation signal, not an
estimate of clinical benefit. Every surface says so.
