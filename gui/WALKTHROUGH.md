# MOSAIC-AML — end-to-end walkthrough: what happens at every step

This traces the **complete journey** of using the decision board: from double-clicking the
launcher, through uploading an scRNA sample, to reading the rendered patient report — naming the
exact file, function, input, and output at each step, with the real numbers from the validated
Grimes AML-CITE-Seq runs. For the broader engine (cohort pipeline, all 9 witnesses, Phases A–D) see
`../OVERVIEW.md`; this document is the **runtime path a user actually exercises**.

> **TL;DR of the spine.** `start_mosaic_board.bat` → SSH tunnel → `gui_server.py` on the cluster →
> browser board → **+ Add patient** → `POST /api/ingest` → `bsub` → `ingest_patient.py` on a compute
> node → (cell-state composition → cohort-trained subtype + optional genetic anchor → arbiter) →
> `patient_report.json` → polled by the browser → rendered as the decision board.

---

## 0. The cast (which machine does what)

| Actor | Role |
|-------|------|
| **Your laptop** (Windows) | Runs the launcher + browser. Holds the SSH tunnel. Renders the board. Does **no** computation. |
| **Cluster login node** (`bmiclusterp2`, reached via the `bmiclusterp-head` alias) | Runs the tiny `gui_server.py` (stdlib only). Receives uploads and **submits LSF jobs**. Does no heavy compute itself. |
| **LSF compute node** (e.g. `bmi-200m5-01`) | Runs `ingest_patient.py` — the actual analysis (composition, CV, arbiter). Has the analysis Python (`/usr/local/anaconda3-2020/bin/python`, 3.12, with anndata/scanpy/sklearn). |
| **Shared `/data` filesystem** | Holds the atlas, the cellHarmony reference, the imputation bundles, the pipeline code, and the `runs/` output dir. Visible to every node. |

Everything heavy runs where the data lives (the cluster). The laptop is a thin viewer + a remote control.

```
  ┌─ YOUR LAPTOP ───────────────┐        ┌─ CLUSTER login node (bmiclusterp2) ─────────────┐
  │ start_mosaic_board.bat      │        │  gui_server.py  (stock /usr/bin/python3, 3.6+)  │
  │   ├─ ssh: start server ──────┼───────▶│   GET /                -> mosaic_board.html     │
  │   ├─ browser  ◀──────────────┼────────│   GET /api/runs|report -> reads runs/*.json     │
  │   └─ ssh -L 8766:…:8766 ─────┼═tunnel═▶│   POST /api/ingest ─┐                           │
  │        (held open)          │        │   GET /api/jobs      │ bsub                       │
  └─────────────────────────────┘        │                      ▼                           │
                                          │             ┌─ LSF COMPUTE node ──────────────┐  │
                                          │             │ ingest_patient.py (py3.12)      │  │
                                          │             │  scRNA → composition → subtype  │  │
                                          │             │  + genetic anchor → arbiter     │  │
                                          │             │  → runs/<id>/patient_report.json│  │
                                          │             └─────────────────────────────────┘  │
                                          │  /data/salomonis-archive/.../AML-multimodal/      │
                                          └──────────────────────────────────────────────────┘
```

---

## PART A — Launching and opening the board

### Step 1 — You double-click `start_mosaic_board.bat`
**File:** `gui/start_mosaic_board.bat` (runs in `cmd.exe` on your laptop).
It sets four variables (`SSH_HOST=bmiclusterp-head`, `PORT=8766`, `GUI_DIR=…/AML-multimodal/gui`,
`URL=http://localhost:8766/`) and then runs three steps.

### Step 2 — Ensure the board server is running on the cluster
```
ssh bmiclusterp-head "pgrep -f '[g]ui_server.py 8766' >/dev/null 2>&1 && echo already-running \
   || (cd <GUI_DIR> && BROWSER=none nohup /usr/bin/python3 gui_server.py 8766 >/tmp/matrixgui.log 2>&1 </dev/null & sleep 2 && echo started)"
```
- **`pgrep -f '[g]ui_server.py 8766'`** checks whether the server is already up. The `[g]` is a
  deliberate trick: the regex `[g]ui_server.py 8766` matches the *server's* command line
  (`…gui_server.py 8766`) but **not the pgrep command's own** command line (which literally contains
  `[g]ui_server.py`), so the check can't match itself and falsely conclude "already running."
- If not running, it launches `gui_server.py` with the **stock `/usr/bin/python3`** under `nohup`
  (survives your SSH session closing) with output to `/tmp/matrixgui.log`.
- **Why stock python:** the server is **stdlib-only** and made **3.6-safe**, so it runs on the
  guaranteed-present `/usr/bin/python3` (3.6.8) on *any* login node. (`/usr/local/anaconda3-2020` is
  per-host and the login alias is a pool — it isn't on every node, which is why we don't depend on it
  for the server.)
- `if errorlevel 1 goto :ssherror` — if SSH itself can't connect, the batch prints a friendly error
  and stops, instead of leaving you with a dead tunnel.

### Step 3 — Schedule the browser to open
```
start "" /min powershell -NoProfile -WindowStyle Hidden -Command "Start-Sleep -Seconds 4; Start-Process 'http://localhost:8766/'"
```
A hidden helper waits 4 s (long enough for the tunnel in Step 4 to bind) and then opens your default
browser at `http://localhost:8766/`.

### Step 4 — Open the SSH tunnel (and keep it open)
```
ssh -N -o ExitOnForwardFailure=yes -o ServerAliveInterval=30 -L 8766:127.0.0.1:8766 bmiclusterp-head
```
- `-L 8766:127.0.0.1:8766` forwards **your laptop's** `localhost:8766` to **the cluster's**
  `127.0.0.1:8766` (where the server listens). So when the browser hits `localhost:8766`, the bytes
  travel encrypted over SSH to the cluster server.
- `-N` = no remote command, just the tunnel. This command **blocks** — that's why the window stays
  open. Closing the window (or Ctrl+C) tears down the tunnel.
- `ExitOnForwardFailure=yes` makes SSH bail immediately if local port 8766 is already taken (instead
  of connecting but silently not forwarding). `ServerAliveInterval=30` keeps the tunnel from idling out.
- **Pool caveat:** the server (Step 2) and the tunnel (Step 4) must land on the *same* login node.
  The alias usually lands consistently; if `localhost:8766` won't connect, run `hostname` in both and
  target the host explicitly.

> **Net effect of Part A:** a server is running on the cluster, an encrypted tunnel maps your
> laptop's 8766 to it, and your browser is open on the board. The cluster server keeps running after
> you close the window; only the tunnel goes away.

---

## PART B — The board loads in your browser

### Step 5 — Browser requests the page
`GET http://localhost:8766/` → through the tunnel → `gui_server.py` `do_GET` matches `/` and returns
`mosaic_board.html` (read fresh from disk each time, so edits need no server restart).
The HTML is one self-contained file: CSS + vanilla JS, no build step, no external assets.

### Step 6 — The page boots (`boot()` in `mosaic_board.html`)
1. `fetch('/api/runs')` — if it succeeds, the page is in **server mode**; if it fails (e.g. the file
   was opened directly with no server) it falls back to **standalone mode** (drag-drop a report).
2. `fetch('/api/capabilities')` → `{ingest, lsf, inbox, python}`. Because `bsub` exists on the
   cluster, this returns `ingest:true, lsf:true`, so the **"+ Add patient"** button is shown and
   `pollJobs()` starts.
3. The roster (left sidebar) is rendered from `/api/runs`, and the first/`?run=`-named report is
   auto-loaded into the board.

### Step 7 — What `/api/runs` returns
`gui_server.scan_runs()` walks `runs/`, and for **every subdirectory containing a
`patient_report.json`** reads it and emits one summary row: `{run, sample_key, annotation, dataset,
leading_hypothesis, overall_confidence, leading_confirmed_by_genetics, has_deliberation,
knowledge_version}`. Directories without that file (logs, diagnostics) are skipped, so the roster
shows only real patients.

---

## PART C — Uploading a patient

### Step 8 — You click "+ Add patient" and submit the form
Fields: **name**, **sample** (a cluster path to a 10x dir / `filtered_feature_bc_matrix.h5` /
`.h5ad`; or pick from the **inbox** dropdown populated by `GET /api/samples`), optional
**mutations**, and **dataset** label. Clicking **Process** runs `submitIngest()`, which does:
```
POST /api/ingest   body: {name, sample, mutations, dataset}
```

### Step 9 — The server dispatches the job (`dispatch_ingest()` in `gui_server.py`)
1. Derives a `run_id = "ingest_" + slug(name)` and writes an initial
   `runs/<run_id>/status.json = {state:"queued", step:"submitting"}` so the UI shows it instantly.
2. Builds the command **as an argv list (never a shell string)** — so the name/sample/mutations you
   typed are passed as literal arguments and **cannot inject shell commands**:
   ```
   [PYTHON, ingest_patient.py, --sample <sample>, --name <name>, --dataset <dataset>,
    --run-id <run_id>, --out-root <RUNS_DIR>, (--mutations <mutations> if given)]
   ```
3. Wraps it in `bsub -q normal -J aml_<run_id> -M 8000 -R rusage[mem=8000] -o …/_job.log …` and runs
   it. **`bsub` returns immediately** with a job id — it does **not** wait for the analysis.
4. **`PYTHON` here is the analysis python** (`/usr/local/anaconda3-2020/bin/python`), chosen because
   the job runs on a **compute node** where that 3.12 interpreter (with anndata/scanpy/sklearn) lives.
   This is deliberately **decoupled** from the stock python that runs the server.
5. Returns `{run_id, job_id, mode:"lsf"}` to the browser.

### Step 10 — The browser starts polling (`pollJobs()`)
Every 4 s it calls `GET /api/jobs` (which reads every `runs/*/status.json`) and renders a
**Processing** panel in the sidebar with the live `state` + `step`. When a job leaves the active set
(done/error), it refreshes the roster; if the finished job is the one you submitted, it auto-opens it.

---

## PART D — The ingest job runs on a compute node (`ingest_patient.py`)

LSF dispatches the job to a compute node. Each sub-step writes `status.json` so the sidebar updates
live. Here is **exactly** what the script computes.

### Step 11 — `status: running / loading_cohort` — build the context
`build_context()` (`amlmm/context.py`) loads the atlas tables (`composition`, `samples`, …) and the
versioned **knowledge base** (`amlmm/knowledge.py`, e.g. `kb-2026.07`). Provenance (python + library
versions + layout + the `deferred` witness list) is recorded for the report.

### Step 12 — `status: running / composition` (a) — read the scRNA sample
`load_query()` opens the sample:
- 10x `filtered_feature_bc_matrix.h5` (incl. **CITE-seq**) → `scanpy.read_10x_h5(gex_only=False)` then
  **keep only `feature_types == "Gene Expression"`** (the Antibody-Capture/ADT rows are dropped — for
  AML-7 that's 36,601 genes kept, 142 ADT dropped); or
- a gene-level `.h5ad` → read directly.

### Step 13 — `composition` (b) — assign every cell to a bone-marrow state
`amlmm/scrna.composition_from_query()` → `assign_cells()`:
1. Load the reference `Hs-MarrowAtlas-L3M.txt` — **2,870 marker genes × 89 populations**, where the
   89 populations are **exactly** the atlas's 89 cell-states (verified 89/89).
2. Intersect the query's gene symbols with the 2,870 markers (AML-7: **2,833 shared**; <50 shared
   aborts loudly — that's how junction-level data is rejected).
3. **CP10k + log1p normalize each cell over all genes**, then restrict to the shared markers. This is
   done **sparse-aware**: total counts are summed on the sparse matrix and only the ~2,833 marker
   columns are densified, so a 10k-cell × 36k-gene matrix uses ~300 MB instead of ~1.5 GB.
4. **Cosine similarity** of each cell vs each of the 89 population columns; `argmax` assigns the cell
   to its best-matching state. (Each cell is ≈ one state — this is why scRNA works where bulk
   deconvolution failed: no unmixing of collinear states.)
5. `composition()` = fraction of cells in each of the 89 states.
- **AML-7 result:** 8,495 cells, mean cosine 0.567, 58/89 states present, dominated by
  Intermediate-Mono / MDP-2 / MPP-MEP (monocyte/myeloid). **BF71-CD34** (CD34-sorted): 10,461 cells,
  cosine 0.619, 80/89 states, dominated by HSC/MPP/GMP/MEP/LMPP — a stem/progenitor profile that
  matches the sort (face-validity that the assignment is real).

### Step 14 — `composition` (c) — predict the subtype (cohort-trained)
The sample is **new** — it has no held-out cross-val prediction of its own. So:
1. `evidence_predictive(ctx, "subtype", "composition", …)` runs the **validated cohort CV**
   (`amlmm/cv.nested_cv_evaluate`, donor-grouped nested CV + permutation baseline) to get the
   **cohort balanced-accuracy** and **permutation p-value**, and it leaves the cohort feature matrix
   `X` (rows = atlas samples, columns = `comp::<state>`, L1-normalized) and labels `y` (canonical
   subtype) in the context. AML-7 cohort balanced-accuracy = **0.6751**.
2. Fit a fresh RandomForest (`models.build(["rf"])`) on the whole cohort `(X, y)`.
3. Build the new sample's feature vector by aligning its 89-state composition to `X`'s `comp::<state>`
   columns (by name; missing → 0) and L1-normalizing — i.e. exactly the cohort's feature space.
4. `predict` + `predict_proba` → subtype + probability. **AML-7 → FLT3, prob 0.4125.**
- The evidence records `held_out:false` and the opinion's caveat says *"external sample: cohort-trained
  classifier, NOT held-out for this sample; out-of-cohort domain shift possible."*
- Grounding/independence: `honest_cv` / `independent`.

### Step 15 — `status: running / genetic` — the genetic witness
`normalize_mutations()` maps what you typed to the arbiter's flag vocabulary (`ANCHOR_MAP`): gene
names (`TP53`, `FLT3`, `NPM1`, …) and cytogenetics/synonyms (`inv16`, `del7`, `t15_17`/APL, `kmt2a`,
`complex`, …). It builds `present` (the observed drivers) and `targetable` (drivers with a KB
therapy). **If you supplied no mutations, `present = []`** and the witness is recorded with low weight
— the run is **anchor-free**. Grounding/independence: `deterministic_fact` / `independent`.

### Step 16 — `status: running / arbiter` — reconcile (`amlmm/arbiter.reconcile_patient`)
Two `AgentResult`s (composition, genetic) are appended to a `Ledger` (`amlmm/ledger.py`); evidence is
recorded **immutably**. Then the **deterministic pre-pass** runs — the LLM does **not** make the
decision:
1. **Effective weight** per witness = `reliability_weight × GROUNDING_FACTOR × INDEPENDENCE_FACTOR`:
   - `GROUNDING_FACTOR = {deterministic_fact:1.0, honest_cv:0.9, classifier_call:0.7, descriptive_aggregate:0.5}`
   - `INDEPENDENCE_FACTOR = {independent:1.0, rna_derived:0.6, imputed_from_RNA:0.5, discovery:0.7}`
   - So composition weighs `0.85×0.9×1.0 = 0.765`; an observed-driver genetic witness weighs
     `0.85×1.0×1.0 = 0.85`.
2. **Genetic anchor:** if a driver is *present*, the first one by `ANCHOR_PRIORITY`
   (`t15_17 > TP53 > FLT3 > NPM1 > …`) becomes the **leading hypothesis** — **by rule, not by weight.**
   A disagreeing prediction is recorded as a *conflict* but **cannot outrank** it (this is the
   "TP53 fix"). With no driver supplied, the leading hypothesis is the composition prediction.
3. **Concordance** = leading's summed weight ÷ total signal weight. (With only one voting witness, as
   in an anchor-free upload, concordance is trivially 1.00 — so it is *not* taken as strength; the
   real uncertainty is conveyed by the "confirm by sequencing" framing + the cohort accuracy on the
   card.)
4. **Therapies** come *only* from the KB keyed on **observed** drivers — never on a predicted subtype.
   So an anchor-free FLT3 prediction yields **no therapy** (FLT3 wasn't observed); instead it yields a
   **validation** recommendation (sequencing). Validations are drawn from KB rules for the present
   claim types.
5. **Confidence:** `anchor & concordance≥0.6 → high`; `anchor or concordance≥0.5 → medium`; else `low`.
6. **LLM narration (optional):** if the gateway is reachable, the LLM writes the prose rationale and
   may set confidence — but **only within constraints**: it cannot change the leading hypothesis or
   therapies, and (honesty guard) it **cannot raise confidence above the deterministic ceiling for a
   call that is not genetically confirmed**. On any LLM error, the deterministic result stands.
- **AML-7 (anchor-free):** leading **FLT3**, `leading_confirmed_by_genetics:false`, confidence
  **medium**, 0 therapies, validations = NGS panel + flow. Supplying `TP53` instead →
  `leading_confirmed_by_genetics:true`, **high**.

### Step 17 — `status: done` — write outputs
The script writes into `runs/<run_id>/`:
- **`patient_report.json`** — the full record (mode, sample_key, annotation, dataset, provenance,
  `panel[]` of witnesses, `consensus{}`, `deliberation` (null for uploads), and an `ingest{}` block
  with the source path, composition quality, supplied mutations, and the **`deferred` witness list**).
- **`ledger.json`** — the immutable evidence ledger + `deterministic_evidence_hash` (reproducibility).
- **`PATIENT.md`** — a human-readable summary.
- **`cv_result.json`** — the cohort CV detail produced in Step 14.
- **`status.json`** — flipped to `{state:"done"}` last, so the UI only refreshes once the report exists.

> **Honest-scope note:** an upload runs **two** witnesses today — composition + genetic. The six
> imputed/descriptive witnesses (GRN/ADT/metabolic/lipid/LSC/UDON/cell-comm) are listed under
> `ingest.deferred` and are simply absent from the panel; the arbiter handles a partial roster. (The
> full 9-witness roster + Phase C deliberation run on the in-cohort `panel.py --patient` path.)

---

## PART E — The report appears and renders

### Step 18 — Polling notices completion
On the next `pollJobs()` tick, the job's `status.json` reads `done`. The front-end calls
`refreshRoster()` (re-fetches `/api/runs`, so the new patient appears in the sidebar) and, since it
was the one you submitted, auto-loads it. (A short delayed re-refresh absorbs any shared-filesystem
lag between the status flip and the report becoming visible.)

### Step 19 — Fetch + render the report
`loadRun(run_id)` → `GET /api/report?run=<run_id>` returns the `patient_report.json` (with the ledger
`deterministic_evidence_hash` merged in). `renderReport()` then draws the board:
- **Hero:** the leading hypothesis (big), the anchor line — green *"⚓ Genetically anchored — observed
  driver"* if confirmed, else amber *"◷ Witness consensus — confirm by sequencing"* — a confidence
  badge, and a concordance gauge.
- **Conflict banner:** green if no conflicts; otherwise it names them and, when there's an anchor,
  explains the disagreeing prediction was *"overruled by the observed driver … (the anchor cannot be
  outvoted)."*
- **Therapies / Validations:** two columns; therapies show biomarker → drug + evidence level + source;
  surface-marker hypotheses (if any) are flagged flow-pending.
- **Evidence panel:** one card per witness with **`VOTES`/`CONTEXT`** tags, the **grounding ×
  independence badges showing the ×factor**, and a literal **effective-weight bar** computed in JS to
  match the arbiter's `_eff_weight` — so imputed witnesses visibly weigh less. Each card expands to
  the raw evidence (prediction + probability + cohort accuracy, drivers, etc.).
- **Deliberation** (only when present) and **descriptive context** (non-voting) sections.
- **Footer:** the decision-support disclaimer + the KB version + the evidence hash.

---

## Appendix 1 — Files written per uploaded patient (`runs/<run_id>/`)
| File | Written by | Contents |
|------|-----------|----------|
| `status.json` | `gui_server` then `ingest_patient` | live state: queued → running(loading_cohort/composition/genetic/arbiter) → done/error |
| `patient_report.json` | `ingest_patient` | the full decision record the GUI renders |
| `ledger.json` | `ingest_patient` | immutable evidence ledger + `deterministic_evidence_hash` |
| `PATIENT.md` | `ingest_patient` | human-readable summary |
| `cv_result.json` | the `classify` step | cohort nested-CV detail (accuracy, permutation p, OOF) |
| `_job.log` | LSF (`bsub -o`) | the job's stdout/stderr |

## Appendix 2 — Weighting reference (why imputed counts for less)
`effective_weight = reliability_weight × GROUNDING_FACTOR × INDEPENDENCE_FACTOR`

| grounding | × | | independence | × |
|-----------|---|---|--------------|---|
| deterministic_fact | 1.0 | | independent | 1.0 |
| honest_cv | 0.9 | | discovery | 0.7 |
| classifier_call | 0.7 | | rna_derived | 0.6 |
| descriptive_aggregate | 0.5 | | imputed_from_RNA | 0.5 |

## Appendix 3 — The honest-design invariants threaded through every step
1. **The decision is deterministic.** A reproducible pre-pass decides; the LLM only narrates within
   fixed conclusions and can be removed (`--no-llm`) without changing the call.
2. **Observed beats predicted, by rule.** A present driver anchors the leading hypothesis; a
   conflicting imputed/predicted call is logged but cannot outrank it.
3. **Imputed evidence is down-weighted and never prescribes.** Therapies key only on observed drivers;
   imputed signals become hypotheses to confirm.
4. **Confidence cannot be inflated.** The LLM may lower confidence but not raise it above the
   deterministic ceiling for an unconfirmed call.
5. **Scope is stated, not hidden.** Deferred witnesses, out-of-cohort prediction, and
   single-witness concordance are all surfaced rather than papered over.

## Appendix 4 — Stopping things
- **Disconnect (laptop):** close the launcher window (ends the tunnel). The cluster server keeps running.
- **Stop the cluster server:** `ssh bmiclusterp-head "pkill -f 'gui_server.py 8766'"`.
- **A stuck upload** shows as `state:"error"` in the Processing panel with the message; the `_job.log`
  in that run dir has the traceback.
