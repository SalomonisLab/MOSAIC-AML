# MOSAIC-AML decision board (GUI)

A zero-dependency viewer **and uploader** for the per-patient reports the engine writes
(`runs/<run>/patient_report.json`). It renders the tumor-board readout, and it can **ingest a new
scRNA sample as a patient** — composition → cohort-trained subtype + (optional) genetic anchor →
arbiter — dispatching the work as an LSF job on the cluster.

Viewing is read-only; ingest runs the pipeline's own `ingest_patient.py` (argv list, no shell).

## Run it

### On the cluster (recommended — enables upload)
The atlas, the cellHarmony reference, and the compute all live on the cluster, so run the server on
a login node and reach it from your laptop over an SSH tunnel:

```bash
# on the cluster (login node) — the stock python3 is fine; the server is STDLIB-ONLY (3.6-safe):
python3 /data/.../AML-multimodal/gui/gui_server.py 8765
# on your laptop:
ssh -L 8765:127.0.0.1:8765 bmiclusterp-head
# then open http://127.0.0.1:8765/
```
It auto-detects `bsub` and submits ingest jobs to LSF (queue `normal`, 8 GB). The ingest *job*
uses the analysis python (`/usr/local/anaconda3-2020/bin/python`, present on compute nodes) — the
server interpreter is independent of it. Override the job python with `AMLMM_PYTHON` if needed.

### Locally (view only)
```bash
python AML-multimodal/gui/gui_server.py        # scans ../runs, serves http://127.0.0.1:8765
```
No `bsub` locally, and the reference lives on the cluster — so **upload is a cluster feature**.

## Add a patient (upload + process)
Click **+ Add patient** and fill in:
- **Patient name** — e.g. `AML-7`.
- **Sample** — a cluster path to a 10x directory, a `filtered_feature_bc_matrix.h5` (CITE-seq is
  fine — the Antibody-Capture features are dropped automatically), or a gene-level `.h5ad`. Drop
  files into the **inbox** (default `AML-multimodal/inbox/`, override via `AMLMM_INBOX`) and pick
  from the dropdown, or paste any path you can read.
- **Observed mutations** (optional) — e.g. `TP53, FLT3`. Supplying a driver engages the
  deterministic **genetic anchor**; leaving it blank runs **anchor-free** (the call is a hypothesis
  to confirm by sequencing).
- **Dataset label** — free text.

Click **Process** → a job goes to LSF; the sidebar shows live status (queued → composition →
arbiter → done). When it finishes the report joins the roster and opens automatically. Each ingest
writes `runs/ingest_<name>/{patient_report.json, ledger.json, PATIENT.md, status.json}`.

**What runs today (v1):** the **composition** witness (89-state cosine assignment + a cohort-trained
subtype classifier, applied out-of-cohort) and the **genetic** witness (your mutations). The imputed
descriptive witnesses (GRN/ADT/metabolic/lipid/LSC/UDON/cell-comm) are not yet wired for uploads and
are listed under `ingest.deferred` in the report.

## Endpoints
| route | returns |
|-------|---------|
| `GET /` | the board HTML |
| `GET /api/runs` | roster summary rows |
| `GET /api/report?run=NAME` | the full report JSON (+ ledger evidence hash) |
| `GET /api/capabilities` | `{ingest, lsf, inbox, python}` |
| `GET /api/jobs` | in-flight + finished ingests (from each run's `status.json`) |
| `GET /api/samples` | candidate inputs in the inbox dir |
| `POST /api/ingest` | `{name, sample, mutations?, dataset?}` → dispatch a new ingest, returns `{run_id, job_id, mode}` |

## Files
- `mosaic_board.html` — the entire front-end (HTML + CSS + JS, no build step). Also works
  standalone: open the file and drag a `patient_report.json` onto it (view only).
- `gui_server.py` — discovery / serve + ingest dispatch (stdlib `http.server`).
- `start_mosaic_board.bat` — Windows double-click launcher: starts the server on the cluster, opens
  the SSH tunnel, and pops your browser.
- `WALKTHROUGH.md` — **step-by-step of what happens at every stage** (launch → tunnel → upload →
  `ingest_patient.py` → arbiter → render), with the real functions, data shapes, and numbers.

## How to read the board
The decision is a **deterministic pre-pass**, not the LLM. An **observed driver mutation anchors**
the leading hypothesis; a disagreeing *predicted* subtype is recorded as a conflict but cannot
override it. Each witness's weight = `reliability × grounding × independence` (shown on every badge),
so imputed-from-RNA witnesses visibly count for less. The LLM narrates within the fixed decision and
may lower confidence but **never raise it above the deterministic ceiling for an unconfirmed call**.
