# Setup — MOSAIC-AML

Getting a functional install: the engine, the Python environment, the (separately distributed) data, and how to run.

## 1. Clone with the engine submodule
The AltAnalyze3 engine (RNA→ADT/GRN/Lipid/Metabolite imputers, cellHarmony cell-state assignment, fastComm cell-communication) is pinned as a submodule at `engine-code/altanalyze3`.

```bash
git clone --recurse-submodules https://github.com/SalomonisLab/MOSAIC-AML.git
# or, if already cloned:
git submodule update --init --recursive
```

## 2. Python environment
```bash
pip install -r pipeline/requirements-lock.txt
```
Core deps: numpy / scipy / scikit-learn / pandas / anndata / h5py. The full-cohort steps are designed to run on an LSF compute node.

## 3. Data (not in git — large `.h5ad` matrices)
The expression matrices are distributed separately (too large for git). Place them under `data/`:

| File | Modality |
|---|---|
| `data/RNA/pseudobulk_counts_hashed.h5ad` | RNA (per sample × cell-state pseudobulks) |
| `data/ADT/pseudobulk_adt_imputed.h5ad` | ADT (surface protein) |
| `data/GRN/imputed_grn_all_pseudobulks.h5ad` | GRN (regulon activity) |
| `data/Lipid/pseudobulk_imputed_lipid_aml.h5ad` | Lipid |
| `data/Metabolite/pseudobulk_imputed_metabolite_aml.h5ad` | Metabolite |
| `data/cell-communication/combined_sample_by_interaction.h5ad` | Cell-communication |

Cell-state / genotype labels and UDON program assignments are already in `labels/`. The per-patient prediction reports are in `runs/`. For the `.h5ad` files, see the data release / contact the Salomonis lab.

## 4. Run
```bash
# Decision board (browse per-patient calls, evidence, add a patient)
python gui/gui_server.py runs 8765          # -> http://127.0.0.1:8765

# Train / deploy the multimodal predictor (LSF compute node)
python pipeline/train_predictor.py

# Ingest a new scRNA-seq sample -> full multimodal mutation panel
python pipeline/ingest_patient.py "PATIENT=/path/to/sample.h5"
```

## Layout
- `pipeline/` — `amlmm` package: context, dataio, cell-state assignment (`scrna.py`), genetics anchor, per-(mutation × modality) classifiers + late fusion (`predictor.py`), multi-agent arbiter/witnesses, control gate, per-patient ingest + reporting.
- `gui/` — decision board (`gui_server.py` + `mosaic_board.html`) and the CEBPA evidence drill-down.
- `scripts/` — per-modality algorithms (LSC RandomForest, fastComm, imputers).
- `engine-code/altanalyze3` — the AltAnalyze3 engine (submodule).

## The two prediction layers

- **CIPHER-AML** — the mutation predictor. Train `pipeline/train_predictor.py` (multimodal) and
  `pipeline/train_bulk_predictor.py` (bulk variant-level, the deployed primary caller).
  Validation pages: `gui/validation.html`, `gui/calibration.html`, `gui/evidence.html`.
- **COMPASS-AML** — the drug predictor. Needs the BeatAML2 inhibitor screen in
  `data/external/beataml/` (fetch from https://github.com/biodev/beataml2.0_data), then:

  ```
  python pipeline/build_state_signatures.py     # lineage signatures from the atlas
  python pipeline/train_drug_model.py           # Model A  -> drug_response_model.pkl
  python pipeline/build_drug_score_refs.py      # cohort-matched score references (required)
  python pipeline/eval_drug_model.py            # the validation battery
  ```

  Validation page: `gui/rx_validation.html`; per-patient view: `gui/therapy.html`.

## Try it without any of the above

`python pipeline/build_static_site.py` freezes the whole board and every validation page into
`deliverables/mosaic_static/` — open `index.html` in a browser, no server needed.
