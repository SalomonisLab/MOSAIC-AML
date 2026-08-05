#!/usr/bin/env python3
"""The COMPASS-AML hook the patient pipeline calls: expression in, drug_report.json out.

Kept in its own module and behind a single try/except at the call site, because the drug layer is an
addition to an existing clinical-decision pipeline: a missing model file, an unbuilt score reference or
an odd input must degrade to "no drug report" and never take down the mutation panel.

  run_for_adata(adata, run_dir, ...)   single-cell -> per-cell-state pseudobulks -> Models A + B + C
  run_for_bulk(series, run_dir, ...)   one bulk expression vector -> Models A + C (no state layer)
"""
from __future__ import annotations
import os, json, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
MODEL = os.path.join(HERE, "drug_response_model.pkl")


def available():
    return os.path.exists(MODEL)


def _state_pseudobulks(adata):
    """Assign every cell to a reference marrow state, then sum RAW counts per state.

    Re-runs the cosine assignment rather than threading labels through the caller: it is a few seconds,
    and it keeps this hook a pure function of the AnnData so the drug layer can be run standalone on
    any sample the pipeline has already ingested."""
    import scipy.sparse as sp
    from amlmm import scrna
    ref = scrna.load_reference()
    labels, scores, n_shared = scrna.assign_cells(adata, ref)
    labels = np.asarray(labels)
    genes = [str(g) for g in adata.var_names]
    X = adata.X
    rows, ncells = {}, {}
    for st in pd.unique(labels):
        m = labels == st
        if m.sum() < 20:                      # too few cells for a usable pseudobulk
            continue
        sub = X[m]
        v = np.asarray(sub.sum(axis=0)).ravel() if sp.issparse(sub) else np.asarray(sub).sum(axis=0).ravel()
        rows[str(st)] = v
        ncells[str(st)] = int(m.sum())
    if not rows:
        return None, None, {"n_shared_markers": int(n_shared), "mean_cosine": float(np.mean(scores))}
    return (pd.DataFrame(rows, index=genes).T, pd.Series(ncells),
            {"n_shared_markers": int(n_shared), "mean_cosine": float(np.mean(scores)),
             "n_states": len(rows), "n_cells_assigned": int(sum(ncells.values()))})


def _mutations_for_drug_layer(mut_preds, observed):
    """Feed the mutation layer's own calls into Model C as EVIDENCE, tagged by provenance.

    Observed genotypes and classifier predictions are both accepted but never conflated: the mechanism
    model records which is which, so a report can say 'pathway-activating lesion (predicted)' rather
    than implying a sequencing result."""
    out = {}
    for p in (mut_preds or []):
        cat = p.get("category") or p.get("mutation")
        if not cat:
            continue
        out[str(cat)] = float(p.get("probability") or 0.0)
    for m in (observed or []):
        out[str(m)] = 1.0
    out["__observed__"] = [str(m) for m in (observed or [])]
    return out or None


def _write(rep, run_dir):
    import predict_drugs as PD
    os.makedirs(run_dir, exist_ok=True)
    j = os.path.join(run_dir, "drug_report.json")
    json.dump(rep, open(j, "w"), indent=1, default=str)
    open(os.path.join(run_dir, "DRUG_REPORT.md"), "w", encoding="utf-8").write(PD.markdown(rep))
    return j


def _summary(rep):
    """The compact block that goes into patient_report.json, so the decision board can show the top
    candidates without loading the full drug report."""
    top = {}
    for tier, blk in (rep.get("ranked") or {}).items():
        top[tier] = [{"inhibitor": r["inhibitor"], "utility": r["utility"],
                      "prob_sensitive": r["components"]["sensitivity"],
                      "n_challenges": len(r.get("challenges") or [])}
                     for r in blk["ranked"][:5]]
    return {"available": True, "n_drugs_modelled": rep.get("n_drugs_modelled"),
            "n_drugs_reported": rep.get("n_drugs_reported"), "n_abstained": rep.get("n_abstained"),
            "patient": rep.get("patient"), "top_by_tier": top,
            "caveat": ("predicted EX-VIVO sensitivity from the BeatAML2 functional screen; a "
                       "prioritisation for trial matching or laboratory validation, not a treatment "
                       "recommendation")}


def run_for_adata(adata, run_dir, mut_preds=None, observed=None, clinical=None):
    import predict_drugs as PD
    if not available():
        return {"available": False, "reason": "drug_response_model.pkl not built "
                                              "(run train_drug_model.py + build_drug_score_refs.py)"}
    counts, ncells, qc = _state_pseudobulks(adata)
    if counts is None:
        return {"available": False, "reason": "no cell state reached the 20-cell floor", "qc": qc}
    rep = PD.run(state_counts=counts, n_cells=ncells,
                 mutations=_mutations_for_drug_layer(mut_preds, observed), clinical=clinical)
    rep["input"] = {"kind": "single-cell", **qc}
    _write(rep, run_dir)
    return _summary(rep)


def run_for_bulk(series, run_dir, mut_preds=None, observed=None, clinical=None):
    import predict_drugs as PD
    if not available():
        return {"available": False, "reason": "drug_response_model.pkl not built"}
    rep = PD.run(bulk=series, mutations=_mutations_for_drug_layer(mut_preds, observed),
                 clinical=clinical)
    rep["input"] = {"kind": "bulk", "n_genes": int(len(series))}
    _write(rep, run_dir)
    s = _summary(rep)
    s["note"] = ("bulk input: cell-state coverage and escape-clone reasoning are unavailable, so the "
                 "utility score is computed from the terms that could be evaluated")
    return s


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Run the COMPASS-AML drug layer for one already-ingested run.")
    ap.add_argument("--sample", required=True, help="10x dir / .h5 / .h5ad, or a bulk expression table")
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--bulk", action="store_true")
    a = ap.parse_args()
    if a.bulk:
        import ingest_patient as IP
        ser, col, scale = IP.parse_bulk_expression(a.sample)
        print(json.dumps(run_for_bulk(ser, a.run_dir), indent=1, default=str))
    else:
        import ingest_patient as IP
        adata, how = IP.load_query(a.sample)
        print(json.dumps(run_for_adata(adata, a.run_dir), indent=1, default=str))
