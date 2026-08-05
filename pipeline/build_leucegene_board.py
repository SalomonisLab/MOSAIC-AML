#!/usr/bin/env python3
"""Add the Leucegene bulk cohort (n=367) to the decision board as external-validation reports.

Leucegene is bulk RNA the caller was never trained on. Each sample is scored with the bulk caller
(ref='leucegene'), gets deterministic therapy/tests panels, and GENE-LEVEL validation against the
dbGAP variant calls (bundle lg_L). No composition/cytogenetics (bulk input, no single cells).

Reports land in runs/leucegene_<sample>/ (mode 'bulk_panel'); gui_server groups the 'leucegene_' prefix.
Run on an LSF compute node (needs numpy):
  bsub -q test -W 60 -M 16000 -R "rusage[mem=16000]" -o lg.log \
    /usr/local/anaconda3-2020/bin/python build_leucegene_board.py
"""
import os, sys, json, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

from ingest_patient import bulk_mutation_result
from amlmm import therapy

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
RUNS = os.environ.get("AMLMM_RUNS") or os.path.join(PROJ, "runs")
BUNDLE = os.path.join(os.path.dirname(PROJ), "aml-bakeoff", "bundle_data.npz")

CYTO = {"inv16", "complex", "del5", "del7", "trisomy8", "kmt2a"}


def unit(cat):
    cl = str(cat).lower()
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl:
        return "inv16"
    if "kmt2a" in cl:
        return "KMT2A"
    if cl in CYTO:
        return cl
    return str(cat).split("_")[0].split("-")[0].upper()


def main():
    d = np.load(BUNDLE, allow_pickle=True)
    genes = [str(g) for g in d["genes"]]
    cats = [str(c) for c in d["drivers"]]
    X = pd.DataFrame(d["lg_X"].astype(float), index=[str(s) for s in d["lg_samples"]], columns=genes)
    L = pd.DataFrame(d["lg_L"].astype(float), index=X.index, columns=cats)

    n = 0
    dist = {}
    for s in X.index:
        preds, caller, _ = bulk_mutation_result(X.loc[s], [], ref="leucegene")
        if preds is None:
            continue
        # per-category truth from dbGAP calls (category-level) -> also drives gene-level validation
        cat_truth = {}
        for c in cats:
            v = L.loc[s, c]
            if pd.notna(v):
                cat_truth[c] = "present" if v == 1 else "absent"
        for p in preds:
            p["true_label"] = cat_truth.get(p["mutation"])       # detail-view checkmark where dbGAP labels it
        # gene-level validation: aggregate model calls + truth to the gene unit
        gt, call = {}, {}
        for c in cats:
            u = unit(c)
            if c in cat_truth:
                # present wins; only mark absent if not already present
                if cat_truth[c] == "present":
                    gt[u] = "present"
                else:
                    gt.setdefault(u, "absent")
        for p in preds:
            u = unit(p["mutation"])
            if p.get("call") == "present":
                call[u] = "present"
            else:
                call.setdefault(u, "absent")
        units = {u: {"truth": t, "call": call.get(u, "absent"), "correct": (call.get(u, "absent") == t)}
                 for u, t in sorted(gt.items())}
        panels = therapy.build_panels(preds, [])                 # deterministic (no agent on 367 samples)

        rep = {
            "mode": "bulk_panel", "sample_key": "Leucegene::" + s, "dataset": "Leucegene",
            "specimen_class": None, "control_gate": None, "validation": True,
            "panel": [], "mutation_predictions": preds, "mutation_caller": caller,
            "treatment_panel": panels["treatments"], "tests_panel": panels["tests"],
            "panels_note": panels["note"],
            "validation_gene": {"n_labeled": len(units),
                                "n_correct": sum(1 for v in units.values() if v["correct"]), "units": units},
            "consensus": {"leading_hypothesis": "external validation (Leucegene bulk RNA)",
                          "overall_confidence": "mutation panel only"},
            "ingest": {"input_kind": "bulk_rna", "source": "Leucegene bundle", "name": s,
                       "bulk_ref": "leucegene", "note": "external bulk cohort, dbGAP variant truth"},
        }
        outd = os.path.join(RUNS, "leucegene_" + "".join(ch if ch.isalnum() else "_" for ch in s))
        os.makedirs(outd, exist_ok=True)
        json.dump(rep, open(os.path.join(outd, "patient_report.json"), "w"), default=str, indent=1)
        n += 1
        dist[len(units)] = dist.get(len(units), 0) + 1

    print("wrote %d Leucegene board reports" % n)
    print("gene-validation-units distribution:", dict(sorted(dist.items())))
    # cohort-level gene accuracy
    tot = cor = 0
    for run in os.listdir(RUNS):
        if not run.startswith("leucegene_"):
            continue
        v = json.load(open(os.path.join(RUNS, run, "patient_report.json"))).get("validation_gene") or {}
        tot += v.get("n_labeled", 0)
        cor += v.get("n_correct", 0)
    print("cohort gene-level accuracy: %d/%d = %.1f%%" % (cor, tot, 100.0 * cor / max(tot, 1)))
    print("LEUCEGENE BOARD OK")


if __name__ == "__main__":
    main()
