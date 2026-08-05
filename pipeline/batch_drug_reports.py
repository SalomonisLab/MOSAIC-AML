#!/usr/bin/env python3
"""Generate COMPASS-AML drug reports for runs the platform has already ingested.

Only samples that exist in the atlas pseudobulk file can get the state-resolved layer, so those are
what this backfills; anything else is skipped with a reason rather than silently produced from a
degraded input.

  python batch_drug_reports.py [--limit N] [--force]
"""
import os, sys, json, time, glob, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, h5py

import predict_drugs as PD
from amlmm.drug.h5rows import obs_column

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "runs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    t0 = time.time()

    with h5py.File(PD.H5, "r") as f:
        keys = set("%s::%s" % (x, y) for x, y in zip(obs_column(f, "Dataset"), obs_column(f, "Sample")))

    todo, skipped = [], {}
    for p in sorted(glob.glob(os.path.join(RUNS, "*", "patient_report.json"))):
        run = os.path.basename(os.path.dirname(p))
        out = os.path.join(os.path.dirname(p), "drug_report.json")
        if os.path.exists(out) and not a.force:
            skipped.setdefault("already_done", 0)
            skipped["already_done"] += 1
            continue
        try:
            sk = json.load(open(p)).get("sample_key")
        except Exception:
            sk = None
        if sk in keys:
            todo.append((run, sk, os.path.dirname(p)))
        else:
            skipped.setdefault("not_in_atlas", 0)
            skipped["not_in_atlas"] += 1
    if a.limit:
        todo = todo[:a.limit]
    print("to generate: %d | skipped: %s" % (len(todo), skipped))

    ok = 0
    for i, (run, sk, d) in enumerate(todo, 1):
        try:
            counts, ncells = PD.atlas_sample(sk)
            rep = PD.run(state_counts=counts, n_cells=ncells)
            rep["input"] = {"kind": "single-cell (atlas pseudobulk)", "sample_key": sk, "run": run}
            json.dump(rep, open(os.path.join(d, "drug_report.json"), "w"), indent=1, default=str)
            open(os.path.join(d, "DRUG_REPORT.md"), "w", encoding="utf-8").write(PD.markdown(rep))
            # mirror the compact summary into patient_report.json so the decision board can show the
            # top candidates in context without fetching the full drug report
            import drug_layer
            pr = os.path.join(d, "patient_report.json")
            j = json.load(open(pr))
            j["drug_response"] = drug_layer._summary(rep)
            json.dump(j, open(pr, "w"), indent=1, default=str)
            ok += 1
        except Exception as e:
            print("  !! %s: %s: %s" % (run, type(e).__name__, e))
        if i % 10 == 0 or i == len(todo):
            print("   %d/%d (%.0fs)" % (i, len(todo), time.time() - t0))
    print("wrote %d drug reports (%.0fs)" % (ok, time.time() - t0))


if __name__ == "__main__":
    main()
