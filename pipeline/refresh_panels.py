#!/usr/bin/env python3
"""Refresh EVERY patient_report.json with the therapy + tests panels (idempotent, no retrain).

The panels are a pure function of a report's own `mutation_predictions`, so every existing report can be
brought up to the current format without re-running the (expensive) predictor. Reports that predate the
mutation caller entirely (no predictions) are reported as stale rather than silently touched.

  python refresh_panels.py [runs_dir]      # default: ../runs
"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amlmm import therapy

RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "runs")

done, stale, failed = [], [], []
for f in sorted(glob.glob(os.path.join(RUNS, "*", "patient_report.json"))):
    run = os.path.basename(os.path.dirname(f))
    try:
        with open(f, encoding="utf-8") as fh:
            r = json.load(fh)
    except Exception as e:
        failed.append((run, "unreadable: %s" % e)); continue
    preds = r.get("mutation_predictions")
    if not isinstance(preds, list) or not preds:
        stale.append((run, r.get("mode")))                 # predates the mutation caller
        continue
    # what counts as SEQUENCED for this report: a known label on the sealed held-out board reports,
    # plus anything the user supplied at ingest.
    supplied = [p.get("mutation") for p in preds if p.get("true_label") == "present"]
    supplied += list((r.get("ingest") or {}).get("mutations_supplied") or [])
    try:
        panels = therapy.build_panels(preds, supplied)
    except Exception as e:
        failed.append((run, "panel build: %s" % e)); continue
    r["treatment_panel"] = panels["treatments"]
    r["tests_panel"] = panels["tests"]
    r["panels_note"] = panels["note"]
    with open(f, "w", encoding="utf-8") as fh:
        json.dump(r, fh, default=str, indent=1)
    done.append((run, len(panels["treatments"]), len(panels["tests"])))

print("refreshed %d report(s):" % len(done))
for run, nt, nx in done[:60]:
    print("   %-46s therapy=%-2d tests=%-2d" % (run, nt, nx))
if stale:
    print("\n%d STALE report(s) with no mutation_predictions (predate the caller — cannot refresh):" % len(stale))
    for run, mode in stale:
        print("   %-46s mode=%s" % (run, mode))
if failed:
    print("\n%d FAILED:" % len(failed))
    for run, why in failed:
        print("   %-46s %s" % (run, why))
print("\nREFRESH OK")
