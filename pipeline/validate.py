#!/usr/bin/env python3
"""MOSAIC-AML retrospective clinical validation CLI (Phase D).

Tests whether the engine's deterministic cohort outputs (the genetic-anchored driver + its
ELN-expected risk, the LSC stemness call) ASSOCIATE with the sparse clinical labels actually
present (ELN_risk, clinical_response, overall_survival). Honest, underpowered-aware.

  python validate.py
Outputs validation_report.json + VALIDATION.md under runs/<run_id>/.
"""
from __future__ import annotations
import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amlmm
from amlmm.retrospective import run_validation


def main(argv=None):
    ap = argparse.ArgumentParser(description="MOSAIC-AML retrospective clinical validation")
    ap.add_argument("--base", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--run-id", default="validation")
    args = ap.parse_args(argv)

    d = amlmm.Config()
    cfg = amlmm.Config(base_dir=args.base or d.base_dir, out_dir=args.out or d.out_dir,
                       run_id=args.run_id)
    ctx = amlmm.build_context(cfg)
    print(f"[validate] layout={ctx.layout}")
    rep = run_validation(ctx)

    cov, ts = rep["coverage"], rep["tests"]
    print(f"\ncoverage: anchored_driver={cov['n_anchored_driver']} eln={cov['n_eln']} "
          f"responder={cov['n_responder']} survival={cov['n_survival']} lsc={cov['n_lsc_confident']}")
    ec = ts["anchored_eln_concordance"]
    print(f"ELN concordance: agreement={ec.get('agreement')} (n={ec.get('n')}, perm p={ec.get('p_perm')})"
          f"{' UNDERPOWERED' if ec.get('underpowered') else ''}")
    for key, label in [("adverse_driver_vs_eln_adverse", "adverse-driver->ELN-adverse"),
                       ("favorable_driver_vs_responder", "favorable-driver->responder"),
                       ("pLSC_vs_eln_adverse", "p-LSC->ELN-adverse")]:
        r = ts[key]
        print(f"{label}: OR={r.get('odds_ratio')} (n={r.get('n')}, Fisher p={r.get('p')})"
              f"{' UNDERPOWERED' if r.get('underpowered') else ''}")
    print(f"\noutputs: {ctx.run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
