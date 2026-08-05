#!/usr/bin/env python3
"""MOSAIC-AML Discovery agent CLI (three-agent rebuild — Discovery first).

Learns per-(metadata field x modality x cell-state) permutation-calibrated weights at pseudobulk
resolution via donor-grouped 5-fold CV, and writes the three Discovery tables + association index +
report under runs/<run_id>/:
  discovery_weights.tsv | discovery_markers.tsv | discovery_associations.tsv (+ _index.json)
  discovery_report.json | DISCOVERY.md  (field predictability + mutation-predictability matrix + optimize-me)

  # smoke (fast, plumbing only — 10 perms can't clear p<0.05 so weights are 0 by construction):
  python discover.py --run-id disc_smoke --fields subtype --modalities composition,ADT \
                     --cell-states-top 5 --screen-permutations 10 --final-permutations 0
  # real sweep:
  python discover.py --run-id discovery --fields subtype,ELN_risk,is_pediatric,disease_category \
                     --modalities composition,RNA,ADT,Metabolite --cell-states-top 40 \
                     --screen-permutations 30 --final-permutations 200
"""
from __future__ import annotations
import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amlmm
from amlmm import step
import amlmm.steps  # noqa: F401  (importing self-registers the `discover` step)


def _csv(s):
    return [x.strip() for x in s.split(",") if x.strip()] if s else None


def main(argv=None):
    ap = argparse.ArgumentParser(description="MOSAIC-AML Discovery agent")
    ap.add_argument("--base", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--run-id", default="discovery")
    ap.add_argument("--fields", default=None, help="comma-separated (default: all candidate fields)")
    ap.add_argument("--mutations-only", action="store_true",
                    help="sweep only the mutation/cyto flags (enumerated at runtime); the headline matrix")
    ap.add_argument("--merge-from", default=None,
                    help="comma-separated run dirs (or 'auto' = glob disc__* under --out) to MERGE into "
                         "--run-id instead of sweeping. Stitches finished fan-out outputs; runs no CV.")
    ap.add_argument("--modalities", default=None, help="comma-separated (default: all modalities)")
    ap.add_argument("--cell-states-top", type=int, default=40)
    ap.add_argument("--no-sample-level", action="store_true",
                    help="skip the sample-level modalities (composition/cell-comm/LSC)")
    ap.add_argument("--screen-permutations", type=int, default=30)
    ap.add_argument("--final-permutations", type=int, default=200)
    ap.add_argument("--screen-promote-alpha", type=float, default=0.10)
    ap.add_argument("--prefilter-features", type=int, default=800)
    ap.add_argument("--max-features", type=int, default=300)
    ap.add_argument("--feature-selector", default="f_classif", choices=["f_classif", "mutual_info"])
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    d = amlmm.Config()
    cfg = amlmm.Config(base_dir=args.base or d.base_dir, out_dir=args.out or d.out_dir,
                       run_id=args.run_id)
    ctx = amlmm.build_context(cfg)
    print(f"[discover] layout={ctx.layout} run_id={args.run_id}")

    # --- merge mode: stitch finished per-field fan-out dirs into one combined sweep (no CV) ---
    if args.merge_from:
        from amlmm import discovery as D
        import glob
        out_root = args.out or d.out_dir
        if args.merge_from.strip() == "auto":
            dirs = sorted(p for p in glob.glob(os.path.join(out_root, "disc__*")) if os.path.isdir(p))
        else:
            dirs = [x.strip() for x in args.merge_from.split(",") if x.strip()]
        print(f"[discover] merging {len(dirs)} run dir(s): {[os.path.basename(x) for x in dirs]}")
        rep = D.merge_discovery(dirs, ctx.run_dir)
        s = rep["summary"]
        print(f"[discover] merged: {s['n_significant']}/{s['n_ok']} significant combos; "
              f"predictable fields: {s['fields_predictable'] or '(none)'}")
        print(f"[discover] outputs: {ctx.run_dir}")
        return 0

    # --- mutations-only: enumerate the mutation/cyto flags at runtime ---
    fields = _csv(args.fields)
    if args.mutations_only:
        from amlmm import discovery as D
        fields = [f for f in D.candidate_fields(ctx) if f.startswith(("mut_", "cyto_"))]
        print(f"[discover] mutations-only: {len(fields)} flags")

    params = {
        "fields": fields, "modalities": _csv(args.modalities),
        "cell_states_top": args.cell_states_top, "sample_level": not args.no_sample_level,
        "screen_permutations": args.screen_permutations, "final_permutations": args.final_permutations,
        "screen_promote_alpha": args.screen_promote_alpha,
        "prefilter_features": args.prefilter_features, "max_features": args.max_features,
        "feature_selector": args.feature_selector, "verbose": not args.quiet,
    }
    res = step.run_step(step.get("discover"), ctx, params)
    if res.status != "ok":
        print(f"[discover] FAILED: {res.error}")
        return 1
    m = res.metrics
    print(f"\n[discover] {m['n_significant']}/{m['n_ok']} significant combos "
          f"(skipped {m['n_skipped']}, error {m['n_error']}) in {res.seconds}s")
    print(f"[discover] predictable fields: {m['fields_predictable'] or '(none)'}")
    print(f"[discover] markers={m['n_markers']} associations={m['n_associations']}")
    print(f"[discover] outputs: {ctx.run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
