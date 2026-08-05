#!/usr/bin/env python3
"""MOSAIC-AML panel CLI.

Cohort mode (default): which witnesses can predict the target across the cohort.
  python panel.py --target subtype
  python panel.py --target subtype --blocks composition,ADT,GRN

Per-patient mode (--patient): a tumor-board report for ONE patient.
  python panel.py --patient "NYU-1::AML-01"
"""
from __future__ import annotations
import os
import re
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amlmm
from amlmm.panel import run_panel, run_patient_panel


def _provenance(ctx, args):
    def _ver(pkg):
        try:
            import importlib.metadata as m
            return m.version(pkg)
        except Exception:
            return "?"
    import platform
    ctx.set("provenance", {"python": platform.python_version(),
                           "scikit-learn": _ver("scikit-learn"), "numpy": _ver("numpy"),
                           "layout": ctx.layout})


def main(argv=None):
    ap = argparse.ArgumentParser(description="MOSAIC-AML per-witness panel")
    ap.add_argument("--target", default="subtype")
    ap.add_argument("--patient", default=None, help="sample_key (Dataset::Sample) -> per-patient report")
    ap.add_argument("--blocks", default="composition",
                    help="comma list of predictive feature blocks (default: composition)")
    ap.add_argument("--no-genetic", action="store_true")
    ap.add_argument("--no-udon", action="store_true")
    ap.add_argument("--strategy", default="donor_kfold",
                    choices=["donor_kfold", "leave_one_cohort_out"])
    ap.add_argument("--permutations", type=int, default=40)
    ap.add_argument("--rounds", type=int, default=0,
                    help="patient mode: feedback-deliberation rounds (0=single pass; >=1 enables Phase C)")
    ap.add_argument("--mode", default="continuous", choices=["continuous", "conflict_triggered"],
                    help="patient feedback mode (only used with --rounds>=1)")
    ap.add_argument("--base", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--run-id", default=None)
    args = ap.parse_args(argv)

    d = amlmm.Config()
    default_id = (("patient_" + re.sub(r"[^A-Za-z0-9]+", "_", args.patient)) if args.patient
                  else f"panel_{args.target}")
    cfg = amlmm.Config(base_dir=args.base or d.base_dir, out_dir=args.out or d.out_dir,
                       run_id=args.run_id or default_id)
    ctx = amlmm.build_context(cfg)
    _provenance(ctx, args)
    blocks = [b.strip() for b in args.blocks.split(",") if b.strip()]

    if args.patient:
        print(f"[patient] {args.patient} layout={ctx.layout} blocks={blocks} "
              f"genetic={not args.no_genetic} udon={not args.no_udon} "
              f"rounds={args.rounds} mode={args.mode}")
        rep = run_patient_panel(ctx, args.patient, blocks=blocks,
                                with_genetic=not args.no_genetic, with_udon=not args.no_udon,
                                strategy=args.strategy, permutations=args.permutations,
                                max_rounds=args.rounds, mode=args.mode)
        print(f"\n=== patient {args.patient} (atlas annotation: {rep.get('annotation')}) ===")
        for it in rep["panel"]:
            o = it["opinion"]
            print(f"  {str(it.get('witness')):18s} ({it.get('independence')}/{it.get('grounding')}) "
                  f"conf={o.get('confidence')} w={o.get('reliability_weight')}: {str(o.get('summary'))[:110]}")
        c = rep["consensus"]
        print(f"\nSUBTYPE CALL [{c.get('overall_confidence')}]: {c.get('subtype_call')} "
              f"(concordance {c.get('concordance')}, confirmed_by_genetics="
              f"{c.get('leading_confirmed_by_genetics')})")
        print(f"CONFLICTS: {c.get('conflicts') or 'none'}")
        ths = c.get("ranked_therapy_hypotheses") or []
        print("THERAPIES: " + ("; ".join(f"{t.get('biomarker')}->{t.get('drug')}"
                                          for t in ths) or "(none observed)"))
        vals = c.get("recommended_validations") or []
        print("VALIDATIONS: " + ("; ".join(f"{v.get('claim')}:{v.get('validation')}"
                                            for v in vals) or "(none)"))
        surf = c.get("surface_therapy_hypotheses") or []
        if surf:
            print("SURFACE (flow-pending): "
                  + "; ".join(str(f"{t.get('biomarker')}->{t.get('drug')}")[:55] for t in surf))
        desc = c.get("descriptive_findings") or []
        if desc:
            who = ", ".join(sorted({str(d.get("witness")) for d in desc}))
            print(f"DESCRIPTIVE FINDINGS: {len(desc)} note(s) from {who}")
        else:
            print("DESCRIPTIVE FINDINGS: none")
        dl = rep.get("deliberation")
        if dl:
            dd = dl.get("drift", {})
            print(f"DELIBERATION [{dl.get('mode')}]: {dl.get('deliberation_rounds')} round(s) "
                  f"({dl.get('stop_reason')}); leading {dd.get('baseline_leading')}->{dd.get('final_leading')} "
                  f"(changed={dd.get('leading_changed')}); concordance "
                  f"{dd.get('baseline_concordance')}->{dd.get('final_concordance')}; "
                  f"groupthink_warning={dd.get('groupthink_warning')}")
        print(f"\noutputs: {ctx.run_dir}")
        return 0

    print(f"[panel] target={args.target} layout={ctx.layout} blocks={blocks} "
          f"genetic={not args.no_genetic} udon={not args.no_udon}")
    rep = run_panel(ctx, args.target, blocks=blocks,
                    with_genetic=not args.no_genetic, with_udon=not args.no_udon,
                    strategy=args.strategy, permutations=args.permutations)
    print("\n=== witnesses ===")
    for it in rep["panel"]:
        e, o = it["evidence"], it["opinion"]
        ba = e.get("balanced_accuracy")
        extra = f"ba={ba} p={e.get('permutation_pvalue')}" if ba is not None else ""
        print(f"  {str(e.get('witness')):16s} kind={str(e.get('kind')):16s} "
              f"conf={o.get('confidence')} w={o.get('reliability_weight')} {extra}")
    c = rep["consensus"]
    print(f"\nCONSENSUS [{c.get('overall_confidence')}]: {str(c.get('consensus'))[:400]}")
    print(f"TARGETABLE: {str(c.get('targetable_summary'))[:300]}")
    print(f"VALIDATIONS: {str(c.get('recommended_validations'))[:300]}")
    print(f"consistency: {c.get('per_witness_consistency')}")
    print(f"\noutputs: {ctx.run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
