#!/usr/bin/env python3
"""Ingest a NEW external sample (single-cell OR bulk RNA) as a patient and emit a decision-board report.

  --sample : gene-level scRNA (10x filtered_feature_bc_matrix.h5 | 10x dir | .h5ad)
        -> cell-state COMPOSITION  (cosine to the 89-state cellHarmony BM reference; amlmm.scrna)
        -> COHORT-TRAINED subtype prediction  (RF on the atlas composition->subtype; predict this sample)
        -> PRIMARY variant-level mutation panel (bulk-equivalent -> bulk_mutation_predictor)
        -> optional GENETIC ANCHOR  (from user-supplied mutations)
        -> arbiter v2  ->  runs/<id>/{patient_report.json, ledger.json, PATIENT.md}

  --bulk   : a BULK RNA expression table (gene x value; tsv/csv/txt)
        -> PRIMARY variant-level mutation panel (bulk_mutation_predictor, --bulk-ref cohort)
        -> optional GENETIC ANCHOR  ->  runs/<id>/patient_report.json  (mode "bulk_panel")
        Cell-state composition / subtype / cytogenetics / control-gate need single cells and are skipped.

The report is byte-shape-identical to what `panel.py --patient` writes, so the GUI renders it the
same way. This is the NEW-PATIENT entry path, with two honest properties:
  * the sample is NOT in the cohort, so the composition witness is a cohort-TRAINED classifier applied
    out-of-cohort (held_out=False) — a hypothesis, weighted by the cohort's honest CV accuracy;
  * with no mutations supplied the arbiter runs ANCHOR-FREE (leading hypothesis = composition's call,
    flagged "confirm by sequencing"); supplying even one driver engages the deterministic anchor.

Imputed descriptive witnesses (GRN / ADT / metabolic / lipid / LSC / UDON / cell-comm) are NOT YET
wired for uploads (each needs its own imputer/baseline) — they are listed under `ingest.deferred`
and simply absent from the panel; the arbiter handles a partial roster.

Run on an LSF COMPUTE node (head node OOM-kills the CV step):
  bsub -q normal -K -M 8000 -R "rusage[mem=8000]" /usr/local/anaconda3-2020/bin/python \
    ingest_patient.py --sample <10x_dir|h5> --name "Patient X" [--mutations "TP53,FLT3"] \
    [--dataset uploaded] [--run-id <id>] [--out-root <runs dir>] [--no-llm]
"""
from __future__ import annotations
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd

from amlmm.context import Config, build_context
from amlmm import scrna, genetics, arbiter
from amlmm import ledger as _led
from amlmm.agent import AgentResult
from amlmm.panel import evidence_predictive, _write_patient_md
from amlmm.bulk_predictor import BulkMutationPredictor
import control_gate as CG

BULK_PKL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bulk_mutation_predictor.pkl")

DEFERRED = ["surfaceome/ADT", "metabolic", "lipid", "GRN-regulon", "LSC",
            "cell-state/UDON", "cell-communication"]
# accepted genetic tokens normalize to the arbiter's de-prefixed flag space (ANCHOR_MAP keys)
SYN = {"NPM1C": "NPM1", "FLT3-ITD": "FLT3", "FLT3ITD": "FLT3", "CKIT": "KIT",
       "INV16": "inv16", "INV(16)": "inv16", "T(8;21)": "t8_21", "T8_21": "t8_21",
       "T(15;17)": "t15_17", "T15_17": "t15_17", "APL": "t15_17", "KMT2A": "kmt2a",
       "MLL": "kmt2a", "DEL7": "del7", "DEL(7)": "del7", "-7": "del7", "MONOSOMY7": "del7",
       "DEL5": "del5", "DEL(5)": "del5", "TRISOMY8": "trisomy8", "+8": "trisomy8",
       "COMPLEX": "complex"}


def _status(run_dir, state, step, message=""):
    """Tiny progress file the GUI polls (state: running|done|error)."""
    try:
        os.makedirs(run_dir, exist_ok=True)
        with open(os.path.join(run_dir, "status.json"), "w", encoding="utf-8") as f:
            json.dump({"state": state, "step": step, "message": message, "ts": time.time()}, f)
    except Exception:
        pass


def load_query(path):
    """gene-level scRNA -> AnnData with Gene-Expression features only (drops ADT). 10x h5/dir via
    scanpy; .h5ad read directly. Returns (adata, note)."""
    p = path
    if os.path.isdir(p):
        h5 = os.path.join(p, "filtered_feature_bc_matrix.h5")
        p = h5 if os.path.exists(h5) else p
    if p.endswith(".h5ad"):
        import anndata as ad
        a = ad.read_h5ad(p)
        a.var_names_make_unique()
        return a, "h5ad"
    # 10x feature-barcode HDF5 (CITE-seq: split off Antibody Capture)
    import scanpy as sc
    a = sc.read_10x_h5(p, gex_only=False)
    if "feature_types" in a.var:
        a = a[:, a.var["feature_types"].astype(str) == "Gene Expression"].copy()
    a.var_names_make_unique()
    return a, "10x_h5"


def normalize_mutations(tokens):
    valid = set(arbiter.ANCHOR_MAP)                      # de-prefixed flags that can anchor
    upper_to_flag = {k.upper(): k for k in valid}
    present, unknown = [], []
    for raw in tokens:
        t = str(raw).strip()
        if not t:
            continue
        u = t.upper()
        if u in SYN:
            present.append(SYN[u])
        elif u in upper_to_flag:
            present.append(upper_to_flag[u])
        else:
            unknown.append(t)                            # recorded but cannot anchor
    # stable de-dup, anchor-priority order so the report is reproducible
    seen, order = set(), []
    for f in present:
        if f not in seen:
            seen.add(f); order.append(f)
    return order, unknown


def composition_result(ctx, adata, permutations):
    """Composition witness for a NEW sample: cosine composition + cohort-trained prediction."""
    res = scrna.composition_from_query(adata)            # {composition, mean_cosine, n_cells, ...}
    comp = res["composition"]                            # Series over the 89 reference states (fractions)
    if res["n_shared_markers"] < 50:
        raise ValueError(f"only {res['n_shared_markers']} shared marker genes — not a gene-level query")

    # cohort balanced-accuracy + permutation p via the validated CV path (also leaves X,y in artifacts)
    ev = evidence_predictive(ctx, "subtype", "composition", strategy="donor_kfold",
                             permutations=permutations)
    X, y = ctx.getart("X"), ctx.getart("y")
    if X is None or ev.get("status") != "ok":
        raise RuntimeError(f"cohort composition model unavailable: {ev.get('reason', ev.get('status'))}")

    # fit a final RF on the whole cohort composition -> subtype, predict THIS sample
    from amlmm import models
    est = models.build(["rf"])["rf"][0]
    est.fit(np.asarray(X.values, float), np.asarray(y).astype(str))
    cols = list(X.columns)                               # 'comp::<state>'
    vec = np.array([float(comp.get(c.split("::", 1)[-1], 0.0)) for c in cols], dtype=float)
    s = vec.sum()
    if s > 0:
        vec = vec / s                                    # L1, matching assemble_features
    pred = str(est.predict(vec.reshape(1, -1))[0])
    prob = (round(float(np.max(est.predict_proba(vec.reshape(1, -1)))), 4)
            if hasattr(est, "predict_proba") else None)

    top = sorted(((str(k), float(v)) for k, v in comp.items()), key=lambda kv: (-kv[1], kv[0]))[:6]
    evidence = {
        "witness": "composition", "kind": "measured", "status": "ok", "held_out": False,
        "patient_prediction": pred, "patient_probability": prob, "true_label": None,
        "cohort_balanced_accuracy": ev.get("balanced_accuracy"),
        "cohort_permutation_p": ev.get("permutation_pvalue"),
        "patient_top_cellstates": {k: round(v, 3) for k, v in top},
        "mean_cosine": round(res["mean_cosine"], 3), "n_cells": res["n_cells"],
        "n_shared_markers": res["n_shared_markers"], "n_states_present": res["n_states_present"],
        "trained_classes": sorted(pd.unique(np.asarray(y).astype(str)).tolist()),
    }
    opinion = {"confidence": 0.5, "reliability_weight": 0.85,
               "summary": f"predicts {pred} (prob {prob}); cohort balanced-acc "
                          f"{ev.get('balanced_accuracy')}",
               "caveats": "external sample: cohort-trained classifier, NOT held-out for this sample; "
                          "out-of-cohort domain shift possible"}
    meta = {"mean_cosine": res["mean_cosine"], "n_cells": res["n_cells"],
            "n_states_present": res["n_states_present"], "n_shared_markers": res["n_shared_markers"]}
    _g = CG.load_gate()                                  # healthy-vs-diseased gate on the composition vector
    meta["control_gate"] = CG.score_gate(_g, comp) if _g is not None else None
    return AgentResult("composition", "predictive", "honest_cv", "independent",
                       status="ok", evidence=evidence, opinion=opinion), meta


def bulk_equiv_from_adata(adata):
    """Collapse a single-cell sample to a whole-sample bulk-equivalent: sum counts across all cells per
    gene -> CP10k (linear), the same representation the atlas 'sc' reference was built from. Returns a
    Series indexed by gene symbol (the predictor maps symbol->ENSG internally)."""
    import scipy.sparse as sp
    X = adata.X
    tot = np.asarray(X.sum(axis=0)).ravel() if sp.issparse(X) else np.asarray(X).sum(axis=0).ravel()
    s = float(tot.sum()) or 1.0
    cp10k = tot / s * 1e4                                    # counts-per-10k; predictor applies log2 + z internally
    return pd.Series(cp10k, index=[str(g) for g in adata.var_names])


def parse_bulk_expression(path, scale="auto", column=None):
    """Parse an uploaded BULK RNA expression file -> linear Series(gene -> value).

    Accepts any delimited table (tsv/csv/txt) with genes in the first column and one or more sample
    columns (ENSG or symbol index; the predictor maps symbols internally). Picks `column` if given, else
    the first numeric column. Linearizes to the scale the predictor expects (log2(x+1) is applied inside):
    log2-normalized (has negatives, e.g. BeatAML norm_exp) -> 2^x; log1p (compressed, non-negative) ->
    expm1; linear counts/TPM/RPKM -> as-is. Auto-detect is overridable via `scale`. Returns (series, col, scale)."""
    df = pd.read_csv(path, sep=None, engine="python", index_col=0)
    num = df.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")
    if num.shape[1] == 0:
        raise ValueError("no numeric expression column found in the bulk file")
    col = column if (column and column in num.columns) else str(num.columns[0])
    ser = num[col].dropna()
    ser.index = ser.index.astype(str).str.split(".").str[0]   # strip ENSG version suffixes
    ser = ser.groupby(ser.index).max()                        # collapse duplicate gene rows
    v = ser.values.astype(float)
    if scale == "auto":
        if np.nanmin(v) < -0.01:
            scale = "log2"
        elif np.nanmax(v) <= 30:
            scale = "log1p"
        else:
            scale = "linear"
    if scale == "log2":
        ser = np.power(2.0, ser)                              # BeatAML norm_exp is log2 -> linearize
    elif scale == "log1p":
        ser = np.expm1(ser)                                   # CP10k+log1p -> linear
    return ser, col, scale


def bulk_mutation_result(expr_linear, present_truth, ref="sc"):
    """PRIMARY mutation caller: predict ~50 variant-level categories from a LINEAR gene-expression Series
    (a single-cell bulk-equivalent with ref='sc', or an uploaded bulk sample with ref='beataml').
    Returns (predictions_list, caller_meta, AgentResult) or (None, None, None) if the model isn't available.
    predictions_list is the GUI-native list (each prediction carries heldout_auc := CV AUROC for the
    reliability/abstain logic); caller_meta holds the panel-level metadata."""
    if not os.path.exists(BULK_PKL):
        return None, None, None
    BP = BulkMutationPredictor.load(BULK_PKL)
    calls = BP.predict(expr_linear, ref=ref)                 # {category: {probability, call, threshold, cv_auroc, confidence}}
    truth = {str(t).upper() for t in present_truth}
    preds = []
    for cat, v in calls.items():
        pr = dict(v); pr["mutation"] = cat
        pr["heldout_auc"] = v.get("cv_auroc")               # GUI keys reliability/abstain off heldout_auc
        pr["source"] = "bulk_rna (BeatAML-trained, variant-level)"
        gene = cat.split("_")[0].upper()                    # coarse gene for a cross-check vs supplied truth
        pr["supplied_truth"] = ("present" if (gene in truth or cat.upper() in truth) else None)
        preds.append(pr)
    n_called = sum(1 for p in preds if p["call"] == "present" and p["confidence"] == "ok")
    _inp = "uploaded bulk RNA" if ref != "sc" else "the sample's single-cell bulk-equivalent"
    caller_meta = {"mode": "bulk_variant_primary", "caller": "BulkMutationPredictor", "input": ref,
                   "predictor": BP.summary(), "n_confident_present": n_called,
                   "note": "PRIMARY mutation caller — ~50 variant-level categories predicted from %s "
                           "(BeatAML-trained, validated cross-cohort). Predicted, not sequenced: high-confidence "
                           "'present' calls are leads to confirm by sequencing. Weak categories (CV AUROC < 0.65) "
                           "are marked abstain." % _inp}
    top = [p["mutation"] for p in preds if p["call"] == "present" and p["confidence"] == "ok"][:8]
    ev = {"witness": "bulk_mutation", "kind": "predictive", "status": "ok",
          "n_confident_present": n_called, "top_predicted_present": top,
          "mean_cv_auroc": BP.summary().get("mean_cv_auroc"), "n_categories": len(BP.categories)}
    op = {"confidence": 0.5, "reliability_weight": 0.6,
          "summary": f"predicts {n_called} driver categories present (top: {top})",
          "caveats": "mutations PREDICTED from bulk-equivalent expression, not sequenced; confirm by DNA"}
    res = AgentResult("bulk_mutation", "predictive", "honest_cv", "independent",
                      status="ok", evidence=ev, opinion=op)
    return preds, caller_meta, res


def genetic_result(present, unknown, eln):
    targetable = {g.replace("mut_", ""): t for g, t in genetics.TARGETABLE.items()
                  if g.replace("mut_", "") in present}
    evidence = {"witness": "genetic", "kind": "independent", "status": "ok",
                "present": present, "targetable": targetable, "eln": eln}
    if unknown:
        evidence["unrecognized_tokens"] = unknown
    opinion = {"confidence": 0.6 if present else 0.2,
               "reliability_weight": 0.85 if present else 0.1,
               "summary": f"present: {present}; targetable: {list(targetable)}"
                          + (f"; unrecognized: {unknown}" if unknown else ""),
               "caveats": "user-supplied mutation calls" if present
                          else "no mutations supplied — arbiter runs anchor-free"}
    return AgentResult("genetic", "genetic", "deterministic_fact", "independent",
                       status="ok", evidence=evidence, opinion=opinion)


def provenance(layout, sample_path, deferred):
    import platform
    out = {"layout": layout, "ingested": True, "source": sample_path,
           "python": platform.python_version(), "deferred_witnesses": deferred}
    for mod in ("scikit-learn", "numpy", "anndata", "scanpy"):
        try:
            import importlib.metadata as im
            out[mod] = im.version(mod)
        except Exception:
            pass
    return out


def main_bulk(args, cfg, run_dir, sample_key):
    """BULK-RNA ingest path: an expression file -> the PRIMARY variant-level mutation panel (+ genetic anchor).
    No cell-state composition / subtype / cytogenetics / control-gate — those require single-cell input — so
    this path skips the atlas load entirely and is fast."""
    try:
        _status(run_dir, "running", "loading_bulk", "parsing bulk expression file")
        ser, col, scale = parse_bulk_expression(args.bulk, scale=args.bulk_scale)

        _status(run_dir, "running", "genetic", "applying genetic anchor")
        present, unknown = normalize_mutations(args.mutations.split(","))
        gen_res = genetic_result(present, unknown, args.eln)

        _status(run_dir, "running", "mutation_calling", "predicting drivers from bulk RNA (primary)")
        mut_preds, mut_caller, _ = bulk_mutation_result(ser, present, ref=args.bulk_ref)
        if mut_preds is None:
            raise RuntimeError("bulk mutation predictor not found (pipeline/bulk_mutation_predictor.pkl)")
        top_present = [p["mutation"] for p in mut_preds
                       if p["call"] == "present" and p["confidence"] == "ok"]

        consensus = {
            "leading_hypothesis": ("predicted drivers: " + ", ".join(top_present[:4])) if top_present
                                  else "no high-confidence drivers predicted",
            "overall_confidence": "mutation panel only — bulk RNA input has no subtype/cytogenetics witness",
            "leading_confirmed_by_genetics": bool(present), "supplied_mutations": present,
        }
        prov = {"layout": "bulk", "ingested": True, "source": args.bulk, "input_kind": "bulk_rna",
                "deferred_witnesses": ["composition/subtype", "cytogenetics", "control-gate"] + DEFERRED}
        report = {
            "mode": "bulk_panel", "sample_key": sample_key, "annotation": None, "dataset": args.dataset,
            "provenance": prov, "specimen_class": None, "control_gate": None,
            "panel": [{"witness": gen_res.name, "grounding": gen_res.grounding,
                       "independence": gen_res.independence, "evidence": gen_res.evidence,
                       "opinion": gen_res.opinion}],
            "mutation_predictions": mut_preds, "mutation_caller": mut_caller,
            "consensus": consensus, "deliberation": None,
            "ingest": {"source": args.bulk, "input_kind": "bulk_rna", "name": args.name,
                       "bulk_ref": args.bulk_ref, "bulk_scale_detected": scale, "bulk_column": col,
                       "n_genes": int(ser.shape[0]), "mutations_supplied": present,
                       "unrecognized_mutations": unknown,
                       "note": "BULK RNA input — variant-level mutation panel only. Cell-state composition, "
                               "subtype prediction, cytogenetics, and the healthy-vs-diseased control gate "
                               "require single-cell input and are not produced for a bulk sample."},
        }
        os.makedirs(run_dir, exist_ok=True)
        json.dump(report, open(os.path.join(run_dir, "patient_report.json"), "w"), default=str, indent=1)
        _status(run_dir, "done", "done",
                f"{sample_key}: {mut_caller['n_confident_present']} drivers predicted (bulk RNA)")
        print("OK %s\n  sample_key: %s\n  input: bulk RNA (%d genes, scale=%s, ref=%s)\n  drivers present: %s"
              % (run_dir, sample_key, ser.shape[0], scale, args.bulk_ref, ", ".join(top_present) or "none"))
        return 0
    except Exception as e:
        import traceback
        traceback.print_exc()
        _status(run_dir, "error", "failed", f"{type(e).__name__}: {e}")
        print(f"ERROR {type(e).__name__}: {e}", file=sys.stderr)
        return 1


def main():
    ap = argparse.ArgumentParser(description="Ingest an scRNA or bulk-RNA sample as a MATRIX-AML patient.")
    ap.add_argument("--sample", help="single-cell input: 10x filtered_feature_bc_matrix.h5, a 10x dir, or .h5ad")
    ap.add_argument("--bulk", help="BULK RNA input: a gene x value table (tsv/csv/txt). Mutually exclusive with --sample")
    ap.add_argument("--bulk-ref", default="beataml", choices=["beataml", "leucegene", "sc"],
                    help="z-score reference cohort for a bulk upload (default beataml = the training cohort)")
    ap.add_argument("--bulk-scale", default="auto", choices=["auto", "linear", "log2", "log1p"],
                    help="expression scale of the bulk file (default auto-detect)")
    ap.add_argument("--name", required=True, help="patient/sample display name")
    ap.add_argument("--mutations", default="", help="comma-separated observed drivers, e.g. 'TP53,FLT3'")
    ap.add_argument("--eln", default=None, help="optional ELN risk label")
    ap.add_argument("--dataset", default="uploaded", help="dataset/source label for the report")
    ap.add_argument("--run-id", default=None, help="runs/<id> subdir (default: slug of --name)")
    ap.add_argument("--out-root", default=None, help="override the runs/ root dir")
    ap.add_argument("--permutations", type=int, default=40)
    ap.add_argument("--no-llm", action="store_true", help="skip LLM narration (deterministic only)")
    args = ap.parse_args()
    if bool(args.sample) == bool(args.bulk):
        ap.error("provide exactly one of --sample (single-cell) or --bulk (bulk RNA)")

    run_id = args.run_id or ("ingest_" + "".join(ch if ch.isalnum() else "_" for ch in args.name).strip("_"))
    cfg = Config(run_id=run_id)
    if args.out_root:
        cfg = Config(out_dir=args.out_root, run_id=run_id)
    run_dir = os.path.join(cfg.out_dir, run_id)
    _status(run_dir, "running", "starting", f"ingesting {args.name}")
    sample_key = f"{args.dataset}::{args.name}"
    if args.bulk:
        return main_bulk(args, cfg, run_dir, sample_key)

    try:
        _status(run_dir, "running", "loading_cohort", "loading atlas + knowledge base")
        ctx = build_context(cfg)
        ctx.set("provenance", provenance(ctx.layout, args.sample, DEFERRED))

        _status(run_dir, "running", "composition", "assigning cells -> composition + subtype")
        adata, how = load_query(args.sample)
        comp_res, meta = composition_result(ctx, adata, args.permutations)

        _status(run_dir, "running", "genetic", "applying genetic anchor")
        present, unknown = normalize_mutations(args.mutations.split(","))
        gen_res = genetic_result(present, unknown, args.eln)

        _status(run_dir, "running", "mutation_calling", "predicting drivers from bulk-equivalent (primary)")
        mut_preds, mut_caller = None, None
        try:
            mut_preds, mut_caller, _ = bulk_mutation_result(bulk_equiv_from_adata(adata), present, ref="sc")
        except Exception as _e:
            print("WARN bulk mutation caller skipped: %s" % _e, file=sys.stderr)

        _status(run_dir, "running", "arbiter", "reconciling witnesses")
        led = _led.Ledger(sample_key, "subtype", provenance=ctx.getart("provenance"))
        ctx.ledger = led
        led.append(comp_res, round=0)
        led.append(gen_res, round=0)
        led.persist(ctx)

        client = None
        if not args.no_llm:
            try:
                from amlmm.llm import LLMClient
                client = LLMClient()
            except Exception:
                client = None
        consensus = arbiter.reconcile_patient(client, ctx, led)
        led.set_arbiter(consensus, round=0)
        led.finalize("single_pass")
        led.persist(ctx)

        gate_call = meta.get("control_gate")
        if gate_call and gate_call.get("call") == "control":     # healthy-vs-diseased gate fires before mutation calling
            consensus["subtype_if_diseased"] = consensus.get("leading_hypothesis")
            consensus["leading_hypothesis"] = "control / no mutation"
            consensus["specimen_class"] = "control"
        report = {
            "mode": "patient_panel", "sample_key": sample_key, "annotation": None,
            "dataset": args.dataset, "provenance": ctx.getart("provenance"),
            "specimen_class": (gate_call or {}).get("call"), "control_gate": gate_call,
            "panel": [{"witness": e["witness"], "grounding": e["grounding"],
                       "independence": e["independence"], "evidence": e["evidence"],
                       "opinion": e["opinion"]} for e in led.current_entries()],
            "mutation_predictions": mut_preds,        # PRIMARY variant-level mutation caller (bulk-equivalent)
            "mutation_caller": mut_caller,
            "consensus": consensus, "deliberation": None,
            "ingest": {"source": args.sample, "input_kind": how, "name": args.name,
                       "mutations_supplied": present, "unrecognized_mutations": unknown,
                       "composition_quality": meta, "deferred_witnesses": DEFERRED},
        }
        ctx.save_json(report, "patient_report.json")
        _write_patient_md(ctx, report)
        _status(run_dir, "done", "done",
                f"{sample_key}: leading {consensus.get('leading_hypothesis')} "
                f"({consensus.get('overall_confidence')})")
        print(f"OK {run_dir}\n  sample_key: {sample_key}\n  specimen_class: {(gate_call or {}).get('call')} "
              f"(gate p_diseased={gate_call.get('p_diseased') if gate_call else 'NA'})\n  "
              f"leading: {consensus.get('leading_hypothesis')} "
              f"(confirmed_by_genetics={consensus.get('leading_confirmed_by_genetics')}, "
              f"{consensus.get('overall_confidence')})")
        return 0
    except Exception as e:
        import traceback
        traceback.print_exc()
        _status(run_dir, "error", "failed", f"{type(e).__name__}: {e}")
        print(f"ERROR {type(e).__name__}: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
