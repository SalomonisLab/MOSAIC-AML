#!/usr/bin/env python3
"""Add GENE-LEVEL validation to each report's roster summary.

The caller is variant-level (58 categories) but the ground truth is mostly gene-level, so a variant-level
"X/N correct" can only score the 17 categories whose exact variant is known. The chosen fix (per the
maintainer) is to validate at the GENE: aggregate a gene's variant categories to one gene-level call
(present iff ANY of its categories is called present) and compare to the gene's known status. This lifts
the denominator to ~40 without asserting variant-level truth we do not have (R882 vs non-R882 etc.).

Truth sources (no numpy -> runs on the login node with the conda python; stdlib csv/json only):
  predict_*  the sc mutation matrix (mutation_matrix_explicit_v2.tsv, gene-level 0/1, keyed by sample_key)
             + the cytogenetic calls already carried on the report (they hold their own true_label)
  trumpp_*   Table S4 lists the DRIVERS ONLY. A gene absent from that table was not necessarily tested
             and found wild-type -- it may simply not have been reported. Those genes are therefore left
             UNLABELLED, not scored as confirmed negatives (see below).
  ingest_*   external upload -> no ground truth -> no validation (left blank, correctly)

WHY THE HEADLINE RATIO IS NOT ACCURACY. Two defects were found in the original "N/M correct":

  1. For Trumpp, every gene the model scores was labelled -- 646 of 688 units across 16 samples came out
     "absent", nearly all of them ASSUMED rather than tested. That is the same never-assayed-counted-as-
     wild-type error already corrected for the bulk cohorts, and it inflates both numerator and
     denominator. Cohorts without confirmed negatives now contribute confirmed POSITIVES only.
  2. Even where negatives are genuinely confirmed (Leucegene, dbGAP), the ratio is dominated by them.
     Leucegene::SRS718881 scored 26/27 -- 26 confirmed-absent genes called absent, and the single real
     driver (IDH2) MISSED. An all-absent predictor scores 26/27 on that sample. A number that a trivial
     predictor maxes out is not a measure of the caller.

So the record now carries the positives separately, and the interface leads with those.

Writes report["validation_gene"] = {n_correct, n_labeled, n_pos_labeled, n_pos_correct,
                                    n_neg_labeled, n_neg_correct, negatives_confirmed,
                                    units:{unit:{truth,call,correct}}}.

  /usr/local/anaconda3-2020/bin/python gene_validation.py [runs_dir]
"""
import os, sys, csv, json, glob

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(HERE), "runs")
MATRIX = os.path.join(HERE, "mutation_matrix_explicit_v2.tsv")

# cytogenetic categories -> a canonical validation unit (folds the aliases together)
CYTO = {"complex": "complex", "del5": "del5", "del7": "del7", "trisomy8": "trisomy8",
        "inv16": "inv16", "inv(16)_cbfb-myh11": "inv16",
        "kmt2a": "KMT2A", "kmt2a-rearrangement": "KMT2A", "kmt2a_fusion": "KMT2A", "kmt2a_ptd": "KMT2A"}


def unit(cat):
    """Map a category to its validation unit: a cytogenetic event, or the SNV gene."""
    cl = str(cat).lower()
    if cl in CYTO:
        return CYTO[cl]
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl:
        return "inv16"
    if "kmt2a" in cl:
        return "KMT2A"
    return str(cat).split("_")[0].split("-")[0].upper()


def load_matrix():
    with open(MATRIX) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        gcol = [(i, c[4:].upper()) for i, c in enumerate(hdr) if c.startswith("mut_")]
        rows = {}
        for r in rd:
            rows[r[0]] = {g: (r[i].strip() not in ("", "0", "0.0", "nan", "NA")) for i, g in gcol}
    return rows


def gene_truth_for(run, rep, matrix):
    """-> {unit: 'present'|'absent'} or {} when no ground truth exists for this sample."""
    key = str(rep.get("sample_key"))
    preds = rep.get("mutation_predictions") or []
    model_units = {unit(p["mutation"]) for p in preds}          # only score units the model actually predicts
    gt = {}
    # cytogenetics: truth already on the report's cyto rows
    for p in preds:
        u = unit(p["mutation"])
        if u in ("complex", "del5", "del7", "trisomy8", "inv16", "KMT2A") and p.get("true_label"):
            gt[u] = p["true_label"]
    if run.startswith("trumpp_"):
        # Table S4 lists the drivers and nothing else, so a gene missing from it is UNKNOWN, not
        # wild-type. Only the listed drivers get a label; everything else stays unlabelled and is
        # excluded from the score. This makes specificity unmeasurable on this cohort, which is the
        # honest position -- previously it was "measured" against 646 assumed negatives.
        for x in (rep.get("known_drivers") or []):
            u = unit(x)
            if u in model_units:
                gt[u] = "present"
        return gt
    if run.startswith("predict_") and key in matrix:
        row = matrix[key]                                        # gene -> mutated? (bool)
        for u in model_units:
            if u in ("complex", "del5", "del7", "trisomy8", "inv16", "KMT2A"):
                continue                                         # cyto handled above
            if u in row:
                gt[u] = "present" if row[u] else "absent"
        return gt
    return gt                                                   # ingest_ / unknown -> no truth


def main():
    matrix = load_matrix()
    updated = 0
    dist = {}
    for f in sorted(glob.glob(os.path.join(RUNS, "*", "patient_report.json"))):
        run = os.path.basename(os.path.dirname(f))
        if not run.startswith(("predict_", "trumpp_", "ingest_")):
            continue
        with open(f, encoding="utf-8") as fh:
            rep = json.load(fh)
        preds = rep.get("mutation_predictions")
        if not isinstance(preds, list) or not preds:
            continue
        gt = gene_truth_for(run, rep, matrix)
        # gene-level call: present iff ANY of the unit's categories is called present
        call = {}
        for p in preds:
            u = unit(p["mutation"])
            if p.get("call") == "present":
                call[u] = "present"
            else:
                call.setdefault(u, "absent")
        units = {}
        for u, t in sorted(gt.items()):
            c = call.get(u, "absent")
            units[u] = {"truth": t, "call": c, "correct": (c == t)}
        pos = [v for v in units.values() if v["truth"] == "present"]
        neg = [v for v in units.values() if v["truth"] == "absent"]
        rep["validation_gene"] = {
            "n_labeled": len(units),
            "n_correct": sum(1 for v in units.values() if v["correct"]),
            # the informative half: confirmed drivers actually recovered. An all-absent predictor scores
            # 0 here and near-perfect on the combined ratio, which is exactly why both are recorded.
            "n_pos_labeled": len(pos), "n_pos_correct": sum(1 for v in pos if v["correct"]),
            "n_neg_labeled": len(neg), "n_neg_correct": sum(1 for v in neg if v["correct"]),
            # No negative is inferred any more, so every labelled negative here is observed. Trumpp is
            # still special: after dropping the assumed wild-types its only negatives are the
            # cytogenetic units that carry an explicit label, so gene-level specificity is not
            # measurable on that cohort and the interface must not imply otherwise.
            "negatives_confirmed": True,
            "negative_scope": "cytogenetics_only" if run.startswith("trumpp_") else "all",
            "units": units}
        with open(f, "w", encoding="utf-8") as fh:
            json.dump(rep, fh, default=str, indent=1)
        updated += 1
        dist.setdefault(len(units), 0)
        dist[len(units)] += 1

    print("updated %d reports with gene-level validation" % updated)

    # Backfill the positive/negative split onto reports written by ANOTHER builder (the 367 Leucegene
    # boards come from build_leucegene_board.py). The new fields are pure derivations of `units`, so
    # recomputing them here is exactly equivalent to rebuilding those boards and costs nothing.
    back = 0
    for f in sorted(glob.glob(os.path.join(RUNS, "*", "patient_report.json"))):
        with open(f, encoding="utf-8") as fh:
            rep = json.load(fh)
        v = rep.get("validation_gene")
        if not isinstance(v, dict) or "n_pos_labeled" in v or not v.get("units"):
            continue
        u = v["units"]
        v["n_pos_labeled"] = sum(1 for d in u.values() if d["truth"] == "present")
        v["n_pos_correct"] = sum(1 for d in u.values() if d["truth"] == "present" and d["correct"])
        v["n_neg_labeled"] = sum(1 for d in u.values() if d["truth"] == "absent")
        v["n_neg_correct"] = sum(1 for d in u.values() if d["truth"] == "absent" and d["correct"])
        v.setdefault("negatives_confirmed", True)
        v.setdefault("negative_scope", "all")
        with open(f, "w", encoding="utf-8") as fh:
            json.dump(rep, fh, default=str, indent=1)
        back += 1
    print("backfilled the positive/negative split onto %d further reports" % back)
    print("n_labeled (genes+cyto units) distribution:", dict(sorted(dist.items())))
    # sample
    for f in sorted(glob.glob(os.path.join(RUNS, "predict_*", "patient_report.json")))[:2]:
        r = json.load(open(f))
        v = r.get("validation_gene") or {}
        print("  %-30s %d/%d correct across %d units"
              % (os.path.basename(os.path.dirname(f)), v.get("n_correct"), v.get("n_labeled"), v.get("n_labeled")))
    print("GENE VALIDATION OK")


if __name__ == "__main__":
    main()
