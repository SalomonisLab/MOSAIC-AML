#!/usr/bin/env python3
"""Benchmark the BULK RNA-ONLY mutation caller across every external cohort we have.

MODEL UNDER TEST: pipeline/bulk_mutation_predictor.pkl — trained on BeatAML2 bulk RNA ONLY.
No imputed modalities, no ADT/Lipid/Metabolite/GRN, single modality. That is exactly the requested
"bulk RNA-seq model, RNA only, no imputation" benchmark. Because it trained on BeatAML alone, every
cohort below is genuinely external to it.

COHORTS
  Leucegene   367  bulk held-out            (dbGAP variant calls, category-level truth)
  sc atlas    387  single-cell held-out     (the caller never saw sc at all; 26 are the sealed set)
  Trumpp       16  single-cell external     (Table S4 truth -> GENE level)
  MDS          64  AML-precursor cohort     (Hs_MDS_UDON; genotype from obs['Group'])

METRICS (per the request)
  AUROC  prevalence-invariant ranking
  AUPRC  reported WITH its baseline (= prevalence). AUPRC is prevalence-sensitive, so on an enriched
         cohort a high AUPRC can be trivial — the lift over baseline is what means something.
  F1     at the DEPLOYED threshold. Note this threshold was F1-max'd on BeatAML's prevalence; applying
         it to a cohort with very different prevalence penalises F1 for reasons that are not the model.

MDS labelling (honest):
  positives = Group contains SRSF2 / U2AF1 / RUNX1 / TET2 (the genotype is in the Group string)
  negatives = ONLY aged / young / Anemia. The 28 empty-Group samples are UNLABELLED and are excluded —
              using them as negatives would invent truth we do not have.
  MISTRG (humanised-mouse xenograft) samples are held out of the primary analysis and reported apart.
"""
import os, sys, json, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score, precision_score, recall_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
from amlmm.bulk_predictor import BulkMutationPredictor

BP = BulkMutationPredictor.load(os.path.join(ROOT, "pipeline", "bulk_mutation_predictor.pkl"))
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
GENES = [str(g) for g in d["genes"]]
CATS = [str(c) for c in d["drivers"]]
HOLD = set(str(x) for x in d["holdout"])
MIN_POS = 3


def metrics(y, p, thr):
    y = np.asarray(y); p = np.asarray(p, float)
    call = (p >= thr).astype(int)
    out = {"n": len(y), "n_pos": int(y.sum()), "prev": float(y.mean())}
    try:
        out["auroc"] = roc_auc_score(y, p)
    except Exception:
        out["auroc"] = np.nan
    out["auprc"] = average_precision_score(y, p)          # baseline == prevalence
    out["lift"] = out["auprc"] / max(out["prev"], 1e-9)   # AUPRC relative to chance
    out["f1"] = f1_score(y, call, zero_division=0)
    out["prec"] = precision_score(y, call, zero_division=0)
    out["rec"] = recall_score(y, call, zero_division=0)
    return out


def score_bundle(Xk, Lk, Sk, ref, restrict=None):
    X = pd.DataFrame(d[Xk].astype(float), index=[str(s) for s in d[Sk]], columns=GENES)
    L = pd.DataFrame(d[Lk].astype(float), index=X.index, columns=CATS)
    if restrict is not None:
        k = [s for s in X.index if s in restrict]
        X, L = X.loc[k], L.loc[k]
    Z = {s: BP._z(BP._clog(BP._align(X.loc[s])), ref) for s in X.index}
    res = {}
    for c in BP.categories:
        if c not in L.columns:
            continue
        y = L[c].dropna()
        if int((y == 1).sum()) < MIN_POS or int((y == 0).sum()) < MIN_POS:
            continue
        ids = list(y.index)
        p = [BP.predict_one(c, Z[s], ref)["probability"] for s in ids]
        res[c] = metrics(y.loc[ids].astype(int).values, p, BP.models[c]["threshold"])
    return res


def score_trumpp():
    gene_of = lambda c: str(c).split("_")[0].split("-")[0].upper()
    rows = []
    for f in sorted(glob.glob(os.path.join(ROOT, "runs", "trumpp_*", "patient_report.json"))):
        r = json.load(open(f))
        kn = set(str(x).upper() for x in (r.get("known_drivers") or []))
        best, thr = {}, {}
        for p in (r.get("mutation_predictions") or []):
            g = gene_of(p["mutation"])
            if p.get("probability") is None:
                continue
            if p["probability"] > best.get(g, -1):
                best[g] = float(p["probability"]); thr[g] = float(p.get("threshold") or 0.5)
        rows.append((kn, best, thr))
    res = {}
    for g in sorted(set().union(*[set(b) for _, b, _ in rows])):
        y = [1 if g in kn else 0 for kn, b, _ in rows if g in b]
        p = [b[g] for _, b, _ in rows if g in b]
        t = np.median([th[g] for _, b, th in rows if g in b])
        if sum(y) >= MIN_POS and (len(y) - sum(y)) >= MIN_POS:
            res[g] = metrics(y, p, t)
    return res


def score_mds(include_mistrg=False):
    m = np.load(os.path.join(ROOT, "labels", "mds_bulk_equiv.npz"), allow_pickle=True)
    X = pd.DataFrame(m["X"].astype(float), index=[str(s) for s in m["samples"]],
                     columns=[str(g) for g in m["genes"]])
    grp = pd.Series([str(x) for x in m["group"]], index=X.index)
    kind = pd.Series([str(x) for x in m["kind"]], index=X.index)
    if not include_mistrg:
        keep = kind == "Primary"
        X, grp = X.loc[keep], grp.loc[keep]
    NEG = {"aged", "young", "Anemia"}                      # the only confident wild-types
    res = {}
    Z = {s: BP._z(BP._clog(BP._align(X.loc[s])), "sc") for s in X.index}   # sc-derived pseudobulk -> sc ref
    for gene, cats in [("SRSF2", ["SRSF2"]), ("U2AF1", ["U2AF1_S34", "U2AF1_Q157/R156"]),
                       ("RUNX1", ["RUNX1_LOF"]), ("TET2", ["TET2"])]:
        ids, y = [], []
        for s in X.index:
            g = grp[s]
            if gene in g:
                ids.append(s); y.append(1)
            elif g in NEG:
                ids.append(s); y.append(0)
            # everything else (empty Group, other genotypes) -> UNLABELLED, excluded
        if sum(y) < MIN_POS or (len(y) - sum(y)) < MIN_POS:
            res[gene] = {"skip": "only %d pos / %d neg" % (sum(y), len(y) - sum(y))}
            continue
        cats = [c for c in cats if c in BP.categories]
        if not cats:
            continue
        p = [max(BP.predict_one(c, Z[s], "sc")["probability"] for c in cats) for s in ids]
        t = float(np.median([BP.models[c]["threshold"] for c in cats]))
        res[gene] = metrics(np.array(y), p, t)
    return res


def table(title, res, note=""):
    print("\n" + "=" * 104)
    print(title + ("   [%s]" % note if note else ""))
    print("=" * 104)
    print("%-24s %5s %5s %6s | %6s | %6s %6s %5s | %5s %5s %5s" %
          ("category", "n", "n+", "prev", "AUROC", "AUPRC", "base", "lift", "F1", "prec", "rec"))
    print("-" * 104)
    ok = {k: v for k, v in res.items() if "skip" not in v}
    for c, m in sorted(ok.items(), key=lambda kv: -(kv[1]["auroc"] if kv[1]["auroc"] == kv[1]["auroc"] else -1)):
        print("%-24s %5d %5d %6.2f | %6.3f | %6.3f %6.3f %5.1fx | %5.2f %5.2f %5.2f" %
              (c, m["n"], m["n_pos"], m["prev"], m["auroc"], m["auprc"], m["prev"], m["lift"],
               m["f1"], m["prec"], m["rec"]))
    for c, m in res.items():
        if "skip" in m:
            print("%-24s  SKIPPED (%s)" % (c, m["skip"]))
    if ok:
        print("-" * 104)
        print("%-24s %5s %5d %6s | %6.3f | %6s %6s %5s | %5.2f %5.2f %5.2f" %
              ("MEAN", "", sum(m["n_pos"] for m in ok.values()), "",
               np.mean([m["auroc"] for m in ok.values()]), "", "", "",
               np.mean([m["f1"] for m in ok.values()]), np.mean([m["prec"] for m in ok.values()]),
               np.mean([m["rec"] for m in ok.values()])))


print("MODEL: %s | trained on: %s | modalities: BULK RNA ONLY (no imputation)"
      % (os.path.basename("bulk_mutation_predictor.pkl"), BP.meta.get("trained_on")))
print("categories: %d | mean BeatAML CV AUROC: %.3f" % (len(BP.categories), BP.summary()["mean_cv_auroc"]))

lg = score_bundle("lg_X", "lg_L", "lg_samples", "leucegene")
table("BULK HELD-OUT: Leucegene (n=367, external)", lg)

sc_ho = score_bundle("sc_X", "sc_L", "sc_samples", "sc", restrict=HOLD)
table("SINGLE-CELL HELD-OUT: sealed set (n=26)", sc_ho)

sc_all = score_bundle("sc_X", "sc_L", "sc_samples", "sc")
table("SINGLE-CELL: full atlas (n=387, all external to the bulk caller)", sc_all)

tp = score_trumpp()
table("SINGLE-CELL EXTERNAL: Trumpp (n=16)", tp, "gene-level truth (Table S4)")

mds = score_mds()
table("MDS (AML-precursor) cohort: Primary only", mds,
      "pos=Group genotype; neg=aged/young/Anemia ONLY; 28 empty-Group samples excluded as unlabelled")
