#!/usr/bin/env python3
"""Expanded evaluation implementing NS's three requested changes, together and separately.

  (2) UNION mutation set: gene-level units with >= MIN_POS (default 6) positives in EITHER the
      single-cell atlas or BeatAML — ~42 units rather than the 28 single-cell-only drivers.
  (3) DROP ungenotyped disease: specimens labelled AML/MDS/T-ALL whose entire mutation row is zero
      carry no genotype evidence, yet were being consumed as negatives for every driver. They are
      removed from training and evaluation. Controls are KEPT — their zeros are true negatives.
  (4) VAF FILTER: a BeatAML positive counts only if the reported t_vaf > VAF_MIN (default 0.10).
      Single-cell VAF is aggregate per mutation, NOT sample-resolved, so the filter cannot be applied
      on that side; this is stated rather than silently skipped.

Arms are cumulative so each change can be attributed:
  baseline      current configuration (>=8 positives, single-cell-only units, no sample filter)
  union         + expanded union mutation set at >=6 positives
  union_geno    + drop ungenotyped disease specimens
  union_geno_vaf + VAF>0.10 filter on BeatAML positives          <- all three changes

  bsub -q test -W 900 -M 72000 -R "rusage[mem=72000]" -o eu.log \
    /usr/local/anaconda3-2020/bin/python expanded_union_model.py
"""
import os, sys, json, csv, math, warnings, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
COMP = os.path.join(ROOT, "engine-code", "altanalyze3", "altanalyze3", "components")
sys.path.insert(0, os.path.join(ROOT, "engine-code", "altanalyze3")); sys.path.insert(0, COMP)
sys.path.insert(0, HERE)
from scipy.optimize import nnls
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from amlmm.context import build_context, Config
from amlmm import discovery as D, dataio, genetics, udon_features as UF
from amlmm.predictor import diff_select, _pct
from amlmm.bulk_predictor import BulkMutationPredictor
def log(m): print(m, flush=True)

OUT = os.path.join(ROOT, "deliverables", "expanded_union_model.json")
MIN_POS = int(os.environ.get("AMLMM_MIN_POS", "6"))
VAF_MIN = float(os.environ.get("AMLMM_VAF_MIN", "0.10"))
ATLAS8 = ["RNA","Composition","ADT","Lipid","Metabolite","GRN","LSC","Cell-comm"]
POOLABLE = ["BulkRNA","ADT","Metabolite","GRN"]

ctx = build_context(Config(run_id="single_modality"))
samples = ctx.tables["samples"]; dg = samples["donor_group"].astype(str); hold = set(ctx.holdout)
import bulk_features as BF

def load_block(mod):
    if mod == "RNA": return UF.canonical_rna(ctx)
    c = ctx.path("_sl_%s.pkl" % mod)
    if os.path.exists(c): return pd.read_pickle(c)
    if mod == "BulkRNA":
        b = BF.bulk_rna_matrix(ctx); b.to_pickle(c); return b
    if mod == "Composition": return D._sample_level_matrix(ctx, "composition", set(samples.index))
    if mod in ("ADT","GRN"): return dataio.sample_modality_matrix(ctx, mod)
    if mod in ("Lipid","Metabolite"): return dataio.sample_modality_matrix(ctx, mod, min_spearman=0.3)
    if mod == "Cell-comm": return dataio.cellcomm_matrix(ctx)
    if mod == "LSC":
        t = ctx.tables.get("lsc_calls")
        cols = [c for c in ["Prob_m-LSC","Prob_p+m-LSC","Prob_p-LSC","MaxProb"] if c in t.columns]
        return t[cols].apply(pd.to_numeric, errors="coerce")

BLK = {}
for m in ATLAS8 + ["BulkRNA"]:
    try:
        b = load_block(m).fillna(0.0)
        BLK[m] = b[~b.index.duplicated(keep="first")].loc[:, lambda d: ~d.columns.duplicated()]
        log("atlas %-11s %s" % (m, BLK[m].shape))
    except Exception as e:
        log("atlas %-11s FAILED %s" % (m, str(e)[:70]))

# ---------------- gene-level unit mapping ----------------
def unit(c):
    cl = str(c).lower(); u = str(c).upper()
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl: return "INV16"
    if "kmt2a" in cl: return "KMT2A"
    if "ITD" in u: return "FLT3-ITD"
    if "TKD" in u or "D835" in u or "I836" in u: return "FLT3-TKD"
    return str(c).split("_")[0].split("-")[0].upper()

# ---------------- single-cell labels, gene-level ----------------
rows = list(csv.reader(open(os.path.join(HERE, "mutation_matrix_explicit.tsv")), delimiter="\t"))
hdr = rows[0]; mcols = [(i, c) for i, c in enumerate(hdr) if c.startswith(("mut_", "cyto_"))]
def is_pos(v): return str(v).strip() not in ("", "0", "0.0", "nan", "NA")
sc_units = {}
for r in rows[1:]:
    d = {}
    for i, c in mcols:
        u = unit(c.replace("mut_", "").replace("cyto_", ""))
        d[u] = d.get(u, 0) or (1 if is_pos(r[i]) else 0)
    sc_units[r[0]] = d
# ungenotyped disease specimens
xl = pd.ExcelFile(os.path.join(ROOT, "labels", "AML_harmonized_metadata_v2_NYU2.xlsx"))
CM = xl.parse("Clinical_Metadata")
CM["key"] = CM["Dataset"].astype(str).str.strip() + "::" + CM["Sample"].astype(str).str.strip()
disease = {k: str(v) for k, v in zip(CM["key"], CM["disease_category"].astype(str))}
UNGENO = {k for k, d in sc_units.items()
          if sum(d.values()) == 0 and str(disease.get(k, "")).upper() in ("AML", "MDS", "T-ALL")}
log("\nungenotyped disease specimens excluded: %d" % len(UNGENO))

# ---------------- BeatAML labels, gene-level, with VAF gate ----------------
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
BP = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
e2s = {v: k for k, v in (BP.sym2ens or {}).items()}
ba_s = [str(s) for s in d["ba_samples"]]; gs = [str(g) for g in d["genes"]]
sy = [e2s.get(g) for g in gs]; kp = [i for i, s in enumerate(sy) if s]
BA_SYM = pd.DataFrame(d["ba_X"][:, kp].astype(float), index=ba_s, columns=[sy[i] for i in kp])
BA_SYM = BA_SYM.T.groupby(level=0).sum().T
cats = [str(c) for c in d["drivers"]]; baL = pd.DataFrame(d["ba_L"].astype(float), index=ba_s, columns=cats)
# per (sample, gene) VAF
vaf = {}
MUTF = os.path.join(ROOT, "data", "external", "beataml", "mutations.txt")
if os.path.exists(MUTF):
    with open(MUTF) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
        si, vi, gi = ci["dbgap_sample_id"], ci["t_vaf"], ci["symbol"]
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) <= max(si, vi, gi): continue
            try: v = float(c[vi])
            except Exception: continue
            k = (c[si], c[gi].upper())
            if v > vaf.get(k, -1): vaf[k] = v
log("BeatAML VAF pairs: %d" % len(vaf))

def ba_unit_labels(u, vaf_gate):
    cc = [c for c in cats if unit(c) == u]
    if not cc: return None
    sub = baL[cc]
    y = (sub.max(axis=1) == 1).astype(float)
    y[sub.isna().all(axis=1)] = np.nan
    if vaf_gate:                       # a positive survives only if its reported VAF exceeds the gate
        for s in y.index[y == 1]:
            sd = (s[:-1] + "D") if s.endswith("R") else s
            v = vaf.get((sd, u))
            if v is not None and v <= VAF_MIN:
                y.loc[s] = np.nan      # neither positive nor a trustworthy negative -> drop
    return y

# ---------------- unit sets ----------------
allidx = [s for s in BLK["RNA"].index]
def sc_pos(u, drop_ungeno):
    return sum(1 for s in allidx if s not in hold and (not drop_ungeno or s not in UNGENO)
               and sc_units.get(s, {}).get(u, 0) == 1)
UNITS_ALL = sorted({u for d0 in sc_units.values() for u in d0} | {unit(c) for c in cats})
def build_units(min_pos, drop_ungeno, vaf_gate):
    keep = []
    for u in UNITS_ALL:
        n_sc = sc_pos(u, drop_ungeno)
        yb = ba_unit_labels(u, vaf_gate)
        n_ba = int(np.nansum(yb.values)) if yb is not None else 0
        if max(n_sc, n_ba) >= min_pos and n_sc >= 3:      # need some sc positives to evaluate on sc
            keep.append((u, n_sc, n_ba))
    return keep

# ---------------- shared feature spaces for augmentation ----------------
from altanalyze3.components.rna2adt.api import Rna2AdtBundle
from altanalyze3.components.rna2grn.api import Rna2GrnBundle
from altanalyze3.components.rna2metabolite.api import load_bundle as load_metab
BAB = {"BulkRNA": BA_SYM}
for nm, ld in [("ADT", lambda: Rna2AdtBundle.load(os.path.join(COMP,"rna2adt","rna2adt_bm_bundle.pkl"))),
               ("GRN", lambda: Rna2GrnBundle.load(os.path.join(COMP,"rna2grn","rna2grn_bundle.pkl.gz"))),
               ("Metabolite", lambda: load_metab(os.path.join(COMP,"rna2metabolite","artifacts",
                                                              "rna2metabolite_aml_bundle.pkl.gz")))]:
    try:
        r = ld().predict_from_dataframe(BA_SYM)
        df = (r.predictions if hasattr(r, "predictions") else r).reindex(BA_SYM.index).fillna(0.0)
        BAB[nm] = df.loc[:, ~df.columns.duplicated()]
    except Exception as e: log("  beataml %s failed: %s" % (nm, str(e)[:60]))
def nkey(s):
    s = str(s)
    if s.startswith("Hu."): s = s[3:]
    return "".join(ch for ch in s.upper() if ch.isalnum())
SHARED, BCOL = {}, {}
for m in POOLABLE:
    if m not in BLK or m not in BAB: continue
    bm = {}
    for c in BAB[m].columns: bm.setdefault(nkey(c), c)
    pr, seen = [], set()
    for c in BLK[m].columns:
        k = nkey(c)
        if k in bm and k not in seen: seen.add(k); pr.append((c, bm[k]))
    if len(pr) >= 20:
        SHARED[m] = [a for a, _ in pr]; BCOL[m] = [b for _, b in pr]
log("poolable: %s" % {k: len(v) for k, v in SHARED.items()})

def auc(y, s):
    s = np.asarray(s, float); ok = ~np.isnan(s)
    return roc_auc_score(y[ok], s[ok]) if (ok.sum() >= 4 and len(set(y[ok])) == 2) else np.nan
def conf(y, call):
    y = np.asarray(y).astype(int); call = np.asarray(call).astype(int)
    tp = int(((call==1)&(y==1)).sum()); fp = int(((call==1)&(y==0)).sum())
    fn = int(((call==0)&(y==1)).sum()); tn = int(((call==0)&(y==0)).sum())
    se = tp/(tp+fn) if tp+fn else 0.0; sp = tn/(tn+fp) if tn+fp else 0.0
    pr = tp/(tp+fp) if tp+fp else 0.0
    return dict(tp=tp,fp=fp,fn=fn,tn=tn,sensitivity=round(se,4),specificity=round(sp,4),
                precision=round(pr,4), f1=round(2*pr*se/(pr+se),4) if pr+se else 0.0)
def nested(comb, y, grp):
    ng = min(3, len(set(grp)))
    if ng < 2: return conf(y, comb >= 0.5)
    call = np.zeros(len(y), int)
    for tri, tei in GroupKFold(ng).split(comb, y, grp):
        bt, bf = 0.5, -1.0
        for t in np.unique(comb[tri]):
            c = conf(y[tri], comb[tri] >= t)
            if c["f1"] > bf: bf, bt = c["f1"], float(t)
        call[tei] = (comb[tei] >= bt).astype(int)
    return conf(y, call)

def run(units, drop_ungeno, vaf_gate, tag):
    out = {}
    for u, n_sc, n_ba in units:
        yba = ba_unit_labels(u, vaf_gate)
        oof_p = {}
        for mod in ATLAS8 + ["BulkRNA"]:
            if mod not in BLK: continue
            poolable = (mod in SHARED) and yba is not None
            cols = SHARED[mod] if poolable else list(BLK[mod].columns)
            A = BLK[mod][cols] if poolable else BLK[mod]
            ids = [s for s in A.index if s not in hold and (not drop_ungeno or s not in UNGENO)]
            y = np.array([sc_units.get(s, {}).get(u, 0) for s in ids])
            if y.sum() < 3 or (y == 0).sum() < 5: continue
            grp = dg.loc[ids].values; Xa = A.loc[ids].values
            extra = None
            if poolable:
                bsel = yba.dropna().index
                yb = yba.loc[bsel].values.astype(int)
                if yb.sum() >= 3 and (yb == 0).sum() >= 5:
                    extra = (BAB[mod].loc[bsel, BCOL[mod]].values, yb)
            oof = np.full(len(ids), np.nan); ng = min(3, len(set(grp)))
            if ng < 2: continue
            for tri, vai in GroupKFold(ng).split(Xa, y, grp):
                if len(set(y[tri])) < 2: continue
                sc_ = StandardScaler().fit(Xa[tri]); Ztr = sc_.transform(Xa[tri]); Zva = sc_.transform(Xa[vai])
                Xf, yf = Ztr, y[tri]
                if extra is not None:
                    Xb, yb2 = extra; keep = Xa[tri].std(0) > 0
                    Xb2 = Xb[:, keep] if Xb.shape[1] == Xa.shape[1] else Xb
                    mb, sb = Xb2.mean(0), Xb2.std(0); sb[sb == 0] = 1.0
                    Zb = (Xb2 - mb)/sb
                    if Zb.shape[1] == Ztr.shape[1]:
                        Xf = np.vstack([Ztr, Zb]); yf = np.concatenate([y[tri], yb2])
                k2 = Xf.std(0) > 0
                if k2.sum() < 2: continue
                sel = diff_select(Xf[:, k2], yf, 500)
                oof[vai] = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000).fit(
                    Xf[:, k2][:, sel], yf).decision_function(Zva[:, k2][:, sel])
            ok = ~np.isnan(oof)
            if ok.sum() < 4: continue
            p = np.full(len(ids), np.nan); p[ok] = _pct(np.sort(oof[ok]), oof[ok])
            a0 = auc(y[ok], p[ok])
            if a0 == a0 and a0 < 0.5: p = 1 - p
            oof_p[mod] = dict(zip(ids, p))
        if not oof_p: continue
        common = sorted(set.intersection(*[set(v) for v in oof_p.values()]))
        yco = np.array([sc_units.get(s, {}).get(u, 0) for s in common])
        if len(common) < 8 or len(set(yco)) != 2 or yco.sum() < 3: continue
        gco = dg.loc[common].values
        cols_ = list(oof_p); O = np.column_stack([[oof_p[m][s] for s in common] for m in cols_])
        Am = np.vstack([O, np.eye(O.shape[1])]); b = np.concatenate([yco.astype(float), np.zeros(O.shape[1])])
        w, _ = nnls(Am, b); w = w/(w.sum() or 1.0); comb = O @ w
        a1 = auc(yco, comb); nm = nested(comb, yco, gco)
        nm.update(unit=u, auroc=round(float(a1),4) if a1==a1 else None,
                  n_pos_sc=int(yco.sum()), n_pos_beataml=n_ba, n_eval=len(common))
        out[u] = nm
    ok = [v for v in out.values() if v["auroc"] is not None]
    summ = {"n_units": len(ok),
            "mean_auroc": round(float(np.mean([v["auroc"] for v in ok])), 4) if ok else None,
            "mean_sensitivity": round(float(np.mean([v["sensitivity"] for v in ok])), 4) if ok else None,
            "mean_specificity": round(float(np.mean([v["specificity"] for v in ok])), 4) if ok else None,
            "mean_f1": round(float(np.mean([v["f1"] for v in ok])), 4) if ok else None}
    tp=sum(v["tp"] for v in ok); fp=sum(v["fp"] for v in ok); fn=sum(v["fn"] for v in ok); tn=sum(v["tn"] for v in ok)
    summ["pooled_sensitivity"]=round(tp/(tp+fn),4) if tp+fn else None
    summ["pooled_specificity"]=round(tn/(tn+fp),4) if tn+fp else None
    log("  %-16s units=%-3d AUROC %.3f sens %.3f spec %.3f (pooled sens %.3f)" %
        (tag, summ["n_units"], summ["mean_auroc"] or 0, summ["mean_sensitivity"] or 0,
         summ["mean_specificity"] or 0, summ["pooled_sensitivity"] or 0))
    return out, summ

ARMS = [("baseline",        build_units(8, False, False), False, False),
        ("union",           build_units(MIN_POS, False, False), False, False),
        ("union_geno",      build_units(MIN_POS, True,  False), True,  False),
        ("union_geno_vaf",  build_units(MIN_POS, True,  True),  True,  True)]
log("\nunit counts per arm: %s" % {t: len(u) for t, u, _, _ in ARMS})
res, summ = {}, {}
for tag, units, dgn, vg in ARMS:
    res[tag], summ[tag] = run(units, dgn, vg, tag)

json.dump({"generated": time.strftime("%Y-%m-%d %H:%M"), "min_pos": MIN_POS, "vaf_min": VAF_MIN,
           "n_ungenotyped_excluded": len(UNGENO),
           "ungenotyped_examples": sorted(UNGENO)[:20],
           "caveat": "Single-cell VAF is aggregate per mutation, not sample-resolved, so the VAF gate "
                     "is applied to BeatAML positives only.",
           "summary": summ, "per_unit": res}, open(OUT, "w"), indent=1)
log("\nwrote %s" % OUT)
log("EXPANDED UNION OK")
