#!/usr/bin/env python3
"""NS's proposal: augment the single-cell multimodal training set with BeatAML, projected into the
same modality space by the rna2* imputers — does more training data capture more mutations?

This is a DIFFERENT question from beataml_impute_experiment.py. That one asked "does imputation add
information WITHIN BeatAML?" (answer: no, it equals a random nonlinear RNA transform). This one asks
"can imputation act as the BRIDGE that lets 707 bulk samples augment 357 single-cell samples?" —
the imputers give both cohorts the SAME ADT/Lipid/Metabolite/GRN feature names, so they can be pooled.

SHARED blocks (both cohorts): BulkRNA gene space + ADT + Lipid + Metabolite + GRN.
(The deployed atlas RNA block is UDON markers + program fractions — not a gene space — so it cannot be
pooled; Composition / LSC / Cell-comm do not exist for bulk at all. Those stay atlas-only.)

ARMS, all evaluated on the SAME donor-grouped CV folds over ATLAS samples (BeatAML only ever joins the
TRAINING side of a fold, never the test side):
  atlas_only            shared blocks trained on atlas alone            <- fair baseline
  atlas_plus_beataml    shared blocks trained on atlas + BeatAML
  atlas_plus_beataml_vaf same, with BeatAML positives weighted by VAF (clonal weighted up)

Each cohort is z-scored against ITS OWN statistics before pooling (the standard domain-shift fix, and
what the deployed bulk caller already does with cohort-matched score references).

  bsub -q test -W 400 -M 56000 -R "rusage[mem=56000]" -o aug.log \
    /usr/local/anaconda3-2020/bin/python beataml_augment_experiment.py
"""
import os, sys, json, warnings, time, csv
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
COMP = os.path.join(ROOT, "engine-code", "altanalyze3", "altanalyze3", "components")
ALT = os.path.join(ROOT, "engine-code", "altanalyze3")
sys.path.insert(0, ALT); sys.path.insert(0, COMP); sys.path.insert(0, HERE)
from scipy.optimize import nnls
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from amlmm.context import build_context, Config
from amlmm import discovery as D, dataio, genetics
from amlmm.predictor import diff_select, _pct
from amlmm.bulk_predictor import BulkMutationPredictor
def log(m): print(m, flush=True)

OUT = os.path.join(ROOT, "deliverables", "beataml_augment_experiment.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
ctx = build_context(Config(run_id="single_modality"))
samples = ctx.tables["samples"]; dg = samples["donor_group"].astype(str); hold = set(ctx.holdout)

# ------------------------------------------------- atlas blocks (shared space only)
import bulk_features as BF
def atlas_block(mod):
    c = ctx.path("_sl_%s.pkl" % mod)
    if mod == "BulkRNA":
        if os.path.exists(c): return pd.read_pickle(c)
        b = BF.bulk_rna_matrix(ctx); b.to_pickle(c); return b
    if os.path.exists(c): return pd.read_pickle(c)
    if mod in ("ADT", "GRN"): return dataio.sample_modality_matrix(ctx, mod)
    return dataio.sample_modality_matrix(ctx, mod, min_spearman=0.3)

ATL = {}
for m in ["BulkRNA", "ADT", "Lipid", "Metabolite", "GRN"]:
    try:
        b = atlas_block(m).fillna(0.0); ATL[m] = b[~b.index.duplicated(keep="first")]
        log("atlas  %-11s %s" % (m, ATL[m].shape))
    except Exception as e:
        log("atlas  %-11s FAILED %s" % (m, str(e)[:90]))

# ------------------------------------------------- BeatAML: RNA + imputed modalities
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
BP = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
ens2sym = {v: k for k, v in (BP.sym2ens or {}).items()}
ba_s = [str(s) for s in d["ba_samples"]]; genes_ens = [str(g) for g in d["genes"]]
BA_ENS = pd.DataFrame(d["ba_X"].astype(float), index=ba_s, columns=genes_ens)
syms = [ens2sym.get(g) for g in genes_ens]; keep = [i for i, s in enumerate(syms) if s]
BA_SYM = pd.DataFrame(d["ba_X"][:, keep].astype(float), index=ba_s, columns=[syms[i] for i in keep])
BA_SYM = BA_SYM.T.groupby(level=0).sum().T
log("BeatAML RNA: %s (ENSG) / %s (symbols)" % (BA_ENS.shape, BA_SYM.shape))

from altanalyze3.components.rna2adt.api import Rna2AdtBundle
from altanalyze3.components.rna2grn.api import Rna2GrnBundle
from altanalyze3.components.rna2lipid.api import Rna2LipidBundle
from altanalyze3.components.rna2metabolite.api import load_bundle as load_metab
BAB = {}
for name, loader in [("ADT", lambda: Rna2AdtBundle.load(os.path.join(COMP, "rna2adt", "rna2adt_bm_bundle.pkl"))),
                     ("GRN", lambda: Rna2GrnBundle.load(os.path.join(COMP, "rna2grn", "rna2grn_bundle.pkl.gz"))),
                     ("Lipid", lambda: Rna2LipidBundle.load(os.path.join(COMP, "rna2lipid", "otherelastic_multitask_try.pkl"))),
                     ("Metabolite", lambda: load_metab(os.path.join(COMP, "rna2metabolite", "artifacts", "rna2metabolite_aml_bundle.pkl.gz")))]:
    try:
        res = loader().predict_from_dataframe(BA_SYM)
        df = (res.predictions if hasattr(res, "predictions") else res).reindex(BA_SYM.index).fillna(0.0)
        BAB[name] = df; log("BeatAML imputed %-11s %s" % (name, df.shape))
    except Exception as e:
        log("BeatAML imputed %-11s FAILED %s" % (name, str(e)[:110]))
# the atlas BulkRNA block is in gene SYMBOL space (var index of the RNA h5ad), so BeatAML must be too
BAB["BulkRNA"] = BA_SYM

# ------------------------------------------------- shared feature spaces
# match on a normalised name (case/punctuation-insensitive): ADT antibody and lipid/metabolite species
# names are formatted differently by the atlas pipeline vs the rna2* bundles, so exact matching drops
# whole modalities that are in fact the same feature set.
def nkey(s): return "".join(ch for ch in str(s).upper() if ch.isalnum())
SHARED, REN = {}, {}
for m in list(ATL):
    if m not in BAB: continue
    ATL[m] = ATL[m].loc[:, ~ATL[m].columns.duplicated()]          # atlas gene symbols contain duplicates
    BAB[m] = BAB[m].loc[:, ~BAB[m].columns.duplicated()]
    bmap = {}
    for c in BAB[m].columns:
        bmap.setdefault(nkey(c), c)
    pairs, _seen = [], set()
    for c in ATL[m].columns:                                       # unique on BOTH sides, order preserved
        k = nkey(c)
        if k in bmap and k not in _seen:
            _seen.add(k); pairs.append((c, bmap[k]))
    if len(pairs) >= 20:
        SHARED[m] = [a for a, _ in pairs]; REN[m] = [b for _, b in pairs]      # positionally paired columns
        log("shared %-11s %d features (atlas %d / beataml %d)" % (m, len(pairs), ATL[m].shape[1], BAB[m].shape[1]))
    else:
        log("shared %-11s only %d overlapping features -> DROPPED" % (m, len(pairs)))
MODS = list(SHARED)

# ------------------------------------------------- labels
M = ctx.tables.get("mutations") or genetics.build_mutation_matrix(ctx)
_m01 = {"present": 1.0, "absent": 0.0}
cats = [str(c) for c in d["drivers"]]; baL = pd.DataFrame(d["ba_L"].astype(float), index=ba_s, columns=cats)
def gene_of(c):
    cl = str(c).lower()
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl: return "INV16"
    if "kmt2a" in cl: return "KMT2A"
    return str(c).split("_")[0].split("-")[0].upper()
def ba_labels_for(short):
    """BeatAML gene-level label for an atlas mutation short-name (present iff any variant category is 1)."""
    s = short.upper().replace("-", "").replace("_", "")
    cc = [c for c in cats if gene_of(c).replace("-", "").replace("_", "") == s]
    if short.upper().startswith("FLT3-ITD"): cc = [c for c in cats if "ITD" in c.upper()]
    if short.upper().startswith("FLT3-TKD"): cc = [c for c in cats if "TKD" in c.upper()]
    if not cc: return None
    sub = baL[cc]
    y = (sub.max(axis=1) == 1).astype(float)
    y[sub.isna().all(axis=1)] = np.nan
    return y

# BeatAML VAF per (sample, gene) for the weighted arm
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
log("VAF pairs parsed: %d" % len(vaf))

def auc(y, s):
    s = np.asarray(s, float); ok = ~np.isnan(s)
    return roc_auc_score(y[ok], s[ok]) if (ok.sum() >= 4 and len(set(y[ok])) == 2) else np.nan
def conf(y, call):
    y = np.asarray(y).astype(int); call = np.asarray(call).astype(int)
    tp = int(((call == 1) & (y == 1)).sum()); fp = int(((call == 1) & (y == 0)).sum())
    fn = int(((call == 0) & (y == 1)).sum()); tn = int(((call == 0) & (y == 0)).sum())
    se = tp/(tp+fn) if tp+fn else 0.0; sp = tn/(tn+fp) if tn+fp else 0.0
    pr = tp/(tp+fp) if tp+fp else 0.0; f1 = 2*pr*se/(pr+se) if pr+se else 0.0
    return dict(tp=tp, fp=fp, tn=tn, fn=fn, sensitivity=round(se,4), specificity=round(sp,4),
                precision=round(pr,4), f1=round(f1,4))
def f1max(comb, y):
    bt, bf = 0.5, -1.0
    for t in np.unique(comb):
        c = conf(y, comb >= t)
        if c["f1"] > bf: bf, bt = c["f1"], float(t)
    return bt

MUTS = []
for f in sorted(c for c in M.columns if str(c).startswith(("mut_", "cyto_"))):
    ym = D.labels_for_field(ctx, f).map(_m01)
    tr = [s for s in ATL["BulkRNA"].index if pd.notna(ym.get(s)) and s not in hold]
    yv = np.array([int(ym[s]) for s in tr])
    if (yv == 1).sum() >= 8 and (yv == 0).sum() >= 8:
        MUTS.append(f)
log("trainable mutations: %d | shared modalities: %s" % (len(MUTS), MODS))

ARMS = ["atlas_only", "atlas_plus_beataml", "atlas_plus_beataml_vaf"]
res = {a: {} for a in ARMS}
for mflag in MUTS:
    short = mflag.replace("mut_", "").replace("cyto_", "")
    yall = D._labels_for_field_raw(ctx, mflag).map(_m01)
    ym = D.labels_for_field(ctx, mflag).map(_m01)
    yba = ba_labels_for(short)
    nba = int(np.nansum(yba.values)) if yba is not None else 0
    for arm in ARMS:
        use_ba = arm != "atlas_only"
        if use_ba and (yba is None or nba < 3): continue
        oof_p, mod_names = {}, []
        for mod in MODS:
            cols = SHARED[mod]
            A = ATL[mod][cols]
            ids = [s for s in A.index if pd.notna(ym.get(s)) and s not in hold]
            y = np.array([int(yall[s]) for s in ids])
            if (y == 1).sum() < 5 or (y == 0).sum() < 5: continue
            grp = dg.loc[ids].values
            Xa = A.loc[ids].values
            if use_ba:
                bsel = yba.dropna().index
                # BeatAML columns are positionally paired to `cols` (SHARED), so slice by position
                Xb = BAB[mod].loc[bsel, REN[mod]].fillna(0.0).values
                yb = yba.loc[bsel].values.astype(int)
                # cohort-wise z-scoring (each cohort standardised on its OWN stats), then pool
                mb, sb = Xb.mean(0), Xb.std(0); sb[sb == 0] = 1.0
                Zb = (Xb - mb) / sb
                wb = np.ones(len(yb))
                if arm.endswith("_vaf"):
                    g = gene_of(short)
                    for i, s in enumerate(bsel):
                        if yb[i] == 1:
                            sd = (s[:-1] + "D") if s.endswith("R") else s
                            v = vaf.get((sd, g))
                            wb[i] = 0.5 + 1.5 * float(v) if v is not None else 1.0   # clonal weighted up
            oof = np.full(len(ids), np.nan)
            ng = min(3, len(set(grp)))
            if ng < 2: continue
            for tri, vai in GroupKFold(ng).split(Xa, y, grp):
                if len(set(y[tri])) < 2: continue
                sc = StandardScaler().fit(Xa[tri]); Za_tr = sc.transform(Xa[tri]); Za_va = sc.transform(Xa[vai])
                Xtr, ytr = Za_tr, y[tri]; w = np.ones(len(ytr))
                if use_ba:
                    Xtr = np.vstack([Za_tr, Zb]); ytr = np.concatenate([y[tri], yb]); w = np.concatenate([np.ones(len(y[tri])), wb])
                keep2 = Xtr.std(0) > 0
                if keep2.sum() < 2: continue
                sel = diff_select(Xtr[:, keep2], ytr, 500)
                est = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000)
                est.fit(Xtr[:, keep2][:, sel], ytr, sample_weight=w)
                oof[vai] = est.decision_function(Za_va[:, keep2][:, sel])
            ok = ~np.isnan(oof)
            if ok.sum() < 4: continue
            p = np.full(len(ids), np.nan); p[ok] = _pct(np.sort(oof[ok]), oof[ok])
            a0 = auc(y[ok], p[ok])
            if a0 == a0 and a0 < 0.5: p = 1 - p
            oof_p[mod] = dict(zip(ids, p)); mod_names.append(mod)
        if not oof_p: continue
        common = sorted(set.intersection(*[set(oof_p[m]) for m in oof_p]))
        yco = np.array([int(yall[s]) for s in common])
        if len(common) < 8 or len(set(yco)) != 2: continue
        cols_ = list(oof_p); O = np.column_stack([[oof_p[m][s] for s in common] for m in cols_])
        if len(cols_) == 1:
            comb = O[:, 0]
        else:
            Amat = np.vstack([O, np.eye(O.shape[1])]); b = np.concatenate([yco.astype(float), np.zeros(O.shape[1])])
            w2, _ = nnls(Amat, b); w2 = w2 / (w2.sum() or 1.0); comb = O @ w2
        a1 = auc(yco, comb); t = f1max(comb, yco); mm = conf(yco, comb >= t)
        mm.update(mutation=short, auroc=round(float(a1), 4) if a1 == a1 else None,
                  n=len(common), n_pos_atlas=int((yco == 1).sum()), n_pos_beataml=nba)
        res[arm][short] = mm
    log("  %-14s atlas n+=%-3d beataml n+=%-4d | %s" % (short,
        res["atlas_only"].get(short, {}).get("n_pos_atlas", 0), nba,
        " ".join("%s=%.3f" % (a.replace("atlas_plus_beataml", "aug").replace("atlas_only", "base"),
                              res[a][short]["auroc"]) for a in ARMS if short in res[a] and res[a][short]["auroc"] is not None)))

shared_m = sorted(set.intersection(*[set(res[a]) for a in ARMS])) if all(res[a] for a in ARMS) else sorted(res["atlas_only"])
summ = {}
for a in ARMS:
    ms = [m for m in shared_m if m in res[a]]
    if not ms: continue
    summ[a] = {k: round(float(np.mean([res[a][m][k] for m in ms if res[a][m].get(k) is not None])), 4)
               for k in ("auroc", "sensitivity", "specificity", "f1")}
    summ[a]["n_mutations"] = len(ms)
payload = {"generated": time.strftime("%Y-%m-%d %H:%M"), "shared_modalities": MODS,
           "shared_feature_counts": {m: len(SHARED[m]) for m in MODS},
           "n_beataml": len(ba_s), "summary": summ, "per_mutation": {m: {a: res[a].get(m) for a in ARMS} for m in shared_m}}
json.dump(payload, open(OUT, "w"), indent=1)
log("\n=== AUGMENTATION RESULT (atlas donor-grouped CV, %d mutations) ===" % len(shared_m))
for a in ARMS:
    if a in summ:
        s = summ[a]; log("  %-24s AUROC %.3f  sens %.3f  spec %.3f  F1 %.3f" % (a, s["auroc"], s["sensitivity"], s["specificity"], s["f1"]))
log("wrote %s" % OUT)
log("AUGMENT EXPERIMENT OK")
