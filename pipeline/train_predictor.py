#!/usr/bin/env python3
"""Train the deployable MutationPredictor on the cohort + emit per-sample board reports.

Self-contained (no external OOF dependency): per (mutation, modality) fit StandardScaler -> differential
top-500 -> LinearSVC on that modality's own non-holdout labelled samples, AND run a donor-grouped 3-fold
CV to get honest OOF percentiles. Per mutation, optimised modality weights = ridge-NNLS of the per-modality
OOF percentiles (aligned on the common train samples) onto truth, floored to the best single modality and
falling back to a uniform mix of the informative (CV-AUC>0.5) modalities so NO mutation trains empty.
Persist -> mutation_predictor.pkl. Then score the sealed held-out (predicted vs known) + control gate and
write one patient_report.json per held-out sample for the board.  Run on an LSF compute node.
"""
import os, sys, json, pickle, warnings, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.optimize import nnls
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from amlmm.context import build_context, Config
from amlmm import discovery as D, dataio, genetics, udon_features as UF, therapy
from amlmm.predictor import MutationPredictor, ModalityModel, diff_select, _pct
try:
    import control_gate as CG
except Exception:
    CG = None

ctx = build_context(Config(run_id="single_modality"))
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS_ROOT = os.path.dirname(os.path.dirname(ctx.path("x")))     # <base>/runs
os.makedirs(RUNS_ROOT, exist_ok=True)
samples = ctx.tables["samples"]; dg = samples["donor_group"].astype(str); hold = set(ctx.holdout)
MODS = ["RNA", "Composition", "ADT", "Lipid", "Metabolite", "GRN", "LSC", "Cell-comm"]
_ENVM = os.environ.get("AMLMM_MODS", "").strip()          # comma-list overrides MODS (e.g. bulk-RNA experiments)
if _ENVM:
    MODS = [m.strip() for m in _ENVM.split(",") if m.strip()]
TAG = os.environ.get("AMLMM_TAG", "").strip()             # output suffix; when set, skip the deployed board reports
SFX = ("_" + TAG) if TAG else ""
def log(m): print(m, flush=True)


def load_block(mod):
    if mod == "RNA":
        return UF.canonical_rna(ctx)
    c = ctx.path("_sl_%s.pkl" % mod)
    if os.path.exists(c):
        return pd.read_pickle(c)
    if mod == "BulkRNA":                                   # whole-sample bulk-equivalent (sum all cell-states -> CP10k+log1p)
        import bulk_features as BF
        b = BF.bulk_rna_matrix(ctx); b.to_pickle(c); return b
    if mod == "Composition":
        return D._sample_level_matrix(ctx, "composition", set(samples.index))
    if mod in ("ADT", "GRN"):
        return dataio.sample_modality_matrix(ctx, mod)
    if mod in ("Lipid", "Metabolite"):
        return dataio.sample_modality_matrix(ctx, mod, min_spearman=0.3)
    if mod == "Cell-comm":
        return dataio.cellcomm_matrix(ctx)
    if mod == "LSC":
        t = ctx.tables.get("lsc_calls"); cols = [c for c in ["Prob_m-LSC", "Prob_p+m-LSC", "Prob_p-LSC", "MaxProb"] if c in t.columns]
        return t[cols].apply(pd.to_numeric, errors="coerce")


BLK = {}
for m in MODS:
    try:
        b = load_block(m).fillna(0.0); BLK[m] = b[~b.index.duplicated(keep="first")]
    except Exception as e:
        log("skip %s: %s" % (m, e))
MODS = [m for m in MODS if m in BLK]
log("modalities: %s" % ", ".join(MODS))

# ---- optional BeatAML augmentation (AMLMM_AUGMENT=1) ------------------------------------------
# The rna2* imputers give BeatAML the SAME ADT/Metabolite/GRN feature names as the atlas, and
# bulk_features puts both in the same gene space — so 707 bulk samples can be pooled into the
# TRAINING side of the poolable blocks. Validated at AUROC 0.889 -> 0.908 (production_fused_model.py).
AUGMENT = os.environ.get("AMLMM_AUGMENT") == "1"
SHARED, BCOL, BAB, BA_LABELS = {}, {}, {}, (lambda s: None)
if AUGMENT:
    try:
        _ROOT = os.path.dirname(HERE)
        _COMP = os.path.join(_ROOT, "engine-code", "altanalyze3", "altanalyze3", "components")
        sys.path.insert(0, os.path.join(_ROOT, "engine-code", "altanalyze3")); sys.path.insert(0, _COMP)
        import bulk_features as _BF                                     # noqa: F401  (atlas BulkRNA block)
        from amlmm.bulk_predictor import BulkMutationPredictor as _BMP
        _d = np.load(os.path.join(os.path.dirname(_ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
        _bp = _BMP.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
        _e2s = {v: k for k, v in (_bp.sym2ens or {}).items()}
        _ba = [str(s) for s in _d["ba_samples"]]; _gs = [str(g) for g in _d["genes"]]
        _sy = [_e2s.get(g) for g in _gs]; _kp = [i for i, s in enumerate(_sy) if s]
        _BASYM = pd.DataFrame(_d["ba_X"][:, _kp].astype(float), index=_ba, columns=[_sy[i] for i in _kp])
        _BASYM = _BASYM.T.groupby(level=0).sum().T
        from altanalyze3.components.rna2adt.api import Rna2AdtBundle as _AdtB
        from altanalyze3.components.rna2grn.api import Rna2GrnBundle as _GrnB
        from altanalyze3.components.rna2metabolite.api import load_bundle as _load_metab
        BAB["BulkRNA"] = _BASYM
        for _n, _ld in [("ADT", lambda: _AdtB.load(os.path.join(_COMP, "rna2adt", "rna2adt_bm_bundle.pkl"))),
                        ("GRN", lambda: _GrnB.load(os.path.join(_COMP, "rna2grn", "rna2grn_bundle.pkl.gz"))),
                        ("Metabolite", lambda: _load_metab(os.path.join(_COMP, "rna2metabolite", "artifacts",
                                                                        "rna2metabolite_aml_bundle.pkl.gz")))]:
            try:
                _r = _ld().predict_from_dataframe(_BASYM)
                _df = (_r.predictions if hasattr(_r, "predictions") else _r).reindex(_BASYM.index).fillna(0.0)
                BAB[_n] = _df.loc[:, ~_df.columns.duplicated()]
            except Exception as _e:
                log("  beataml %s failed: %s" % (_n, str(_e)[:80]))

        def _nkey(s):
            s = str(s)
            if s.startswith("Hu."):                                     # rna2adt prefixes every antibody
                s = s[3:]
            return "".join(ch for ch in s.upper() if ch.isalnum())
        for _m in ["BulkRNA", "ADT", "Metabolite", "GRN"]:
            if _m not in BAB:
                continue
            if _m not in BLK:                                           # BulkRNA is not a deployed modality
                try:
                    _b = load_block(_m).fillna(0.0); BLK[_m] = _b[~_b.index.duplicated(keep="first")]
                    MODS.append(_m)
                except Exception:
                    continue
            BLK[_m] = BLK[_m].loc[:, ~BLK[_m].columns.duplicated()]
            _bm = {}
            for _c in BAB[_m].columns:
                _bm.setdefault(_nkey(_c), _c)
            _pairs, _seen = [], set()
            for _c in BLK[_m].columns:
                _k = _nkey(_c)
                if _k in _bm and _k not in _seen:
                    _seen.add(_k); _pairs.append((_c, _bm[_k]))
            if len(_pairs) >= 20:
                SHARED[_m] = [a for a, _ in _pairs]; BCOL[_m] = [b for _, b in _pairs]
                BLK[_m] = BLK[_m][SHARED[_m]]                           # train/deploy on the shared space
                log("  augmentable %-11s %d shared features" % (_m, len(_pairs)))
        _cats = [str(c) for c in _d["drivers"]]
        _baL = pd.DataFrame(_d["ba_L"].astype(float), index=_ba, columns=_cats)
        def _gene_of(c):
            cl = str(c).lower()
            if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl: return "INV16"
            if "kmt2a" in cl: return "KMT2A"
            return str(c).split("_")[0].split("-")[0].upper()
        def BA_LABELS(short):                                           # noqa: F811
            s = short.upper().replace("-", "").replace("_", "")
            cc = [c for c in _cats if _gene_of(c).replace("-", "").replace("_", "") == s]
            if short.upper().startswith("FLT3-ITD"): cc = [c for c in _cats if "ITD" in c.upper()]
            if short.upper().startswith("FLT3-TKD"): cc = [c for c in _cats if "TKD" in c.upper()]
            if not cc: return None
            sub = _baL[cc]; y = (sub.max(axis=1) == 1).astype(float); y[sub.isna().all(axis=1)] = np.nan
            return y
        log("AUGMENT ON: BeatAML n=%d, poolable blocks %s" % (len(_ba), list(SHARED)))
    except Exception as _e:
        log("AUGMENT requested but unavailable (%s) — training atlas-only" % str(_e)[:110])
        AUGMENT = False

M = ctx.tables.get("mutations") or genetics.build_mutation_matrix(ctx)
_m01 = {"present": 1.0, "absent": 0.0}
def auc(y, s):
    s = np.asarray(s, float); ok = ~np.isnan(s)
    return roc_auc_score(y[ok], s[ok]) if (ok.sum() >= 4 and len(set(y[ok])) == 2) else np.nan

allidx = set(BLK["RNA"].index) if "RNA" in BLK else set().union(*[set(BLK[m].index) for m in MODS])
MUTS = []
for f in sorted(c for c in M.columns if str(c).startswith(("mut_", "cyto_"))):
    ym = D.labels_for_field(ctx, f).map(_m01)
    tr = [s for s in allidx if pd.notna(ym.get(s)) and s not in hold]
    yv = np.array([int(ym[s]) for s in tr])
    if (yv == 1).sum() >= 8 and (yv == 0).sum() >= 8:        # well-powered only
        MUTS.append(f)
log("trainable mutations: %d" % len(MUTS))

P = MutationPredictor(); P.modalities = MODS
P.mutations = [m.replace("mut_", "").replace("cyto_", "") for m in MUTS]
NAME = {m.replace("mut_", "").replace("cyto_", ""): m for m in MUTS}
heldout_pred = {}                                              # sample -> {mut_short: prob}


def _aug_rows(extra, keep):
    """BeatAML rows for one modality, z-scored on BeatAML's OWN statistics (cohort-matched
    normalisation) and masked to the same columns the atlas fold kept."""
    Xb, yb = extra
    Xb = Xb[:, keep]
    mb, sb = Xb.mean(0), Xb.std(0); sb[sb == 0] = 1.0
    return (Xb - mb) / sb, yb


def cv_oof(B, ids, y, grp, extra=None):
    """donor-grouped 3-fold OOF percentile scores for linSVM on this modality block.
    `extra` = (BeatAML matrix aligned to B.columns, labels); it joins only the TRAINING side of a fold."""
    X = B.loc[ids].values; keep = X.std(0) > 0
    if keep.sum() < 2 or len(set(grp)) < 2:
        return None
    oof = np.full(len(ids), np.nan)
    for tri, vai in GroupKFold(min(3, len(set(grp)))).split(X, y, grp):
        if len(set(y[tri])) < 2:
            continue
        sc = StandardScaler().fit(X[tri][:, keep]); Ztr = sc.transform(X[tri][:, keep]); Zva = sc.transform(X[vai][:, keep])
        Xf, yf = Ztr, y[tri]
        if extra is not None:
            Zb, yb = _aug_rows(extra, keep)
            Xf = np.vstack([Ztr, Zb]); yf = np.concatenate([y[tri], yb])
        sel = diff_select(Xf, yf, 500)
        d = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000).fit(Xf[:, sel], yf).decision_function(Zva[:, sel])
        oof[vai] = d
    ok = ~np.isnan(oof)
    if ok.sum() < 4:
        return None
    p = np.full(len(ids), np.nan); p[ok] = _pct(np.sort(oof[ok]), oof[ok])
    a = auc(y[ok], p[ok])
    if a == a and a < 0.5:
        p = 1 - p
    return p, (max(a, 1 - a) if a == a else 0.5)


for mflag in MUTS:
    mshort = mflag.replace("mut_", "").replace("cyto_", "")
    yall = D._labels_for_field_raw(ctx, mflag).map(_m01); ym = D.labels_for_field(ctx, mflag).map(_m01)
    oof_p, cvauc, te_p, mod_names = {}, {}, {}, []
    yba = BA_LABELS(mshort) if AUGMENT else None                 # BeatAML gene-level labels, or None
    for mod in MODS:
        B = BLK[mod]
        tr = [s for s in B.index if pd.notna(ym.get(s)) and s not in hold]
        te = [s for s in B.index if s in hold and pd.notna(yall.get(s))]
        ytr = np.array([int(yall[s]) for s in tr])
        if (ytr == 1).sum() < 5 or (ytr == 0).sum() < 5:
            continue
        # BeatAML augmentation for the poolable blocks (bulk rows added to TRAINING only)
        extra = None
        if AUGMENT and yba is not None and mod in SHARED:
            bsel = yba.dropna().index
            yb = yba.loc[bsel].values.astype(int)
            if yb.sum() >= 5 and (yb == 0).sum() >= 5:
                extra = (BAB[mod].loc[bsel, BCOL[mod]].values, yb)
        Xtr = B.loc[tr].values; keep = Xtr.std(0) > 0
        sc = StandardScaler().fit(Xtr[:, keep]); Ztr = sc.transform(Xtr[:, keep])
        Xf, yf = Ztr, ytr
        if extra is not None:
            Zb, ybb = _aug_rows(extra, keep)
            Xf = np.vstack([Ztr, Zb]); yf = np.concatenate([ytr, ybb])
        sel = diff_select(Xf, yf, 500)
        svm = LinearSVC(C=0.02, class_weight="balanced", max_iter=5000).fit(Xf[:, sel], yf)
        dtr = svm.decision_function(Ztr[:, sel]); ssorted = np.sort(dtr)
        a = auc(ytr, _pct(ssorted, dtr)); sign = 1.0 if (a != a or a >= 0.5) else -1.0
        cv = cv_oof(B, tr, ytr, dg.loc[tr].values, extra=extra)
        P.models[(mshort, mod)] = ModalityModel(list(B.columns), sc, keep, sel, svm, ssorted, sign,
                                                round(cv[1], 3) if cv else None)
        mod_names.append(mod)
        if cv:
            oof_p[mod] = dict(zip(tr, cv[0])); cvauc[mod] = cv[1]
        if te:
            dte = svm.decision_function(sc.transform(B.loc[te].values[:, keep])[:, sel])
            q = _pct(ssorted, dte); te_p[mod] = dict(zip(te, q if sign > 0 else 1 - q))
    # optimised weights: NNLS of per-modality OOF percentiles (common train) onto truth ; floor + uniform fallback
    wm = {m: 0.0 for m in mod_names}
    common = sorted(set.intersection(*[set(oof_p[m]) for m in oof_p])) if oof_p else []
    yco = np.array([int(yall[s]) for s in common]) if common else np.array([])
    if len(common) >= 8 and len(set(yco)) == 2:
        cols = [m for m in mod_names if m in oof_p]
        O = np.column_stack([[oof_p[m].get(s, 0.5) for s in common] for m in cols])
        A = np.vstack([O, np.eye(O.shape[1])]); b = np.concatenate([yco.astype(float), np.zeros(O.shape[1])])
        w, _ = nnls(A, b)
        single = {m: auc(yco, np.array([oof_p[m][s] for s in common])) for m in cols}
        bs = max(single, key=lambda k: (single[k] if single[k] == single[k] else -1))
        comb = auc(yco, O @ w / (w.sum() or 1)) if w.sum() > 0 else np.nan
        if not (comb == comb and comb >= (single[bs] if single[bs] == single[bs] else -1) - 1e-9):
            w = np.zeros(len(cols)); w[cols.index(bs)] = 1.0
        tot = w.sum() or 1.0
        for m, wi in zip(cols, w):
            wm[m] = round(float(wi / tot), 3)
    if sum(wm.values()) == 0:                                  # fallback: uniform over informative modalities
        info = [m for m in mod_names if cvauc.get(m, 0) > 0.52] or mod_names
        for m in info:
            wm[m] = round(1.0 / len(info), 3)
    P.weights[mshort] = wm
    # OOF-calibrated present/absent threshold — maximize F1 on the deployed OOF blend (was Youden's J).
    # Youden (TPR-FPR) ignores base rate: for rare drivers a low threshold keeps FPR small yet admits many
    # false positives (over-calling, low precision). F1 penalizes false positives via precision, so it
    # raises the threshold and stops over-calling. Out-of-fold percentiles generalize (no in-sample
    # optimism). Falls back to 0.5 when the OOF is too small.
    P.thresholds[mshort] = 0.5
    P.train_auc[mshort] = None
    if len(common) >= 8 and len(set(yco)) == 2:
        wc = [m for m in mod_names if m in oof_p and wm.get(m, 0) > 0]
        if wc:
            Oc = np.column_stack([[oof_p[m].get(s, 0.5) for s in common] for m in wc])
            wv = np.array([wm[m] for m in wc]); comb = Oc @ wv / (wv.sum() or 1.0)
            _ta = auc(yco, comb); P.train_auc[mshort] = round(float(_ta), 3) if _ta == _ta else None  # fused CV-OOF AUROC on train
            bt, bf1 = 0.5, -1.0
            for t in np.unique(comb):
                pr = comb >= t
                tp = int((pr & (yco == 1)).sum()); fn = int(((~pr) & (yco == 1)).sum()); fp = int((pr & (yco == 0)).sum())
                prec = tp / (tp + fp) if (tp + fp) else 0.0
                rec = tp / (tp + fn) if (tp + fn) else 0.0
                f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
                if f1 > bf1:
                    bf1, bt = f1, float(t)
            P.thresholds[mshort] = round(bt, 3)
    P.prevalence[mshort] = round(float(np.mean([int(yall[s]) for s in allidx if pd.notna(ym.get(s)) and s not in hold])), 3)
    # held-out AUC of the deployed combiner
    te_all = sorted(set().union(*[set(te_p[m]) for m in te_p])) if te_p else []
    if te_all:
        yte = np.array([int(yall[s]) for s in te_all]); agg = np.full(len(te_all), np.nan)
        num = np.zeros(len(te_all)); den = np.zeros(len(te_all))
        for i, s in enumerate(te_all):
            for m in mod_names:
                if s in te_p.get(m, {}) and wm.get(m, 0) > 0:
                    num[i] += wm[m] * te_p[m][s]; den[i] += wm[m]
            if den[i] > 0:
                agg[i] = num[i] / den[i]; heldout_pred.setdefault(s, {})[mshort] = float(agg[i])
        ha = auc(yte, agg)
        P.heldout_auc[mshort] = round(float(ha), 3) if ha == ha else None
    log("  %-13s w=%s auc=%s" % (mshort, {k: v for k, v in wm.items() if v > 0}, P.heldout_auc.get(mshort)))

valid = [v for v in P.heldout_auc.values() if v is not None]
P.meta = {"rna": "raw+prog (UDON markers + program fractions)", "base_model": "linSVM",
          "n_train_samples": len([s for s in allidx if s not in hold]),
          "trained": time.strftime("%Y-%m-%d %H:%M")}
P.save(os.path.join(HERE, "mutation_predictor%s.pkl" % SFX))
log("\nsaved mutation_predictor%s.pkl | modalities=[%s] | %d mutations, mean held-out AUC %.3f over %d"
    % (SFX, ",".join(MODS), len(P.mutations), float(np.mean(valid)) if valid else float("nan"), len(valid)))

# ---- per-held-out-sample board reports (ALL sealed holdout, calibrated calls + abstention) ----
# skipped for tagged experiment runs (AMLMM_TAG) so they don't clobber the deployed board reports.
gate = CG.load_gate() if (CG and not TAG) else None
NPOS = {m: int(sum(1 for s in allidx if s not in hold
                   and D._labels_for_field_raw(ctx, NAME[m]).map(_m01).get(s) == 1)) for m in P.mutations}
nrep = 0; no_data = []
for s in (sorted(hold) if not TAG else []):
    sample = {mod: BLK[mod].loc[s] for mod in MODS if s in BLK[mod].index}
    preds = []
    for mshort in P.mutations:
        if not any((mshort, mod) in P.models for mod in MODS):
            continue
        pr = P.predict_one(mshort, sample); pr["mutation"] = mshort
        tl = D._labels_for_field_raw(ctx, NAME[mshort]).map(_m01).get(s)
        pr["true_label"] = ("present" if tl == 1 else "absent") if pd.notna(tl) else None
        # abstention: a call is low-confidence if underpowered (few positives) OR the held-out reliability
        # is weak (AUC < 0.65) / unknown. Thresholds alone can't fix weak drivers (KIT/PTPN11/WT1/IDH2) — the
        # separable signal isn't there — so flag them as leads-to-sequence rather than hard calls.
        ha = P.heldout_auc.get(mshort)
        if NPOS.get(mshort, 0) < 3:
            pr["confidence"] = "abstain: underpowered (n+<3)"
        elif ha is None:
            pr["confidence"] = "abstain: no held-out positives (reliability unknown)"
        elif ha < 0.65:
            pr["confidence"] = "abstain: weak separability (held-out AUC %.2f < 0.65)" % ha
        else:
            pr["confidence"] = "ok"
        preds.append(pr)
    preds.sort(key=lambda p: -(p["probability"] or 0))
    spec = None
    if gate is not None and "Composition" in sample:
        comp = sample["Composition"]
        spec = CG.score_gate(gate, pd.Series(comp.values, index=comp.index))
    # closes the loop: calls -> therapy hypotheses + confirmatory tests. On these sealed held-out reports a
    # KNOWN label (true_label == present) is what "sequenced" means, so those drivers read as confirmed.
    _pan = therapy.build_panels(preds, [p["mutation"] for p in preds if p.get("true_label") == "present"])
    rep = {"mode": "mutation_panel", "sample_key": s, "dataset": s.split("::")[0] if "::" in s else None,
           "specimen_class": (spec["call"] if spec else None), "control_gate": spec,
           "mutation_predictions": preds, "predictor": P.summary(), "validation": True,
           "treatment_panel": _pan["treatments"], "tests_panel": _pan["tests"], "panels_note": _pan["note"],
           "modalities_available": sorted(sample.keys()),
           "note": "sealed held-out sample — predicted vs known labels (OOF-calibrated thresholds)"}
    if not sample:
        rep["note"] = "held-out sample has no scRNA/atlas modality data — cannot be scored"; no_data.append(s)
    d = os.path.join(RUNS_ROOT, "predict_" + "".join(ch if ch.isalnum() else "_" for ch in s))
    os.makedirs(d, exist_ok=True)
    json.dump(rep, open(os.path.join(d, "patient_report.json"), "w"), default=str, indent=1)
    nrep += 1
log("wrote %d held-out reports (all %d holdout); no-data: %s" % (nrep, len(hold), no_data or "none"))
log("thresholds: median %.2f range %.2f-%.2f" % (
    float(np.median(list(P.thresholds.values()))) if P.thresholds else 0.5,
    min(P.thresholds.values()) if P.thresholds else 0.5, max(P.thresholds.values()) if P.thresholds else 0.5))

# ---- model card: per-mutation TRAIN (donor-grouped CV-OOF) vs sealed HELD-OUT AUROC + operating point ----
card = {"trained": P.meta.get("trained"), "n_train_samples": P.meta.get("n_train_samples"),
        "mean_heldout_auc": round(float(np.mean(valid)), 3) if valid else None,
        "note": "train_auc = donor-grouped cross-validated out-of-fold AUROC of the deployed fused combiner on "
                "the training cohort; heldout_auc = same combiner on the sealed 29-sample held-out set.",
        "mutations": {}}
for m in P.mutations:
    tw = sorted([(mm, w) for mm, w in P.weights.get(m, {}).items() if w > 0], key=lambda x: -x[1])[:3]
    ta, ha = P.train_auc.get(m), P.heldout_auc.get(m)
    card["mutations"][m] = {"train_auc": ta, "heldout_auc": ha, "n_pos": NPOS.get(m),
                            "gap": (round(ta - ha, 3) if (ta is not None and ha is not None) else None),
                            "threshold": P.thresholds.get(m), "top_modalities": [mm for mm, _ in tw]}
card["meta"] = {"modalities": MODS, "tag": TAG or "deployed"}
json.dump(card, open(os.path.join(HERE, "model_card%s.json" % SFX), "w"), indent=1)
log("wrote model_card%s.json (%d mutations; train CV-OOF vs sealed held-out AUROC; modalities=%s)"
    % (SFX, len(card["mutations"]), MODS))
log("TRAIN PREDICTOR OK")
