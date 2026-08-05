#!/usr/bin/env python3
"""Cell-state localization map — WHERE in the 89-state marrow each driver's signal lives.

For every board driver x cell-state, score how much that single cell-state carries the driver's signal,
using two orthogonal, non-imputed axes:

  1. EXPRESSION localization — donor-grouped 3-fold OOF AUROC of a linSVM (same recipe as the deployed
     predictor: StandardScaler -> diff_select -> LinearSVC C=0.02 balanced) trained on the state's OWN
     measured-RNA pseudobulk, mutant vs non-mutant. High AUROC in state S = the transcriptional signal the
     classifier keys on is readable in compartment S. Leakage-aware: holdout samples excluded, donor-grouped
     folds, per-fold feature selection.
  2. COMPOSITION shift — Mann-Whitney U on state abundance (cell fraction) mutant vs non-mutant. Tells whether
     the driver reshapes the marrow toward/away from compartment S, independent of expression.

Both axes are MEASURED / COMPUTED (RNA + composition), NOT imputed from RNA, so "where the signal lives" is
not a circular RNA->RNA claim. Descriptive interpretation layer, consistent with (not overriding) the
sample-level call. Runs locally (no LSF): reads data/RNA/pseudobulk_counts_hashed.h5ad directly.

Usage:  python cellstate_localize.py [--drivers NPM1,TP53,...] [--out gui/cellstate_localization.json]
"""
import os, sys, json, time, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from amlmm.context import build_context, Config
from amlmm import discovery as D, genetics
from amlmm.predictor import diff_select, _pct

MIN_POS, MIN_NEG = 5, 8          # per-state minimum mutant / non-mutant samples to score a cell
MIN_CELLS = 10                    # a sample "has" state S only if >= this many cells there (h5ad floor already 10)
TOPK = 300                        # diff_select genes per fold
GENE_PREVALENCE = 0.15            # keep genes nonzero in >= this fraction of a state's rows
CELLSTATE_COL = "Hs-BM-titrated-reference-centroid"


def log(m): print(m, flush=True)


def load_h5ad_perstate(path):
    """Return (X csr rows=(sample,state), sample_key[], state[], n_cells[], genes[])."""
    with h5py.File(path, "r") as f:
        def col(n):
            g = f["obs"][n]
            if isinstance(g, h5py.Group):
                cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
                return cats[g["codes"][:]]
            v = g[:]
            return np.array([x.decode() if isinstance(x, bytes) else x for x in v]) if v.dtype.kind == "S" else v
        ds, sm, st = col("Dataset"), col("Sample"), col(CELLSTATE_COL)
        nc = np.asarray(col("n_cells"), float)
        genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]])
        Xg = f["X"]
        shape = tuple(Xg.attrs["shape"]) if "shape" in Xg.attrs else (len(sm), len(genes))
        X = csr_matrix((Xg["data"][:], Xg["indices"][:], Xg["indptr"][:]), shape=shape)
    skey = np.array(["%s::%s" % (d, s) for d, s in zip(ds, sm)])
    return X, skey, st, nc, genes


def cp10k_log1p(M):
    tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1.0
    return np.log1p(M / tot[:, None] * 1e4)


def oof_auroc(Xd, y, groups):
    """donor-grouped 3-fold OOF AUROC, deployed recipe (scale -> diff_select -> linSVM)."""
    ng = len(set(groups))
    if ng < 2 or (y == 1).sum() < 2 or (y == 0).sum() < 2:
        return None
    oof = np.full(len(y), np.nan)
    for tri, vai in GroupKFold(min(3, ng)).split(Xd, y, groups):
        if len(set(y[tri])) < 2:
            continue
        keep = Xd[tri].std(0) > 0
        if keep.sum() < 2:
            continue
        sc = StandardScaler().fit(Xd[tri][:, keep])
        Ztr, Zva = sc.transform(Xd[tri][:, keep]), sc.transform(Xd[vai][:, keep])
        sel = diff_select(Ztr, y[tri], min(TOPK, Ztr.shape[1]))
        d = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000).fit(Ztr[:, sel], y[tri]).decision_function(Zva[:, sel])
        oof[vai] = d
    ok = ~np.isnan(oof)
    if ok.sum() < (MIN_POS + MIN_NEG) or len(set(y[ok])) < 2:
        return None
    a = roc_auc_score(y[ok], oof[ok])
    return max(a, 1 - a)         # orientation-free: "separability", direction reported via composition


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drivers", default="")
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-states", type=int, default=0)
    args = ap.parse_args()

    t0 = time.time()
    ctx = build_context(Config(run_id="cellstate_localize"))
    samples = ctx.tables["samples"]
    # donor grouping (same as train_predictor: donor_group; fall back to Donor_ID/dataset)
    if "donor_group" in samples.columns:
        dg = samples["donor_group"].astype(str)
    elif "Donor_ID" in samples.columns:
        dg = samples["Donor_ID"].astype(str)
    else:
        dg = samples["dataset"].astype(str)
    hold = set(ctx.holdout)
    comp = ctx.tables.get("composition")             # per-sample x state fractions
    M = ctx.tables.get("mutations")
    if M is None:
        M = genetics.build_mutation_matrix(ctx)

    all_flags = sorted(c for c in M.columns if str(c).startswith(("mut_", "cyto_")))
    want = set(x.strip() for x in args.drivers.split(",") if x.strip())
    def short(f): return f.replace("mut_", "").replace("cyto_", "")
    flags = [f for f in all_flags if (not want or short(f) in want or f in want)]

    log("loading RNA per-state pseudobulks ...")
    X, skey, state, ncell, genes = load_h5ad_perstate(ctx._modality_paths["RNA"])
    states = sorted(set(state))
    if args.max_states:
        states = states[:args.max_states]
    log("  %d rows, %d samples, %d states, %d genes  (%.1fs)" % (X.shape[0], len(set(skey)), len(states), len(genes), time.time() - t0))
    # index rows by state -> row positions
    rows_by_state = {s: np.where((state == s) & (ncell >= MIN_CELLS))[0] for s in states}

    out = {"meta": {"metric": "per-state donor-grouped OOF AUROC (mutant vs non-mutant), measured RNA + composition shift",
                    "min_pos": MIN_POS, "min_neg": MIN_NEG, "topk": TOPK, "n_states": len(states),
                    "built": time.strftime("%Y-%m-%d %H:%M")}, "drivers": {}}

    for f in flags:
        d0 = time.time()
        sh = short(f)
        raw = D._labels_for_field_raw(ctx, f)
        lab = {s: (1 if raw.get(s) == "present" else 0 if raw.get(s) == "absent" else None) for s in samples.index}
        mut = set(s for s in samples.index if lab[s] == 1 and s not in hold)
        non = set(s for s in samples.index if lab[s] == 0 and s not in hold)
        if len(mut) < MIN_POS:
            log("  %-14s skip (only %d mutant)" % (sh, len(mut)))
            continue
        per_state = {}
        for s in states:
            rows = rows_by_state[s]
            if len(rows) == 0:
                continue
            sk = skey[rows]
            ymask = np.array([1 if k in mut else 0 if k in non else -1 for k in sk])
            keep = ymask >= 0
            if keep.sum() == 0:
                continue
            rr, yy, gg = rows[keep], ymask[keep], dg.reindex(sk[keep]).values
            npos, nneg = int((yy == 1).sum()), int((yy == 0).sum())
            if npos < MIN_POS or nneg < MIN_NEG:
                continue
            sub = cp10k_log1p(X[rr].toarray())
            prev = (sub > 0).mean(0)
            gkeep = prev >= GENE_PREVALENCE
            if gkeep.sum() < 20:
                continue
            a = oof_auroc(sub[:, gkeep], yy, gg)
            if a is None:
                continue
            # composition shift for this state (mutant vs non-mutant abundance)
            cshift, cp = None, None
            if comp is not None and s in comp.columns:
                mv = comp.loc[[k for k in mut if k in comp.index], s].astype(float).values
                nv = comp.loc[[k for k in non if k in comp.index], s].astype(float).values
                if len(mv) >= 3 and len(nv) >= 3:
                    try:
                        u, cp = mannwhitneyu(mv, nv, alternative="two-sided")
                        cshift = float(np.median(mv) - np.median(nv))
                    except Exception:
                        pass
            per_state[s] = {"auroc": round(float(a), 3), "n_mut": npos, "n_non": nneg,
                            "comp_shift": (round(cshift, 4) if cshift is not None else None),
                            "comp_p": (round(float(cp), 4) if cp is not None else None)}
        if not per_state:
            log("  %-14s no scorable states" % sh)
            continue
        top = sorted(per_state.items(), key=lambda kv: -kv[1]["auroc"])[:5]
        out["drivers"][sh] = {"n_mut": len(mut), "n_non": len(non), "states": per_state,
                              "top_states": [t[0] for t in top]}
        log("  %-14s %2d states scored | top: %s  (%.1fs)"
            % (sh, len(per_state), ", ".join("%s %.2f" % (k, v["auroc"]) for k, v in top[:3]), time.time() - d0))

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        json.dump(out, open(args.out, "w"), indent=1)
        log("wrote %s (%d drivers)" % (args.out, len(out["drivers"])))
    log("DONE %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
