"""A5 — do the RNA-imputed modalities add anything to drug-response prediction?

The imputers (rna2adt / rna2grn / rna2lipid) only unpickle on Linux, so they were run on the cluster
against BeatAML bulk and the outputs shipped back. Coverage of their expected input genes was 94-96%.
No metabolite bundle exists in the cluster checkout, so that modality is absent, not failed.

Everything is imputed FROM the RNA that is already in the model, so the honest prior is that they add
little: they cannot carry information RNA does not already contain. What they can do is re-express it
in a form the linear model finds easier -- which is exactly what happened for several drivers in the
mutation layer, where metabolite and lipid were the best single blocks.

Scored on the interaction target (patient main effect removed), the honest measure.
"""
import os, sys, json, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from eval_drug_model import load_all
from train_drug_model import build_space, fit_predict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMP = os.path.join(ROOT, "data", "external", "beataml", "beataml_imputed_blocks.npz")
OUT = os.path.join(ROOT, "deliverables", "exp_compass_imputed.json")

A = load_all(); long, row_of, drugs = A["long"], A["row_of"], A["drugs"]
d = np.load(IMP, allow_pickle=True)
imp_samples = [str(s) for s in d["samples"]]
pos = {s: i for i, s in enumerate(imp_samples)}
order = [pos.get(s) for s in A["ba"]]
def blk(name, k=60):
    X = d[name][[i if i is not None else 0 for i in order]].astype(float)
    X[[i is None for i in order]] = 0.0
    if X.shape[1] > k:                       # GRN is 7,486 wide; compress before it swamps the fusion
        from sklearn.decomposition import PCA
        X = PCA(n_components=k, random_state=0).fit_transform((X - X.mean(0)) / (X.std(0) + 1e-9))
    return X
IMPB = {"ADT": blk("ADT"), "GRN": blk("GRN"), "Lipid": blk("Lipid")}
print("imputed blocks aligned:", {k: v.shape for k, v in IMPB.items()})

specs = sorted(long["specimen"].unique())
subj = long.drop_duplicates("specimen").set_index("specimen")["subject"].to_dict()
ARMS = {"deployed": [], "+ADT": ["ADT"], "+GRN": ["GRN"], "+Lipid": ["Lipid"],
        "+all imputed": ["ADT", "GRN", "Lipid"]}
pred = {a: {} for a in ARMS}
t0 = time.time()
for k, (tri, tei) in enumerate(GroupKFold(n_splits=5).split(specs, groups=[subj[s] for s in specs])):
    trs = {specs[i] for i in tri}; tes = {specs[i] for i in tei}
    f_tr = long[long["specimen"].isin(trs)]; f_te = long[long["specimen"].isin(tes)]
    rows_tr = sorted({row_of[s] for s in trs})
    fs, Z = build_space(A["X"], rows_tr, A["meta"], A["mut"], 100, 4000, A["genes"])
    for arm, extra in ARMS.items():
        Xp, _, sl = fs.transform(Z, meta=A["meta"], mut=A["mut"], blocks=["rna", "state", "mut", "clin"])
        for e in extra:                                    # append the imputed block + its slice
            sl = dict(sl); sl[e] = slice(Xp.shape[1], Xp.shape[1] + IMPB[e].shape[1])
            Xp = np.hstack([Xp, IMPB[e]])
        from amlmm.drug import model as M
        norm = M.fit_norm(f_tr, np.ones(len(f_tr), bool), tail=0.20)
        ltr, lte = M.apply_norm(f_tr, norm), M.apply_norm(f_te, norm)
        mod = M.DrugResponseModel(fs).prepare_targets(Z, A["sym2ens"]); mod.sym2ens = A["sym2ens"]
        mod.fit(ltr, Xp, row_of, drugs, sl)
        p = mod.predict_rows(lte, Xp, row_of)
        for s, dr, v in zip(lte["specimen"], lte["inhibitor"], p.loc[lte.index]):
            if v == v:
                pred[arm][(s, dr)] = v
    print("   fold %d/5 (%.0fs)" % (k + 1, time.time() - t0), flush=True)

# interaction target
Mx = long.pivot_table(index="specimen", columns="inhibitor", values="auc", aggfunc="mean")
Zt = -(Mx - Mx.median()) / (1.4826 * (Mx - Mx.median()).abs().median())
Zt = Zt.sub(Zt.median(axis=1), axis=0)
res = {}
print("\n== A5 · imputed modalities, scored on the interaction target ==")
print("  %-16s %10s %10s" % ("arm", "Spearman", "AUROC"))
for arm in ARMS:
    sp, au = [], []
    for dr in drugs:
        xs = [(Zt.at[s, dr], v) for (s, d2), v in pred[arm].items()
              if d2 == dr and s in Zt.index and np.isfinite(Zt.at[s, dr])]
        if len(xs) < 50:
            continue
        y = np.array([a for a, _ in xs]); p = np.array([b for _, b in xs])
        sp.append(float(spearmanr(y, p).statistic))
        lo, hi = np.quantile(y, [.2, .8]); yy = np.where(y >= hi, 1, np.where(y <= lo, 0, -1)); kk = yy >= 0
        if kk.sum() >= 20 and len(set(yy[kk])) == 2:
            au.append(float(roc_auc_score(yy[kk], p[kk])))
    res[arm] = {"n_drugs": len(sp), "spearman": round(float(np.mean(sp)), 4),
                "auroc": round(float(np.mean(au)), 4)}
    print("  %-16s %10.3f %10.3f" % (arm, res[arm]["spearman"], res[arm]["auroc"]))
json.dump(res, open(OUT, "w"), indent=1)
print("\nwrote %s" % OUT)
