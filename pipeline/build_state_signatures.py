#!/usr/bin/env python3
"""Derive differentiation-state signatures from our own marrow atlas, in the shared ENSG feature space.

BeatAML2's own analysis found differentiation state to be a broad determinant of ex-vivo drug response,
so the drug model needs a state axis for BULK specimens -- which have no cell-state structure at all.
Full 89-state deconvolution of bulk was already shown to be too collinear to trust here, so instead we
take the robust half of that idea: collapse the 89 atlas states into 10 lineage groups, learn a marker
signature per group from the atlas itself, and score any expression vector against them.

The signatures live in the bundle's 14,237-ENSG space, which means the identical score is computable on
(a) a BeatAML bulk specimen, (b) a single-cell sample's bulk-equivalent, and (c) one cell state's
pseudobulk -- the last being what makes the state-resolved response model possible.

  python build_state_signatures.py  ->  labels/cellstate_signatures.json
"""
import os, sys, json, re
import numpy as np, h5py
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
H5 = os.path.join(ROOT, "data", "RNA", "pseudobulk_counts_hashed.h5ad")
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")
NORM_EXP = os.path.join(ROOT, "data", "external", "beataml", "norm_exp.txt")
OUT = os.path.join(ROOT, "labels", "cellstate_signatures.json")
STATE_COL = "Hs-BM-titrated-reference-centroid"
TOPK = 60

# 89 atlas states -> 10 lineage groups. Ordered: the first pattern that matches wins, so the specific
# patterns (MEP-Eryth, MultiLin-GMP) must precede the general ones (MEP, Multilin).
GROUPS = [
    ("stroma",      r"^(MSC|Osteoblast|Stromal)"),
    ("T_NK",        r"^(CD4 |CD8 |MAIT|NK-|T CD)"),
    ("B_lineage",   r"^(Pro-B|pre-B|Transitional-B|B Memory|Plasma Cell|BMCP)"),
    ("DC",          r"^(cDC|pDC|ASDC|pre-DC)"),
    ("erythroid",   r"^(ERP-|Erythroblast|MEP-Eryth)"),
    ("MEP_Mk",      r"^(MEP-|MKP-|MK-Platelet|MPP-MEP)"),
    ("monocytic",   r"^(Classical-Mono|Intermediate Mono|Mono-|Non-Classical Mono|Mac$|MDP-)"),
    ("granulocytic", r"^(cMOP|preNeu|immNeu|Myeloid intermediate)"),
    ("GMP",         r"^(MultiLin-GMP)"),
    ("LMPP_CLP",    r"^(LMPP|CLP)"),
    ("HSC_MPP",     r"^(HSC-|MPP-|Multilin-)"),
]
# a coarse primitive->mature ordering used to derive a single differentiation axis
PRIMITIVE = ["HSC_MPP", "LMPP_CLP", "GMP"]
MATURE = ["monocytic", "granulocytic", "DC"]


def group_of(state):
    for name, pat in GROUPS:
        if re.match(pat, state):
            return name
    return None


def main():
    d = np.load(BUNDLE, allow_pickle=True)
    keep_ens = [str(g) for g in d["genes"]]
    ens_set = set(keep_ens)

    # gene symbol -> ENSG, from the BeatAML expression header (two columns only, cheap)
    import pandas as pd
    hdr = pd.read_csv(NORM_EXP, sep="\t", usecols=["stable_id", "display_label"])
    sym2ens = {}
    for e, s in zip(hdr["stable_id"].astype(str), hdr["display_label"].astype(str)):
        if e in ens_set:
            sym2ens.setdefault(s, e)

    with h5py.File(H5, "r") as f:
        def col(n):
            g = f["obs"][n]
            if isinstance(g, h5py.Group):
                cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
                return cats[g["codes"][:]]
            return g[:]
        states = col(STATE_COL)
        genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]])
        Xg = f["X"]; shape = tuple(Xg.attrs["shape"])
        X = csr_matrix((Xg["data"][:], Xg["indices"][:], Xg["indptr"][:]), shape=shape)

    print("atlas pseudobulks %s | genes look like %s" % (str(shape), list(genes[:3])))
    # map atlas gene symbols onto the shared ENSG space
    gi, ens_for_col = [], []
    seen = set()
    for j, g in enumerate(genes):
        e = g if g in ens_set else sym2ens.get(g)
        if e in ens_set and e not in seen:
            seen.add(e); gi.append(j); ens_for_col.append(e)
    gi = np.array(gi)
    print("atlas genes mapped into the shared space: %d / %d" % (len(gi), len(keep_ens)))

    grp = np.array([group_of(s) or "" for s in states])
    names = [n for n, _ in GROUPS][::-1]
    unmapped = sorted(set(states[grp == ""]))
    if unmapped:
        print("NOT assigned to a lineage group (dropped):", unmapped)

    # per-group aggregate profile: sum RAW counts across all its pseudobulk rows, THEN CP10k + log1p
    prof = {}
    for name in names:
        m = grp == name
        if not m.sum():
            continue
        v = np.asarray(X[m][:, gi].sum(0)).ravel().astype(np.float64)
        tot = v.sum() or 1.0
        prof[name] = np.log1p(v / tot * 1e4)
        print("  %-13s %5d pseudobulks  %3d states" % (name, int(m.sum()), len(set(states[m]))))

    P = np.vstack([prof[n] for n in names])
    sigs = {}
    for i, name in enumerate(names):
        others = np.delete(P, i, axis=0).max(0)
        spec = P[i] - others                      # specificity: over the best competing lineage
        idx = np.argsort(spec)[::-1][:TOPK]
        idx = [j for j in idx if spec[j] > 0][:TOPK]
        sigs[name] = [ens_for_col[j] for j in idx]

    out = {"space": "ENSG (bundle 14,237)", "topk": TOPK, "state_column": STATE_COL,
           "groups": {n: [p for g, p in GROUPS if g == n][0] for n in names},
           "signatures": sigs,
           "axes": {"primitive": PRIMITIVE, "mature": MATURE},
           "n_states_used": int((grp != "").sum())}
    json.dump(out, open(OUT, "w"), indent=1)
    print("\nwrote %s" % OUT)
    ens2sym = {v: k for k, v in sym2ens.items()}
    for n in names:
        print("  %-13s %s" % (n, ", ".join(ens2sym.get(e, e) for e in sigs[n][:10])))


if __name__ == "__main__":
    main()
