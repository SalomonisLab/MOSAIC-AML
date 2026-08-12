#!/usr/bin/env python3
"""Does ANY cheap statistic flag a single specimen the mutation caller's reference does not describe?

The n=1 gap is real: `predict_cohort` re-references a batch of >=8 against itself, but an upload is one
patient, so it silently inherits whatever offset the nearest reference carries. GSE281087 is the worked
example -- 54 of 65 positive calls wrong, with nothing on the report indicating anything was unusual.

Two candidates have already FAILED and are recorded so nobody re-tries them:

  present-call count   GSE281087 specimens make 4-7 present calls against 3.7 expected for a reference
                       sample, so a "too many calls" rule fired on 0/15 out-of-distribution specimens
                       and 13/20 in-distribution ones. Removed from bulk_predictor.py.
  cellHarmony cosine   the three SORTED NORMAL populations (CD34-1 0.693, CD34-2 0.727, GMP-1 0.787)
                       have the HIGHEST mapping cosines in the cohort, not the lowest -- clean
                       homogeneous populations map crisply and messy AML maps poorly, so the signal
                       runs backwards.

This tests the expression-space statistics, which is where the offset actually lives. After z-scoring
against reference R, an in-distribution specimen's z should look like N(0,1) over R's own genes. So:

  mean|z|         E|N(0,1)| = 0.798 for a specimen the reference describes
  frac |z| > 3    ~0.003 for a matched specimen
  mean z          a pure location shift, which is what a sorted population should produce
  sd z            a scale shift

Each is scored the only way that matters: does it separate GSE281087 from the specimens the reference
WAS built from, with a threshold that does not fire on the latter? A statistic that cannot do that is
reported as a failure and not deployed.

  python exp_ood_guard.py  ->  deliverables/exp_ood_guard.json
"""
import os, sys, json, glob, time, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_ood_guard.json")
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")
GSE_H5 = "/data/salomonis-archive/FASTQs/PublicDatasets/DATASET_8_GSE281087/counts_h5"


def stats_for(z, sel):
    """The four candidate statistics on one specimen's z-vector, over the model's selected genes."""
    v = np.asarray(z, float)[sel]
    v = v[np.isfinite(v)]
    if not len(v):
        return None
    return {"mean_abs_z": float(np.mean(np.abs(v))), "frac_z_gt3": float(np.mean(np.abs(v) > 3)),
            "mean_z": float(np.mean(v)), "sd_z": float(np.std(v))}


def bulk_equiv_from_h5(path, genes, sym2ens):
    """Pseudobulk a 10x h5 into the model's gene space (sum of counts per gene, CPM, log2)."""
    import h5py
    from scipy.sparse import csc_matrix
    with h5py.File(path, "r") as f:
        grp = f["matrix"] if "matrix" in f else f[list(f.keys())[0]]
        names = [x.decode() if isinstance(x, bytes) else str(x)
                 for x in (grp["features/name"][:] if "features" in grp else grp["gene_names"][:])]
        data, idx, ptr = grp["data"][:], grp["indices"][:], grp["indptr"][:]
        shape = tuple(grp["shape"][:])
        M = csc_matrix((data, idx, ptr), shape=shape)
    tot = np.asarray(M.sum(axis=1)).ravel()                      # per-gene total across cells
    gsum = {}
    for nm, v in zip(names, tot):
        gsum[nm] = gsum.get(nm, 0.0) + float(v)
    s = sum(gsum.values()) or 1.0
    gpos = {g: i for i, g in enumerate(genes)}
    out = np.zeros(len(genes))
    hit = 0
    for nm, v in gsum.items():
        j = gpos.get(nm if nm in gpos else sym2ens.get(nm))
        if j is not None:
            out[j] = 1e6 * v / s; hit += 1
    return out, hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse-dir", default=GSE_H5)
    a = ap.parse_args()
    t0 = time.time()
    import pickle
    from amlmm.bulk_predictor import BulkMutationPredictor as B
    m = pickle.load(open(os.path.join(HERE, "bulk_mutation_predictor.pkl"), "rb"))
    m.__class__ = B
    genes, sym2ens = m.genes, (m.sym2ens or {})
    sel = np.unique(np.concatenate([mm["sel"] for mm in m.models.values()]))
    print("model: %d genes, %d selected across %d categories" % (len(genes), len(sel), len(m.models)),
          flush=True)

    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "n_selected_genes": int(len(sel)),
           "failed_candidates": {
               "present_call_count": "fired on 0/15 OOD and 13/20 in-distribution; removed",
               "cellharmony_cosine": "runs backwards -- sorted normals have the HIGHEST cosines"},
           "groups": {}}

    # ---- in-distribution: the BeatAML specimens the `beataml` reference was built from -------------
    d = np.load(BUNDLE, allow_pickle=True)
    X = d["ba_X"].astype(np.float64)
    ba = [str(s) for s in d["ba_samples"]]
    bgenes = [str(x) for x in d["genes"]]
    # map the bundle's gene order onto the model's
    gpos = {g: i for i, g in enumerate(genes)}
    cols = np.array([gpos.get(g, -1) for g in bgenes])
    ok = cols >= 0
    A = np.zeros((X.shape[0], len(genes)))
    A[:, cols[ok]] = X[:, ok]
    for ref in ("beataml", "sc"):
        if ref not in m.refs:
            continue
        Z = m._z(np.vstack([m._clog(r) for r in A]), ref)
        rows = [stats_for(z, sel) for z in Z]
        rows = [r for r in rows if r]
        res["groups"]["beataml_specimens__ref_%s" % ref] = {
            "n": len(rows),
            **{k: {"median": round(float(np.median([r[k] for r in rows])), 4),
                   "p95": round(float(np.percentile([r[k] for r in rows], 95)), 4)}
               for k in ("mean_abs_z", "frac_z_gt3", "mean_z", "sd_z")}}
        print("  beataml specimens vs ref=%s: median mean|z| %.3f"
              % (ref, np.median([r["mean_abs_z"] for r in rows])), flush=True)

    # ---- out-of-distribution: GSE281087 -----------------------------------------------------------
    files = sorted(glob.glob(os.path.join(a.gse_dir, "*.h5")))
    print("GSE281087 h5 files: %d" % len(files), flush=True)
    gse = []
    for f in files:
        try:
            vec, hit = bulk_equiv_from_h5(f, genes, sym2ens)
            for ref in ("sc", "beataml"):
                if ref not in m.refs:
                    continue
                z = m._z(m._clog(vec), ref)
                st = stats_for(z, sel)
                if st:
                    st.update({"sample": os.path.basename(f), "ref": ref, "genes_hit": hit})
                    gse.append(st)
        except Exception as ex:
            print("   %s failed: %s" % (os.path.basename(f), str(ex)[:90]), flush=True)
    for ref in ("sc", "beataml"):
        rows = [r for r in gse if r["ref"] == ref]
        if not rows:
            continue
        res["groups"]["gse281087__ref_%s" % ref] = {
            "n": len(rows),
            **{k: {"median": round(float(np.median([r[k] for r in rows])), 4),
                   "min": round(float(np.min([r[k] for r in rows])), 4),
                   "max": round(float(np.max([r[k] for r in rows])), 4)}
               for k in ("mean_abs_z", "frac_z_gt3", "mean_z", "sd_z")},
            "per_sample": sorted(rows, key=lambda r: -r["mean_abs_z"])}

    # ---- does any statistic separate, with a threshold that spares the in-distribution set? --------
    verdict = {}
    for ref in ("sc", "beataml"):
        idk, ook = "beataml_specimens__ref_%s" % ref, "gse281087__ref_%s" % ref
        if idk not in res["groups"] or ook not in res["groups"]:
            continue
        for k in ("mean_abs_z", "frac_z_gt3", "sd_z"):
            thr = res["groups"][idk][k]["p95"]
            rows = [r[k] for r in gse if r["ref"] == ref]
            if not rows:
                continue
            caught = float(np.mean([x > thr for x in rows]))
            verdict["%s @ ref=%s" % (k, ref)] = {
                "threshold_p95_of_in_distribution": round(float(thr), 4),
                "fraction_of_gse281087_flagged": round(caught, 3),
                "false_alarm_rate_by_construction": 0.05,
                "usable": bool(caught >= 0.80)}
    res["verdict"] = verdict
    print("\n== verdict ==", flush=True)
    for k, v in verdict.items():
        print("  %-26s flags %.0f%% of GSE281087 at a 5%% false-alarm threshold -> %s"
              % (k, 100 * v["fraction_of_gse281087_flagged"], "USABLE" if v["usable"] else "not usable"),
              flush=True)
    if not any(v["usable"] for v in verdict.values()):
        res["conclusion"] = ("no cheap expression-space statistic separates this cohort at a 5% false "
                             "alarm rate; the n=1 out-of-distribution gap stays OPEN and documented")

    json.dump(res, open(OUT, "w"), indent=1)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
