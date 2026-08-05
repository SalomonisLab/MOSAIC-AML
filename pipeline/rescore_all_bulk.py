#!/usr/bin/env python3
"""Re-score EVERY board report through the PRIMARY bulk variant-level caller (50 categories).

Why: existing reports were produced by the old 26-mutation sc multimodal predictor, which never even
EVALUATED ~25 of the bulk caller's categories (SRSF2, SF3B1_*, U2AF1_*, JAK2_V617F, ZRSR2, and every
variant split: DNMT3A_R882 vs nonR882, FLT3_ITD vs TKD, TP53 hotspot vs LOF, ...). Those weren't
"predicted absent" — they were absent from the model.

MERGE, don't replace: the bulk caller has NO cytogenetics, and inv16/del5/del7/complex/trisomy8/KMT2A
are exactly what the sc multimodal system uniquely does well (inv16 1.00, del7 0.99). So each report
becomes  bulk 50 variant-level calls  +  the sc system's cytogenetic calls (with their known labels).

Each sample's bulk-equivalent is recomputed from its ORIGINAL source:
  A) atlas samples (predict_*) -> bulk_features.bulk_rna_matrix(ctx)          [CP10k+log1p -> expm1]
  B) Trumpp cohort (trumpp_*)  -> sum raw counts per obs['Library'] from the cellHarmony h5ad
  C) uploads (ingest_*)        -> the 10x recorded in report['ingest']['source']
Samples present in no source are left UNTOUCHED and listed as orphans.

Run on an LSF compute node (the Trumpp h5ad is large):
  bsub -q test -W 240 -M 40000 -R "rusage[mem=40000]" -o rescore.log \
    /usr/local/anaconda3-2020/bin/python rescore_all_bulk.py
"""
import os, sys, json, glob, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

from amlmm import therapy
from ingest_patient import bulk_mutation_result, bulk_equiv_from_adata, load_query, BULK_PKL

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
RUNS = os.environ.get("AMLMM_RUNS") or os.path.join(PROJ, "runs")
TRUMPP_H5AD = "/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_AML_UDON/output/combined_with_umap_and_markers.h5ad"
LIBCOL = "Library"
DRY = os.environ.get("DRY_RUN") == "1"

# the sc system's exclusive territory — keep these calls, the bulk caller has no cytogenetics
SC_CYTO = {"complex", "del5", "del7", "inv16", "inv(16)_CBFB-MYH11", "kmt2a",
           "KMT2A-rearrangement", "trisomy8"}
# old sc name -> new bulk category, so a known label carries over
TRUTH_ALIAS = {"FLT3-ITD": "FLT3_ITD", "FLT3-TKD": "FLT3_TKD_D835/I836", "NPM1": "NPM1_exon12_frameshift",
               "IDH1": "IDH1_R132", "RUNX1": "RUNX1_LOF", "WT1": "WT1_LOF"}


def log(m): print(m, flush=True)


# ---------------- sources ----------------
def source_atlas():
    """{sample_key: linear CP10k Series over gene symbols} for every atlas sample."""
    try:
        from amlmm.context import build_context, Config
        import bulk_features as BF
        ctx = build_context(Config(run_id="rescore_all"))
        B = BF.bulk_rna_matrix(ctx)                      # CP10k + log1p
        lin = np.expm1(B)                                # -> linear CP10k (what the 'sc' ref expects)
        log("  A) atlas: %d samples x %d genes" % lin.shape)
        return {str(i): lin.loc[i] for i in lin.index}
    except Exception as e:
        log("  A) atlas source unavailable: %s" % e)
        return {}


def source_trumpp():
    """{'Trumpp::<Library>': linear CP10k Series} by summing raw counts per Library."""
    if not os.path.exists(TRUMPP_H5AD):
        log("  B) trumpp h5ad not found (%s) — skipping" % TRUMPP_H5AD)
        return {}
    try:
        import anndata as ad
        import scipy.sparse as sp
        a = ad.read_h5ad(TRUMPP_H5AD)
        X = a.layers["counts"] if ("counts" in getattr(a, "layers", {})) else a.X
        libs = a.obs[LIBCOL].astype(str).values
        genes = [str(g) for g in a.var_names]
        out = {}
        for lib in pd.unique(libs):
            m = (libs == lib)
            sub = X[m]
            tot = np.asarray(sub.sum(axis=0)).ravel() if sp.issparse(sub) else np.asarray(sub).sum(axis=0).ravel()
            s = float(tot.sum()) or 1.0
            out["Trumpp::" + lib] = pd.Series(tot / s * 1e4, index=genes)
        log("  B) trumpp: %d libraries from %d cells" % (len(out), len(libs)))
        return out
    except Exception as e:
        log("  B) trumpp source failed: %s" % e)
        return {}


def bulk_equiv_from_10x_h5(path):
    """Scanpy-FREE bulk-equivalent from a 10x filtered_feature_bc_matrix.h5 -> linear CP10k Series.

    The cluster python has NumPy 2.4, which numba (and therefore scanpy) refuses to load, so the normal
    scanpy reader dies here. We only need per-gene totals, so read the CSC arrays with h5py and bincount
    the gene indices — no scanpy, no numba. CITE-seq files carry Antibody Capture rows too, so keep only
    'Gene Expression' features (matching load_query).
    """
    import h5py
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        data = g["data"][:]
        indices = g["indices"][:]
        shape = g["shape"][:]
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in g["features/name"][:]]
        try:
            ftype = [x.decode() if isinstance(x, bytes) else str(x) for x in g["features/feature_type"][:]]
        except Exception:
            ftype = ["Gene Expression"] * len(names)
    n_genes = int(shape[0])
    tot = np.bincount(indices.astype(np.int64), weights=data.astype(float), minlength=n_genes)
    keep = [i for i, t in enumerate(ftype) if t == "Gene Expression"]
    tot = tot[keep]; names = [names[i] for i in keep]
    s = float(tot.sum()) or 1.0
    ser = pd.Series(tot / s * 1e4, index=names)
    return ser.groupby(ser.index).sum()          # collapse duplicate gene symbols


def source_upload(rep):
    """bulk-equivalent for an uploaded sample, from the 10x path recorded in its report."""
    src = ((rep.get("ingest") or {}).get("source")) or ((rep.get("provenance") or {}).get("source"))
    if not src or not os.path.exists(src):
        return None
    p = src
    if os.path.isdir(p):
        h5 = os.path.join(p, "filtered_feature_bc_matrix.h5")
        if os.path.exists(h5):
            p = h5
    if str(p).endswith(".h5"):
        return bulk_equiv_from_10x_h5(p)         # preferred: no scanpy/numba
    if str(p).endswith(".h5ad"):
        import anndata as ad                     # anndata alone is fine on this python
        return bulk_equiv_from_adata(ad.read_h5ad(p))
    adata, _ = load_query(src)                   # last resort (needs scanpy)
    return bulk_equiv_from_adata(adata)


# ---------------- rescore ----------------
def main():
    if not os.path.exists(BULK_PKL):
        log("FATAL: bulk_mutation_predictor.pkl not found at %s" % BULK_PKL); return 1
    log("re-scoring reports under %s" % RUNS)
    A = source_atlas()
    B = source_trumpp()
    lookup = {}; lookup.update(A); lookup.update(B)

    only = os.environ.get("ONLY")                # e.g. ONLY=ingest_  -> re-run just the uploads
    done, orphan, failed = [], [], []
    for f in sorted(glob.glob(os.path.join(RUNS, "*", "patient_report.json"))):
        run = os.path.basename(os.path.dirname(f))
        if only and not run.startswith(only):
            continue
        try:
            with open(f, encoding="utf-8") as fh:
                rep = json.load(fh)
        except Exception as e:
            failed.append((run, "unreadable: %s" % e)); continue
        old = rep.get("mutation_predictions")
        if not isinstance(old, list) or not old:
            continue                                     # predates the caller entirely; leave alone
        key = str(rep.get("sample_key"))

        ser = lookup.get(key)
        if ser is None and run.startswith("ingest_"):
            try:
                ser = source_upload(rep)
            except Exception as e:
                failed.append((run, "upload load: %s" % e)); continue
        if ser is None:
            orphan.append((run, key)); continue

        # known labels from the OLD report carry over where the category is the same thing
        truth = {}
        for p in old:
            t = p.get("true_label")
            if not t:
                continue
            m = p.get("mutation")
            truth[m] = t
            if m in TRUTH_ALIAS:
                truth[TRUTH_ALIAS[m]] = t
        supplied = [m for m, t in truth.items() if t == "present"]

        try:
            preds, caller, _ = bulk_mutation_result(ser, supplied, ref="sc")
        except Exception as e:
            failed.append((run, "bulk score: %s" % e)); continue
        if preds is None:
            failed.append((run, "caller returned nothing")); continue

        for p in preds:                                  # attach truth where we actually have it
            p["true_label"] = truth.get(p["mutation"])
        # keep the sc system's cytogenetics (bulk has none) — they carry their own known labels
        cyto = [p for p in old if p.get("mutation") in SC_CYTO]
        for p in cyto:
            p["source"] = "sc multimodal (cytogenetics)"
        merged = preds + cyto
        merged.sort(key=lambda p: -(p.get("probability") or 0))

        panels = therapy.build_panels(merged, supplied)
        rep["mutation_predictions"] = merged
        rep["mutation_caller"] = dict(caller, note=caller.get("note", "") +
                                      " Cytogenetic calls (%d) retained from the sc multimodal predictor, "
                                      "which the bulk caller does not cover." % len(cyto))
        rep["treatment_panel"] = panels["treatments"]
        rep["tests_panel"] = panels["tests"]
        rep["panels_note"] = panels["note"]
        if not DRY:
            with open(f, "w", encoding="utf-8") as fh:
                json.dump(rep, fh, default=str, indent=1)
        npres = sum(1 for p in merged if p.get("call") == "present" and p.get("confidence") == "ok")
        done.append((run, len(preds), len(cyto), npres))

    log("\nre-scored %d report(s)%s:" % (len(done), " [DRY RUN — nothing written]" if DRY else ""))
    for run, nb, nc, npres in done:
        log("   %-44s bulk=%-3d +cyto=%-2d  present=%d" % (run, nb, nc, npres))
    if orphan:
        log("\n%d ORPHAN report(s) — sample not in the atlas / Trumpp / its upload source (left untouched):" % len(orphan))
        for run, key in orphan:
            log("   %-44s %s" % (run, key))
    if failed:
        log("\n%d FAILED:" % len(failed))
        for run, why in failed:
            log("   %-44s %s" % (run, why))
    log("\nRESCORE OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
