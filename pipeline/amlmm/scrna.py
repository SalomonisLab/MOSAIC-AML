"""scRNA new-patient entry path (Phase D) — assign each query cell to a reference bone-marrow
state by cosine similarity to the cellHarmony reference, yielding cell-state COMPOSITION with no
deconvolution. The reference's 89 populations are exactly the atlas's 89 cell-states, so the
composition is a drop-in for the panel. This is the route bulk couldn't take: per-cell
classification (each cell is ~one state) instead of unmixing 89 collinear states from a bulk
average — which is why scRNA works and bulk doesn't.

Reference: cellHarmony 'hs_bm_reference' -> Hs-MarrowAtlas-L3M.txt (marker x population). We do a
focused cosine assignment (cellHarmony's core step) rather than driving its full pipeline, since
the panel only needs composition.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import sparse

import os as _os

# The cellHarmony marrow reference. A single hardcoded cluster path meant an scRNA upload failed on
# every machine that was not the cluster -- including the laptop someone would naturally try the tool
# on first. Resolve in order: explicit env override, the copy that ships inside the vendored
# altanalyze3 checkout, then the cluster archive.
_REF_NAME = "Hs-MarrowAtlas-L3M.txt"
_REF_CANDIDATES = [
    _os.environ.get("AMLMM_CELLHARMONY_REF"),
    _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))),
                  "engine-code", "altanalyze3", "altanalyze3", "components", "cellHarmony", "flask",
                  "references", "Human", "BoneMarrow", "Zhang-2024", _REF_NAME),
    ("/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/"
     "components/cellHarmony/flask/references/Human/BoneMarrow/Zhang-2024/" + _REF_NAME),
]
DEFAULT_REFERENCE = next((p for p in _REF_CANDIDATES if p and _os.path.exists(p)),
                         _REF_CANDIDATES[-1])


def load_reference(path=DEFAULT_REFERENCE):
    """marker-gene x population reference matrix (rows = genes, cols = the 89 cell-states)."""
    ref = pd.read_csv(path, sep="\t", index_col=0)
    ref.index = ref.index.astype(str)
    ref.columns = [str(c) for c in ref.columns]
    return ref


def _unit(M, axis):
    n = np.linalg.norm(M, axis=axis, keepdims=True)
    n[n == 0] = 1.0
    return M / n


def assign_cells(adata, reference, min_score=0.0):
    """Cosine-assign each query cell to a reference population.
    Returns (labels: np.ndarray[str], scores: np.ndarray[float], n_shared_markers: int).
    Query is CP10k+log1p normalized over ALL genes, then restricted to shared marker genes."""
    genes = [str(g) for g in adata.var_names]
    gset = set(genes)
    shared = [g for g in reference.index if g in gset]
    if len(shared) < 50:
        raise ValueError(f"only {len(shared)} shared marker genes between query and reference — "
                         "check that the query var_names are gene symbols")
    gpos = {g: i for i, g in enumerate(genes)}
    cols = np.array([gpos[g] for g in shared])
    X = adata.X
    # CP10k+log1p over ALL genes, then restrict to shared markers. Keep it sparse-safe: a real
    # 10x matrix (10k cells x 36k genes) densified whole is GBs and OOMs the node — instead sum
    # over all genes on the sparse matrix and densify ONLY the marker columns (cells x ~few-k).
    if sparse.issparse(X):
        tot = np.asarray(X.sum(axis=1), dtype=np.float64).reshape(-1, 1)
        sub = np.asarray(X[:, cols].todense(), dtype=np.float32)
    else:
        Xd = np.asarray(X, dtype=np.float32)
        tot = Xd.sum(axis=1, keepdims=True)
        sub = Xd[:, cols]
    tot[tot == 0] = 1.0
    Q = np.log1p(sub / tot * 1e4).astype(np.float32)  # cells x shared-markers, CP10k+log1p
    R = reference.loc[shared].to_numpy(dtype=np.float32)   # shared-markers x populations
    cos = _unit(Q, 1) @ _unit(R, 0)                   # cells x populations cosine
    j = cos.argmax(axis=1)
    scores = cos.max(axis=1)
    pops = list(reference.columns)
    labels = np.array([pops[k] for k in j], dtype=object)
    if min_score > 0:
        labels = np.where(scores >= min_score, labels, "unassigned")
    return labels, scores, len(shared)


def composition(labels, reference):
    """Cell-state composition (fraction per state) over the full 89-population vocabulary."""
    pops = [str(c) for c in reference.columns]
    s = pd.Series(labels)
    s = s[s != "unassigned"]
    freq = s.value_counts(normalize=True) if len(s) else pd.Series(dtype=float)
    return freq.reindex(pops).fillna(0.0)


def composition_from_query(adata, reference_path=DEFAULT_REFERENCE, min_score=0.0):
    """One call: query AnnData -> (composition Series over the 89 atlas states, mean cosine,
    n_cells, n_shared_markers). The composition row is a drop-in for the panel's `composition`."""
    ref = load_reference(reference_path)
    labels, scores, n_shared = assign_cells(adata, ref, min_score=min_score)
    comp = composition(labels, ref)
    return {"composition": comp, "mean_cosine": float(np.mean(scores)),
            "n_cells": int(len(labels)), "n_shared_markers": int(n_shared),
            "n_states_present": int((comp > 0).sum())}
