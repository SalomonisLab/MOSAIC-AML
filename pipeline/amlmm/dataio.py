"""Data layer — load the local OR cluster AML multimodal atlas into canonical tables.

Auto-detects layout from base_dir so the same code runs in both places:
  * 'local'   — reorganized working copy: modalities under data/, labels flat under labels/
  * 'cluster' — original deposit: modalities at base/<modality>/, labels scattered
                (Metadata/, RNA/clusters/, LSC-prediction/algorithm/)

Produces (see Context):
  ctx.tables['pseudobulks'|'samples'|'composition'|'udon_programs'|'lsc_calls'|'coverage']
  ctx._modality_paths   name -> .h5ad path (opened lazily, backed mode)
Modality matrices are NOT loaded into memory here.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

# modality .h5ad relative to the modality root (same names in both layouts)
MODALITY_FILES = {
    "RNA":               "RNA/pseudobulk_counts_hashed.h5ad",
    "GRN":               "GRN/imputed_grn_all_pseudobulks.h5ad",
    "Metabolite":        "Metabolite/pseudobulk_imputed_metabolite_aml.h5ad",
    "Lipid":             "Lipid/pseudobulk_imputed_lipid_aml.h5ad",
    "ADT":               "ADT/pseudobulk_adt_imputed.h5ad",
    "cell-communication":"cell-communication/combined_sample_by_interaction.h5ad",
}
IMPUTED_MODALITIES = ("GRN", "Metabolite", "Lipid", "ADT")
CELLSTATE_COL = "Hs-BM-titrated-reference-centroid"

LAYOUTS = {
    "local": {
        "modality_root": "data",
        "metadata":      "labels/AML_clinical_metadata_CLEAN.tsv",
        "composition":   "labels/sample_cellstate_counts.tsv",
        "lsc":           "labels/LSC_prediction_subtype_RF.tsv",
        "udon_programs": "labels/RNA_UDON_final_program_assignments.tsv",
    },
    "cluster": {
        "modality_root": "",
        "metadata":      "Metadata/AML_clinical_metadata_CLEAN.tsv",
        "composition":   "LSC-prediction/algorithm/sample_cellstate_counts.tsv",
        "lsc":           "LSC-prediction/prediction_subtype_classifications_RF.tsv",
        "udon_programs": "RNA/clusters/UDON_final_program_assignments.tsv",
    },
}


def _key(dataset, sample) -> str:
    return f"{dataset}::{sample}"


def _mode_str(s) -> str:
    """Most-common non-null value as a plain str. Robust to (a) pandas-2.x pyarrow-backed
    string columns whose value_counts().index is an Arrow array, and (b) all-null groups
    where value_counts() is empty (index[0] would raise) -> returns 'nan'."""
    vc = s.astype("object").dropna().astype(str).value_counts()
    return vc.index[0] if len(vc) else "nan"


def _read_xlsx(path) -> pd.DataFrame:
    """Minimal stdlib .xlsx reader (no openpyxl): first worksheet -> all-string DataFrame.
    Handles shared strings, inline strings, numeric cells, and ragged rows."""
    import zipfile
    import re
    import xml.etree.ElementTree as ET
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        shared = []
        if "xl/sharedStrings.xml" in names:
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(f"{ns}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{ns}t")))
        sheet = "xl/worksheets/sheet1.xml"
        if sheet not in names:
            sheet = next(n for n in names
                         if n.startswith("xl/worksheets/sheet") and n.endswith(".xml"))
        data = ET.fromstring(z.read(sheet)).find(f"{ns}sheetData")

    def col_idx(ref):
        m = re.match(r"([A-Z]+)", ref or "A")
        n = 0
        for ch in (m.group(1) if m else "A"):
            n = n * 26 + (ord(ch) - 64)
        return n - 1

    rows = []
    for r in data.findall(f"{ns}row"):
        cells, maxi = {}, -1
        for c in r.findall(f"{ns}c"):
            i = col_idx(c.get("r", "A1")); maxi = max(maxi, i)
            t, v = c.get("t"), c.find(f"{ns}v")
            if t == "s":
                val = shared[int(v.text)] if (v is not None and v.text) else ""
            elif t == "inlineStr":
                isn = c.find(f"{ns}is")
                val = "".join(n.text or "" for n in isn.iter(f"{ns}t")) if isn is not None else ""
            else:
                val = v.text if v is not None else ""
            cells[i] = val
        rows.append([cells.get(i, "") for i in range(maxi + 1)])
    if not rows:
        return pd.DataFrame()
    header = [str(h) for h in rows[0]]
    w = len(header)
    body = [(r + [""] * (w - len(r)))[:w] for r in rows[1:]]
    return pd.DataFrame(body, columns=header)


def _read_table(path) -> pd.DataFrame:
    if str(path).lower().endswith(".xlsx"):
        return _read_xlsx(path)
    return pd.read_csv(path, sep="\t", dtype=str)


def _find_meta(dirpath, default_name) -> str:
    """Find the clinical-metadata file in dirpath regardless of exact name/format
    (upstream has swapped .tsv <-> .xlsx and renamed it). Prefer tsv/csv, then 'clinical'."""
    try:
        cands = [f for f in os.listdir(dirpath)
                 if "metadata" in f.lower() and f.lower().endswith((".tsv", ".csv", ".xlsx"))]
    except OSError:
        cands = []
    if not cands:
        return os.path.join(dirpath, default_name)
    cands.sort(key=lambda f: (0 if f.lower().endswith((".tsv", ".csv")) else 1,
                              0 if "clinical" in f.lower() else 1, f))
    return os.path.join(dirpath, cands[0])


def resolve_paths(base_dir) -> dict:
    """Detect layout and return absolute paths for every input we read."""
    for name, lyt in LAYOUTS.items():
        rna = os.path.join(base_dir, lyt["modality_root"], MODALITY_FILES["RNA"])
        if os.path.exists(rna):
            mods = {m: os.path.join(base_dir, lyt["modality_root"], rel)
                    for m, rel in MODALITY_FILES.items()
                    if os.path.exists(os.path.join(base_dir, lyt["modality_root"], rel))}
            md = os.path.join(base_dir, lyt["metadata"])
            return {
                "layout": name,
                "modalities": mods,
                "metadata": _find_meta(os.path.dirname(md), os.path.basename(md)),
                "composition": os.path.join(base_dir, lyt["composition"]),
                "lsc": os.path.join(base_dir, lyt["lsc"]),
                "udon_programs": os.path.join(base_dir, lyt["udon_programs"]),
            }
    raise FileNotFoundError(
        f"Could not find the RNA matrix under {base_dir} for any known layout "
        f"(tried {[l['modality_root'] for l in LAYOUTS.values()]}).")


def load_into(ctx) -> None:
    import anndata as ad
    p = resolve_paths(ctx.config.base_dir)
    ctx._modality_paths = dict(p["modalities"])
    ctx.layout = p["layout"]

    # ---- pseudobulks from RNA obs ----
    obs = ad.read_h5ad(p["modalities"]["RNA"], backed="r").obs.copy()
    pb = pd.DataFrame({
        "sample": obs["Sample"].astype(str).values,
        "dataset": obs["Dataset"].astype(str).values,
        "cell_state": obs[CELLSTATE_COL].astype(str).values,
        "annotation": obs["Annotation"].astype(str).values,
        "n_cells": pd.to_numeric(obs["n_cells"], errors="coerce").fillna(0).astype(int).values,
    }, index=obs.index.astype(str))
    pb["sample_key"] = [_key(d, s) for d, s in zip(pb["dataset"], pb["sample"])]
    ctx.tables["pseudobulks"] = pb

    # ---- samples (aggregate pseudobulks) ----
    g = pb.groupby("sample_key")
    samples = pd.DataFrame({
        "dataset": g["dataset"].first(),
        "sample": g["sample"].first(),
        "annotation": g["annotation"].agg(_mode_str),
        "n_states": g.size(),
        "n_cells_total": g["n_cells"].sum(),
    })

    # ---- clinical metadata join ----
    meta = _read_table(p["metadata"]).fillna("")
    meta["sample_key"] = [_key(d, s) for d, s in zip(meta["Dataset"], meta["Sample"])]
    meta = meta.replace("", np.nan).drop_duplicates("sample_key").set_index("sample_key")
    keep = ["Donor_ID", "disease_category", "is_pediatric", "timepoint", "therapy",
            "drug_treatment", "age", "sex", "race",
            "FAB", "blast_pct", "WHO_classification", "ELN_risk", "clinical_response",
            "karyotype"]
    samples = samples.join(meta[[c for c in keep if c in meta.columns]], how="left")
    samples["donor_id"] = samples.get("Donor_ID")
    samples["donor_group"] = samples["donor_id"].where(
        samples["donor_id"].notna() & (samples["donor_id"].astype(str) != "nan"),
        samples.index.to_series(),
    ).astype(str)
    for c in ("age", "blast_pct"):
        if c in samples.columns:
            samples[c] = pd.to_numeric(samples[c], errors="coerce")
    ctx.tables["samples"] = samples

    # ---- LSC calls + SampleID crosswalk ----
    lsc = pd.read_csv(p["lsc"], sep="\t", dtype=str)
    lsc["sample_key"] = [_key(d, s) for d, s in zip(lsc["Dataset"], lsc["Sample"])]
    xwalk = dict(zip(lsc["SampleID"].astype(str), lsc["sample_key"]))
    ctx.tables["lsc_calls"] = lsc.drop_duplicates("sample_key").set_index("sample_key")

    # ---- composition (cell frequency) ----
    comp = pd.read_csv(p["composition"], sep="\t")
    sid = comp.iloc[:, 0].astype(str)
    comp = comp.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    comp.index = [xwalk.get(i, i) for i in sid]
    comp = comp[~comp.index.duplicated(keep="first")]
    ctx.tables["composition"] = comp

    # ---- UDON RNA programs ----
    ctx.tables["udon_programs"] = pd.read_csv(p["udon_programs"], sep="\t", dtype=str)

    # ---- coverage ----
    ctx.tables["coverage"] = {
        "layout": p["layout"],
        "n_samples": int(len(samples)),
        "n_pseudobulks": int(len(pb)),
        "n_datasets": int(samples["dataset"].nunique()),
        "samples_with_metadata": int(samples.get("disease_category").notna().sum())
            if "disease_category" in samples else 0,
        "samples_with_composition": int(samples.index.isin(comp.index).sum()),
        "samples_with_donor_id": int(samples["donor_id"].notna().sum()),
        "modalities": list(ctx._modality_paths),
    }


def feature_fidelity(ctx, modality) -> "pd.Series | None":
    """Per-feature held-out Spearman (var['heldout_spearman']) for imputed modalities."""
    a = ctx.open_modality(modality)
    if "heldout_spearman" in a.var.columns:
        return pd.to_numeric(a.var["heldout_spearman"], errors="coerce")
    return None


def sample_modality_matrix(ctx, modality, sample_keys=None, min_spearman=None) -> pd.DataFrame:
    """Aggregate a modality's pseudobulks to one vector per sample (n_cells-weighted
    mean over the sample's cell-states). Optionally keep only features with held-out
    Spearman >= min_spearman. Loads only the needed rows from the backed matrix."""
    if modality not in ctx._modality_paths:     # modality not deposited -> empty (witness skips)
        return pd.DataFrame(index=[])
    pb = ctx.tables["pseudobulks"]
    # NON-backed cached load + in-memory slice. Backed SPARSE fancy-indexing (`a[idx].to_memory()`) is
    # broken on the upgraded cluster stack (anndata 0.10.8 + scipy 1.15 / numpy 2.4 raise `isintlike`
    # inhomogeneous-shape in the backed CSR reconstruction). Reuses pseudobulk_io's workaround (read once,
    # slice every call); local import avoids a circular module load (pseudobulk_io imports dataio).
    from . import pseudobulk_io as _pio
    a = _pio._mem_modality(ctx, modality)
    keys = list(sample_keys) if sample_keys is not None else list(ctx.tables["samples"].index)
    keyset = set(keys)
    rows = pb.index[pb["sample_key"].isin(keyset)]
    pos = {name: i for i, name in enumerate(a.obs_names.astype(str))}
    rows = [r for r in rows if r in pos]
    idx = np.asarray([pos[r] for r in rows], dtype=np.intp)
    if len(idx) == 0:
        return pd.DataFrame(index=[])
    X = a.X[idx]                                            # in-memory row slice (sparse or dense)
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    df = pd.DataFrame(X, index=rows, columns=list(a.var_names.astype(str)))
    df["_w"] = pb.loc[rows, "n_cells"].clip(lower=1).values.astype(float)
    df["_k"] = pb.loc[rows, "sample_key"].values

    def wmean(block):
        # pandas 2.x excludes the grouping column (`_k`) from the per-group block; older pandas keeps it.
        # drop with errors="ignore" so the weighted mean is computed over the feature columns either way.
        ww = block["_w"].values
        vals = block.drop(columns=["_w", "_k"], errors="ignore")
        return pd.Series(np.average(vals.values, axis=0, weights=ww), index=vals.columns)

    out = df.groupby("_k", sort=False).apply(wmean)
    out = out.reindex([k for k in keys if k in set(df["_k"])])
    if min_spearman is not None and modality in IMPUTED_MODALITIES:
        fid = feature_fidelity(ctx, modality)
        if fid is not None:
            keep = set(fid[fid >= min_spearman].index)
            out = out[[c for c in out.columns if c in keep]]
    return out


def cohort_modality_matrix(ctx, modality, min_spearman=None) -> pd.DataFrame:
    """All-sample aggregated matrix (the z-score / robust-z baseline a descriptive
    witness needs), cached on ctx so a batch loop pays the aggregation once.
    Cache key includes min_spearman so a filtered + unfiltered view never collide."""
    key = f"_cohort::{modality}::{min_spearman}"
    cached = ctx.tables.get(key)
    if cached is not None:
        return cached
    M = sample_modality_matrix(ctx, modality, min_spearman=min_spearman)
    ctx.tables[key] = M
    return M


def cellcomm_matrix(ctx, sample_keys=None) -> pd.DataFrame:
    """Direct reader for the ALREADY sample-level cell-communication modality
    (obs = one row per sample, var = ligand-receptor interactions). NOT the
    pseudobulk-aggregation path: this h5ad has no pseudobulk rows.

    Loaded NON-backed (full) and cached on ctx: the matrix is small (~383 x 141k,
    ~4% dense CSR, ~25 MB) and BACKED fancy-row indexing of this file SEGFAULTS on
    the cluster's anndata 0.10.8 / HDF5 stack. Full load + in-memory slicing is
    stable in both environments. Keys rows by sample_key = Dataset::<obs idx (Sample)>."""
    a = ctx.tables.get("_cellcomm_ad")
    if a is None:
        path = ctx._modality_paths.get("cell-communication")
        if not path:
            return pd.DataFrame(index=[])
        import anndata as ad
        a = ad.read_h5ad(path)
        ctx.tables["_cellcomm_ad"] = a
    obs = a.obs
    allkeys = [_key(d, s) for d, s in zip(obs["Dataset"].astype(str), obs.index.astype(str))]
    pos = {k: i for i, k in enumerate(allkeys)}
    keys = list(sample_keys) if sample_keys is not None else allkeys
    idx = [pos[k] for k in keys if k in pos]
    rows = [k for k in keys if k in pos]
    if not idx:
        return pd.DataFrame(index=[])
    X = a.X[idx]
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    return pd.DataFrame(X, index=rows, columns=list(a.var_names.astype(str)))
