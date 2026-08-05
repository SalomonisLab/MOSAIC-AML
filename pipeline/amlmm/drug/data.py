"""BeatAML2 ex-vivo inhibitor response: load, QC, normalise within drug, define response classes.

The released `beataml_probit_curve_fits_v4_dbgap.txt` holds one fitted dose-response curve per
(specimen x inhibitor). The response variable is the area under the fitted viability curve (`auc`,
here on a 0-300 scale): LOWER auc = MORE sensitive. That number is not comparable across inhibitors -
concentration ranges, potency and cohort distributions all differ - so every downstream use is of a
WITHIN-DRUG standardised value, exactly as the Beat AML heat maps did it.

Three things this module is careful about, because each is a way to get a good-looking model that is
measuring nothing:

1. **Curve quality.** A probit fit can converge onto noise. Rows are flagged (not silently dropped)
   for non-convergence, an off-panel concentration range, an implausible *increasing* viability curve,
   and a within-drug extreme deviance. `qc_pass` is the conjunction; the flags stay on the row so any
   analysis can report what it excluded.
2. **Within-drug normalisation.** Robust (median / MAD) z plus a rank percentile. Robust rather than
   mean/SD because the AUC distribution for a potent drug is strongly left-skewed and a handful of
   exquisitely sensitive specimens would otherwise set the scale.
3. **Tail classes, not a global cutoff.** sensitive = bottom `tail` of within-drug AUC, resistant =
   top `tail`, and the middle is 'indeterminate': retained for regression, excluded from the binary
   task. A single universal AUC cutoff across 165 inhibitors would mostly encode drug potency.

  load()               -> tidy long frame, one row per (specimen x inhibitor), QC'd and normalised
  drug_summary(df)     -> one row per inhibitor: n, dynamic range, wave shift, eligibility
  eligible_drugs(...)  -> the inhibitors that pass the initial-implementation inclusion criteria
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))          # .../AML-multimodal
BA_DIR = os.path.join(ROOT, "data", "external", "beataml")

CURVES = "beataml_probit_curve_fits_v4_dbgap.txt"
CLIN = "beataml_wv1to4_clinical.xlsx"

SENS, RES, IND = "sensitive", "resistant", "indeterminate"


# ---------------------------------------------------------------- loading ----
def _read_curves(ba_dir):
    p = os.path.join(ba_dir, CURVES)
    if not os.path.exists(p):
        raise FileNotFoundError(
            "%s not found. Fetch it from the public BeatAML2 data repository:\n"
            "  https://github.com/biodev/beataml2.0_data/raw/main/%s" % (p, CURVES))
    d = pd.read_csv(p, sep="\t", low_memory=False)
    d["inhibitor"] = d["inhibitor"].astype(str).str.strip()
    d["specimen"] = d["dbgap_rnaseq_sample"].astype(str)
    d["subject"] = d["dbgap_subject_id"].astype(str)
    return d


def _read_clinical(ba_dir):
    cl = pd.read_excel(os.path.join(ba_dir, CLIN))
    cl["specimen"] = cl["dbgap_rnaseq_sample"].astype(str)
    cl["subject"] = cl["dbgap_subject_id"].astype(str)
    return cl


# -------------------------------------------------------------------- QC ----
def qc_annotate(d):
    """Attach per-row curve-quality flags and the conjunction `qc_pass`.

    Flags are additive and kept on the frame so that any report can state exactly what it dropped.
    `increasing` curves (viability rising with dose) are kept - they are genuine non-responders - but
    marked, because their AUC is an upper bound rather than a fitted potency.
    """
    d = d.copy()
    d["flag_not_converged"] = ~d["converged"].astype(bool)
    d["flag_increasing"] = d["curve_type"].astype(str).eq("increasing")
    d["flag_auc_missing"] = ~np.isfinite(pd.to_numeric(d["auc"], errors="coerce"))

    # off-panel concentration range: the assay panel changed for a minority of rows; comparing an AUC
    # integrated over a different concentration window to the drug's modal window is not like-for-like
    key = d["inhibitor"].astype(str)
    modal = d.groupby(key)[["min_conc", "max_conc"]].agg(lambda s: s.round(6).mode().iloc[0]
                                                         if len(s.mode()) else np.nan)
    mn = key.map(modal["min_conc"]); mx = key.map(modal["max_conc"])
    d["flag_offpanel_conc"] = (~np.isclose(d["min_conc"], mn, rtol=1e-3, equal_nan=True) |
                               ~np.isclose(d["max_conc"], mx, rtol=1e-3, equal_nan=True))

    # within-drug extreme lack of fit: deviance above the drug's own 99th percentile
    dev = pd.to_numeric(d["deviance"], errors="coerce")
    cut = dev.groupby(key).transform(lambda s: s.quantile(0.99) if s.notna().sum() >= 20 else np.inf)
    d["flag_poor_fit"] = dev > cut

    d["qc_pass"] = ~(d["flag_not_converged"] | d["flag_auc_missing"] |
                     d["flag_offpanel_conc"] | d["flag_poor_fit"])
    return d


# ------------------------------------------------- within-drug normalising ----
def normalize_within_drug(d, min_n=20):
    """Add `auc_z` (robust z, sign kept: higher = more resistant), `auc_pct` (rank in [0,1]) and
    `sens_score` = -auc_z so that larger is always 'better response'."""
    d = d.copy()
    a = pd.to_numeric(d["auc"], errors="coerce")
    g = d.groupby("inhibitor")["auc"]
    med = g.transform("median")
    mad = g.transform(lambda s: np.nanmedian(np.abs(s - np.nanmedian(s))))
    scale = 1.4826 * mad
    scale = scale.where(scale > 1e-9, g.transform("std").replace(0, np.nan))     # degenerate MAD fallback
    d["auc_z"] = (a - med) / scale
    d["auc_pct"] = g.rank(pct=True, method="average")
    n = d.groupby("inhibitor")["auc"].transform("size")
    d.loc[n < min_n, ["auc_z", "auc_pct"]] = np.nan                             # too few to normalise against
    d["sens_score"] = -d["auc_z"]
    return d


def assign_classes(d, tail=0.20):
    """sensitive / resistant / indeterminate by within-drug AUC tails.

    An `increasing`-curve row can never be called sensitive: viability that rises with dose is not a
    response, and letting it into the positive class would teach the classifier assay noise.
    """
    d = d.copy()
    lo, hi = tail, 1.0 - tail
    cls = pd.Series(IND, index=d.index, dtype=object)
    cls[d["auc_pct"] <= lo] = SENS
    cls[d["auc_pct"] >= hi] = RES
    cls[d["auc_pct"].isna()] = None
    cls[(cls == SENS) & d["flag_increasing"]] = IND
    d["resp_class"] = cls
    d["y_sensitive"] = np.where(cls == SENS, 1.0, np.where(cls == RES, 0.0, np.nan))
    return d


# ------------------------------------------------------------- public API ----
def load(ba_dir=BA_DIR, specimens=None, tail=0.20, qc_only=True, single_agent_only=True):
    """Tidy (specimen x inhibitor) response table joined to specimen-level metadata.

    specimens : optional iterable restricting to specimens we hold expression for.
    """
    d = _read_curves(ba_dir)
    if single_agent_only:
        d = d[d["type"].astype(str).eq("single-agent")]
    if specimens is not None:
        d = d[d["specimen"].isin(set(map(str, specimens)))]
    d = qc_annotate(d)
    n_all = len(d)
    dropped = {f: int(d.loc[~d["qc_pass"], f].sum()) for f in d.columns if f.startswith("flag_")}
    if qc_only:
        d = d[d["qc_pass"]]
    d = normalize_within_drug(d)
    d = assign_classes(d, tail=tail)

    cl = _read_clinical(ba_dir)
    meta_cols = ["specimen", "cohort", "centerID", "ageAtSpecimenAcquisition", "consensus_sex",
                 "specimenType", "diseaseStageAtSpecimenCollection", "dxAtSpecimenAcquisition",
                 "isRelapse", "isDenovo", "isTransformed", "priorMDS", "ELN2017",
                 "fabBlastMorphology", "%.Blasts.in.BM", "%.Blasts.in.PB", "wbcCount"]
    meta = cl[[c for c in meta_cols if c in cl.columns]].drop_duplicates("specimen")
    d = d.merge(meta, on="specimen", how="left")
    d.attrs["n_rows_before_qc"] = n_all
    d.attrs["n_rows_after_qc"] = len(d)
    d.attrs["qc_dropped_by_flag"] = dropped       # counted BEFORE filtering, else every count is 0
    return d


def drug_summary(d):
    """One row per inhibitor: how much data, how much dynamic range, and is the response wave-stable.

    `wave_shift` is the standardised difference in median AUC between the Waves1+2 and Waves3+4
    acquisitions. A drug whose whole distribution moved between waves is measuring the assay batch as
    much as the biology, and is excluded from the initial deployable panel.
    """
    rows = []
    for drug, g in d.groupby("inhibitor"):
        a = pd.to_numeric(g["auc"], errors="coerce")
        w1 = a[g["cohort"].astype(str).str.contains("1", na=False)]
        w2 = a[g["cohort"].astype(str).str.contains("3", na=False)]
        sd = a.std(ddof=1)
        shift = (np.nan if (len(w1) < 15 or len(w2) < 15 or not sd or not np.isfinite(sd))
                 else float((w2.median() - w1.median()) / sd))
        rows.append({
            "inhibitor": drug,
            "n_obs": int(len(g)),
            "n_specimens": int(g["specimen"].nunique()),
            "n_subjects": int(g["subject"].nunique()),
            "auc_median": float(a.median()),
            "auc_mad": float(np.nanmedian(np.abs(a - np.nanmedian(a)))),
            "auc_iqr": float(a.quantile(0.75) - a.quantile(0.25)),
            "frac_increasing": float(g["flag_increasing"].mean()),
            "n_wave12": int(len(w1)), "n_wave34": int(len(w2)),
            "wave_shift": shift,
            "n_sensitive": int((g["resp_class"] == SENS).sum()),
            "n_resistant": int((g["resp_class"] == RES).sum()),
        })
    return pd.DataFrame(rows).sort_values("n_specimens", ascending=False).reset_index(drop=True)


def eligible_drugs(summary, min_specimens=120, min_subjects=110, min_iqr=20.0,
                   max_frac_increasing=0.60, max_wave_shift=0.6, min_class=25):
    """Initial-implementation inclusion criteria, applied to `drug_summary` output.

    Deliberately conservative: a model that reports on 165 inhibitors of which half are unlearnable is
    worse than one that reports on the subset where the experiment actually produced signal.
    """
    s = summary
    keep = ((s["n_specimens"] >= min_specimens) & (s["n_subjects"] >= min_subjects) &
            (s["auc_iqr"] >= min_iqr) & (s["frac_increasing"] <= max_frac_increasing) &
            (s["n_sensitive"] >= min_class) & (s["n_resistant"] >= min_class) &
            (s["wave_shift"].abs().fillna(0.0) <= max_wave_shift))
    out = s.copy()
    out["eligible"] = keep
    out["exclusion"] = ""
    out.loc[s["n_specimens"] < min_specimens, "exclusion"] += "n_specimens;"
    out.loc[s["n_subjects"] < min_subjects, "exclusion"] += "n_subjects;"
    out.loc[s["auc_iqr"] < min_iqr, "exclusion"] += "no_dynamic_range;"
    out.loc[s["frac_increasing"] > max_frac_increasing, "exclusion"] += "mostly_nonresponsive_curves;"
    out.loc[(s["n_sensitive"] < min_class) | (s["n_resistant"] < min_class), "exclusion"] += "small_tails;"
    out.loc[s["wave_shift"].abs().fillna(0.0) > max_wave_shift, "exclusion"] += "wave_shift;"

    # A drug knocked out ONLY by wave_shift still has ample data and dynamic range - its AUC
    # distribution simply moved between acquisition waves. Rather than deleting clinically central
    # agents (Cytarabine is the induction backbone) we give them their own tier: modelled, reported,
    # but only ever validated leave-wave-out, and always shown with the caveat attached.
    out["tier"] = np.where(keep, "primary",
                           np.where(out["exclusion"] == "wave_shift;", "wave_conditional", "excluded"))
    return out


if __name__ == "__main__":
    import sys
    d = load()
    print("rows %d (of %d before QC)  specimens %d  subjects %d  inhibitors %d"
          % (len(d), d.attrs["n_rows_before_qc"], d["specimen"].nunique(),
             d["subject"].nunique(), d["inhibitor"].nunique()))
    for f in [c for c in d.columns if c.startswith("flag_")]:
        print("  %-22s %d" % (f, int(d[f].sum())))
    s = drug_summary(d)
    e = eligible_drugs(s)
    print("\neligible inhibitors: %d / %d" % (int(e["eligible"].sum()), len(e)))
    print(e[e["eligible"]].head(12)[["inhibitor", "n_specimens", "n_subjects", "auc_iqr",
                                     "n_sensitive", "n_resistant", "wave_shift"]].to_string(index=False))
    print("\nexcluded reasons:", e.loc[~e["eligible"], "exclusion"].value_counts().to_dict())
