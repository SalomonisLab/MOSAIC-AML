"""Item 1 — rebuild the drug-response LABEL from the raw wells, where the real ceiling lives.

Two changes, both from `raw_inhibitor.txt` (555,583 wells) rather than the released probit fits:

  replicate-averaged AUC   Averaging n replicate measurements raises reliability by Spearman-Brown,
                           2r/(1+r): at the measured median r = 0.529 two replicates give 0.692, lifting
                           the achievable Spearman ceiling from 0.727 to 0.832. This is the only change
                           available that raises the CEILING rather than chasing it.
  dose-anchored AUC        The released AUC integrates viability across the whole tested range, up to
                           10 uM -- concentrations no patient reaches. Re-integrating over <= 1 uM asks
                           the clinically meaningful question instead of rewarding cytotoxicity at
                           unreachable doses.

AUC is the trapezoidal integral of viability over log10(concentration), rescaled to the released 0-300
convention so everything downstream is comparable.

  python build_better_targets.py -> data/external/beataml/better_targets.tsv
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "external", "beataml", "raw_inhibitor.txt")
OUT = os.path.join(ROOT, "data", "external", "beataml", "better_targets.tsv")
CLINICAL_MAX_UM = 1.0          # anchor: plausible free plasma exposure for most of this panel


def auc_from_series(c, v):
    """Trapezoid over log10 concentration, scaled to the released 0-300 AUC convention."""
    o = np.argsort(c)
    c, v = np.asarray(c)[o], np.asarray(v)[o]
    if len(c) < 3 or c[0] <= 0:
        return np.nan
    x = np.log10(c)
    span = x[-1] - x[0]
    if span <= 0:
        return np.nan
    return float(np.trapezoid(np.clip(v, 0, 200), x) / span * (300.0 / 100.0))


t0 = time.time()
r = pd.read_csv(RAW, sep="\t", low_memory=False)
r = r[r["passed_qc"].astype(str).str.upper().eq("TRUE")]
r["v"] = pd.to_numeric(r["normalized_viability"], errors="coerce")
r["c"] = pd.to_numeric(r["well_concentration"], errors="coerce")
r = r.dropna(subset=["v", "c", "dbgap_rnaseq_sample"])
r["specimen"] = r["dbgap_rnaseq_sample"].astype(str)
print("usable wells: %d | specimens %d | inhibitors %d"
      % (len(r), r["specimen"].nunique(), r["inhibitor"].nunique()), flush=True)

# average replicate wells at each concentration BEFORE integrating -- that is where the noise reduction
# happens; averaging two AUCs afterwards is equivalent only if the dose grids match, which they do not
g = r.groupby(["specimen", "inhibitor", "c"])
avg = g["v"].mean().rename("v_mean").reset_index()
nrep = g["v"].size().rename("n_wells").reset_index()
avg = avg.merge(nrep, on=["specimen", "inhibitor", "c"])

rows = []
for (sp, dr), s in avg.groupby(["specimen", "inhibitor"]):
    full = auc_from_series(s["c"].values, s["v_mean"].values)
    lowm = s["c"] <= CLINICAL_MAX_UM
    anch = auc_from_series(s.loc[lowm, "c"].values, s.loc[lowm, "v_mean"].values) if lowm.sum() >= 3 else np.nan
    rows.append((sp, dr, full, anch, float(s["n_wells"].mean()), int(len(s))))
t = pd.DataFrame(rows, columns=["specimen", "inhibitor", "auc_repavg", "auc_dose_anchored",
                                "mean_wells_per_conc", "n_concentrations"])
t.to_csv(OUT, sep="\t", index=False)
print("wrote %s: %d (specimen x inhibitor) rows | %.1f%% have >1 well per concentration | (%.0fs)"
      % (OUT, len(t), 100 * (t["mean_wells_per_conc"] > 1).mean(), time.time() - t0))
print("  dose-anchored computable for %.0f%% of rows" % (100 * t["auc_dose_anchored"].notna().mean()))
