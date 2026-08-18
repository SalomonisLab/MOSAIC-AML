#!/usr/bin/env python3
"""Items 6 and 7: the parts of the supplied literature the first pass used only crudely, or not at all.

The first ELN 2022 benchmark applied the Pollyea TP53/FLT3-ITD/NRAS/KRAS rule as a FLAT three-tier
score. Loewenberg/Doehner (Blood 2024, less-intensive therapies) report something richer: signalling
mutations modify prognosis WITHIN a genotype, roughly multiplicatively --

    NPM1mut       39.0 months without a signalling mutation vs  9.9 with
    IDH2mut       36.9                                      vs 12.2
    RUNX1mut      32.5                                      vs  9.3
    MR genes      22.9                                      vs 12.9

and DDX41 is reported as a particularly favourable group under HMA therapy. Neither the interaction
structure nor DDX41 has been tested here.

  E1  interaction    does a genotype x signalling interaction beat the flat 4-gene tier?
  E2  per-genotype   the published table, reproduced (or not) genotype by genotype
  E3  DDX41          is it favourable in this cohort, and is it even present?
  E4  age bands      the ASH 2025 older-adult guidelines target >=75s specifically; does the model hold
                     there, and is the signalling effect age-dependent?

  python exp_signaling_interactions.py  ->  deliverables/exp_signaling_interactions.json
"""
import os, sys, json, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold

from amlmm.survival import data as SD, coxph as CX
from train_survival_model import BUNDLE
from eln2022 import load_mutations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_signaling_interactions.json")
PUBLISHED = {"NPM1": (39.0, 9.9), "IDH2": (36.9, 12.2), "RUNX1": (32.5, 9.3), "MR_genes": (22.9, 12.9)}
MR = ["ASXL1", "BCOR", "EZH2", "RUNX1", "SF3B1", "SRSF2", "STAG2", "U2AF1", "ZRSR2"]


def km_median_months(t, e):
    kt, ks = CX.km(np.asarray(t, float), np.asarray(e, int))
    for x, y in zip(kt, ks):
        if y <= 0.5:
            return round(12 * float(x), 1)
    return None


def main():
    t0 = time.time()
    d = np.load(BUNDLE, allow_pickle=True)
    ba = [str(s) for s in d["ba_samples"]]
    cl = SD.load_cohort(specimens=ba)
    t = cl["time_years"].values.astype(float)
    e = cl["event"].values.astype(int)
    spec = cl["specimen"].astype(str).values
    dna = (cl["dbgap_dnaseq_sample"].astype(str).values
           if "dbgap_dnaseq_sample" in cl.columns else spec)
    MUT = load_mutations(0.0)

    def gene(g, vaf=0.10):
        return np.array([any(v >= vaf for v in ((MUT.get(str(a)) or MUT.get(str(b)) or {}).get(g, [])))
                         for a, b in zip(spec, dna)])

    itd = cl["FLT3-ITD"].astype(str).str.lower().str.startswith("pos").values
    sig = itd | gene("NRAS") | gene("KRAS")          # the signalling / venetoclax-resistance axis
    mr = np.zeros(len(spec), bool)
    for g in MR:
        mr |= gene(g)
    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "n": int(len(t)), "deaths": int(e.sum()),
           "signalling_definition": "FLT3-ITD or NRAS or KRAS mutated",
           "n_signalling": int(sig.sum()),
           "source": "Loewenberg/Doehner, Blood 2024 (ELN in less-intensive therapies)"}

    # ---- E2: the published table, genotype by genotype -------------------------------------------
    groups = {"NPM1": cl["NPM1"].astype(str).str.lower().str.startswith("pos").values,
              "IDH2": gene("IDH2"), "RUNX1": gene("RUNX1"), "MR_genes": mr}
    e2 = {}
    print("== E2: does the published signalling interaction reproduce? ==")
    for name, m in groups.items():
        a, b = m & ~sig, m & sig
        if a.sum() < 8 or b.sum() < 8:
            e2[name] = {"n_without": int(a.sum()), "n_with": int(b.sum()), "note": "too few"}
            print("  %-10s too few (%d without / %d with)" % (name, a.sum(), b.sum()))
            continue
        ma, mb = km_median_months(t[a], e[a]), km_median_months(t[b], e[b])
        chi, p = CX.logrank(t[m], e[m], sig[m].astype(int))
        e2[name] = {"n_without": int(a.sum()), "median_without_months": ma,
                    "n_with": int(b.sum()), "median_with_months": mb,
                    "logrank_chi2": round(float(chi), 3), "logrank_p": float(p),
                    "published_without": PUBLISHED[name][0], "published_with": PUBLISHED[name][1],
                    "same_direction": bool(ma is not None and mb is not None and ma > mb)}
        print("  %-10s n=%3d median %6s mo  |  n=%3d median %6s mo  |  p=%.3g  (published %.1f vs %.1f)"
              % (name, a.sum(), ma, b.sum(), mb, p, PUBLISHED[name][0], PUBLISHED[name][1]))
    res["E2_per_genotype"] = e2

    # ---- E1: interaction model vs the flat tier ---------------------------------------------------
    tp53 = gene("TP53")
    flat = np.where(tp53, 2.0, np.where(sig, 1.0, 0.0))
    npm1 = groups["NPM1"]
    Zi = np.vstack([tp53.astype(float), sig.astype(float), npm1.astype(float), mr.astype(float),
                    (npm1 & sig).astype(float), (mr & sig).astype(float)]).T
    Zf = flat.reshape(-1, 1)
    g = cl["subject"].values
    oof_f, oof_i = np.full(len(t), np.nan), np.full(len(t), np.nan)
    for tr, te in GroupKFold(n_splits=5).split(np.zeros(len(t)), groups=g):
        try:
            oof_f[te] = CX.CoxPH(alpha=0.05).fit(Zf[tr], t[tr], e[tr]).risk(Zf[te])
            oof_i[te] = CX.CoxPH(alpha=0.05).fit(Zi[tr], t[tr], e[tr]).risk(Zi[te])
        except Exception:
            pass
    mask = np.isfinite(oof_f) & np.isfinite(oof_i)
    res["E1_interaction_vs_flat"] = {
        "n_scored": int(mask.sum()),
        "flat_4gene_c_index": round(float(CX.c_index(t[mask], e[mask], oof_f[mask])), 4),
        "genotype_x_signalling_c_index": round(float(CX.c_index(t[mask], e[mask], oof_i[mask])), 4)}
    print("\n== E1: flat 4-gene tier %.4f  vs  genotype x signalling %.4f =="
          % (res["E1_interaction_vs_flat"]["flat_4gene_c_index"],
             res["E1_interaction_vs_flat"]["genotype_x_signalling_c_index"]))

    # ---- E3: DDX41 --------------------------------------------------------------------------------
    dd = gene("DDX41")
    res["E3_DDX41"] = {"n_mutated": int(dd.sum())}
    print("\n== E3: DDX41 mutated in %d specimens ==" % dd.sum())
    if dd.sum() >= 5:
        chi, p = CX.logrank(t, e, dd.astype(int))
        res["E3_DDX41"].update({"median_mutated_months": km_median_months(t[dd], e[dd]),
                                "median_wildtype_months": km_median_months(t[~dd], e[~dd]),
                                "logrank_p": float(p)})
        print("   %s vs %s months, p=%.3g" % (res["E3_DDX41"]["median_mutated_months"],
                                              res["E3_DDX41"]["median_wildtype_months"], p))
    else:
        res["E3_DDX41"]["note"] = ("too few DDX41-mutated specimens to assess; the published favourable "
                                   "effect is reported under HMA therapy")
        print("   too few to assess")

    # ---- E4: the older-adult population the ASH 2025 guidelines target ----------------------------
    age = pd.to_numeric(cl["ageAtDiagnosis"], errors="coerce").values
    e4 = {}
    print("\n== E4: by age band (ASH 2025 targets >=75) ==")
    for nm, m in (("<60", age < 60), ("60-74", (age >= 60) & (age < 75)), (">=75", age >= 75)):
        if m.sum() < 25 or e[m].sum() < 8:
            e4[nm] = {"n": int(m.sum()), "note": "too small"}
            print("  %-7s n=%d  too small" % (nm, m.sum()))
            continue
        mm = m & mask
        e4[nm] = {"n": int(m.sum()), "deaths": int(e[m].sum()),
                  "signalling_prevalence": round(float(sig[m].mean()), 3),
                  "flat_4gene_c_index": round(float(CX.c_index(t[mm], e[mm], oof_f[mm])), 4),
                  "interaction_c_index": round(float(CX.c_index(t[mm], e[mm], oof_i[mm])), 4)}
        v = e4[nm]
        print("  %-7s n=%3d deaths=%3d  signalling %.0f%%  flat %.3f  interaction %.3f"
              % (nm, v["n"], v["deaths"], 100 * v["signalling_prevalence"],
                 v["flat_4gene_c_index"], v["interaction_c_index"]))

    json.dump(res, open(OUT, "w"), indent=1)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
