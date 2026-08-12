#!/usr/bin/env python3
"""Re-baseline risk prediction against ELN 2022, and test where ELN 2022 is known to fail.

Everything in this platform has been benchmarked against BeatAML's shipped ELN2017 column. That is the
wrong bar: 2022 is current practice, it differs from 2017 in exactly the patients we get wrong, and it
has to be INFERRED here (see eln2022.py). Four questions, three of them taken straight from the
literature Nathan supplied:

  Q1  the bar          how much survival discrimination do ELN2017, ELN2022@10% VAF and ELN2022@40%
                       VAF actually carry, and does the VAF threshold move it?

  Q2  adverse is not   Rollig/Bill (Leukemia 2025): patients whose only adverse feature is a
      one thing        myelodysplasia-related gene mutation live markedly longer than other adverse
                       patients (14.7 vs 8.3 months). If that reproduces here, ELN 2022's adverse
                       category is heterogeneous and splitting it is free discrimination.

  Q3  it fails under   Pollyea/Dohner (Blood 2024): "ELN prognostic classifiers did not provide
      less-intensive   clinically meaningful risk stratification" under venetoclax-azacitidine;
      therapy          TP53/FLT3-ITD/NRAS/KRAS split those patients into 26.5 / 12.1 / 5.5 month
                       groups instead. We independently measured C-index 0.554 for our model and
                       0.481 -- below chance -- for age+ELN in non-intensively-treated patients. This
                       tests whether the published 4-gene rule rescues that stratum.

  Q4  does MOSAIC add  the deployed model against the ELN 2022 bar rather than the 2017 one.
      anything

  python exp_eln2022_benchmark.py  ->  deliverables/exp_eln2022_benchmark.json
"""
import os, sys, json, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold

from amlmm.survival import data as SD, coxph as CX
from train_survival_model import build_blocks, fit_arm, ARM_BLOCKS, BUNDLE

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_eln2022_benchmark.json")
ORDER = {"Favorable": 0.0, "Intermediate": 1.0, "Adverse": 2.0}
MR_GENES = ["ASXL1", "BCOR", "EZH2", "RUNX1", "SF3B1", "SRSF2", "STAG2", "U2AF1", "ZRSR2"]


def km_median(t, e):
    kt, ks = CX.km(np.asarray(t, float), np.asarray(e, int))
    for x, y in zip(kt, ks):
        if y <= 0.5:
            return float(x)
    return None


def group_table(t, e, g, label):
    rows = {}
    for lv in pd.unique(g):
        m = g == lv
        if m.sum() < 5:
            continue
        med = km_median(t[m], e[m])
        rows[str(lv)] = {"n": int(m.sum()), "deaths": int(e[m].sum()),
                         "median_months": None if med is None else round(12 * med, 1)}
    return rows


def ordinal_c(t, e, g):
    """C-index of an ordinal risk label (Favorable<Intermediate<Adverse) against survival."""
    s = np.array([ORDER.get(str(x), np.nan) for x in g], float)
    m = np.isfinite(s)
    if m.sum() < 20 or e[m].sum() < 5:
        return None
    return round(float(CX.c_index(t[m], e[m], s[m])), 4)


def main():
    t0 = time.time()
    d = np.load(BUNDLE, allow_pickle=True)
    genes = [str(x) for x in d["genes"]]
    ba = [str(s) for s in d["ba_samples"]]
    X = d["ba_X"].astype(np.float64)
    mut_all = pd.DataFrame(d["ba_L"].astype(float), index=ba, columns=[str(c) for c in d["drivers"]])
    row_of = {s: i for i, s in enumerate(ba)}
    cl = SD.load_cohort(specimens=ba)
    rows = np.array([row_of[s] for s in cl["specimen"]])
    t = cl["time_years"].values.astype(float)
    e = cl["event"].values.astype(int)
    g = cl["subject"].values
    spec = cl["specimen"].astype(str).values

    res = {"generated": time.strftime("%Y-%m-%d %H:%M"),
           "cohort": {"n": int(len(t)), "deaths": int(e.sum())},
           "sources": {
               "ELN2022": "Dohner et al., Blood 140(12):1345-1377 (2022), Table 6",
               "Q2": "Rollig/Bill, Leukemia 2025 - MDS-related mutations in European AML",
               "Q3": "Pollyea/Dohner, Blood 144(21):2211 (2024) - ELN 2022 does not stratify ven/aza"}}

    # ---- load the two inferred ELN 2022 labels -------------------------------------------
    lab = {}
    for v in (10, 40):
        p = os.path.join(ROOT, "labels", "eln2022_beataml_vaf%02d.tsv" % v)
        if not os.path.exists(p):
            continue
        df = pd.read_csv(p, sep="\t").drop_duplicates("specimen").set_index("specimen")
        lab[v] = df
    eln17 = cl["ELN2017"].astype(str).values

    # ---------------- Q1: the bar --------------------------------------------------------
    q1 = {"ELN2017 (shipped)": {"c_index": ordinal_c(t, e, eln17),
                                "groups": group_table(t, e, eln17, "ELN2017")}}
    covered = {}
    for v, df in lab.items():
        s = pd.Series(spec).map(df["ELN2022"]).values
        covered[v] = s
        q1["ELN2022 @ VAF>=%d%%" % v] = {"c_index": ordinal_c(t, e, s),
                                         "n_labelled": int(pd.notna(s).sum()),
                                         "groups": group_table(t, e, pd.Series(s).fillna("NA").values, "x")}
    res["Q1_the_bar"] = q1
    print("== Q1: how much does the risk label carry ==")
    for k, v in q1.items():
        print("  %-24s C-index %s   %s" % (k, v["c_index"],
              {kk: vv["median_months"] for kk, vv in v["groups"].items() if kk in ORDER}))

    # ---------------- Q2: is `adverse` one thing? ----------------------------------------
    # Gene status comes from mutations.txt, not the bundle's driver matrix: the bundle is missing
    # RUNX1, SF3B1 and U2AF1, which are three of the nine ELN 2022 myelodysplasia-related genes, and
    # sourcing from the variant file is also what lets the VAF threshold apply consistently.
    from eln2022 import load_mutations
    MUT = load_mutations(0.0)
    dna = (cl["dbgap_dnaseq_sample"].astype(str).values if "dbgap_dnaseq_sample" in cl.columns
           else spec)

    def gene_calls(gene, vaf=0.10):
        out = np.zeros(len(spec), bool)
        for i, (s1, s2) in enumerate(zip(spec, dna)):
            d = MUT.get(str(s1)) or MUT.get(str(s2)) or {}
            out[i] = any(v >= vaf for v in d.get(gene, []))
        return out

    if 10 in covered:
        s10 = covered[10]
        mr = np.zeros(len(spec), bool)
        for gname in MR_GENES:
            mr |= gene_calls(gname)
        adv = pd.Series(s10).eq("Adverse").values
        sub = np.where(adv & mr, "adverse_MR_gene", np.where(adv & ~mr, "adverse_other", "not_adverse"))
        tab = group_table(t, e, sub, "adverse split")
        a1 = adv & mr; a2 = adv & ~mr
        chi = p = None
        if a1.sum() >= 5 and a2.sum() >= 5:
            lh = a1 | a2
            chi, p = CX.logrank(t[lh], e[lh], (a1[lh]).astype(int))
        res["Q2_adverse_heterogeneity"] = {
            "groups": tab, "logrank_chi2": None if chi is None else round(float(chi), 3),
            "logrank_p": None if p is None else float(p),
            "published": "Rollig/Bill: MR-gene adverse 14.7 vs non-MR adverse 8.3 months (p<0.001)"}
        print("\n== Q2: splitting ELN 2022 adverse by MDS-related gene status ==")
        for k, v in tab.items():
            print("  %-18s n=%3d deaths=%3d  median %s months" % (k, v["n"], v["deaths"], v["median_months"]))
        if p is not None:
            print("  log-rank MR vs non-MR within adverse: chi2 %.2f, p = %.3g" % (chi, p))

    # ---------------- Q3: the less-intensive stratum --------------------------------------
    ty = cl["typeInductionTx"].astype(str)
    intensive = ty.str.contains("Standard Chemo", na=False).values
    q3 = {}
    for nm, mask in (("intensive", intensive), ("non_intensive_or_unknown", ~intensive)):
        if mask.sum() < 25:
            continue
        row = {"n": int(mask.sum()), "deaths": int(e[mask].sum()),
               "ELN2017_c_index": ordinal_c(t[mask], e[mask], eln17[mask])}
        if 10 in covered:
            row["ELN2022_c_index"] = ordinal_c(t[mask], e[mask], covered[10][mask])
        # the published 4-gene rule: TP53 / FLT3-ITD / NRAS / KRAS
        has = gene_calls
        tp53 = has("TP53")
        itd = cl["FLT3-ITD"].astype(str).str.lower().str.startswith("pos").values if "FLT3-ITD" in cl.columns else np.zeros(len(t), bool)
        ras = has("NRAS") | has("KRAS")
        tier = np.where(tp53, 2.0, np.where(itd | ras, 1.0, 0.0))     # lower = higher benefit
        mm = mask
        if e[mm].sum() >= 5:
            row["pollyea_4gene_c_index"] = round(float(CX.c_index(t[mm], e[mm], tier[mm])), 4)
            lv = np.where(tier == 0, "lower_risk(no TP53/ITD/RAS)", np.where(tier == 1, "intermediate(ITD/RAS)", "higher_risk(TP53)"))
            row["pollyea_groups"] = group_table(t[mm], e[mm], lv[mm], "pollyea")
        q3[nm] = row
    res["Q3_treatment_strata"] = q3
    res["Q3_note"] = ("BeatAML records induction type, not venetoclax-azacitidine specifically, so the "
                      "non-intensive stratum here is a proxy for the population Pollyea studied, not a "
                      "replication of it.")
    print("\n== Q3: does the risk label work outside intensive induction ==")
    for nm, r in q3.items():
        print("  %-26s n=%3d deaths=%3d | ELN2017 %s | ELN2022 %s | Pollyea 4-gene %s"
              % (nm, r["n"], r["deaths"], r.get("ELN2017_c_index"), r.get("ELN2022_c_index"),
                 r.get("pollyea_4gene_c_index")))

    # ---------------- Q4: does MOSAIC add over the ELN 2022 bar? --------------------------
    if 10 in covered:
        s10 = pd.Series(covered[10])
        ok = s10.notna().values
        idx = np.where(ok)[0]
        oof = np.full(len(t), np.nan)
        gk = GroupKFold(n_splits=5)
        for i_in, i_out in gk.split(np.zeros(len(idx)), groups=g[idx]):
            tr, te = idx[i_in], idx[i_out]
            try:
                _, B = build_blocks(cl, X, mut_all, rows, sorted(rows[tr]), genes, 60, 4000, train_idx=tr)
                r, _ = fit_arm("deployed", {b: B[b][tr] for b in ARM_BLOCKS["deployed"]},
                               {b: B[b][te] for b in ARM_BLOCKS["deployed"]}, t[tr], e[tr], g[tr])
                oof[te] = r
            except Exception as ex:
                print("   fold failed: %s" % str(ex)[:110])
        m = np.isfinite(oof)
        c_mosaic = round(float(CX.c_index(t[m], e[m], oof[m])), 4) if m.sum() > 20 else None
        c_eln = ordinal_c(t[m], e[m], covered[10][m])
        # combined: ELN 2022 ordinal + the model score, stacked by a 2-covariate Cox
        comb = None
        if c_mosaic is not None:
            Z = np.vstack([oof[m], [ORDER[str(x)] for x in covered[10][m]]]).T
            cm = CX.CoxPH(alpha=0.01).fit(Z, t[m], e[m])
            comb = round(float(CX.c_index(t[m], e[m], cm.risk(Z))), 4)
        res["Q4_mosaic_vs_eln2022"] = {"n_scored": int(m.sum()), "ELN2022_c_index": c_eln,
                                       "MOSAIC_deployed_c_index": c_mosaic,
                                       "combined_c_index": comb,
                                       "note": ("combined is fitted in-sample on two covariates and is "
                                                "indicative only; the honest comparison is the two "
                                                "standalone columns")}
        print("\n== Q4: MOSAIC against the ELN 2022 bar (n=%d) ==" % m.sum())
        print("  ELN2022 %s | MOSAIC deployed %s | both %s" % (c_eln, c_mosaic, comb))

    json.dump(res, open(OUT, "w"), indent=1)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
