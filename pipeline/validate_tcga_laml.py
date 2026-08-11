#!/usr/bin/env python3
"""B2 — external validation of the MOSAIC-AML survival layer on TCGA-LAML.

Everything reported so far comes from BeatAML2: one cohort, one sequencing pipeline, one set of
clinical conventions. Patient-grouped CV and a sealed hold-out protect against learning a patient, but
neither protects against learning a COHORT. This script asks the only question that can settle it —
does the frozen model rank survival in patients it has never seen, measured on a different platform,
treated at different institutions, a decade earlier?

  TCGA-LAML   n=173 with expression (Illumina HiSeq RSEM, gene symbols), 186 with follow-up
  BeatAML2    n=444 training patients (245 deaths)

WHAT IS AND IS NOT REFIT. The Cox coefficients, the PCA rotation, the variable-gene selection and the
NNLS fusion weights are all loaded frozen from survival_model.pkl and never touched. The ONE thing
refit on TCGA is the per-gene z-reference (median/MAD), because RSEM log2 counts and BeatAML's units
are not the same scale and a model applied across platforms without cohort-matched normalisation is
measuring the platform. That is the same cohort-matched-reference mechanism the drug layer already
uses for single-cell input, not a refit of the predictor.

HONEST LIMITATIONS, stated up front rather than buried:
  * TCGA mutation calls are not in this download, so the `mut` block enters as unreported — exactly
    the degraded path a real upload without a mutation caller takes, with the trailing
    fraction-missing column set to 1.
  * TCGA-LAML is essentially treatment-homogeneous (standard induction), so the +0.035 C-index that
    baseline induction type contributes cannot be externally validated here. The `deployed` arm is
    reported, but its treatment block carries no information in this cohort.
  * CALGB cytogenetic risk is mapped onto ELN 2017 (Favorable/Intermediate/Poor -> Fav/Int/Adverse).
    Close, not identical; the age+ELN bar is therefore approximate.

  python validate_tcga_laml.py  ->  deliverables/validation_tcga_laml.json
                                    deliverables/fig_Sv4_tcga_km.png
"""
import os, sys, json, gzip, time, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

from amlmm.survival import data as SD, coxph as CX

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TCGA = os.path.join(ROOT, "data", "external", "tcga_laml")
BUNDLE = os.path.join(HERE, "survival_model.pkl")
OUT = os.path.join(ROOT, "deliverables", "validation_tcga_laml.json")
FIG = os.path.join(ROOT, "deliverables", "figures", "Sv4_tcga_external_validation.png")

HORIZONS = [1.0, 2.0, 5.0]
RISK_MAP = {"Favorable": "Favorable", "Intermediate/Normal": "Intermediate", "Poor": "Adverse"}

# `data/` is gitignored (the BeatAML expression matrix alone is 281 MB), so fetch on first run rather
# than leaving a reproduction step that only works if someone already has the files.
XENA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
FILES = {"HiSeqV2.gz": "TCGA.LAML.sampleMap%2FHiSeqV2.gz",
         "LAML_clinicalMatrix": "TCGA.LAML.sampleMap%2FLAML_clinicalMatrix",
         "LAML_survival.txt": "survival%2FLAML_survival.txt"}


def ensure_data():
    import urllib.request
    os.makedirs(TCGA, exist_ok=True)
    for name, path in FILES.items():
        dst = os.path.join(TCGA, name)
        if os.path.exists(dst) and os.path.getsize(dst) > 1000:
            continue
        print("fetching %s from the UCSC Xena hub ..." % name)
        urllib.request.urlretrieve(XENA + path, dst)


def load_tcga():
    """Expression, clinical and survival for one primary specimen per patient."""
    with gzip.open(os.path.join(TCGA, "HiSeqV2.gz"), "rt") as f:
        E = pd.read_csv(f, sep="\t", index_col=0)
    surv = pd.read_csv(os.path.join(TCGA, "LAML_survival.txt"), sep="\t")
    cl = pd.read_csv(os.path.join(TCGA, "LAML_clinicalMatrix"), sep="\t", low_memory=False)
    cl = cl.rename(columns={cl.columns[0]: "sample"})

    surv = surv[surv["OS.time"].notna() & surv["OS"].notna()]
    keep = [s for s in E.columns if s in set(surv["sample"])]
    E = E[keep]
    surv = surv.set_index("sample").loc[keep]
    cl = cl.set_index("sample").reindex(keep)

    # one specimen per patient (all TCGA-LAML are -03 peripheral blood; guard anyway)
    pat = surv["_PATIENT"].values
    seen, sel = set(), []
    for i, p in enumerate(pat):
        if p not in seen:
            seen.add(p); sel.append(i)
    E, surv, cl = E.iloc[:, sel], surv.iloc[sel], cl.iloc[sel]

    t = surv["OS.time"].astype(float).values / 365.25          # days -> years
    e = surv["OS"].astype(float).values
    ok = np.where(t > 0)[0]                                    # a zero follow-up time is uninformative
    return E.iloc[:, ok], cl.iloc[ok], t[ok], e[ok]


def tcga_clinical_frame(cl):
    """Map TCGA field names onto the BeatAML clinical schema the blocks were fitted on."""
    risk = cl["acute_myeloid_leukemia_calgb_cytogenetics_risk_category"].map(RISK_MAP)
    return pd.DataFrame({
        "ageAtDiagnosis": pd.to_numeric(cl["age_at_initial_pathologic_diagnosis"], errors="coerce"),
        "ELN2017": risk.values,
        "consensus_sex": cl["gender"].str.capitalize().values,
        "%.Blasts.in.PB": pd.to_numeric(cl["lab_procedure_blast_cell_outcome_percentage_value"],
                                        errors="coerce").values,
        "%.Blasts.in.BM": pd.to_numeric(
            cl["lab_procedure_bone_marrow_blast_cell_outcome_percent_value"], errors="coerce").values,
        "wbcCount": pd.to_numeric(cl["lab_procedure_leukocyte_result_unspecified_value"],
                                  errors="coerce").values,
        # TCGA-LAML is initial-diagnosis de-novo AML by construction of the cohort
        "isRelapse": 0.0, "isDenovo": 1.0, "isTransformed": 0.0, "priorMDS": 0.0,
        # essentially uniform standard induction -> this block cannot be validated here
        "typeInductionTx": "Standard Chemo",
    }, index=cl.index)


def build_blocks(bundle, E, clf):
    """TCGA feature blocks in the model's own space. Only the z-reference is cohort-matched."""
    fs = bundle["feature_space"]
    s2e = bundle.get("sym2ens") or {}
    gpos = {g: i for i, g in enumerate(fs.genes)}
    n_s = E.shape[1]
    A = np.zeros((n_s, len(fs.genes)))
    matched = 0
    for sym, row in zip(E.index, E.values):
        j = gpos.get(str(sym) if str(sym) in gpos else s2e.get(str(sym)))
        if j is not None:
            A[:, j] = row; matched += 1
    frac = matched / len(fs.genes)

    # Xena HiSeqV2 is log2(RSEM+1); add_reference expects linear space, as for BeatAML
    fs.add_reference("tcga", np.power(2.0, A) - 1.0)
    Z = fs.z(np.power(2.0, A) - 1.0, "tcga")

    P = fs.pca.transform(Z[:, fs.sel])
    S, _ = fs._state_block(Z)
    cols = bundle.get("mut_columns") or []
    M = np.zeros((n_s, len(cols)))
    M[:, -1] = 1.0                                    # whole panel unreported
    C, _ = SD.clinical_block(clf)
    AE, _ = SD.age_eln_block(clf)
    age = AE[:, [0]]
    kn = bundle.get("age_knots") or [45.0, 60.0, 71.0]
    AGS = np.hstack([age] + [np.clip(age - q, 0, None) ** 3 for q in kn])
    TXB = np.tile(np.array([[1.0, 0.0, 0.0]]), (n_s, 1))
    return {"rna": P, "state": S, "mut": M, "clin": C, "age_eln": AE,
            "age_spline": AGS, "txbase": TXB}, frac, matched


def boot_c(t, e, r, B=2000, seed=0):
    rng = np.random.RandomState(seed)
    n, out = len(t), []
    for _ in range(B):
        i = rng.randint(0, n, n)
        if e[i].sum() < 5:
            continue
        out.append(CX.c_index(t[i], e[i], r[i]))
    q = np.percentile(out, [2.5, 97.5])
    return round(float(q[0]), 4), round(float(q[1]), 4)


def boot_delta(t, e, r1, r0, B=2000, seed=0):
    rng = np.random.RandomState(seed)
    n, d = len(t), []
    for _ in range(B):
        i = rng.randint(0, n, n)
        if e[i].sum() < 5:
            continue
        d.append(CX.c_index(t[i], e[i], r1[i]) - CX.c_index(t[i], e[i], r0[i]))
    d = np.asarray(d)
    return round(float(d.mean()), 4), round(float(np.percentile(d, 2.5)), 4), \
        round(float(np.percentile(d, 97.5)), 4), round(float((d <= 0).mean()), 4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fig", action="store_true")
    a = ap.parse_args()
    t0 = time.time()

    import pickle
    ensure_data()
    bundle = pickle.load(open(BUNDLE, "rb"))
    E, cl, t, e = load_tcga()
    clf = tcga_clinical_frame(cl)
    blocks, frac, matched = build_blocks(bundle, E, clf)
    print("TCGA-LAML: %d patients, %d deaths, %d censored | median follow-up %.2f y"
          % (len(t), int(e.sum()), int((1 - e).sum()), float(np.median(t))))
    print("gene overlap: %d / %d of the model's space (%.1f%%)" % (matched, len(bundle["feature_space"].genes), 100 * frac))
    print("cytogenetic risk: %s" % clf["ELN2017"].value_counts(dropna=False).to_dict())

    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "cohort": "TCGA-LAML (UCSC Xena, HiSeqV2)",
           "n": int(len(t)), "deaths": int(e.sum()), "censored": int(len(t) - e.sum()),
           "median_followup_y": round(float(np.median(t)), 3),
           "max_followup_y": round(float(t.max()), 3),
           "gene_overlap": {"matched": matched, "of": len(bundle["feature_space"].genes),
                            "fraction": round(frac, 4)},
           "frozen": ["cox coefficients", "PCA rotation", "variable-gene selection", "NNLS fusion"],
           "refit_on_tcga": ["per-gene z-reference (median/MAD)"],
           "arms": {}}

    risks = {}
    for arm, m in bundle["models"].items():
        names = bundle["arm_blocks"][arm]
        if not all(b in blocks for b in names):
            continue
        sub = {b: blocks[b] for b in names}
        try:
            r = (m.risk(sub) if hasattr(m, "stack_blocks") else m.risk(sub[names[0]])).ravel()
        except Exception as ex:
            res["arms"][arm] = {"error": str(ex)}; continue
        ok = np.isfinite(r)
        if ok.sum() < 30:
            res["arms"][arm] = {"error": "only %d scored" % ok.sum()}; continue
        risks[arm] = r
        ci = CX.c_index(t[ok], e[ok], r[ok])
        lo, hi = boot_c(t[ok], e[ok], r[ok])
        row = {"blocks": names, "n_scored": int(ok.sum()), "c_index": round(ci, 4),
               "c_index_ci95": [lo, hi]}
        for h in HORIZONS:
            au, n = CX.td_auc(t[ok], e[ok], r[ok], h)
            row["auc_%gy" % h] = None if au is None else round(au, 4)
            row["auc_%gy_n" % h] = n
        res["arms"][arm] = row

    print("\n== transferred arms (coefficients frozen on BeatAML) ==")
    print("  %-12s %8s %-16s %8s %8s %8s" % ("arm", "C-index", "95% CI", "AUC 1y", "AUC 2y", "AUC 5y"))
    for arm, row in res["arms"].items():
        if "c_index" not in row:
            print("  %-12s  %s" % (arm, row.get("error"))); continue
        print("  %-12s %8.3f  [%.3f, %.3f] %8s %8s %8s"
              % (arm, row["c_index"], row["c_index_ci95"][0], row["c_index_ci95"][1],
                 row.get("auc_1y"), row.get("auc_2y"), row.get("auc_5y")))

    # ---- the comparison that matters: molecular gain over age + cytogenetics, in TCGA ------------
    if "age_eln" in risks:
        res["gain_over_age_eln"] = {}
        print("\n== gain over age + cytogenetic risk (bootstrap over patients) ==")
        for arm in ("molecular", "full", "deployed"):
            if arm not in risks:
                continue
            m = np.isfinite(risks[arm]) & np.isfinite(risks["age_eln"])
            d, lo, hi, p = boot_delta(t[m], e[m], risks[arm][m], risks["age_eln"][m])
            res["gain_over_age_eln"][arm] = {"delta_c": d, "ci95": [lo, hi], "p_one_sided": p}
            print("  %-12s dC = %+.3f  [%+.3f, %+.3f]  P(dC<=0) = %.3f" % (arm, d, lo, hi, p))

    # ---- risk-group separation, on the arm the layer actually deploys ---------------------------
    best = "deployed" if "deployed" in risks else ("full" if "full" in risks else None)
    if best:
        r = risks[best]
        ok = np.isfinite(r)
        tt, ee, rr = t[ok], e[ok], r[ok]
        cut = np.quantile(rr, [1 / 3, 2 / 3])
        grp = np.digitize(rr, cut)
        km = {}
        print("\n== risk tertiles, %s arm ==" % best)
        for g, lab in enumerate(("low", "intermediate", "high")):
            s = grp == g
            kt, ks = CX.km(tt[s], ee[s])
            two = float(np.interp(2.0, kt, ks, left=1.0, right=ks[-1] if len(ks) else 1.0))
            med = next((float(x) for x, y in zip(kt, ks) if y <= 0.5), None)
            km[lab] = {"n": int(s.sum()), "deaths": int(ee[s].sum()),
                       "surv_2y": round(two, 4), "median_y": None if med is None else round(med, 3)}
            print("  %-13s n=%3d  deaths=%3d  2-y survival %5.1f%%  median %s"
                  % (lab, s.sum(), ee[s].sum(), 100 * two,
                     "n.r." if med is None else "%.2f y" % med))
        lh = grp != 1                                          # low vs high, dropping the middle
        chi, p = CX.logrank(tt[lh], ee[lh], (grp[lh] == 2).astype(int))
        km["logrank_low_vs_high"] = {"chi2": round(float(chi), 3), "p": float(p)}
        print("  log-rank low vs high: chi2 = %.1f, p = %.3g" % (chi, p))
        res["tertiles"] = {"arm": best, **km}

        if not a.no_fig:
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(6.2, 4.4))
                for g, (lab, col) in enumerate(zip(("low", "intermediate", "high"),
                                                   ("#2c6e49", "#b08900", "#9b2226"))):
                    s = grp == g
                    kt, ks = CX.km(tt[s], ee[s])
                    ax.step(np.concatenate([[0], kt]), np.concatenate([[1], ks]), where="post",
                            color=col, lw=1.9, label="%s risk (n=%d)" % (lab, s.sum()))
                ax.set_xlabel("years from diagnosis"); ax.set_ylabel("overall survival")
                ax.set_title("TCGA-LAML external validation — MOSAIC-AML risk tertiles\n"
                             "coefficients frozen on BeatAML2 · C-index %.3f · log-rank p = %.2g"
                             % (res["arms"][best]["c_index"], p), fontsize=9)
                ax.set_ylim(0, 1.02); ax.set_xlim(0, min(8, tt.max()))
                ax.grid(alpha=0.25, lw=0.6); ax.legend(frameon=False, fontsize=8)
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
                fig.tight_layout(); fig.savefig(FIG, dpi=170); plt.close(fig)
                print("\nwrote %s" % FIG)
            except Exception as ex:
                print("figure skipped: %s" % ex)

    json.dump(res, open(OUT, "w"), indent=1)
    print("wrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
